"""Extract fitted regression state and map it to spatial sites."""

from collections.abc import Callable, Sequence
from dataclasses import dataclass
from typing import Any

import numpy as np
import pandas as pd
from numpy.typing import NDArray
from statsmodels.regression.linear_model import OLS

from . import _native

FloatMatrix = NDArray[np.float64]
EARTH_CIRCUMFERENCE_M = _native.build_info()["earth_circumference_m"]
RANK_TOLERANCE = 1e-10


@dataclass(frozen=True)
class Fit:
    names: tuple[str, ...]
    estimates: FloatMatrix
    design: FloatMatrix
    instruments: FloatMatrix | None
    scores: FloatMatrix
    bread: FloatMatrix
    labels: pd.Index | None
    frame: pd.DataFrame | None
    host: str
    fe_id: NDArray[np.int64] | None = None
    residualize: Callable[[FloatMatrix], FloatMatrix] | None = None


@dataclass(frozen=True)
class Sites:
    coordinates: FloatMatrix
    ids: pd.Index
    membership: NDArray[np.int64] | None
    geographic: bool
    units: str


def numeric_matrix(values: Any, name: str, *, width: int | None = None) -> FloatMatrix:
    matrix = np.asarray(values, dtype=float)
    if matrix.ndim != 2 or not matrix.shape[1] or (width is not None and matrix.shape[1] != width):
        raise ValueError(f"{name} must be a matrix of the required width")
    if not np.isfinite(matrix).all():
        raise ValueError(f"{name} contains missing or non-finite values")
    return matrix


def factor_codes(values: Sequence, name: str) -> tuple[NDArray[np.int64], pd.Index]:
    if np.ndim(values) != 1:
        raise ValueError(f"{name} must be a one-dimensional vector of site labels")
    codes, levels = pd.factorize(values if isinstance(values, pd.Series) else pd.Index(values))
    if (codes < 0).any():
        raise ValueError(f"{name} contains missing labels")
    return codes.astype(np.int64, copy=False), pd.Index(levels)


def scaled_svd(x: FloatMatrix):
    scale = np.max(np.abs(x), axis=0)
    safe = np.where(scale > 0, scale, 1)
    u, s, vt = np.linalg.svd(x / safe, full_matrices=x.shape[0] < x.shape[1])
    rank = np.count_nonzero(s > s[0] * RANK_TOLERANCE)
    return u, s, vt, safe, rank


def coefficient_bread(x: FloatMatrix, names: tuple[str, ...]) -> FloatMatrix:
    _, s, vt, scale, rank = scaled_svd(x)
    if rank < x.shape[1]:
        nonidentified = np.any(np.abs(vt[rank:]) > 1e-7, axis=0)
        aliases = ", ".join(name for name, bad in zip(names, nonidentified, strict=True) if bad)
        raise ValueError(f"Fitted coefficients are not separately identified: {aliases}")
    factor = vt.T / s / scale[:, None]
    return factor @ factor.T


def fitted_directions(x: FloatMatrix, z: FloatMatrix | None) -> FloatMatrix:
    if z is None:
        return x
    u, _, _, _, rank = scaled_svd(z)
    qz = u[:, :rank]
    return qz @ (qz.T @ x)


def _extract_statsmodels(model) -> Fit:
    base = model.model
    names = tuple(base.exog_names)
    x = numeric_matrix(base.exog, "Fitted design", width=len(names))
    residuals = np.asarray(model.resid, dtype=float)
    data = base.data
    frame = getattr(data, "frame", None)
    labels = pd.Index(data.row_labels) if data.row_labels is not None else None
    if isinstance(frame, pd.DataFrame) and labels is not None and frame.index.is_unique:
        frame = frame.loc[labels]
    else:
        frame = None
    return Fit(
        names,
        np.asarray(model.params, dtype=float),
        x,
        None,
        x * residuals[:, None],
        coefficient_bread(x, names),
        labels,
        frame,
        "statsmodels",
    )


def _extract_pyfixest(model) -> Fit:
    # PyFixest stores the exact within-transformed second-stage matrices. Its
    # separately exposed first-stage fit can be evaluated on a different sample.
    if getattr(model, "_method", None) != "feols" or type(model).__name__ not in {"Feols", "Feiv"}:
        raise ValueError("SCPC accepts one unweighted PyFixest feols() or ordinary 2SLS fit")
    if getattr(model, "_weights_name", None) is not None:
        raise ValueError("SCPC requires unweighted regression")
    if getattr(model, "_lean", False) or not getattr(model, "_store_data", False):
        raise ValueError("Fit PyFixest with lean=False and store_data=True")
    names = tuple(map(str, model._coefnames))
    x = numeric_matrix(model._X, "Stored PyFixest design", width=len(names))
    z = numeric_matrix(model._Z, "Stored PyFixest instruments") if model._is_iv else None
    direction = fitted_directions(x, z)
    bread = coefficient_bread(direction, names)
    residuals = np.asarray(model._u_hat, dtype=float)
    frame = model._data
    if not isinstance(frame, pd.DataFrame) or len(frame) != len(x):
        raise ValueError("Stored PyFixest estimation rows do not match the fit")
    fe_id = None
    residualize = None
    if model._has_fixef:
        groups = model._fe
        if groups.shape[1] == 1:
            fe_id, _ = factor_codes(groups.iloc[:, 0], "Fitted fixed effect")
        else:
            fe = np.ascontiguousarray(groups.to_numpy(), dtype=np.uintp)
            weights = np.ones(len(x))
            demean = model._demean_func
            tol = min(1e-10, model._fixef_tol)
            maxiter = model._fixef_maxiter

            def residualize(values):
                output, success = demean(
                    x=values,
                    flist=fe,
                    weights=weights,
                    tol=tol,
                    maxiter=maxiter,
                )
                if not success:
                    raise ValueError(
                        "PyFixest did not converge while projecting spatial directions"
                    )
                return output

    return Fit(
        names,
        np.asarray(model._beta_hat).ravel(),
        x,
        z,
        direction * residuals[:, None],
        bread,
        None,
        frame,
        "pyfixest",
        fe_id,
        residualize,
    )


def extract_fit(model) -> Fit:
    if isinstance(getattr(model, "model", None), OLS) and hasattr(model, "resid"):
        return _extract_statsmodels(model)
    if type(model).__module__.startswith("pyfixest."):
        return _extract_pyfixest(model)
    raise ValueError("SCPC supports unweighted Statsmodels OLS and PyFixest OLS or ordinary 2SLS")


def _match_keys(fitted: pd.Index, supplied: pd.Index, name: str) -> NDArray[np.intp]:
    if not fitted.is_unique or not supplied.is_unique or fitted.hasnans or supplied.hasnans:
        raise ValueError(f"{name} must uniquely identify fitted and supplied rows")
    positions = supplied.get_indexer(fitted)
    if (positions < 0).any():
        raise ValueError(f"{name} does not identify every fitted row in the supplied data")
    return positions


def _same_rows(a: pd.DataFrame, b: pd.DataFrame) -> bool:
    return a.reset_index(drop=True).equals(b.reset_index(drop=True))


def align_frame(
    fit: Fit,
    data: pd.DataFrame | None,
    row_id: str | None,
    spatial_columns: Sequence[str],
) -> tuple[pd.DataFrame, pd.Index]:
    if data is None:
        if fit.frame is None:
            raise ValueError(
                "Supply a frame with fitted row identifiers or coordinates in fitted order"
            )
        labels = fit.labels if fit.labels is not None else fit.frame.index
        if row_id is not None:
            labels = pd.Index(fit.frame[row_id])
            _match_keys(labels, labels, row_id)
        return fit.frame, labels

    if not isinstance(data, pd.DataFrame) or not data.columns.is_unique:
        raise ValueError("data must be a pandas DataFrame with unique column names")
    stored = fit.frame
    common = (
        []
        if stored is None
        else [c for c in stored.columns if c in data and c not in spatial_columns]
    )
    if row_id is not None:
        if stored is None or row_id not in stored or row_id not in data:
            raise ValueError(
                f"Fitted and supplied data must both contain row identifier {row_id!r}"
            )
        labels = pd.Index(stored[row_id])
        positions = _match_keys(labels, pd.Index(data[row_id]), row_id)
    elif fit.labels is not None:
        labels = fit.labels
        positions = _match_keys(labels, data.index, "Data index")
    else:
        # PyFixest replaces the input index with positional labels. Compare raw
        # fitted rows by value or use an explicit stable identifier column.
        if not common or stored is None:
            raise ValueError("PyFixest needs shared fitted columns or a stable row_id column")
        keys = [
            pd.Index(pd.util.hash_pandas_object(frame[common], index=False))
            for frame in (stored, data)
        ]
        positions = _match_keys(*keys, "Shared PyFixest rows (supply row_id if ambiguous)")
        labels = data.index.take(positions)
        if not labels.is_unique:
            raise ValueError("Supply unique data index labels or row_id for the fitted rows")
    sample = data.iloc[positions]
    if common and stored is not None and not _same_rows(stored[common], sample[common]):
        raise ValueError("Supplied data values disagree with the fitted observations")
    return sample, labels


def coordinate_columns(lon, lat, coords_euclidean) -> tuple[str, ...]:
    if coords_euclidean is not None:
        if lon is not None or lat is not None:
            raise ValueError("Choose geographic or Euclidean coordinates")
        if isinstance(coords_euclidean, str) or len(coords_euclidean) == 0:
            raise ValueError("coords_euclidean must be a sequence of column names")
        return tuple(coords_euclidean)
    if lon is None or lat is None or lon == lat:
        raise ValueError("Supply both lon and lat, or coords_euclidean")
    return (lon, lat)


def coordinates_from_frame(frame, columns) -> FloatMatrix:
    axes = []
    for column in columns:
        if column not in frame:
            raise ValueError(f"Coordinate column {column!r} is missing")
        if not pd.api.types.is_numeric_dtype(frame[column]):
            raise ValueError(f"Coordinate column {column!r} must be numeric")
        values = frame[column].to_numpy(dtype=float, na_value=np.nan)
        if not np.isfinite(values).all():
            raise ValueError(f"Coordinate column {column!r} contains missing or non-finite values")
        axes.append(values)
    return np.column_stack(axes)


def canonical_coordinates(array: FloatMatrix, geographic: bool) -> FloatMatrix:
    if not geographic:
        return array
    coordinates = numeric_matrix(array, "Geographic coordinates", width=2).copy()
    lon, lat = coordinates.T
    if np.any(np.abs(lat) > 90):
        raise ValueError("Latitude must be between -90 and 90 degrees")
    lon[:] = np.fmod(lon, 360)
    lon[lon >= 180] -= 360
    lon[lon < -180] += 360
    lon[np.abs(lat) == 90] = 0
    return coordinates


def spatial_sites(
    coordinates: FloatMatrix,
    obs: pd.Index,
    cluster: Sequence | None,
    geographic: bool,
    units: str | None = None,
) -> Sites:
    if len(coordinates) != len(obs):
        raise ValueError("Coordinates must have one row per fitted observation")
    coordinates = canonical_coordinates(coordinates, geographic)
    units = units or ("fraction of Earth circumference" if geographic else "coordinate units")
    if cluster is None:
        return Sites(coordinates, obs, None, geographic, units)
    if len(cluster) != len(obs):
        raise ValueError("Cluster labels must have one entry per fitted observation")
    codes, levels = factor_codes(cluster, "Spatial cluster")
    _, first = np.unique(codes, return_index=True)
    sites = coordinates[first]
    if np.any(coordinates != sites[codes]):
        raise ValueError("Coordinates must be constant within each spatial cluster")
    return Sites(sites, levels, codes, geographic, units)
