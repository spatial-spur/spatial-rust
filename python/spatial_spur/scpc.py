"""Inference on a fitted OLS or ordinary 2SLS regression."""

from collections.abc import Sequence
from dataclasses import dataclass
from operator import index
from typing import Any

import numpy as np
import pandas as pd
from numpy.typing import ArrayLike

from . import _native
from ._model import (
    FloatMatrix,
    align_frame,
    coordinate_columns,
    coordinates_from_frame,
    extract_fit,
    numeric_matrix,
    spatial_sites,
)

STAT_COLUMNS = ("coef", "se", "t", "p", "ci_lower", "ci_upper")
CV_COLUMNS = ("32%", "10%", "5%", "1%")
CV_LEVELS = {0.68: 0, 0.90: 1, 0.95: 2, 0.99: 3}


def integer(value: Any, name: str, minimum: int, maximum: int | None = None) -> int:
    try:
        if isinstance(value, (bool, np.bool_)):
            raise TypeError
        number = index(value)
    except TypeError as exc:
        raise ValueError(f"{name} must be an integer") from exc
    if number < minimum or (maximum is not None and number > maximum):
        raise ValueError(f"{name} is outside its supported range")
    return number


def real_float(
    value: Any,
    name: str,
    *,
    gt: float | None = None,
    lt: float | None = None,
) -> float:
    if isinstance(value, (bool, np.bool_, complex, np.complexfloating, str, bytes)):
        raise ValueError(f"{name} must be a finite real number")
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be a finite real number") from exc
    if not np.isfinite(number):
        raise ValueError(f"{name} must be finite")
    if gt is not None and not (number > gt):
        raise ValueError(f"{name} must be greater than {gt}")
    if lt is not None and not (number < lt):
        raise ValueError(f"{name} must be less than {lt}")
    return number


def _resolve_parm_indices(
    coef_names: tuple[str, ...],
    parm: str | int | Sequence[str] | Sequence[int] | None,
) -> list[int]:
    if parm is None:
        return list(range(len(coef_names)))
    if isinstance(parm, str):
        if parm not in coef_names:
            raise ValueError(f"Unknown coefficient {parm!r}")
        return [coef_names.index(parm)]
    if isinstance(parm, (int, np.integer)) and not isinstance(parm, (bool, np.bool_)):
        idx = int(parm)
        if not 0 <= idx < len(coef_names):
            raise ValueError(f"Coefficient index {idx} out of range")
        return [idx]
    values = list(parm)
    if not values:
        return []
    if all(isinstance(v, str) for v in values):
        out = []
        for name in values:
            if name not in coef_names:
                raise ValueError(f"Unknown coefficient {name!r}")
            out.append(coef_names.index(name))
        return out
    out = []
    for v in values:
        idx = integer(v, "parm", 0, len(coef_names) - 1)
        out.append(idx)
    return out


def native_info(backend: str | None, nthreads: int | None) -> dict:
    if backend not in (None, "auto", "rust"):
        raise ValueError("backend must be 'rust' or 'auto'")
    info = dict(_native.build_info())
    info["nthreads_requested"] = None if nthreads is None else integer(nthreads, "nthreads", 1)
    return info


@dataclass(frozen=True)
class SCPCResult:
    scpcstats: FloatMatrix
    scpccvs: FloatMatrix | None
    w: FloatMatrix
    coef_names: tuple[str, ...]
    critical_95: FloatMatrix
    cv: float
    q: int
    c0: float
    cmax: float
    avc: float
    inference: str
    estimator: str
    method: str
    large_n_seed: int | None
    landmarks: int | None
    landmark_samples: int | None
    cv_tolerance: float
    observations: pd.Index
    sites: pd.Index
    site_membership: tuple[int, ...]
    coordinates: FloatMatrix
    coordinate_units: str
    distance_scale: float
    tied: bool
    execution: dict

    @property
    def backend(self) -> str:
        return "rust"

    @property
    def n_obs(self) -> int:
        return len(self.observations)

    @property
    def n_sites(self) -> int:
        return len(self.sites)

    @property
    def diagnostics(self) -> dict[str, bool]:
        return {"tied_eigenvalues": self.tied}

    @property
    def table(self) -> pd.DataFrame:
        return pd.DataFrame(self.scpcstats, index=self.coef_names, columns=STAT_COLUMNS)

    def __repr__(self) -> str:
        return (
            f"SCPCResult(ncoef={len(self.coef_names)}, q={self.q}, "
            f"avc={self.avc!r}, method={self.method!r}, inference={self.inference!r})"
        )

    def coef(self) -> pd.Series:
        return pd.Series(self.scpcstats[:, 0], index=self.coef_names, name="coef")

    def confint(
        self,
        parm: str | int | Sequence[str] | Sequence[int] | None = None,
        level: float = 0.95,
    ) -> pd.DataFrame:
        level = real_float(level, "level", gt=0.0, lt=1.0)
        idx = _resolve_parm_indices(self.coef_names, parm)
        names = [self.coef_names[i] for i in idx]
        if level == 0.95:
            values = self.scpcstats[idx, 4:6]
        else:
            if self.scpccvs is None or level not in CV_LEVELS:
                raise ValueError(f"Confidence level {level} is not available (pass cvs=True for 68%/90%/99%)")
            cv_vals = self.scpccvs[idx, CV_LEVELS[level]]
            coef_vals = self.scpcstats[idx, 0]
            se_vals = self.scpcstats[idx, 1]
            values = np.column_stack((coef_vals - cv_vals * se_vals, coef_vals + cv_vals * se_vals))
        lower = 100 * (1 - level) / 2
        upper = 100 * (1 + level) / 2
        return pd.DataFrame(values, index=names, columns=[f"{lower:g} %", f"{upper:g} %"])

    def summary(self) -> str:
        table = self.table.assign(critical_95=self.critical_95)
        provider = self.execution.get("linalg_provider", "unknown")
        details = [
            f"SCPC {self.inference}, {self.method}",
            f"q={self.q}",
            f"avc={self.avc:g}",
            f"{self.n_obs} observations at {self.n_sites} sites",
            f"backend=rust ({provider})",
        ]
        if self.method == "nystrom":
            details.append(
                f"seed={self.large_n_seed}, landmarks={self.landmarks}, samples={self.landmark_samples}"
            )
        lines = [", ".join(details), table.to_string()]
        if self.scpccvs is not None:
            cvs = pd.DataFrame(self.scpccvs, index=self.coef_names, columns=CV_COLUMNS)
            lines.extend(["", "Two-sided critical values:", cvs.to_string()])
        if self.tied:
            lines.extend(["", "Note: Leading spatial eigenvalues exhibited a numerical tie."])
        return "\n".join(lines)

    def __str__(self) -> str:
        return self.summary()


def scpc(
    model,
    data: pd.DataFrame | None = None,
    *,
    lon: str | None = None,
    lat: str | None = None,
    coords_euclidean: Sequence[str] | None = None,
    coords: ArrayLike | None = None,
    geographic: bool = False,
    units: str | None = None,
    cluster: str | Sequence | None = None,
    row_id: str | None = None,
    ncoef: int | None = None,
    avc: float = 0.03,
    inference: str | None = None,
    uncond: bool = False,
    method: str = "auto",
    large_n_seed: int = 1,
    landmarks: int = 1000,
    landmark_samples: int = 20,
    cv_tolerance: float = 0.001,
    cvs: bool = False,
    backend: str | None = None,
    nthreads: int | None = None,
) -> SCPCResult:
    """Run SCPC with coordinates joined to fitted rows or supplied in fitted order.

    Set cluster to one actual spatial site identifier when clustering observations.
    The returned weight matrix is in result.sites order and the zero-based
    site_membership maps each fitted observation to a row of that matrix.
    """
    execution = native_info(backend, nthreads)
    avc = real_float(avc, "avc", gt=0.001, lt=0.99)
    cv_tolerance = real_float(cv_tolerance, "cv_tolerance", gt=0.0)
    if (
        not isinstance(uncond, (bool, np.bool_))
        or not isinstance(cvs, (bool, np.bool_))
        or not isinstance(geographic, (bool, np.bool_))
    ):
        raise ValueError("uncond, cvs, and geographic must be booleans")
    uncond = bool(uncond)
    cvs = bool(cvs)
    geographic = bool(geographic)
    landmarks = integer(landmarks, "landmarks", 2)
    landmark_samples = integer(landmark_samples, "landmark_samples", 1)
    fit = extract_fit(model)
    estimator = "ols" if fit.instruments is None else "2sls"
    procedure = inference or ("unconditional" if uncond else f"conditional_{estimator}")
    if procedure not in {"unconditional", f"conditional_{estimator}"} or (
        uncond and procedure != "unconditional"
    ):
        raise ValueError(
            f"Inference {procedure!r} is incompatible with the fitted {estimator} model"
        )

    if coords is None:
        if coords_euclidean is not None and geographic:
            raise ValueError("Do not set geographic=True with Euclidean coordinates")
        columns = coordinate_columns(lon, lat, coords_euclidean)
        keys = (*columns, cluster) if isinstance(cluster, str) else columns
        sample, obs = align_frame(fit, data, row_id, keys)
        coordinates = coordinates_from_frame(sample, columns)
        geographic = lon is not None
        if isinstance(cluster, str):
            if cluster not in sample:
                raise ValueError(f"Spatial cluster column {cluster!r} is missing")
            cluster = sample[cluster]
    else:
        if any(value is not None for value in (lon, lat, coords_euclidean, row_id, data)):
            raise ValueError(
                "With explicit coords, supply coordinates and labels in fitted row order"
            )
        if isinstance(cluster, str):
            raise ValueError("With explicit coords, supply cluster as fitted-order labels")
        coordinates = numeric_matrix(coords, "Coordinates")
        obs = fit.labels if fit.labels is not None else pd.RangeIndex(len(fit.design))
    sites = spatial_sites(coordinates, obs, cluster, geographic, units)
    if method == "approx":
        method = "nystrom"
    if method not in {"auto", "exact", "nystrom"}:
        raise ValueError("method must be 'auto', 'exact', or 'nystrom'")
    if method == "auto":
        method = "nystrom" if len(sites.ids) >= 4500 else "exact"
    if method == "nystrom":
        seed = integer(large_n_seed, "large_n_seed", 0, 2**32 - 1)
    else:
        seed = None
    conditional = procedure != "unconditional"
    k = len(fit.names) if ncoef is None else integer(ncoef, "ncoef", 1, len(fit.names))
    out = _native.run_model_scpc(
        sites.coordinates,
        fit.estimates[:k],
        fit.bread[:k],
        fit.scores,
        avc,
        inference=procedure,
        method=method,
        design=fit.design if conditional else None,
        instruments=fit.instruments if conditional else None,
        cluster_id=sites.membership,
        fe_id=fit.fe_id if conditional else None,
        residualizer=fit.residualize if conditional else None,
        latlong=sites.geographic,
        cvs=cvs,
        nthreads=execution["nthreads_requested"],
        seed=seed,
        landmarks=landmarks,
        samples=landmark_samples,
        cv_tolerance=cv_tolerance,
    )
    return SCPCResult(
        out["scpcstats"],
        out["scpccvs"],
        out["wfin"],
        fit.names[:k],
        out["critical_95"],
        out["cv"],
        out["q"],
        out["c0"],
        out["cmax"],
        avc,
        procedure,
        estimator,
        method,
        seed,
        min(landmarks, len(sites.ids)) if seed is not None else None,
        landmark_samples if seed is not None else None,
        cv_tolerance,
        obs.copy(),
        sites.ids.copy(),
        tuple(range(len(obs))) if sites.membership is None else tuple(map(int, sites.membership)),
        sites.coordinates.copy(),
        sites.units,
        out["distance_scale"],
        out["tied"],
        execution,
    )
