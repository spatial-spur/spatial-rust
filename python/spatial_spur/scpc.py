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

    def summary(self) -> str:
        table = self.table.assign(critical_95=self.critical_95)
        return (
            f"SCPC {self.inference}, {self.method}, q={self.q}, "
            f"{self.n_obs} observations at {self.n_sites} sites\n{table.to_string()}"
        )

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
        columns = coordinate_columns(lon, lat, coords_euclidean)
        keys = (*columns, cluster) if isinstance(cluster, str) else columns
        sample, obs = align_frame(fit, data, row_id, keys)
        coordinates = coordinates_from_frame(sample, columns)
        geographic = geographic or lon is not None
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
        landmarks = integer(landmarks, "landmarks", 2)
        landmark_samples = integer(landmark_samples, "landmark_samples", 1)
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
