"""Spatial persistence diagnostics, confidence sets, and transformations."""

from collections.abc import Sequence
from dataclasses import dataclass
from importlib.resources import files

import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from patsy import EvalEnvironment
from numpy.typing import ArrayLike

from . import _native
from ._model import (
    EARTH_CIRCUMFERENCE_M,
    FloatMatrix,
    canonical_coordinates,
    coefficient_bread,
    coordinate_columns,
    coordinates_from_frame,
    factor_codes,
    numeric_matrix,
)
from .scpc import SCPCResult, integer, native_info, scpc

TESTS = ("i0", "i1", "i0resid", "i1resid")


@dataclass(frozen=True)
class Simulation:
    q: int
    requested_draws: int
    inference_draws: int
    tuning_draws: int
    seed: int | None


@dataclass(frozen=True)
class TestResult:
    test: str
    lr: float
    pvalue: float
    cvalue: FloatMatrix
    ha_param: float
    best_tuning_power: float | None
    tied: bool
    all_directions: bool
    simulation: Simulation
    observations: pd.Index
    coordinate_units: str
    normalization: str
    execution: dict
    null_grid: FloatMatrix | None = None
    null_pvalues: FloatMatrix | None = None

    @property
    def backend(self) -> str:
        return "rust"


@dataclass(frozen=True)
class HalfLifeResult:
    ci_lower: float
    ci_upper: float
    max_dist: float
    level: float
    normdist: bool
    distance_units: str
    pv_mat: FloatMatrix
    hl_grid: FloatMatrix
    components: FloatMatrix
    tied: bool
    all_directions: bool
    simulation: Simulation
    observations: pd.Index
    execution: dict

    @property
    def backend(self) -> str:
        return "rust"

    def __str__(self) -> str:
        units = "sample diameters" if self.normdist else self.distance_units
        return f"{self.level:g}% spatial half-life confidence set ({units})\n{self.components}"


@dataclass(frozen=True)
class PipelineResult:
    tests: dict[str, TestResult]
    levels: SCPCResult
    transformed: SCPCResult
    transformed_data: pd.DataFrame

    def summary(self) -> str:
        diagnostics = pd.Series({name: result.pvalue for name, result in self.tests.items()})
        return (
            f"SPUR diagnostic p-values\n{diagnostics.to_string()}\n\n"
            f"Levels\n{self.levels.summary()}\n\nTransformed\n{self.transformed.summary()}"
        )


def _draw_blocks(q, nrep, seed, draws, tuning_draws, *, tune):
    q = integer(q, "q", 2)
    if draws is None:
        if tuning_draws is not None:
            raise ValueError("Supply both inference and tuning draws for a diagnostic")
        nrep = 100_000 if nrep is None else integer(nrep, "nrep", 1)
        seed = integer(seed, "seed", 0)
        rng = np.random.default_rng(seed)
        block = rng.standard_normal((q * (2 if tune else 1), nrep))
    else:
        inference = numeric_matrix(draws, "Inference draws")
        if inference.shape[0] != q:
            raise ValueError("Inference draws must have shape (q, nrep)")
        if nrep is not None and integer(nrep, "nrep", 1) != inference.shape[1]:
            raise ValueError("nrep does not match the supplied draws")
        nrep = inference.shape[1]
        seed = None
        if tune:
            if tuning_draws is None:
                raise ValueError("Diagnostics require a separate tuning draw block")
            tuning = numeric_matrix(tuning_draws, "Tuning draws")
            if tuning.shape != inference.shape:
                raise ValueError("Inference and tuning blocks must have the same shape")
            block = np.vstack((inference, tuning))
        elif tuning_draws is not None:
            raise ValueError("Half-life uses one supplied draw block")
        else:
            block = inference
    return block, Simulation(q, nrep, nrep, nrep if tune else 0, seed)


def _coordinates(frame, lon, lat, coords_euclidean):
    columns = coordinate_columns(lon, lat, coords_euclidean)
    geographic = lon is not None
    return canonical_coordinates(coordinates_from_frame(frame, columns), geographic), geographic


def _diagnostic(test, native, sim, rows, geographic, execution):
    return TestResult(
        test,
        float(native["lr"][0]),
        float(native["pvalue"][0]),
        native["cvalue"],
        native["ha_param"],
        native["best_tuning_power"],
        native["tied"],
        native["all_directions"],
        sim,
        rows.copy(),
        "meters" if geographic else "coordinate units",
        "sample diameter",
        execution,
        native.get("rho_grid"),
        native.get("pvalue_mat"),
    )


def _response(data, column):
    if column not in data or not pd.api.types.is_numeric_dtype(data[column]):
        raise ValueError(f"Response {column!r} must name a numeric column")
    values = data[column].to_numpy(dtype=float, na_value=np.nan)
    if not np.isfinite(values).all():
        raise ValueError(f"Response {column!r} contains missing or non-finite values")
    return values


def _formula_fit(formula, data, environment):
    if not isinstance(data, pd.DataFrame) or not data.columns.is_unique:
        raise ValueError("data must be a pandas DataFrame with unique columns")
    if not data.index.is_unique:
        raise ValueError("Use a unique index to identify fitted regression rows")
    base = smf.ols(formula, data, eval_env=environment)
    sample = data.loc[base.data.row_labels]
    y = np.asarray(base.endog, dtype=float)
    x = numeric_matrix(base.exog, "Fitted SPUR design")
    if y.ndim != 1 or not np.isfinite(y).all():
        raise ValueError("SPUR requires one finite numeric response")
    coefficient_bread(x, tuple(base.exog_names))
    return base, sample, y, x


def spurtest(
    formula: str,
    data: pd.DataFrame,
    *,
    test: str,
    lon: str | None = None,
    lat: str | None = None,
    coords_euclidean: Sequence[str] | None = None,
    q: int = 15,
    nrep: int | None = None,
    seed: int = 42,
    draws: ArrayLike | None = None,
    tuning_draws: ArrayLike | None = None,
    backend: str | None = None,
    nthreads: int | None = None,
    _environment=None,
) -> TestResult:
    if test not in TESTS:
        raise ValueError("test must be i0, i1, i0resid, or i1resid")
    execution = native_info(backend, nthreads)
    if test.endswith("resid"):
        environment = _environment or EvalEnvironment.capture(1)
        _, sample, y, x = _formula_fit(formula, data, environment)
        model_inputs = (y, x)
    else:
        sample = data
        column, _, rhs = formula.partition("~")
        if rhs.strip() not in ("", "1"):
            raise ValueError("Levels diagnostics take one response column")
        model_inputs = (_response(data, column.strip()),)
    coords, geographic = _coordinates(sample, lon, lat, coords_euclidean)
    if integer(q, "q", 2) >= len(coords):
        raise ValueError("q must be smaller than the fitted sample size")
    block, sim = _draw_blocks(q, nrep, seed, draws, tuning_draws, tune=True)
    compute = {
        "i0": _native.run_spur_i0,
        "i1": _native.run_spur_i1,
        "i0resid": _native.run_spur_i0_resid,
        "i1resid": _native.run_spur_i1_resid,
    }[test]
    result = compute(
        coords,
        *model_inputs,
        block,
        sim.q,
        sim.inference_draws,
        latlong=geographic,
        nthreads=execution["nthreads_requested"],
    )
    return _diagnostic(test, result, sim, sample.index, geographic, execution)


def spurtest_i0(var: str, data: pd.DataFrame, **kwargs) -> TestResult:
    return spurtest(var, data, test="i0", **kwargs)


def spurtest_i1(var: str, data: pd.DataFrame, **kwargs) -> TestResult:
    return spurtest(var, data, test="i1", **kwargs)


def spurtest_i0resid(formula: str, data: pd.DataFrame, **kwargs) -> TestResult:
    return spurtest(
        formula, data, test="i0resid", _environment=EvalEnvironment.capture(1), **kwargs
    )


def spurtest_i1resid(formula: str, data: pd.DataFrame, **kwargs) -> TestResult:
    return spurtest(
        formula, data, test="i1resid", _environment=EvalEnvironment.capture(1), **kwargs
    )


def spurhalflife(
    var: str,
    data: pd.DataFrame,
    *,
    lon: str | None = None,
    lat: str | None = None,
    coords_euclidean: Sequence[str] | None = None,
    q: int = 15,
    nrep: int | None = None,
    level: float = 95,
    normdist: bool = False,
    seed: int = 42,
    draws: ArrayLike | None = None,
    backend: str | None = None,
    nthreads: int | None = None,
) -> HalfLifeResult:
    execution = native_info(backend, nthreads)
    q = integer(q, "q", 2, 100)
    if not 0 < level < 100:
        raise ValueError("level must lie strictly between 0 and 100")
    coords, geographic = _coordinates(data, lon, lat, coords_euclidean)
    if q >= len(coords):
        raise ValueError("q must be smaller than the fitted sample size")
    block, sim = _draw_blocks(q, nrep, seed, draws, None, tune=False)
    out = _native.run_spur_halflife(
        coords,
        _response(data, var),
        block,
        q,
        sim.inference_draws,
        level=level / 100,
        latlong=geographic,
        nthreads=execution["nthreads_requested"],
    )
    diameter = out["max_dist"] * (EARTH_CIRCUMFERENCE_M if geographic else 1)
    factor = 1 if normdist else diameter
    return HalfLifeResult(
        out["ci_lower"] * factor,
        out["ci_upper"] * factor,
        diameter,
        level,
        normdist,
        "meters" if geographic else "coordinate units",
        out["pv_mat"],
        out["hl_grid_ho"] * factor,
        np.asarray(out["components"], dtype=float).reshape(-1, 2) * factor,
        out["tied"],
        out["all_directions"],
        sim,
        data.index.copy(),
        execution,
    )


def spurtransform(
    formula: str | Sequence[str],
    data: pd.DataFrame,
    *,
    lon: str | None = None,
    lat: str | None = None,
    coords_euclidean: Sequence[str] | None = None,
    prefix: str = "h_",
    transformation: str = "lbmgls",
    radius: float | None = None,
    clustvar: str | None = None,
    replace: bool = False,
    backend: str | None = None,
    nthreads: int | None = None,
) -> pd.DataFrame:
    execution = native_info(backend, nthreads)
    columns = (
        tuple(
            term.strip()
            for term in formula.replace("~", "+").split("+")
            if term.strip() not in ("", "1", "0")
        )
        if isinstance(formula, str)
        else tuple(formula)
    )
    if not columns or len(set(columns)) != len(columns) or not prefix:
        raise ValueError("Supply distinct response columns and a nonempty output prefix")
    outputs = [prefix + name for name in columns]
    if not replace and any(name in data for name in outputs):
        raise ValueError("The requested output columns already exist")
    if transformation not in {"lbmgls", "nn", "iso", "cluster"}:
        raise ValueError("transformation must be lbmgls, nn, iso, or cluster")
    if (transformation == "iso") != (radius is not None):
        raise ValueError("Supply radius only for the ISO transformation")
    if (transformation == "cluster") != (clustvar is not None):
        raise ValueError("Supply clustvar only for the cluster transformation")
    response = np.column_stack([_response(data, name) for name in columns])
    if clustvar is None:
        coords, geographic = _coordinates(data, lon, lat, coords_euclidean)
        groups = None
    else:
        coords, geographic = None, False
        groups, _ = factor_codes(data[clustvar], "Spatial cluster")
    if radius is not None and geographic:
        radius /= EARTH_CIRCUMFERENCE_M
    result = _native.run_spur_transform(
        transformation,
        len(data),
        response,
        coords=coords,
        latlong=geographic,
        radius=radius,
        cluster_id=groups,
        nthreads=execution["nthreads_requested"],
    )["transformed"]
    return data.assign(**dict(zip(outputs, result.T, strict=True)))


def spur(
    formula: str,
    data: pd.DataFrame,
    *,
    lon: str | None = None,
    lat: str | None = None,
    coords_euclidean: Sequence[str] | None = None,
    q: int = 15,
    nrep: int | None = None,
    seed: int = 42,
    draws: ArrayLike | None = None,
    tuning_draws: ArrayLike | None = None,
    avc: float = 0.03,
    uncond: bool = False,
    cvs: bool = False,
    method: str = "auto",
    large_n_seed: int = 1,
    backend: str | None = None,
    nthreads: int | None = None,
) -> PipelineResult:
    execution = native_info(backend, nthreads)
    base, sample, y, x = _formula_fit(formula, data, EvalEnvironment.capture(1))
    terms = base.data.design_info.term_name_slices
    if "Intercept" not in terms:
        raise ValueError("The SPUR pipeline requires a fitted intercept")
    intercept = terms["Intercept"].start
    regressors = [i for i in range(x.shape[1]) if i != intercept]
    coords, geographic = _coordinates(sample, lon, lat, coords_euclidean)
    if integer(q, "q", 2) > len(coords) - x.shape[1]:
        raise ValueError("q exceeds the residual degrees of freedom")
    block, sim = _draw_blocks(q, nrep, seed, draws, tuning_draws, tune=True)
    out = _native.run_spur_pipeline(
        coords,
        y,
        x,
        block,
        sim.q,
        sim.inference_draws,
        np.column_stack((y, x[:, regressors])),
        latlong=geographic,
        nthreads=execution["nthreads_requested"],
    )
    tests = {
        test: _diagnostic(test, out[test], sim, sample.index, geographic, execution)
        for test in TESTS
    }
    raw = out["transformed_vars"]
    transformed_design = x.copy()
    transformed_design[:, regressors] = raw[:, 1:]
    column_names = base.exog_names
    transformed_fit = sm.OLS(
        pd.Series(raw[:, 0], index=sample.index),
        pd.DataFrame(transformed_design, index=sample.index, columns=column_names),
    ).fit()
    options = dict(
        coords=coords,
        geographic=geographic,
        avc=avc,
        uncond=uncond,
        cvs=cvs,
        method=method,
        large_n_seed=large_n_seed,
        backend=backend,
        nthreads=nthreads,
    )
    levels = scpc(base.fit(), **options)
    transformed = scpc(transformed_fit, **options)
    names = [base.endog_names, *(column_names[i] for i in regressors)]
    transformed_data = sample.copy()
    for name, values in zip(names, raw.T, strict=True):
        column = f"h_{name}"
        while column in transformed_data:
            column = "h_" + column
        transformed_data[column] = values
    return PipelineResult(tests, levels, transformed, transformed_data)


def load_chetty_data() -> pd.DataFrame:
    """Load the public Chetty et al. commuting-zone example."""
    with files("spatial_spur").joinpath("data/chetty.csv").open("rb") as source:
        return pd.read_csv(source)


def standardize(data: pd.DataFrame, columns: Sequence[str]) -> pd.DataFrame:
    """Center and scale columns using their sample standard deviations."""
    values = data[list(columns)]
    return data.assign(**((values - values.mean()) / values.std()).to_dict("series"))
