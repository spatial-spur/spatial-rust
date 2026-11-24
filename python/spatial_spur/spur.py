"""Spatial persistence diagnostics, confidence sets, and transformations."""

from collections.abc import Sequence
from dataclasses import dataclass
from importlib.resources import files
import re
import warnings

import numpy as np
import pandas as pd
import statsmodels.api as sm
from numpy.typing import ArrayLike
from patsy import INTERCEPT, EvalEnvironment, ModelDesc, PatsyError, dmatrices, dmatrix

from . import _native
from ._model import (
    EARTH_CIRCUMFERENCE_M,
    FloatMatrix,
    canonical_coordinates,
    coefficient_bread,
    coordinate_columns,
    coordinates_from_frame,
    factor_codes,
    is_real_numeric_dtype,
    numeric_matrix,
)
from .scpc import SCPCResult, integer, native_info, real_float, scpc

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
    cvalue_mat: FloatMatrix | None = None

    def __repr__(self) -> str:
        return (
            f"TestResult(test={self.test!r}, lr={self.lr:.4f}, "
            f"pvalue={self.pvalue:.4g}, q={self.simulation.q})"
        )

    def summary(self) -> str:
        stat_name = "LFUR" if self.test.startswith("i1") else "LFST"
        lines = [
            f"Spatial {self.test.upper()} Test Results (q={self.simulation.q}, backend=rust)",
            "-" * 45,
            f"Test Statistic ({stat_name}):  {self.lr:9.4f}",
            f"P-value:                {self.pvalue:9.4f}",
            f"CV 1%:                  {self.cvalue[0]:9.4f}",
            f"CV 5%:                  {self.cvalue[1]:9.4f}",
            f"CV 10%:                 {self.cvalue[2]:9.4f}",
            "-" * 45,
        ]
        return "\n".join(lines)

    def __str__(self) -> str:
        return self.summary()

    @property
    def backend(self) -> str:
        return "rust"

    @property
    def n_obs(self) -> int:
        return len(self.observations)

    @property
    def n_sites(self) -> int:
        return len(self.observations)


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

    @property
    def n_obs(self) -> int:
        return len(self.observations)

    @property
    def n_sites(self) -> int:
        return len(self.observations)

    def __repr__(self) -> str:
        return (
            f"HalfLifeResult(level={self.level:g}, ci=({self.ci_lower:.4g}, {self.ci_upper:.4g}), "
            f"q={self.simulation.q}, normdist={self.normdist})"
        )

    def summary(self) -> str:
        units = "sample diameters" if self.normdist else self.distance_units
        lines = [f"{self.level:g}% spatial half-life confidence set ({units})"]
        if self.components.size:
            lines.extend(f"[{lo:.6g}, {hi:.6g}]" for lo, hi in self.components)
        else:
            lines.append("Empty confidence set")
        lines.append(f"Max distance in sample {self.max_dist:.6g} {self.distance_units}")
        return "\n".join(lines)

    def __str__(self) -> str:
        return self.summary()


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


def _require_dataframe(data: pd.DataFrame) -> None:
    if not isinstance(data, pd.DataFrame) or not data.columns.is_unique:
        raise ValueError("data must be a pandas DataFrame with unique columns")


def _coerce_masked_columns(data: pd.DataFrame) -> pd.DataFrame:
    out = data.copy(deep=False)
    for col in out.columns:
        series = out[col]
        if isinstance(series.dtype, pd.api.extensions.ExtensionDtype) and is_real_numeric_dtype(
            series
        ):
            out[col] = series.to_numpy(dtype=float, na_value=np.nan)
    return out


def _warn_spur_inference(
    q: int,
    nrep: int,
    tests: Sequence[TestResult],
    *,
    level: float | None = None,
    all_directions: bool = False,
) -> None:
    notes = []
    all_directions = all_directions or any(t.all_directions for t in tests)
    if q >= 40 or all_directions:
        scope = " uses all positive directions and" if all_directions else ""
        notes.append(
            f"q={q}{scope} retains higher spatial frequencies. "
            "Results can be sensitive to local noise."
        )
    if level is None:
        unattainable = [f"{alpha}%" for alpha in (1, 5, 10) if nrep + 1 < 100 // alpha]
        if unattainable:
            notes.append(
                f"With nrep={nrep}, the minimum p-value is {1 / (nrep + 1):.4g}. "
                f"Rejection is unattainable at {', '.join(unattainable)}."
            )
        for test in tests:
            if test.best_tuning_power is not None and test.best_tuning_power < 0.45:
                notes.append(
                    f"{test.test} did not reach the 50% tuning target. "
                    f"The best estimated power was {test.best_tuning_power:.1%}."
                )
    else:
        if q <= 5:
            notes.append(f"q={q} provides few spatial contrasts for estimating the half-life.")
        if nrep / (nrep + 1) < level / 100:
            notes.append(
                f"With nrep={nrep}, the minimum p-value is {1 / (nrep + 1):.4g}. "
                f"The {level:g}% confidence set necessarily covers [0, infinity]."
            )
    if notes:
        warnings.warn(" ".join(notes), UserWarning, stacklevel=3)


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
        native.get("cvalue_mat"),
    )


def _response(data, column):
    if column not in data or not is_real_numeric_dtype(data[column]):
        raise ValueError(f"Response {column!r} must name a numeric column")
    values = data[column].to_numpy(dtype=float, na_value=np.nan)
    if not np.isfinite(values).all():
        raise ValueError(f"Response {column!r} contains missing or non-finite values")
    return values


def _formula_design(formula, data, environment):
    _require_dataframe(data)
    if not data.index.is_unique:
        raise ValueError("Use a unique index to identify fitted regression rows")
    if isinstance(formula, str) and re.search(r"\boffset\s*\(", formula):
        raise ValueError("SPUR does not support regression offsets")
    prepared = _coerce_masked_columns(data)
    # Keep the design and its term slices together across Statsmodels versions.
    try:
        response, design = dmatrices(
            formula, prepared, eval_env=environment, return_type="dataframe"
        )
    except PatsyError as exc:
        raise ValueError(f"Invalid SPUR formula: {exc}") from exc
    if len(design) == 0:
        raise ValueError("SPUR formula has no complete finite observations")
    if response.shape[1] != 1 or not np.isfinite(response.to_numpy()).all():
        raise ValueError("SPUR requires one finite numeric response")
    x = numeric_matrix(design, "Fitted SPUR design")
    coefficient_bread(x, tuple(design.columns))
    return response.iloc[:, 0], design, data.loc[design.index]


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
    _require_dataframe(data)
    execution = native_info(backend, nthreads)
    q = integer(q, "q", 2)
    if test.endswith("resid"):
        environment = _environment or EvalEnvironment.capture(1)
        response, design, sample = _formula_design(formula, data, environment)
        model_inputs = (response.to_numpy(), design.to_numpy())
        max_q = len(sample) - design.shape[1]
    else:
        sample = data
        column, _, rhs = formula.partition("~")
        if rhs.strip() not in ("", "1"):
            raise ValueError("Levels diagnostics take one response column")
        model_inputs = (_response(data, column.strip()),)
        max_q = len(sample) - 1
    coords, geographic = _coordinates(sample, lon, lat, coords_euclidean)
    if q > max_q:
        raise ValueError("q exceeds the available spatial or residual degrees of freedom")
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
    out = _diagnostic(test, result, sim, sample.index, geographic, execution)
    _warn_spur_inference(sim.q, sim.inference_draws, [out])
    return out


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
    _require_dataframe(data)
    if not isinstance(normdist, (bool, np.bool_)):
        raise ValueError("normdist must be a boolean")
    normdist = bool(normdist)
    execution = native_info(backend, nthreads)
    q = integer(q, "q", 2, 100)
    level = real_float(level, "level", gt=0.0, lt=100.0)
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
    _warn_spur_inference(
        sim.q, sim.inference_draws, [], level=level, all_directions=out["all_directions"]
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


def _transform_inputs(
    formula: str | Sequence[str], data: pd.DataFrame, environment: EvalEnvironment
) -> tuple[tuple[str, ...], pd.DataFrame]:
    work = _coerce_masked_columns(data).reset_index(drop=True)
    if isinstance(formula, str):
        stripped = formula.strip()
        if stripped in work:
            columns = (stripped,)
        elif "~" not in formula and all(
            part.strip() in work for part in formula.split("+") if part.strip()
        ):
            columns = tuple(part.strip() for part in formula.split("+") if part.strip())
        else:
            if re.search(r"\boffset\s*\(", formula):
                raise ValueError("SPUR transformation does not support offsets")
            try:
                desc = ModelDesc.from_formula(formula)
                terms = [
                    term
                    for term in desc.lhs_termlist + desc.rhs_termlist
                    if term != INTERCEPT
                ]
                if not terms:
                    raise ValueError("Supply at least one response variable to transform")
                pos_work = work.reset_index(drop=True)
                design = dmatrix(
                    ModelDesc([], terms),
                    pos_work,
                    eval_env=environment,
                    NA_action="drop",
                    return_type="dataframe",
                )
                design = design.reindex(pos_work.index)
                design.index = work.index
            except PatsyError as exc:
                raise ValueError(f"Invalid SPUR transformation formula: {exc}") from exc
            return tuple(design.columns), design
    else:
        columns = tuple(formula)
    if not columns or not all(isinstance(col, str) for col in columns):
        raise ValueError("Supply distinct response columns")
    missing = [col for col in columns if col not in work]
    if missing:
        raise ValueError(f"Response {missing[0]!r} must name a numeric column")
    selected = pd.DataFrame(index=work.index)
    for col in columns:
        if not is_real_numeric_dtype(work[col]):
            raise ValueError(f"Response {col!r} must name a numeric column")
        selected[col] = work[col].to_numpy(dtype=float, na_value=np.nan)
    return columns, selected


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
    _require_dataframe(data)
    if not isinstance(replace, (bool, np.bool_)):
        raise ValueError("replace must be a boolean")
    replace = bool(replace)
    if not isinstance(prefix, str) or not prefix:
        raise ValueError("Supply a nonempty string output prefix")
    execution = native_info(backend, nthreads)
    if transformation not in {"lbmgls", "nn", "iso", "cluster"}:
        raise ValueError("transformation must be lbmgls, nn, iso, or cluster")
    if (transformation == "iso") != (radius is not None):
        raise ValueError("Supply radius only for the ISO transformation")
    if (transformation == "cluster") != (clustvar is not None):
        raise ValueError("Supply clustvar only for the cluster transformation")
    if transformation == "cluster" and any(
        val is not None for val in (lon, lat, coords_euclidean)
    ):
        raise ValueError("Do not supply spatial coordinates for the cluster transformation")
    if radius is not None:
        radius = real_float(radius, "radius", gt=0.0)

    columns, selected = _transform_inputs(formula, data, EvalEnvironment.capture(1))
    if not columns or len(set(columns)) != len(columns):
        raise ValueError("Supply distinct response columns and a nonempty output prefix")
    outputs = [prefix + name for name in columns]
    if not replace and any(name in data for name in outputs):
        raise ValueError("The requested output columns already exist")

    resp_raw = selected.to_numpy(dtype=float, na_value=np.nan)
    if np.isinf(resp_raw).any():
        raise ValueError("Response contains non-finite values")
    valid_mask = ~np.isnan(resp_raw).any(axis=1)

    if clustvar is None:
        coord_cols = coordinate_columns(lon, lat, coords_euclidean)
        geographic = lon is not None
        work = data.reset_index(drop=True)
        for col in coord_cols:
            if col not in work:
                raise ValueError(f"Coordinate column {col!r} is missing")
            if not is_real_numeric_dtype(work[col]):
                raise ValueError(f"Coordinate column {col!r} must be numeric")
            vals = work[col].to_numpy(dtype=float, na_value=np.nan)
            if np.isinf(vals).any():
                raise ValueError(f"Coordinate column {col!r} contains non-finite values")
            valid_mask &= ~np.isnan(vals)
        valid_rows = np.flatnonzero(valid_mask)
        min_rows = 2
        if len(valid_rows) < min_rows:
            raise ValueError("SPUR transformation requires at least 2 complete observations")
        coords = canonical_coordinates(
            coordinates_from_frame(work.iloc[valid_rows], coord_cols), geographic
        )
        groups = None
    else:
        if not isinstance(clustvar, str) or clustvar not in data:
            raise ValueError(f"Spatial cluster column {clustvar!r} is missing")
        geographic = False
        cluster_series = data[clustvar].reset_index(drop=True)
        valid_mask &= ~cluster_series.isna().to_numpy()
        valid_rows = np.flatnonzero(valid_mask)
        if len(valid_rows) < 1:
            raise ValueError("Cluster transformation requires at least 1 complete observation")
        coords = None
        groups, _ = factor_codes(cluster_series.iloc[valid_rows], "Spatial cluster")

    if radius is not None and geographic:
        radius /= EARTH_CIRCUMFERENCE_M
    result = _native.run_spur_transform(
        transformation,
        len(valid_rows),
        np.ascontiguousarray(resp_raw[valid_rows]),
        coords=coords,
        latlong=geographic,
        radius=radius,
        cluster_id=groups,
        nthreads=execution["nthreads_requested"],
    )["transformed"]
    padded = np.full((len(data), len(outputs)), np.nan, dtype=float)
    padded[valid_rows] = result
    out = data.copy()
    for name, col in zip(outputs, padded.T, strict=True):
        out[name] = col
    return out


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
    avc = real_float(avc, "avc", gt=0.001, lt=0.99)
    response, design, sample = _formula_design(formula, data, EvalEnvironment.capture(1))
    intercept = design.design_info.term_slices.get(INTERCEPT)
    if intercept is None:
        raise ValueError("The SPUR pipeline requires a fitted intercept")
    y, x = response.to_numpy(), design.to_numpy()
    regressors = [i for i in range(x.shape[1]) if i != intercept.start]
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
    _warn_spur_inference(sim.q, sim.inference_draws, list(tests.values()))
    raw = out["transformed_vars"]
    column_names = design.columns
    slope_names = [column_names[i] for i in regressors]
    transformed_fit = sm.OLS(
        pd.Series(raw[:, 0], index=sample.index),
        pd.DataFrame(raw[:, 1:], index=sample.index, columns=slope_names),
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
    levels = scpc(sm.OLS(response, design).fit(), **options)
    transformed = scpc(transformed_fit, **options)
    names = [response.name, *slope_names]
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
    _require_dataframe(data)
    cols = list(columns)
    missing = [col for col in cols if col not in data]
    if missing:
        raise ValueError(f"Column {missing[0]!r} is missing from data")
    out = data.copy()
    for col in cols:
        values = out[col]
        out[col] = (values - values.mean()) / values.std()
    return out
