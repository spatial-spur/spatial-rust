import warnings

import numpy as np
import pandas as pd
import pytest

from reference.spatial import brownian_gls, distances
from spatial_spur import (
    load_chetty_data,
    spur,
    spurhalflife,
    spurtest_i0,
    spurtest_i1resid,
    spurtransform,
)


def test_geographic_nearest_neighbors_and_iso_radius():
    d = pd.DataFrame(dict(lon=-120 + np.arange(17) / 4, lat=0, y=np.arange(17)))
    nn = spurtransform("~y", d, transformation="nn", lon="lon", lat="lat", nthreads=1)
    np.testing.assert_array_equal(nn.h_y, np.r_[-1, np.zeros(15), 1])
    coords = pd.DataFrame(dict(s=[0.0, 0.0, 1.0, 2.1], y=[2, 6, 9, 10]))
    iso = spurtransform(
        "y", coords, coords_euclidean=["s"], transformation="iso", radius=1, nthreads=1
    )
    np.testing.assert_array_equal(iso.h_y, [-4, 4, 0, 0])
    nn = spurtransform("y", coords, coords_euclidean=["s"], transformation="nn", nthreads=1)
    np.testing.assert_allclose(nn.h_y, [-7, -3, 5, 1])


def test_stable_neighbors_under_common_offsets():
    d = pd.DataFrame(dict(s=np.arange(1000), y=(np.arange(1000) % 13) / 8))
    args = dict(coords_euclidean=["s"], transformation="iso", radius=2000, nthreads=2)
    ref = spurtransform("y", d, **args)
    np.testing.assert_array_equal(ref.h_y, spurtransform("y", d.assign(y=d.y + 1e12), **args).h_y)
    np.testing.assert_array_equal(spurtransform("y", d.assign(y=1e12), **args).h_y, 0)
    d = pd.DataFrame({"s": [-1.0, 0.0, 1.0], "y": [1.0, 4.0, 3.0]})
    assert spurtransform("y", d, coords_euclidean=["s"], transformation="nn").h_y[1] == 2.0


def test_lbmgls_and_pipeline_use_fitted_sample_and_named_inputs():
    rng = np.random.default_rng(19)
    d = pd.DataFrame(rng.normal(size=(48, 4)), columns=["a", "b", "y", "x"])
    d.loc[[4, 5], "y"] = np.nan
    clean = d.dropna()
    expected = brownian_gls(clean[["y", "x"]].to_numpy(), distances(clean[["a", "b"]].to_numpy()))
    options = dict(coords_euclidean=["a", "b"], backend="rust", nthreads=1)
    separate = spurtransform("y~x", clean, **options)
    np.testing.assert_allclose(separate[["h_y", "h_x"]], expected, atol=1e-11)
    draws = rng.standard_normal((4, 240))
    tuning = rng.standard_normal((4, 240))
    result = spur("y~x", d, q=4, draws=draws, tuning_draws=tuning, **options)
    np.testing.assert_allclose(result.transformed_data[["h_y", "h_x"]], expected, atol=1e-11)
    assert result.levels.n_obs == result.transformed.n_obs == len(clean)
    single = spurtest_i0("y", clean, q=4, draws=draws, tuning_draws=tuning, **options)
    assert single.pvalue == result.tests["i0"].pvalue
    assert single.simulation.inference_draws == single.simulation.tuning_draws == 240
    assert single.simulation.seed is None
    with pytest.raises(ValueError, match="separate tuning draw"):
        spurtest_i0("y", clean, q=4, draws=draws, **options)


def test_pipeline_and_residual_diagnostics_preserve_formula_terms():
    rng = np.random.default_rng(27)
    data = pd.DataFrame(rng.normal(size=(60, 4)), columns=["a", "b", "y", "x"])
    data["group"] = pd.Categorical(np.tile(["a", "b", "c"], 20))
    data.loc[14, "x"] = np.nan

    def squared(values):
        return values**2

    draws = rng.standard_normal((4, 240))
    options = dict(
        coords_euclidean=["a", "b"],
        q=4,
        draws=draws,
        tuning_draws=rng.standard_normal(draws.shape),
        nthreads=1,
    )
    formula = "y ~ C(group) + squared(x)"
    result = spur(formula, data, **options)
    assert result.levels.table.index.tolist() == [
        "Intercept",
        "C(group)[T.b]",
        "C(group)[T.c]",
        "squared(x)",
    ]
    assert result.levels.observations.equals(data.index.delete(14))
    assert "h_squared(x)" in result.transformed_data
    direct = spurtest_i1resid(formula, data, **options)
    assert direct.pvalue == result.tests["i1resid"].pvalue
    np.testing.assert_allclose(direct.lr, result.tests["i1resid"].lr)
    with pytest.raises(ValueError, match="requires a fitted intercept"):
        spur("y ~ 0 + Intercept + x", data.assign(Intercept=1.0), **options)


def test_half_life_saved_units_and_components():
    d = load_chetty_data().loc[:, ["lat", "lon", "am"]].dropna().iloc[:39]
    draws = np.random.default_rng(14).standard_normal((4, 256))
    opts = dict(lon="lon", lat="lat", q=4, draws=draws, backend="rust", nthreads=1)
    raw = spurhalflife("am", d, **opts)
    normalized = spurhalflife("am", d, normdist=True, **opts)
    assert raw.distance_units == "meters"
    assert str(normalized).startswith("95% spatial half-life confidence set (sample diameters)")
    np.testing.assert_allclose(raw.components, normalized.components * raw.max_dist)
    np.testing.assert_allclose(raw.hl_grid, normalized.hl_grid * raw.max_dist)
    assert raw.simulation.tuning_draws == 0


def test_i1_gaussian_calibration_against_independent_angular_ratio():
    from reference.spatial import i1_calibration
    from spatial_spur import spurtest_i1

    rng = np.random.default_rng(713)
    d = pd.DataFrame(rng.normal(size=(36, 3)), columns=["a", "b", "y"])
    normals = rng.standard_normal((8, 253))
    result = spurtest_i1(
        "y",
        d,
        coords_euclidean=["a", "b"],
        q=4,
        nthreads=1,
        draws=normals[:4],
        tuning_draws=normals[4:],
    )
    dmat = distances(d[["a", "b"]].to_numpy())
    lr, pvalue, cutoffs = i1_calibration(d.y.to_numpy(), dmat, 4, result.ha_param, normals[:4])
    assert result.pvalue == pvalue
    np.testing.assert_allclose(result.lr, lr, rtol=1e-11)
    np.testing.assert_allclose(result.cvalue, cutoffs, rtol=1e-11)


def test_nn_tolerance_preserves_ties_and_distinguishes_resolved_neighbors():
    matrix = np.array([[0, 0], [1, 0], [0, 1], [-1, 0]], dtype=float)
    theta = 0.381
    rotation = np.array([[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]])
    options = dict(coords_euclidean=["a", "b"], nthreads=1, transformation="nn")
    origin = np.array([520_000, 4_510_000])

    def transformed(points):
        data = pd.DataFrame(points @ rotation + origin, columns=["a", "b"])
        return spurtransform("y", data.assign(y=[1, -1, 2, -2]), **options).h_y.iloc[0]

    assert transformed(matrix) == pytest.approx(4 / 3)
    matrix[2:] *= 1.000001
    assert transformed(matrix) == 2


def test_spur_booleans_nullable_formulas_self_col_and_summaries():
    rng = np.random.default_rng(77)
    d = pd.DataFrame(
        {
            "a": rng.normal(size=30),
            "b": rng.normal(size=30),
            "y": pd.Series(rng.normal(size=30), dtype="Float64"),
            "x": pd.Series(rng.normal(size=30), dtype="Float64"),
            "self": rng.normal(size=30),
        }
    )
    d.loc[2, "y"] = pd.NA
    d.loc[5, "x"] = pd.NA
    opts = dict(coords_euclidean=["a", "b"], nthreads=1)

    with pytest.raises((TypeError, ValueError)):
        spurhalflife("y", d, level=True, q=4, nrep=50, **opts)
    with pytest.raises((TypeError, ValueError)):
        spurtransform("y", d, transformation="iso", radius=True, **opts)
    with pytest.raises((TypeError, ValueError)):
        spur("y ~ x", d, q=4, nrep=50, avc=True, **opts)

    tr = spurtransform("y ~ I(x ** 2) + self", d, transformation="lbmgls", **opts)
    assert "h_y" in tr.columns and "h_I(x ** 2)" in tr.columns and "h_self" in tr.columns
    assert np.isnan(tr.loc[2, "h_y"]) and np.isnan(tr.loc[5, "h_y"])
    assert tr["h_y"].notna().sum() == 28

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        t_res = spurtest_i0("y", d.dropna(), q=4, nrep=80, **opts)
        hl_res = spurhalflife("y", d.dropna(), q=4, nrep=10, level=95, **opts)
        s_res = spur("y ~ x", d, q=4, nrep=80, **opts)
    assert any("unattainable" in str(w.message) for w in caught)
    assert any("few spatial contrasts" in str(w.message) for w in caught)
    assert t_res.cvalue_mat.shape == (len(t_res.null_grid), 3) and len(t_res.null_grid) >= 75
    assert "P-value:" in t_res.summary() and repr(t_res).startswith("TestResult(")
    assert (
        "half-life confidence set" in hl_res.summary()
        and repr(hl_res).startswith("HalfLifeResult(")
    )
    assert "SPUR diagnostic p-values" in s_res.summary()
    assert s_res.transformed.table.index.tolist() == ["x"]


def test_spur_rejects_zero_draws_reports_25_draw_power_and_is_row_permutation_invariant():
    from spatial_spur import spurtest_i1

    d = pd.DataFrame(
        {
            "s": np.linspace(0.0, 1.0, 24),
            "y": np.sin(np.linspace(0.2, 2.7, 24)),
        }
    )
    opts = dict(coords_euclidean=["s"], q=4, nthreads=1)
    zeros = np.zeros((4, 25))
    valid = np.random.default_rng(11).standard_normal((4, 25))
    for fn in (spurtest_i0, spurtest_i1):
        with pytest.raises(ValueError, match="positive column norms"):
            fn("y", d, draws=zeros, tuning_draws=valid, **opts)
        with pytest.raises(ValueError, match="positive column norms"):
            fn("y", d, draws=valid, tuning_draws=zeros, **opts)

    tuned = spurtest_i0("y", d, nrep=25, seed=9, **opts)
    assert tuned.best_tuning_power is not None and np.isfinite(tuned.best_tuning_power)

    rng = np.random.default_rng(31)
    draws = rng.standard_normal((4, 80))
    tuning = rng.standard_normal((4, 80))
    d_rev = d.iloc[::-1].reset_index(drop=True)
    for fn in (spurtest_i0, spurtest_i1):
        a = fn("y", d, draws=draws, tuning_draws=tuning, **opts)
        b = fn("y", d_rev, draws=draws, tuning_draws=tuning, **opts)
        assert np.isclose(a.lr, b.lr, rtol=0, atol=1e-11)
        assert np.isclose(a.pvalue, b.pvalue, rtol=0, atol=1e-12)
        assert np.isclose(a.ha_param, b.ha_param, rtol=0, atol=1e-11)
