import numpy as np
import pandas as pd
import pytest
import statsmodels.api as sm
import statsmodels.formula.api as smf

from reference.spatial import distances, regression, scpc_at_dimension
from spatial_spur import TiedEigenspaceError, __version__, build_info, scpc


def fixture(n=72):
    rng = np.random.default_rng(19)
    d = pd.DataFrame(rng.normal(size=(n, 6)), columns=list("yxuzwa"))
    d.index = [f"row{i}" for i in range(n)]
    d["b"] = rng.normal(size=n)
    d["u"] = 2 * d.z + d.u
    d["y"] = d.y + d.x + d.u
    d["g"] = np.arange(n) % 4
    d["h"] = np.arange(n) // 6
    d["id"] = np.arange(n) + 1000
    return d


def test_fitted_order_rank_and_nullable_columns():
    d = fixture()
    d["a"] = d.a.astype("Float64")
    d.loc[d.index[5], "y"] = np.nan
    m = smf.ols("y ~ x", d).fit()
    opt = dict(coords_euclidean=["a", "b"], backend="rust", nthreads=1, cvs=True)
    expected = scpc(m, d, **opt)
    actual = scpc(m, d.iloc[::-1], **opt)
    np.testing.assert_array_equal(expected.scpcstats, actual.scpcstats)
    np.testing.assert_allclose(actual.critical_95, actual.scpccvs[:, 2])
    np.testing.assert_allclose(actual.w[:, 1:].T @ actual.w[:, 1:], np.eye(actual.q), atol=2e-14)
    assert actual.observations.tolist() == m.model.data.row_labels.tolist()
    assert actual.n_obs == actual.n_sites == 71
    bad = d.copy()
    bad.loc[d.index[10], "a"] = pd.NA
    with pytest.raises(ValueError, match="Coordinate column 'a'.*missing"):
        scpc(m, bad, **opt)
    d["alias"] = 2 * d.x
    rank_deficient = smf.ols("y ~ x + alias", d).fit()
    for uncond in (True, False):
        with pytest.raises(ValueError, match="not separately identified.*x, alias"):
            scpc(rank_deficient, d, uncond=uncond, **opt)


@pytest.mark.parametrize("method", ["exact", "nystrom"])
def test_clustered_public_weights_reconstruct_standard_errors(method):
    d = fixture(72)
    d["site"] = np.arange(len(d)) // 2
    d[["a", "b"]] = np.repeat(d[["a", "b"]].to_numpy()[::2], 2, axis=0)
    d = d.sample(frac=1, random_state=42)
    m = smf.ols("y ~ x", d).fit()
    out = scpc(
        m,
        d.iloc[::-1],
        coords_euclidean=["a", "b"],
        cluster="site",
        method=method,
        large_n_seed=17,
        landmarks=22,
        landmark_samples=4,
        nthreads=1,
    )
    eig = out.w[:, 1:]
    np.testing.assert_allclose(out.w.T @ out.w, np.eye(out.q + 1), atol=2e-12)
    score = m.model.exog * m.resid.to_numpy()[:, None] @ m.normalized_cov_params.to_numpy()
    by_site = np.zeros((out.n_sites, score.shape[1]))
    np.add.at(by_site, np.asarray(out.site_membership), score)
    se = np.sqrt(out.n_sites / out.q) * np.linalg.norm(eig.T @ by_site, axis=0)
    np.testing.assert_allclose(out.scpcstats[:, 1], se, rtol=1e-12)
    if method == "nystrom":
        assert (out.large_n_seed, out.landmarks, out.landmark_samples) == (17, 22, 4)


def test_complete_tied_blocks_and_geographic_colocation():
    d = pd.DataFrame({"x": [0, 1, 2, 4, 8, 5], "y": [2, 5, 1, 2, 8, 3], "a": np.eye(6).tolist()})
    coordinates = np.asarray(d.pop("a").tolist())
    model = smf.ols("y ~ x", d).fit()
    out = scpc(model, coords=coordinates, nthreads=1)
    assert out.q == 5 and out.tied
    order = [3, 0, 2, 5, 1, 4]
    again = scpc(smf.ols("y ~ x", d.iloc[order]).fit(), coords=coordinates[order], nthreads=1)
    np.testing.assert_allclose(out.scpcstats, again.scpcstats, atol=2e-13)
    coordinates = np.column_stack(([-180, 180, 0, 120, 22, 138], [0, 0, 90, 90, 20, 20]))
    with pytest.raises(ValueError, match="co-located-pair fraction.*4 distinct"):
        scpc(model, coords=coordinates, geographic=True, nthreads=1)


def test_dense_reference_for_clustered_conditional_ols():
    d = fixture(18)
    d["s"] = np.repeat(np.arange(9), 2)
    d[["a", "b"]] = np.repeat(d[["a", "b"]].to_numpy()[::2], 2, axis=0)
    m = smf.ols("y ~ x", d).fit()
    out = scpc(m, d, coords_euclidean=["a", "b"], cluster="s", nthreads=1, cv_tolerance=1e-6)
    membership = np.eye(out.n_sites)[np.asarray(out.site_membership)]
    ref = scpc_at_dimension(
        d.y.to_numpy(),
        m.model.exog,
        distances(out.coordinates),
        out.q,
        out.avc,
        membership=membership,
    )
    columns = ["coef", "se", "t", "p"]
    expected = np.column_stack([ref[col] for col in columns])
    np.testing.assert_allclose(out.scpcstats[:, :4], expected, atol=2e-8)
    np.testing.assert_allclose(out.critical_95, ref["critical"], atol=1.1e-6)
    assert np.all(out.critical_95 + 1e-6 >= ref["endpoint"])


def test_pyfixest_fitted_design_and_rows():
    pf = pytest.importorskip("pyfixest")
    d = fixture(90)
    d.loc[d.index[::4], "y"] = np.nan
    d["clipw"] = d.w.clip(*d.w.quantile([0.1, 0.9]))
    opt = dict(coords_euclidean=["a", "b"], nthreads=1, cvs=True)
    clip = "np.clip(w, np.quantile(w, .1), np.quantile(w, .9))"
    for extra in (" | g", " | g + h", " | g + h | u ~ z"):
        a = pf.feols(f"y ~ x + {clip}{extra}", d)
        b = pf.feols(f"y ~ x + clipw{extra}", d)
        result = scpc(a, d.iloc[::-1], row_id="id", **opt)
        np.testing.assert_allclose(result.scpcstats, scpc(b, d, **opt).scpcstats, atol=2e-12)
        assert result.estimator == ("2sls" if "u ~ z" in extra else "ols")
        assert result.n_obs == d.y.notna().sum()
        with pytest.raises(ValueError, match="row_id"):
            scpc(a, d[["a", "b"]], **opt)


def test_pyfixest_spline_only_instruments_match_fitted_matrix_oracle():
    pf = pytest.importorskip("pyfixest")
    d = fixture(84)
    d.loc[d.index[::5], "y"] = np.nan
    m = pf.feols("y ~ x | g | u ~ bs(z, df=4)", d)
    out = scpc(m, d, coords_euclidean=["a", "b"], backend="rust", nthreads=1)
    beta, u, influence, _ = regression(m._Y.ravel(), m._X, m._Z)
    np.testing.assert_allclose(out.scpcstats[:, 0], beta, atol=1e-12)
    target = np.sqrt(out.n_sites / out.q) * np.linalg.norm(
        out.w[:, 1:].T @ (influence * u[:, None]), axis=0
    )
    np.testing.assert_allclose(out.scpcstats[:, 1], target, rtol=1e-11)


def test_chetty_fracfor_independence_limit_sets_cutoff():
    from reference.spatial import conditional_weights, critical
    from spatial_spur import load_chetty_data, standardize

    d = load_chetty_data().query("state not in ['AK', 'HI']")
    d = standardize(d[["am", "fracfor", "lon", "lat"]].dropna(), ["am", "fracfor"])
    model = smf.ols("am ~ 0 + fracfor", d).fit()
    fits = [
        scpc(
            model,
            d,
            lon="lon",
            lat="lat",
            method=method,
            landmarks=len(d),
            cv_tolerance=1e-6,
            nthreads=1,
        )
        for method in ("exact", "nystrom")
    ]
    _, _, directions, projection = regression(model.model.endog, model.model.exog)
    for out in fits:
        assert out.q == 8
        wc = conditional_weights(out.w, directions[:, 0], projection, np.eye(len(d)))
        bound = critical([out.w.T @ out.w, wc.T @ wc])
        assert abs(bound - 3.441575) < 2e-6
        assert abs(out.critical_95[0] - bound) < 2e-6


def test_statsmodels_panel_index_aligns_by_observation():
    d = fixture()
    d.index = pd.MultiIndex.from_arrays([d.g, d.id], names=["unit", "row"])
    model = smf.ols("y ~ x", d).fit()
    options = dict(coords_euclidean=["a", "b"], nthreads=1)
    aligned = scpc(model, d, **options)
    shuffled = scpc(model, d.iloc[::-1], **options)
    assert aligned.observations.equals(d.index)
    np.testing.assert_array_equal(aligned.scpcstats, shuffled.scpcstats)


def test_scpc_validation_methods_and_metadata():
    assert TiedEigenspaceError.__module__ == "spatial_spur._native"
    info = build_info()
    assert info["package_version"] == __version__

    d = fixture(32)
    m = smf.ols("y ~ x", d).fit()
    base = dict(coords_euclidean=["a", "b"], nthreads=1)
    for bad in (
        dict(avc=True),
        dict(cv_tolerance=True),
        dict(method="nystrom", landmarks=True),
        dict(method="nystrom", landmark_samples=False),
        dict(geographic=True),
    ):
        with pytest.raises((TypeError, ValueError)):
            scpc(m, d, **{**base, **bad})
    with pytest.raises(ValueError, match="duplicate"):
        scpc(m, d, coords_euclidean=["a", "a"], nthreads=1)

    off = smf.ols("y ~ x", d, offset=d.z).fit()
    with pytest.raises(ValueError, match="offsets"):
        scpc(off, d, **base)
    wls = smf.wls("y ~ x", d, weights=np.arange(1.0, len(d) + 1.0)).fit()
    with pytest.raises(ValueError, match="unweighted"):
        scpc(wls, d, **base)

    res = scpc(m, d, cvs=True, **base)
    np.testing.assert_allclose(res.coef(), m.params.to_numpy())
    ci95 = res.confint()
    ci90 = res.confint(level=0.90)
    assert ci95.shape == (2, 2)
    assert np.all(ci90.iloc[:, 1] - ci90.iloc[:, 0] < ci95.iloc[:, 1] - ci95.iloc[:, 0])
    assert "SCPC" in repr(res)
    summary_text = res.summary()
    assert "avc=" in summary_text and "Two-sided critical values:" in summary_text

    d_unhash = d.copy()
    d_unhash["tags"] = [[i] for i in range(len(d))]
    m_unhash = smf.ols("y ~ x", d_unhash).fit()
    aligned = scpc(m_unhash, d_unhash, **base)
    np.testing.assert_allclose(aligned.scpcstats, res.scpcstats)


def test_pyfixest_weighted_and_residualizer_exception_preservation():
    from spatial_spur import _native

    d = fixture(40)
    coords = d[["a", "b"]].to_numpy()
    coef = np.array([1.0])
    bread = np.array([[1.0]])
    scores = d[["y"]].to_numpy()

    class CustomResidualizerError(RuntimeError):
        pass

    def raise_interrupt(_matrix):
        raise KeyboardInterrupt("user stop")

    def raise_custom(_matrix):
        raise CustomResidualizerError("boom from callback")

    for err_cls, fn in (
        (KeyboardInterrupt, raise_interrupt),
        (CustomResidualizerError, raise_custom),
    ):
        with pytest.raises(err_cls):
            _native.run_model_scpc(
                coords,
                coef,
                bread,
                scores,
                0.03,
                inference="conditional_ols",
                method="exact",
                design=d[["x"]].to_numpy(),
                direction_design=d[["x"]].to_numpy(),
                residualizer=fn,
                nthreads=1,
            )

    pf = pytest.importorskip("pyfixest")
    wfit = pf.feols("y ~ x", d.assign(w=np.arange(1.0, len(d) + 1.0)), weights="w")
    with pytest.raises(ValueError, match="unweighted"):
        scpc(wfit, d, coords_euclidean=["a", "b"], nthreads=1)
    unhash = pf.feols("y ~ x", d)
    d_list = d.copy()
    d_list["meta"] = [[i] for i in range(len(d))]
    out = scpc(unhash, d_list.iloc[::-1], coords_euclidean=["a", "b"], nthreads=1)
    assert out.n_obs == len(d)


def test_native_build_info_matches_package_version():
    from spatial_spur import _native

    assert isinstance(sm, object)
    assert _native.build_info()["package_version"] == __version__


def test_pyfixest_constant_shared_column_rejected_on_external_frame():
    pf = pytest.importorskip("pyfixest")
    d = fixture(36).assign(k=1)
    m = pf.feols("y ~ x", d)
    external = d[["a", "b", "k"]].iloc[::-1].reset_index(drop=True)
    with pytest.raises(ValueError, match="row_id"):
        scpc(m, external, coords_euclidean=["a", "b"], nthreads=1)


def test_nystrom_single_sample_does_not_stop_at_initial_limit():
    rng = np.random.default_rng(123)
    d = pd.DataFrame(
        {
            "a": rng.uniform(0, 10, 48),
            "b": rng.uniform(0, 10, 48),
            "x": rng.normal(size=48),
        }
    )
    d["y"] = 1.0 + 0.5 * d["x"] + rng.normal(scale=0.3, size=48)
    m = smf.ols("y ~ x", d).fit()
    out = scpc(
        m,
        d,
        coords_euclidean=["a", "b"],
        avc=0.005,
        method="nystrom",
        landmarks=27,
        landmark_samples=1,
        large_n_seed=1,
        nthreads=1,
    )
    assert out.q > 10
