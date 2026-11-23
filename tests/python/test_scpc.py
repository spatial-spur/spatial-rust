import numpy as np
import pandas as pd
import pytest
import statsmodels.formula.api as smf

from reference.spatial import distances, regression, scpc_at_dimension
from spatial_spur import scpc


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
