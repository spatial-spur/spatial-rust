import numpy as np
import pandas as pd
import pytest

from reference.spatial import brownian_gls, distances
from spatial_spur import load_chetty_data, spur, spurhalflife, spurtest_i0, spurtransform


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
