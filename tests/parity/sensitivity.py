"""Measure approximation sensitivity and memory on designs near the automatic switch."""

import argparse
import json
import resource
import sys
import time

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

from spatial_spur import scpc


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("method", choices=["auto", "exact", "nystrom"])
    parser.add_argument("--n", type=int, default=4500)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--nthreads", type=int, default=4)
    args = parser.parse_args()
    rng = np.random.default_rng(73)
    data = pd.DataFrame(rng.uniform(size=(args.n, 2)), columns=["a", "b"])
    data["treat"] = (data.a + data.b).rank() <= 30
    data["treat"] = data.treat.astype(float)
    data["x"] = rng.normal(size=args.n)
    data["y"] = 0.6 * data.x + 0.9 * data.treat + np.sin(4 * data.a) + rng.normal(size=args.n)
    model = smf.ols("y ~ x + treat", data).fit()
    start = time.perf_counter()
    out = scpc(
        model,
        data,
        coords_euclidean=["a", "b"],
        method=args.method,
        large_n_seed=args.seed,
        nthreads=args.nthreads,
        backend="rust",
    )
    elapsed = time.perf_counter() - start
    peak_rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    if sys.platform != "darwin":
        peak_rss *= 1024
    assert out.method == (
        "nystrom"
        if args.method == "auto" and args.n >= 4500
        else "exact"
        if args.method == "auto"
        else args.method
    )
    assert np.all(np.isfinite(out.scpcstats))
    print(
        json.dumps(
            {
                "requested_method": args.method,
                "method": out.method,
                "seed": out.large_n_seed,
                "n": out.n_obs,
                "q": out.q,
                "seconds": elapsed,
                "peak_rss_mib": peak_rss / 2**20,
                "nthreads": args.nthreads,
                "execution": out.execution,
                "uncertainty": out.table[["se", "ci_lower", "ci_upper"]].to_dict("index"),
            }
        )
    )


if __name__ == "__main__":
    main()
