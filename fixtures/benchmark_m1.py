#!/usr/bin/env python3
"""Benchmark Milestone 1 exact conditional OLS SCPC across Python, R, and Stata.

Measures cold (first invocation in a fresh process) and warm (median of 5 timed
calls after 1 warm-up) elapsed milliseconds on bench_200.csv and bench_1000.csv
for y ~ x1 + x2, avc = 0.03, uncond = False, cvs = False.
"""

from __future__ import annotations

import json
import os
import statistics
import subprocess
import time
from pathlib import Path

ROOT = Path("/Users/pauldavidboll/Github/spatial-spur")
DATA_DIR = ROOT / "verification" / "data"
PY_BIN = ROOT / "scpc-python" / ".venv" / "bin" / "python"
R_LIB = ROOT / ".Rlib"
STATA_BIN = "/Applications/StataNow/StataSE.app/Contents/MacOS/stata-se"
OUT_JSON = ROOT / "spatial-rust" / "fixtures" / "benchmark_m1_results.json"

SIZES = [200, 1000]
WARM_RUNS = 5


def bench_python(n: int, backend: str, nthreads: int | None) -> dict:
    code = f"""
import json, os, time, statistics
for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[k] = "1"
import pandas as pd
import statsmodels.formula.api as smf
from scpc import scpc

df = pd.read_csv("{DATA_DIR / f'bench_{n}.csv'}")
fit = smf.ols("y ~ x1 + x2", data=df).fit()
t0 = time.perf_counter()
res0 = scpc(fit, df, coords_euclidean=["s_1", "s_2"], avc=0.03, backend="{backend}", nthreads={repr(nthreads)})
cold_ms = (time.perf_counter() - t0) * 1000.0
assert res0.backend == "{backend}", f"expected {backend}, got {{res0.backend}}"

times = []
for _ in range({WARM_RUNS}):
    t0 = time.perf_counter()
    res = scpc(fit, df, coords_euclidean=["s_1", "s_2"], avc=0.03, backend="{backend}", nthreads={repr(nthreads)})
    times.append((time.perf_counter() - t0) * 1000.0)
print(json.dumps({{
    "cold_ms": round(cold_ms, 2),
    "warm_median_ms": round(statistics.median(times), 2),
    "warm_min_ms": round(min(times), 2),
    "q": int(res0.q),
    "backend": res0.backend,
}}))
"""
    proc = subprocess.run([str(PY_BIN), "-c", code], capture_output=True, text=True, check=True)
    return json.loads(proc.stdout.strip().splitlines()[-1])


def bench_r(n: int, backend: str, nthreads: int | None) -> dict:
    nt_arg = "NULL" if nthreads is None else f"{nthreads}L"
    script = f"""
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
.libPaths(c("{R_LIB}", .libPaths()))
suppressPackageStartupMessages({{ library(scpcR); library(jsonlite) }})
df <- read.csv("{DATA_DIR / f'bench_{n}.csv'}")
fit <- lm(y ~ x1 + x2, data = df)
t0 <- proc.time()[["elapsed"]]
res0 <- scpc(fit, data = df, coords_euclidean = c("s_1", "s_2"), avc = 0.03, backend = "{backend}", nthreads = {nt_arg})
cold_ms <- (proc.time()[["elapsed"]] - t0) * 1000
stopifnot(identical(res0$backend, "{backend}"))
times <- numeric({WARM_RUNS})
for (i in seq_len({WARM_RUNS})) {{
  t0 <- proc.time()[["elapsed"]]
  res <- scpc(fit, data = df, coords_euclidean = c("s_1", "s_2"), avc = 0.03, backend = "{backend}", nthreads = {nt_arg})
  times[i] <- (proc.time()[["elapsed"]] - t0) * 1000
}}
cat(toJSON(list(
  cold_ms = round(cold_ms, 2),
  warm_median_ms = round(median(times), 2),
  warm_min_ms = round(min(times), 2),
  q = as.integer(res0$q),
  backend = res0$backend
), auto_unbox = TRUE))
"""
    proc = subprocess.run(["Rscript", "-e", script], capture_output=True, text=True, check=True)
    return json.loads(proc.stdout.strip().splitlines()[-1])


def bench_stata(n: int, backend: str, nthreads: int | None) -> dict:
    tmp_dir = Path("/tmp/scpc_m1_bench")
    tmp_dir.mkdir(parents=True, exist_ok=True)
    do_path = tmp_dir / f"bench_{n}_{backend}_{nthreads}.do"
    out_csv = tmp_dir / f"bench_{n}_{backend}_{nthreads}.csv"
    nt_opt = "" if nthreads is None else f"nthreads({nthreads})"
    do_text = f"""
clear all
set more off
adopath ++ "{ROOT / 'scpc-stata' / 'src'}"
import delimited "{DATA_DIR / f'bench_{n}.csv'}", clear
drop s_1_clustered s_2_clustered
regress y x1 x2, robust

timer clear 1
timer on 1
quietly scpc, avc(0.03) backend({backend}) {nt_opt}
timer off 1
assert "`e(backend)'" == "{backend}"
local qval = e(q)
quietly timer list 1
local cold_ms = r(t1) * 1000

forvalues i = 1/{WARM_RUNS} {{
    timer clear 2
    timer on 2
    quietly scpc, avc(0.03) backend({backend}) {nt_opt}
    timer off 2
    quietly timer list 2
    local w_`i' = r(t2) * 1000
}}

mata:
fh = fopen("{out_csv}", "w")
fput(fh, "cold_ms,w1,w2,w3,w4,w5,q")
fput(fh, sprintf("%12.3f,%12.3f,%12.3f,%12.3f,%12.3f,%12.3f,%6.0f",
    strtoreal(st_local("cold_ms")),
    strtoreal(st_local("w_1")),
    strtoreal(st_local("w_2")),
    strtoreal(st_local("w_3")),
    strtoreal(st_local("w_4")),
    strtoreal(st_local("w_5")),
    strtoreal(st_local("qval"))))
fclose(fh)
end
"""
    do_path.write_text(do_text)
    subprocess.run([STATA_BIN, "-b", "do", str(do_path)], cwd=str(tmp_dir), capture_output=True, text=True, check=True)
    lines = out_csv.read_text().strip().splitlines()
    vals = [float(x.strip()) for x in lines[1].split(",")]
    warm = vals[1:6]
    return {
        "cold_ms": round(vals[0], 2),
        "warm_median_ms": round(statistics.median(warm), 2),
        "warm_min_ms": round(min(warm), 2),
        "q": int(vals[6]),
        "backend": backend,
    }


def main() -> None:
    results = {}
    configs = [
        ("rust_4t", "rust", 4),
        ("rust_1t", "rust", 1),
        ("fallback", "fallback", None),
    ]
    for n in SIZES:
        results[str(n)] = {"python": {}, "r": {}, "stata": {}}
        for label, backend, nt in configs:
            print(f"Benchmarking N={n} {label}...", flush=True)
            results[str(n)]["python"][label] = bench_python(n, backend, nt)
            results[str(n)]["r"][label] = bench_r(n, backend, nt)
            results[str(n)]["stata"][label] = bench_stata(n, backend, nt)

    OUT_JSON.write_text(json.dumps(results, indent=2) + "\n")
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
