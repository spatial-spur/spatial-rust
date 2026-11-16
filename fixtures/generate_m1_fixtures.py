#!/usr/bin/env python3
"""Generate frozen Milestone 1 exact conditional & unconditional SCPC fixtures
from pinned main/master releases across Python, R, and Stata.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path("/Users/pauldavidboll/Github/spatial-spur")
sys.path.insert(0, str(ROOT / "scpc-python" / "src"))

import numpy as np
import pandas as pd
import pyfixest as pf
import statsmodels.formula.api as smf
from scpc import scpc

DATA_DIR = ROOT / "verification" / "data"
OUT_PATH = ROOT / "spatial-rust" / "fixtures" / "m1_exact_scpc.json"
R_LIB = str(ROOT / ".Rlib")
STATA_BIN = "/Applications/StataNow/StataSE.app/Contents/MacOS/stata-se"

CASES = [
    ("bench_200_euclid_cond", "bench_200.csv", False, False, True),
    ("bench_200_euclid_uncond", "bench_200.csv", False, True, False),
    ("bench_200_latlong_cond", "bench_200.csv", True, False, True),
    ("bench_1000_euclid_cond", "bench_1000.csv", False, False, True),
    ("bench_1000_euclid_uncond", "bench_1000.csv", False, True, False),
    ("bench_1000_latlong_cond", "bench_1000.csv", True, False, True),
]


def run_python_cases() -> dict:
    results = {}
    for name, csv_name, latlong, uncond, cvs in CASES:
        df = pd.read_csv(DATA_DIR / csv_name)
        coord_kwargs = (
            {"lon": "s_1", "lat": "s_2"}
            if latlong
            else {"coords_euclidean": ["s_1", "s_2"]}
        )

        fit_sm = smf.ols("y ~ x1 + x2", data=df).fit()
        res_sm = scpc(fit_sm, df, avc=0.03, uncond=uncond, cvs=cvs, **coord_kwargs)

        fit_pf = pf.feols("y ~ x1 + x2", data=df, vcov="hetero")
        res_pf = scpc(fit_pf, df, avc=0.03, uncond=uncond, cvs=cvs, **coord_kwargs)

        results[name] = {
            "statsmodels": {
                "c0": float(res_sm.c0),
                "cv": float(res_sm.cv),
                "q": int(res_sm.q),
                "coef_names": list(res_sm.coef_names),
                "scpcstats": res_sm.scpcstats.tolist(),
                "scpccvs": res_sm.scpccvs.tolist() if res_sm.scpccvs is not None else None,
            },
            "pyfixest": {
                "c0": float(res_pf.c0),
                "cv": float(res_pf.cv),
                "q": int(res_pf.q),
                "coef_names": list(res_pf.coef_names),
                "scpcstats": res_pf.scpcstats.tolist(),
                "scpccvs": res_pf.scpccvs.tolist() if res_pf.scpccvs is not None else None,
            },
        }
    return results


def run_r_cases() -> dict:
    case_lines = []
    for name, csv_name, latlong, uncond, cvs in CASES:
        coord_arg = (
            'lon = "s_1", lat = "s_2"'
            if latlong
            else 'coords_euclidean = c("s_1", "s_2")'
        )
        uncond_r = "TRUE" if uncond else "FALSE"
        cvs_r = "TRUE" if cvs else "FALSE"
        case_lines.append(
            f'''
            df <- read.csv("{DATA_DIR / csv_name}")
            fit <- lm(y ~ x1 + x2, data = df)
            res <- scpc(fit, data = df, {coord_arg}, avc = 0.03, uncond = {uncond_r}, cvs = {cvs_r})
            out[["{name}"]] <- list(
              c0 = as.numeric(res$c0),
              cv = as.numeric(res$cv),
              q = as.integer(res$q),
              coef_names = rownames(res$scpcstats),
              scpcstats = unname(res$scpcstats),
              scpccvs = if (!is.null(res$scpccvs)) unname(res$scpccvs) else NULL
            )
            '''
        )

    script = (
        f'.libPaths(c("{R_LIB}", .libPaths()))\n'
        'suppressPackageStartupMessages({ library(scpcR); library(jsonlite) })\n'
        'out <- list()\n'
        + "\n".join(case_lines)
        + '\ncat(toJSON(out, auto_unbox = TRUE, digits = 16, null = "null"))\n'
    )
    proc = subprocess.run(["Rscript", "-e", script], capture_output=True, text=True, check=True)
    return json.loads(proc.stdout.strip())


def run_stata_cases() -> dict:
    results = {}
    tmp_dir = Path("/tmp/scpc_m1_fixtures")
    tmp_dir.mkdir(parents=True, exist_ok=True)
    for name, csv_name, latlong, uncond, cvs in CASES:
        do_path = tmp_dir / f"{name}.do"
        json_path = tmp_dir / f"{name}.json"
        log_path = Path(f"{name}.log")
        # In Stata scpc, latlong treats s_1 as latitude and s_2 as longitude; CSV has s_1=lon, s_2=lat
        coord_setup = (
            "drop s_1_clustered s_2_clustered\nrename s_1 tmp_lon\nrename s_2 s_1\nrename tmp_lon s_2"
            if latlong
            else "drop s_1_clustered s_2_clustered"
        )
        opts = []
        if latlong:
            opts.append("latlong")
        if uncond:
            opts.append("uncond")
        if cvs:
            opts.append("cvs")
        opts_str = ", avc(0.03) " + " ".join(opts)

        # Reorder Stata rows from [x1, x2, _cons] to [(Intercept), x1, x2] for easy comparison
        cvs_write = ""
        if cvs:
            cvs_write = """
    cvs = st_matrix("e(scpccvs)")
    fput(fh, sprintf(`""scpccvs": [[%22.14e, %22.14e, %22.14e, %22.14e], [%22.14e, %22.14e, %22.14e, %22.14e], [%22.14e, %22.14e, %22.14e, %22.14e]]"',
        cvs[3,1], cvs[3,2], cvs[3,3], cvs[3,4],
        cvs[1,1], cvs[1,2], cvs[1,3], cvs[1,4],
        cvs[2,1], cvs[2,2], cvs[2,3], cvs[2,4]))
"""
        else:
            cvs_write = '    fput(fh, `""scpccvs": null"\')'

        do_text = f"""
clear all
set more off
adopath ++ "{ROOT / 'scpc-stata' / 'src'}"
import delimited "{DATA_DIR / csv_name}", clear
{coord_setup}
regress y x1 x2, robust
scpc{opts_str}

mata:
b = st_matrix("e(scpcstats)")
(void) _unlink("{json_path}")
fh = fopen("{json_path}", "w")
fput(fh, "{{")
fput(fh, sprintf(`""c0": %22.14e,"', *findexternal("c0")))
fput(fh, sprintf(`""cv": %22.14e,"', *findexternal("cvfin")))
fput(fh, sprintf(`""q": %10.0f,"', cols(*findexternal("Wfin")) - 1))
fput(fh, sprintf(`""scpcstats": ["'))
fput(fh, sprintf(`"[%22.14e, %22.14e, %22.14e, %22.14e, %22.14e, %22.14e],"', b[3,1], b[3,2], b[3,3], b[3,4], b[3,5], b[3,6]))
fput(fh, sprintf(`"[%22.14e, %22.14e, %22.14e, %22.14e, %22.14e, %22.14e],"', b[1,1], b[1,2], b[1,3], b[1,4], b[1,5], b[1,6]))
fput(fh, sprintf(`"[%22.14e, %22.14e, %22.14e, %22.14e, %22.14e, %22.14e]],"', b[2,1], b[2,2], b[2,3], b[2,4], b[2,5], b[2,6]))
{cvs_write}
fput(fh, "}}")
fclose(fh)
end
"""
        do_path.write_text(do_text)
        subprocess.run([STATA_BIN, "-b", "do", str(do_path)], cwd=str(tmp_dir), capture_output=True, text=True, check=True)
        with open(json_path) as f:
            results[name] = json.load(f)
    return results


def main() -> None:
    print("Generating Python fixtures...")
    py_res = run_python_cases()
    print("Generating R fixtures...")
    r_res = run_r_cases()
    print("Generating Stata fixtures...")
    stata_res = run_stata_cases()

    payload = {
        "python_commit": "83716a8",
        "r_commit": "e6e451c",
        "stata_commit": "064a622+matastrict",
        "python": py_res,
        "r": r_res,
        "stata": stata_res,
    }
    OUT_PATH.write_text(json.dumps(payload, indent=2) + "\n")
    print(f"Saved fixtures to {OUT_PATH}")

    # Cross-check Python statsmodels vs R vs Stata immediately
    for name, _, _, _, _ in CASES:
        p = py_res[name]["statsmodels"]
        r = r_res[name]
        s = stata_res[name]
        p_stats = np.array(p["scpcstats"])
        r_stats = np.array(r["scpcstats"])
        s_stats = np.array(s["scpcstats"])
        print(
            f"{name}: q=(py:{p['q']}, r:{r['q']}, st:{s['q']}), "
            f"c0_diff_py_r={abs(p['c0']-r['c0']):.2e}, "
            f"max_stat_diff_py_r={np.max(np.abs(p_stats-r_stats)):.2e}, "
            f"max_stat_diff_py_st={np.max(np.abs(p_stats-s_stats)):.2e}"
        )


if __name__ == "__main__":
    main()
