"""Run the installed Python, R, and Stata interfaces on common inputs."""

import argparse
import json
import os
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm

from spatial_spur import load_chetty_data, scpc, spur, spurhalflife, spurtransform, standardize

HERE = Path(__file__).resolve().parent
TERMS = """
fracblack racseg segpov25 fraccom15 hipc gini incsh1 tsr tsperc hsdrop scind fracrel crimer
fracsm fracdiv fracmar loctr colpc coltui colgrad manshare chimp tlfpr migirate migorate fracfor
""".split()
GEO = dict(lon="lon", lat="lat", nthreads=1, backend="rust")


def summary(result, case, stage):
    table = result.table.rename(columns={"ci_lower": "lo", "ci_upper": "hi"})
    table["cv"] = result.critical_95
    table["term"] = ["_cons" if name == "Intercept" else name for name in table.index]
    table["case"], table["stage"], table["q"], table["n"] = case, stage, result.q, result.n_obs
    return table.reset_index(drop=True)


def common_inputs(directory, full_chetty):
    chetty = load_chetty_data().query("state not in ['AK', 'HI']")
    common = standardize(
        chetty[["am", "fracblack", "gini", "lon", "lat"]].dropna().iloc[::6],
        ["am", "fracblack", "gini"],
    ).rename(columns={"am": "y", "fracblack": "x", "gini": "x2"})
    common.to_csv(directory / "common.csv", index=False, float_format="%.17g")
    common = pd.read_csv(directory / "common.csv", float_precision="round_trip")
    rng = np.random.default_rng(54)
    normals = rng.standard_normal((12, 400)).astype("<f8")
    half_normals = rng.standard_normal((6, 400)).astype("<f8")
    normals.tofile(directory / "diagnostics.bin")
    half_normals.tofile(directory / "halflife.bin")
    p = spur("y ~ x + x2", common, q=6, draws=normals[:6], tuning_draws=normals[6:], **GEO)
    half = spurhalflife("y", common, q=6, draws=half_normals, normdist=True, **GEO)
    diagnostic = pd.DataFrame(
        [(name, t.lr, t.pvalue, t.ha_param, *t.cvalue) for name, t in p.tests.items()],
        columns=["test", "statistic", "p", "ha", "cv01", "cv05", "cv10"],
    )
    fits = [summary(p.levels, "common", "levels"), summary(p.transformed, "common", "adjusted")]
    if full_chetty:
        cases = []
        for name in TERMS:
            d = standardize(chetty[["am", name, "lon", "lat"]].dropna(), ["am", name])
            cases.append(d.rename(columns={"am": "y", name: "x"}).assign(term=name))
        pd.concat(cases).to_csv(directory / "chetty-cases.csv", index=False, float_format="%.17g")
        for term, data in pd.read_csv(
            directory / "chetty-cases.csv", float_precision="round_trip"
        ).groupby("term", sort=False):
            a = scpc(sm.OLS(data.y, data[["x"]]).fit(), data, **GEO)
            transformed = spurtransform("y~x", data, **GEO)
            b = scpc(
                sm.OLS(transformed.h_y, transformed[["h_x"]].rename(columns={"h_x": "x"})).fit(),
                transformed,
                **GEO,
            )
            fits.extend([summary(a, term, "levels"), summary(b, term, "adjusted")])
    return pd.concat(fits), diagnostic, half.components


def compare(expected, result, keys, label):
    a = expected.set_index(keys).sort_index().sort_index(axis=1)
    b = result.set_index(keys).sort_index().sort_index(axis=1)
    assert a.index.equals(b.index) and a.columns.equals(b.columns), label
    np.testing.assert_allclose(a.to_numpy(), b.to_numpy(), rtol=3e-9, atol=5e-10, err_msg=label)
    return float(np.max(np.abs(a.to_numpy() - b.to_numpy())))


def historical(fits, path):
    archive = json.loads(path.read_text())
    historical = {}
    for stage, name in (("levels", "ols"), ("adjusted", "fgls")):
        actual = fits.loc[(fits.stage == stage) & (fits.case != "common")].set_index("case")
        actual = actual.loc[[name.lower() for name in archive["reg_names"]]]
        original = archive[name]
        np.testing.assert_allclose(actual.coef, original["beta"], atol=2e-10)
        np.testing.assert_allclose(actual.se, original["se"], atol=2e-9)
        low = np.abs(actual.lo - original["lo"])
        high = np.abs(actual.hi - original["hi"])
        historical[stage] = {
            "largest_endpoint_difference": float(max(low.max(), high.max())),
            "case": low.idxmax() if low.max() >= high.max() else high.idxmax(),
            "different_5_percent_decisions": actual.index[
                ((actual.lo > 0) | (actual.hi < 0))
                != ((np.asarray(original["lo"]) > 0) | (np.asarray(original["hi"]) < 0))
            ].tolist(),
        }
    return historical


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--stata-package", type=Path, required=True)
    parser.add_argument("--stata", default="stata")
    parser.add_argument("--rscript", default="Rscript")
    parser.add_argument("--chetty", action="store_true")
    parser.add_argument("--historical-reference", type=Path)
    args = parser.parse_args()
    work = args.work.resolve()
    work.mkdir(parents=True, exist_ok=False)
    (work / "plus").mkdir()
    fits, diagnostic, half = common_inputs(work, args.chetty)
    fits.to_csv(work / "scpc-python.csv", index=False)

    env = os.environ.copy()
    with (work / "r.log").open("w") as log:
        subprocess.run(
            [args.rscript, "--vanilla", HERE / "host.R"],
            cwd=work,
            env=env,
            stderr=subprocess.STDOUT,
            stdout=log,
            check=True,
        )
    subprocess.run(
        [
            args.stata,
            "-b",
            "do",
            str(HERE / "host.do"),
            str(args.stata_package.resolve()),
            str(work / "plus"),
        ],
        cwd=work,
        check=True,
    )
    stata_log = work / "plus.log"
    if not stata_log.exists():
        stata_log = next(path for path in work.glob("*.log") if path.name != "r.log")
    assert "Stata parity completed\n" in stata_log.read_text(), stata_log
    result = {"chetty_fits_per_host": 52 if args.chetty else 0, "common_draws": 400}
    for lang, read, extension in (("r", pd.read_csv, "csv"), ("stata", pd.read_stata, "dta")):
        result[lang] = {
            "scpc_max_abs_difference": compare(
                fits, read(work / f"scpc-{lang}.{extension}"), ["case", "stage", "term"], lang
            ),
            "diagnostic_max_abs_difference": compare(
                diagnostic, read(work / f"diagnostics-{lang}.{extension}"), ["test"], lang
            ),
        }
        matrix = read(work / f"halflife-{lang}.{extension}").to_numpy()
        if lang == "stata":
            matrix[:, 1] = np.nan_to_num(matrix[:, 1], nan=np.inf, posinf=np.inf)
        np.testing.assert_allclose(half, matrix, rtol=2e-9, atol=1e-10, err_msg=lang)
    if args.historical_reference and args.chetty:
        result["historical_matlab"] = historical(fits, args.historical_reference)
    (work / "summary.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
