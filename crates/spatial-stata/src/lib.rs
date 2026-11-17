//! Stata SPI v3 plugin entry points backed by `spatial-core`.

pub mod stata_spi;

use spatial_core::{run_exact_scpc, ExactScpcRequest, Metric};
use stata_spi::{
    active_obs, error, plugin_var_count, read_matrix, read_var_matrix, save_scalar, write_matrix,
    STATA_ERR, STATA_OK,
};
use std::ffi::CStr;
use std::os::raw::{c_char, c_int};
use std::panic::{catch_unwind, AssertUnwindSafe};

const DEFAULT_CV_TOL: f64 = 0.001;

struct PreparedOlsInputs {
    wj: Vec<f64>,
    xj: Option<Vec<f64>>,
    xjs: Option<Vec<f64>>,
    proj_q: Option<Vec<f64>>,
    proj_rank: usize,
}

fn parse_bool_flag(raw: &str, name: &str) -> Result<bool, String> {
    match raw {
        "0" => Ok(false),
        "1" => Ok(true),
        _ => Err(format!("invalid boolean flag {name}='{raw}' (expected 0 or 1)")),
    }
}

fn parse_usize(raw: &str, name: &str) -> Result<usize, String> {
    raw.parse::<usize>()
        .map_err(|_| format!("invalid integer argument {name}='{raw}'"))
}

fn parse_nthreads(raw: &str) -> Result<Option<usize>, String> {
    let val = raw
        .parse::<isize>()
        .map_err(|_| format!("invalid nthreads argument '{raw}'"))?;
    if val == -1 {
        Ok(None)
    } else if val >= 1 {
        Ok(Some(val as usize))
    } else {
        Err(format!("nthreads must be -1 or >= 1, got {val}"))
    }
}

/// Two-pass Modified Gram-Schmidt on `[x_raw, 1]` (or `[x_raw]` if `!has_cons`),
/// returning row-major orthonormal `Q` (`n x rank`) and `rank`.
fn orthonormal_basis_mgs(x_raw: &[f64], n: usize, nx: usize, has_cons: bool) -> (Vec<f64>, usize) {
    let p = nx + usize::from(has_cons);
    if p == 0 || n == 0 {
        return (Vec::new(), 0);
    }

    let mut cols: Vec<Vec<f64>> = Vec::with_capacity(p);
    let mut v = vec![0.0; n];

    for col_idx in 0..p {
        if has_cons && col_idx == p - 1 {
            v.fill(1.0);
        } else {
            for i in 0..n {
                v[i] = x_raw[i * nx + col_idx];
            }
        }
        let orig_norm = v.iter().map(|x| x * x).sum::<f64>().sqrt();
        if !orig_norm.is_finite() || orig_norm <= 1e-14 {
            continue;
        }
        // Two passes ensure near machine-precision orthogonality for ill-conditioned designs.
        for _ in 0..2 {
            for q_col in &cols {
                let dot: f64 = q_col.iter().zip(v.iter()).map(|(a, b)| a * b).sum();
                for (vi, &qi) in v.iter_mut().zip(q_col.iter()) {
                    *vi -= dot * qi;
                }
            }
        }
        let norm = v.iter().map(|x| x * x).sum::<f64>().sqrt();
        if norm.is_finite() && norm > 1e-12 * orig_norm && norm > 1e-14 {
            let inv = 1.0 / norm;
            cols.push(v.iter().map(|&x| x * inv).collect());
        }
    }

    let rank = cols.len();
    let mut q_row_major = vec![0.0; n * rank];
    for (r, q_col) in cols.iter().enumerate() {
        for i in 0..n {
            q_row_major[i * rank + r] = q_col[i];
        }
    }
    (q_row_major, rank)
}

fn prepare_ols_inputs(
    n: usize,
    p: usize,
    k_use: usize,
    has_cons: bool,
    uncond: bool,
    b: &[f64],
    bread: &[f64],
    scores: &[f64],
    x_raw: &[f64],
) -> Result<PreparedOlsInputs, String> {
    let nx = p - usize::from(has_cons);
    if b.len() != p || bread.len() != p * p || scores.len() != n * p {
        return Err("internal dimension mismatch in OLS score inputs".to_string());
    }
    if !uncond && x_raw.len() != n * nx {
        return Err("internal dimension mismatch in conditional regressor inputs".to_string());
    }

    let nf = n as f64;
    let mut wj = vec![0.0; n * k_use];
    for i in 0..n {
        let score_row = &scores[i * p..(i + 1) * p];
        for j in 0..k_use {
            let bread_row = &bread[j * p..(j + 1) * p];
            let dot: f64 = bread_row
                .iter()
                .zip(score_row.iter())
                .map(|(a, s)| a * s)
                .sum();
            wj[i * k_use + j] = nf * dot + b[j];
        }
    }

    if uncond {
        return Ok(PreparedOlsInputs {
            wj,
            xj: None,
            xjs: None,
            proj_q: None,
            proj_rank: 0,
        });
    }

    let mut xj = vec![0.0; n * k_use];
    let mut xjs = vec![0.0; n * k_use];
    for i in 0..n {
        let x_row = &x_raw[i * nx..(i + 1) * nx];
        for j in 0..k_use {
            let bread_row = &bread[j * p..(j + 1) * p];
            let mut dot: f64 = bread_row[..nx]
                .iter()
                .zip(x_row.iter())
                .map(|(a, x)| a * x)
                .sum();
            if has_cons {
                dot += bread_row[p - 1];
            }
            let x_dir = nf * dot;
            xj[i * k_use + j] = x_dir;
            xjs[i * k_use + j] = if x_dir > 0.0 {
                1.0
            } else if x_dir < 0.0 {
                -1.0
            } else {
                0.0
            };
        }
    }

    let (proj_q, proj_rank) = orthonormal_basis_mgs(x_raw, n, nx, has_cons);
    Ok(PreparedOlsInputs {
        wj,
        xj: Some(xj),
        xjs: Some(xjs),
        proj_q: (proj_rank > 0).then_some(proj_q),
        proj_rank,
    })
}

fn handle_exact_scpc(args: &[&str]) -> Result<(), String> {
    if args.len() != 16 {
        return Err(format!(
            "exact_scpc expects 16 arguments after subcommand, got {}",
            args.len()
        ));
    }

    let d = parse_usize(args[0], "d")?;
    let p = parse_usize(args[1], "p")?;
    let k_use = parse_usize(args[2], "k_use")?;
    let has_cons = parse_bool_flag(args[3], "has_cons")?;
    let avc = args[4]
        .parse::<f64>()
        .map_err(|_| format!("invalid avc '{}'", args[4]))?;
    let latlong = parse_bool_flag(args[5], "latlong")?;
    let uncond = parse_bool_flag(args[6], "uncond")?;
    let cvs = parse_bool_flag(args[7], "cvs")?;
    let nthreads = parse_nthreads(args[8])?;

    let b_mat_name = args[9];
    let bread_mat_name = args[10];
    let stats_mat_name = args[11];
    let cvs_mat_name = args[12];
    let c0_scal_name = args[13];
    let cv_scal_name = args[14];
    let q_scal_name = args[15];

    if d == 0 || p == 0 || k_use == 0 || k_use > p {
        return Err(format!(
            "invalid dimensions: d={d}, p={p}, k_use={k_use}"
        ));
    }
    if latlong && d != 2 {
        return Err(format!("latlong requires d=2 coordinates, got d={d}"));
    }
    if !avc.is_finite() || avc <= 0.001 || avc >= 0.99 {
        return Err(format!("avc must lie strictly in (0.001, 0.99), got {avc}"));
    }

    let nx = p - usize::from(has_cons);
    let expected_vars = d + p + if uncond { 0 } else { nx };
    let actual_vars = plugin_var_count();
    if actual_vars != expected_vars {
        return Err(format!(
            "exact_scpc expected {expected_vars} variables in plugin varlist, got {actual_vars}"
        ));
    }

    let obs = active_obs();
    let n = obs.len();
    if n < 2 {
        return Err(format!("too few active observations ({n}) for SCPC"));
    }

    // Read all Stata dataset variables and matrices on the calling thread before Rayon.
    let coords = read_var_matrix(&obs, 1, d, "coordinates (s_*)")?;
    let scores = read_var_matrix(&obs, 1 + d, p, "scores (scpc_score*)")?;
    let x_raw = if !uncond && nx > 0 {
        read_var_matrix(&obs, 1 + d + p, nx, "conditional regressors")?
    } else {
        Vec::new()
    };

    let b = read_matrix(b_mat_name, 1, p)?;
    let bread = read_matrix(bread_mat_name, p, p)?;

    let inputs = prepare_ols_inputs(n, p, k_use, has_cons, uncond, &b, &bread, &scores, &x_raw)?;
    let metric = if latlong {
        Metric::GreatCircleFractionLatLon
    } else {
        Metric::Euclidean
    };

    let out = run_exact_scpc(&ExactScpcRequest {
        coords: &coords,
        n,
        d,
        metric,
        avc,
        uncond,
        cvs,
        coef: &b[..k_use],
        wj: &inputs.wj,
        xj: inputs.xj.as_deref(),
        xjs: inputs.xjs.as_deref(),
        proj_u: inputs.proj_q.as_deref(),
        proj_v: inputs.proj_q.as_deref(),
        proj_rank: inputs.proj_rank,
        cv_tolerance: DEFAULT_CV_TOL,
        nthreads,
    })
    .map_err(|e| e.to_string())?;

    write_matrix(stats_mat_name, k_use, 6, &out.scpcstats)?;
    if cvs {
        let cvs_mat = out
            .scpccvs
            .as_deref()
            .ok_or_else(|| "missing scpccvs output when cvs=1".to_string())?;
        write_matrix(cvs_mat_name, k_use, 4, cvs_mat)?;
    }
    save_scalar(c0_scal_name, out.c0)?;
    save_scalar(cv_scal_name, out.cvfin)?;
    save_scalar(q_scal_name, out.q as f64)?;

    Ok(())
}

fn dispatch(args: &[&str]) -> Result<(), String> {
    let subcmd = *args
        .first()
        .ok_or_else(|| "missing plugin subcommand".to_string())?;
    match subcmd {
        "ping" | "ping_v2" => Ok(()),
        "exact_scpc" => handle_exact_scpc(&args[1..]),
        other => Err(format!("unknown spatial-stata subcommand '{other}'")),
    }
}

/// Entry point invoked from `c/stata_bridge.c`.
#[no_mangle]
pub extern "C" fn rs_stata_call(argc: c_int, argv: *const *const c_char) -> c_int {
    let res = catch_unwind(AssertUnwindSafe(|| {
        if argc <= 0 || argv.is_null() {
            return Err("plugin called without subcommand arguments".to_string());
        }
        let mut args = Vec::with_capacity(argc as usize);
        for idx in 0..(argc as usize) {
            let ptr = unsafe { *argv.add(idx) };
            if ptr.is_null() {
                return Err(format!("null argv[{idx}] pointer"));
            }
            let arg = unsafe { CStr::from_ptr(ptr) }
                .to_str()
                .map_err(|_| format!("argv[{idx}] is not valid UTF-8"))?;
            args.push(arg);
        }
        dispatch(&args)
    }));

    match res {
        Ok(Ok(())) => STATA_OK,
        Ok(Err(msg)) => {
            error(&format!("scpc_rust plugin error: {msg}"));
            STATA_ERR
        }
        Err(_) => {
            error("scpc_rust plugin panicked unexpectedly");
            STATA_ERR
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn load_bench_200_csv() -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>) {
        let path = format!(
            "{}/../../fixtures/../../verification/data/bench_200.csv",
            env!("CARGO_MANIFEST_DIR")
        );
        let text = std::fs::read_to_string(path).expect("bench_200.csv");
        let mut y = Vec::new();
        let mut x1 = Vec::new();
        let mut x2 = Vec::new();
        let mut coords_lon_lat = Vec::new();
        for line in text.lines().skip(1) {
            if line.trim().is_empty() {
                continue;
            }
            let cols: Vec<f64> = line
                .split(',')
                .map(|s| s.trim().parse::<f64>().unwrap())
                .collect();
            y.push(cols[1]);
            x1.push(cols[2]);
            x2.push(cols[3]);
            coords_lon_lat.extend_from_slice(&[cols[6], cols[7]]);
        }
        (y, x1, x2, coords_lon_lat)
    }

    fn invert_symmetric(mat: &[f64], p: usize) -> Vec<f64> {
        let mut a = mat.to_vec();
        let mut inv = vec![0.0; p * p];
        for i in 0..p {
            inv[i * p + i] = 1.0;
        }
        for col in 0..p {
            let mut pivot = col;
            for row in (col + 1)..p {
                if a[row * p + col].abs() > a[pivot * p + col].abs() {
                    pivot = row;
                }
            }
            if pivot != col {
                for k in 0..p {
                    a.swap(col * p + k, pivot * p + k);
                    inv.swap(col * p + k, pivot * p + k);
                }
            }
            let diag = a[col * p + col];
            for k in 0..p {
                a[col * p + k] /= diag;
                inv[col * p + k] /= diag;
            }
            for row in 0..p {
                if row == col {
                    continue;
                }
                let factor = a[row * p + col];
                for k in 0..p {
                    a[row * p + k] -= factor * a[col * p + k];
                    inv[row * p + k] -= factor * inv[col * p + k];
                }
            }
        }
        inv
    }

    #[test]
    fn orthonormal_basis_handles_cons_and_nocons() {
        let n = 4;
        let x_raw = vec![1.0, 2.0, 3.0, 4.0];
        let (q_cons, rank_cons) = orthonormal_basis_mgs(&x_raw, n, 1, true);
        assert_eq!(rank_cons, 2);
        // Check Q' Q == I_2
        for r1 in 0..2 {
            for r2 in 0..2 {
                let dot: f64 = (0..n)
                    .map(|i| q_cons[i * 2 + r1] * q_cons[i * 2 + r2])
                    .sum();
                let exp = if r1 == r2 { 1.0 } else { 0.0 };
                assert!((dot - exp).abs() < 1e-14);
            }
        }
        let (q_nocons, rank_nocons) = orthonormal_basis_mgs(&x_raw, n, 1, false);
        assert_eq!(rank_nocons, 1);
        let norm2: f64 = q_nocons.iter().map(|v| v * v).sum();
        assert!((norm2 - 1.0).abs() < 1e-14);

        // Intercept-only (nx = 0, has_cons = true)
        let (q_int, rank_int) = orthonormal_basis_mgs(&[], n, 0, true);
        assert_eq!(rank_int, 1);
        for val in q_int {
            assert!((val - 0.5).abs() < 1e-14);
        }
    }

    #[test]
    fn ols_preparation_matches_bench_200_euclid_and_latlong() {
        let (y, x1, x2, coords_lon_lat) = load_bench_200_csv();
        let n = y.len();
        let p = 3;
        let nx = 2;
        let mut x_raw = Vec::with_capacity(n * nx);
        let mut x_design = Vec::with_capacity(n * p);
        for i in 0..n {
            x_raw.extend_from_slice(&[x1[i], x2[i]]);
            x_design.extend_from_slice(&[x1[i], x2[i], 1.0]);
        }

        let mut xtx = vec![0.0; p * p];
        let mut xty = vec![0.0; p];
        for i in 0..n {
            for r in 0..p {
                xty[r] += x_design[i * p + r] * y[i];
                for c in 0..p {
                    xtx[r * p + c] += x_design[i * p + r] * x_design[i * p + c];
                }
            }
        }
        let bread = invert_symmetric(&xtx, p);
        let mut b = vec![0.0; p];
        for r in 0..p {
            b[r] = (0..p).map(|c| bread[r * p + c] * xty[c]).sum();
        }
        let mut scores = vec![0.0; n * p];
        for i in 0..n {
            let yhat: f64 = (0..p).map(|c| x_design[i * p + c] * b[c]).sum();
            let resid = y[i] - yhat;
            for c in 0..p {
                scores[i * p + c] = x_design[i * p + c] * resid;
            }
        }

        let inputs = prepare_ols_inputs(n, p, p, true, false, &b, &bread, &scores, &x_raw).unwrap();
        let out_euclid = run_exact_scpc(&ExactScpcRequest {
            coords: &coords_lon_lat,
            n,
            d: 2,
            metric: Metric::Euclidean,
            avc: 0.03,
            uncond: false,
            cvs: true,
            coef: &b,
            wj: &inputs.wj,
            xj: inputs.xj.as_deref(),
            xjs: inputs.xjs.as_deref(),
            proj_u: inputs.proj_q.as_deref(),
            proj_v: inputs.proj_q.as_deref(),
            proj_rank: inputs.proj_rank,
            cv_tolerance: DEFAULT_CV_TOL,
            nthreads: Some(2),
        })
        .unwrap();
        assert_eq!(out_euclid.q, 19);
        assert!((out_euclid.c0 - 3.500409719601247).abs() < 1e-12);
        // Row 0 is x1, row 1 is x2, row 2 is _cons
        assert!((out_euclid.scpcstats[0] - 1.5473280763431876).abs() < 1e-10);
        assert!((out_euclid.scpcstats[12] - 1.4834167407175964).abs() < 1e-10);

        // Stata latlong: s_1=lat, s_2=lon
        let mut coords_lat_lon = Vec::with_capacity(n * 2);
        for pair in coords_lon_lat.chunks_exact(2) {
            coords_lat_lon.extend_from_slice(&[pair[1], pair[0]]);
        }
        let out_geo = run_exact_scpc(&ExactScpcRequest {
            coords: &coords_lat_lon,
            n,
            d: 2,
            metric: Metric::GreatCircleFractionLatLon,
            avc: 0.03,
            uncond: false,
            cvs: true,
            coef: &b,
            wj: &inputs.wj,
            xj: inputs.xj.as_deref(),
            xjs: inputs.xjs.as_deref(),
            proj_u: inputs.proj_q.as_deref(),
            proj_v: inputs.proj_q.as_deref(),
            proj_rank: inputs.proj_rank,
            cv_tolerance: DEFAULT_CV_TOL,
            nthreads: Some(2),
        })
        .unwrap();
        assert_eq!(out_geo.q, 18);
        assert!((out_geo.c0 - 1396.4193247376622).abs() < 1e-10);
        assert!((out_geo.cvfin - 2.1009220402410382).abs() < 1e-10);
    }

    #[test]
    fn arg_parsers_validate_inputs() {
        assert_eq!(parse_nthreads("-1").unwrap(), None);
        assert_eq!(parse_nthreads("1").unwrap(), Some(1));
        assert_eq!(parse_nthreads("4").unwrap(), Some(4));
        assert!(parse_nthreads("0").is_err());
        assert!(parse_nthreads("-2").is_err());
        assert!(!parse_bool_flag("0", "flag").unwrap());
        assert!(parse_bool_flag("1", "flag").unwrap());
        assert!(parse_bool_flag("2", "flag").is_err());
    }
}
