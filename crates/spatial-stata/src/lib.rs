//! Stata SPI v3 plugin entry points backed by `spatial-core`.

mod spur_stata;
mod stata_spi;

use spatial_core::{
    nystrom::normalize_coords, run_approx_scpc, run_exact_scpc, ApproxScpcRequest, CoreError,
    ErrorKind, ExactScpcOutput, ExactScpcRequest, Metric,
};
use stata_spi::{
    active_obs, error, plugin_var_count, read_matrix, read_var_matrix, save_scalar, write_matrix,
    STATA_ERR, STATA_OK,
};
use std::ffi::CStr;
use std::os::raw::{c_char, c_int};
use std::panic::{catch_unwind, AssertUnwindSafe};

const DEFAULT_CV_TOL: f64 = 0.001;
const MGS_RANK_REL_TOL: f64 = 1e-12;

fn error_kind_label(kind: ErrorKind) -> &'static str {
    match kind {
        ErrorKind::InvalidInput => "invalid input",
        ErrorKind::UnsupportedCase => "unsupported case",
        ErrorKind::TiedEigenspace => "tied eigenspace",
        ErrorKind::CalibrationFailed => "calibration failed",
        ErrorKind::DegenerateVariance => "degenerate variance",
        ErrorKind::InternalFailure => "internal failure",
    }
}

pub(crate) fn format_core_error(err: CoreError) -> String {
    format!("{}: {}", error_kind_label(err.kind()), err.message())
}

struct PreparedInputs {
    wj: Vec<f64>,
    xj: Option<Vec<f64>>,
    xjs: Option<Vec<f64>>,
    proj_u: Option<Vec<f64>>,
    proj_v: Option<Vec<f64>>,
    proj_rank: usize,
}

pub(crate) fn parse_bool_flag(raw: &str, name: &str) -> Result<bool, String> {
    match raw {
        "0" => Ok(false),
        "1" => Ok(true),
        _ => Err(format!("invalid boolean flag {name}='{raw}' (expected 0 or 1)")),
    }
}

pub(crate) fn parse_usize(raw: &str, name: &str) -> Result<usize, String> {
    raw.parse::<usize>()
        .map_err(|_| format!("invalid integer argument {name}='{raw}'"))
}

fn parse_seed(raw: &str) -> Result<u32, String> {
    let val = raw
        .parse::<u64>()
        .map_err(|_| format!("invalid seed argument '{raw}' (expected integer in [0, 2^32))"))?;
    u32::try_from(val).map_err(|_| format!("seed argument '{raw}' must be in [0, 2^32)"))
}

fn permute_rows(src: &[f64], perm: &[usize], cols: usize) -> Vec<f64> {
    let mut out = Vec::with_capacity(perm.len() * cols);
    for &src_idx in perm {
        out.extend_from_slice(&src[src_idx * cols..(src_idx + 1) * cols]);
    }
    out
}

pub(crate) fn parse_nthreads(raw: &str) -> Result<Option<usize>, String> {
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

fn normalize_column(v: &mut [f64]) -> bool {
    let orig_norm = v.iter().map(|x| x * x).sum::<f64>().sqrt();
    if orig_norm.is_finite() && orig_norm > 0.0 {
        let inv = 1.0 / orig_norm;
        for vi in v.iter_mut() {
            *vi *= inv;
        }
        return true;
    }
    let mut max_abs = 0.0_f64;
    for &x in v.iter() {
        if !x.is_finite() {
            return false;
        }
        max_abs = max_abs.max(x.abs());
    }
    if max_abs == 0.0 {
        return false;
    }
    let inv_max = 1.0 / max_abs;
    for vi in v.iter_mut() {
        *vi *= inv_max;
    }
    let scaled_norm = v.iter().map(|x| x * x).sum::<f64>().sqrt();
    if !scaled_norm.is_finite() || scaled_norm == 0.0 {
        return false;
    }
    let inv = 1.0 / scaled_norm;
    for vi in v.iter_mut() {
        *vi *= inv;
    }
    true
}

/// Two-pass Modified Gram-Schmidt on arbitrary column vectors, returning row-major `Q` (`n x rank`).
fn orthonormal_cols_mgs(cols_in: &[Vec<f64>], n: usize) -> (Vec<f64>, usize) {
    if cols_in.is_empty() || n == 0 {
        return (Vec::new(), 0);
    }

    let mut cols: Vec<Vec<f64>> = Vec::with_capacity(cols_in.len());
    let mut v = vec![0.0; n];

    for src in cols_in {
        v.copy_from_slice(src);
        if !normalize_column(&mut v) {
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
        if norm.is_finite() && norm > MGS_RANK_REL_TOL {
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

/// Two-pass Modified Gram-Schmidt on `[x_raw, 1]` (or `[x_raw]` if `!has_cons`),
/// returning row-major orthonormal `Q` (`n x rank`) and `rank`.
fn orthonormal_basis_mgs(x_raw: &[f64], n: usize, nx: usize, has_cons: bool) -> (Vec<f64>, usize) {
    let p = nx + usize::from(has_cons);
    if p == 0 || n == 0 {
        return (Vec::new(), 0);
    }
    let mut cols_in = Vec::with_capacity(p);
    for col_idx in 0..nx {
        cols_in.push((0..n).map(|i| x_raw[i * nx + col_idx]).collect());
    }
    if has_cons {
        cols_in.push(vec![1.0; n]);
    }
    orthonormal_cols_mgs(&cols_in, n)
}

fn compute_wj(
    n_sites: usize,
    p: usize,
    k_use: usize,
    b: &[f64],
    bread: &[f64],
    scores: &[f64],
) -> Vec<f64> {
    let nf = n_sites as f64;
    let mut wj = vec![0.0; n_sites * k_use];
    for i in 0..n_sites {
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
    wj
}

fn compute_xjs(
    n_sites: usize,
    m_obs: usize,
    k_use: usize,
    xj: &[f64],
    cluster_id: Option<&[usize]>,
) -> Vec<f64> {
    let mut xjs = vec![0.0; m_obs * k_use];
    if let Some(cl) = cluster_id {
        let mut sq_sum = vec![0.0_f64; n_sites * k_use];
        for (i, &g) in cl.iter().enumerate() {
            for j in 0..k_use {
                let v = xj[i * k_use + j];
                sq_sum[g * k_use + j] += v * v;
            }
        }
        for (i, &g) in cl.iter().enumerate() {
            for j in 0..k_use {
                let norm = sq_sum[g * k_use + j].sqrt();
                xjs[i * k_use + j] = if norm > 0.0 {
                    xj[i * k_use + j] / norm
                } else {
                    0.0
                };
            }
        }
    } else {
        for idx in 0..xj.len() {
            let x_dir = xj[idx];
            xjs[idx] = if x_dir > 0.0 {
                1.0
            } else if x_dir < 0.0 {
                -1.0
            } else {
                0.0
            };
        }
    }
    xjs
}

fn solve_multi_rhs(a_in: &[f64], b_in: &[f64], p: usize, ncols: usize) -> Result<Vec<f64>, String> {
    let mut a = a_in.to_vec();
    let mut b = b_in.to_vec();
    let max_entry = a.iter().map(|v| v.abs()).fold(0.0_f64, f64::max);
    if !max_entry.is_finite() || max_entry == 0.0 {
        return Err("IV projection matrix X' P_Z X is zero or non-finite".to_string());
    }
    let tol = 1e-12 * max_entry;

    for col in 0..p {
        let mut pivot = col;
        let mut max_val = a[col * p + col].abs();
        for row in (col + 1)..p {
            let cand = a[row * p + col].abs();
            if cand > max_val {
                max_val = cand;
                pivot = row;
            }
        }
        if !max_val.is_finite() || max_val <= tol {
            return Err("IV residualizer is rank deficient: X' P_Z X is singular".to_string());
        }
        if pivot != col {
            for k in 0..p {
                a.swap(col * p + k, pivot * p + k);
            }
            for k in 0..ncols {
                b.swap(col * ncols + k, pivot * ncols + k);
            }
        }
        let diag = a[col * p + col];
        let inv_diag = 1.0 / diag;
        for k in col..p {
            a[col * p + k] *= inv_diag;
        }
        for k in 0..ncols {
            b[col * ncols + k] *= inv_diag;
        }
        for row in 0..p {
            if row == col {
                continue;
            }
            let factor = a[row * p + col];
            if factor == 0.0 {
                continue;
            }
            for k in col..p {
                a[row * p + k] -= factor * a[col * p + k];
            }
            for k in 0..ncols {
                b[row * ncols + k] -= factor * b[col * ncols + k];
            }
        }
    }
    Ok(b)
}

#[allow(clippy::too_many_arguments)]
fn prepare_ols_inputs(
    n_sites: usize,
    m_obs: usize,
    p: usize,
    k_use: usize,
    has_cons: bool,
    uncond: bool,
    b: &[f64],
    bread: &[f64],
    scores: &[f64],
    x_raw: &[f64],
    cluster_id: Option<&[usize]>,
) -> Result<PreparedInputs, String> {
    let nx = p - usize::from(has_cons);
    if b.len() != p || bread.len() != p * p || scores.len() != n_sites * p {
        return Err("internal dimension mismatch in OLS score inputs".to_string());
    }
    if !uncond && x_raw.len() != m_obs * nx {
        return Err("internal dimension mismatch in conditional regressor inputs".to_string());
    }
    if let Some(cl) = cluster_id {
        if cl.len() != m_obs {
            return Err("internal dimension mismatch in cluster_id inputs".to_string());
        }
    } else if m_obs != n_sites {
        return Err("m_obs must equal n_sites when unclustered".to_string());
    }

    let wj = compute_wj(n_sites, p, k_use, b, bread, scores);
    if uncond {
        return Ok(PreparedInputs {
            wj,
            xj: None,
            xjs: None,
            proj_u: None,
            proj_v: None,
            proj_rank: 0,
        });
    }

    let nf = n_sites as f64;
    let mut xj = vec![0.0; m_obs * k_use];
    for i in 0..m_obs {
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
            xj[i * k_use + j] = nf * dot;
        }
    }

    let xjs = compute_xjs(n_sites, m_obs, k_use, &xj, cluster_id);
    let (proj_q, proj_rank) = orthonormal_basis_mgs(x_raw, m_obs, nx, has_cons);
    let proj_opt = (proj_rank > 0).then_some(proj_q);
    Ok(PreparedInputs {
        wj,
        xj: Some(xj),
        xjs: Some(xjs),
        proj_u: proj_opt.clone(),
        proj_v: proj_opt,
        proj_rank,
    })
}

#[allow(clippy::too_many_arguments)]
fn prepare_iv_inputs(
    n_sites: usize,
    m_obs: usize,
    p: usize,
    k_use: usize,
    has_cons: bool,
    uncond: bool,
    n_endog: usize,
    n_z_excl: usize,
    b: &[f64],
    bread: &[f64],
    scores: &[f64],
    x_raw: &[f64],
    z_excl: &[f64],
    cluster_id: Option<&[usize]>,
) -> Result<PreparedInputs, String> {
    let nx = p - usize::from(has_cons);
    if n_endog == 0 || n_endog > nx {
        return Err(format!("invalid n_endog={n_endog} for nx={nx}"));
    }
    if b.len() != p || bread.len() != p * p || scores.len() != n_sites * p {
        return Err("internal dimension mismatch in IV score inputs".to_string());
    }
    if !uncond && (x_raw.len() != m_obs * nx || z_excl.len() != m_obs * n_z_excl) {
        return Err("internal dimension mismatch in IV conditional inputs".to_string());
    }
    if let Some(cl) = cluster_id {
        if cl.len() != m_obs {
            return Err("internal dimension mismatch in cluster_id inputs".to_string());
        }
    } else if m_obs != n_sites {
        return Err("m_obs must equal n_sites when unclustered".to_string());
    }

    let wj = compute_wj(n_sites, p, k_use, b, bread, scores);
    if uncond {
        return Ok(PreparedInputs {
            wj,
            xj: None,
            xjs: None,
            proj_u: None,
            proj_v: None,
            proj_rank: 0,
        });
    }

    let mut x_full = vec![0.0; m_obs * p];
    for i in 0..m_obs {
        x_full[i * p..i * p + nx].copy_from_slice(&x_raw[i * nx..(i + 1) * nx]);
        if has_cons {
            x_full[i * p + nx] = 1.0;
        }
    }

    let n_exog = nx - n_endog;
    let mut z_cols = Vec::with_capacity(n_exog + n_z_excl + usize::from(has_cons));
    for col in 0..n_exog {
        let raw_col = n_endog + col;
        z_cols.push((0..m_obs).map(|i| x_raw[i * nx + raw_col]).collect());
    }
    for col in 0..n_z_excl {
        z_cols.push((0..m_obs).map(|i| z_excl[i * n_z_excl + col]).collect());
    }
    if has_cons {
        z_cols.push(vec![1.0; m_obs]);
    }
    let (q_z, rank_z) = orthonormal_cols_mgs(&z_cols, m_obs);
    if rank_z < p {
        return Err("IV residualizer failed: instrument space has insufficient rank".to_string());
    }

    // Q_z' * X_full (rank_z x p), then P_Z X = Q_z * (Q_z' * X_full) (m_obs x p)
    let mut qzt_x = vec![0.0; rank_z * p];
    for i in 0..m_obs {
        let q_row = &q_z[i * rank_z..(i + 1) * rank_z];
        let x_row = &x_full[i * p..(i + 1) * p];
        for r in 0..rank_z {
            let qr = q_row[r];
            for c in 0..p {
                qzt_x[r * p + c] += qr * x_row[c];
            }
        }
    }
    let mut pzx = vec![0.0; m_obs * p];
    for i in 0..m_obs {
        let q_row = &q_z[i * rank_z..(i + 1) * rank_z];
        let p_row = &mut pzx[i * p..(i + 1) * p];
        for r in 0..rank_z {
            let qr = q_row[r];
            for c in 0..p {
                p_row[c] += qr * qzt_x[r * p + c];
            }
        }
    }

    // A = X' P_Z X (p x p)
    let mut gram_a = vec![0.0; p * p];
    for i in 0..m_obs {
        let x_row = &x_full[i * p..(i + 1) * p];
        let p_row = &pzx[i * p..(i + 1) * p];
        for r in 0..p {
            let xr = x_row[r];
            for c in 0..p {
                gram_a[r * p + c] += xr * p_row[c];
            }
        }
    }

    let mut pzx_t = vec![0.0; p * m_obs];
    for i in 0..m_obs {
        for c in 0..p {
            pzx_t[c * m_obs + i] = pzx[i * p + c];
        }
    }
    let v_t = solve_multi_rhs(&gram_a, &pzx_t, p, m_obs)?;
    let mut proj_v = vec![0.0; m_obs * p];
    for i in 0..m_obs {
        for c in 0..p {
            proj_v[i * p + c] = v_t[c * m_obs + i];
        }
    }

    let nf = n_sites as f64;
    let mut xj_base = vec![0.0; m_obs * k_use];
    for i in 0..m_obs {
        let p_row = &pzx[i * p..(i + 1) * p];
        for j in 0..k_use {
            let bread_row = &bread[j * p..(j + 1) * p];
            let dot: f64 = bread_row.iter().zip(p_row.iter()).map(|(a, b)| a * b).sum();
            xj_base[i * k_use + j] = nf * dot;
        }
    }
    let xjs_base = compute_xjs(n_sites, m_obs, k_use, &xj_base, cluster_id);

    // IV conditional denominator weights are Xi_s' M' Xi R where M' = I - V U'.
    Ok(PreparedInputs {
        wj,
        xj: Some(xj_base),
        xjs: Some(xjs_base),
        proj_u: Some(proj_v),
        proj_v: Some(x_full),
        proj_rank: p,
    })
}

fn aggregate_cluster_inputs(
    cl_raw: &[f64],
    coords_obs: &[f64],
    scores_obs: &[f64],
    d: usize,
    p: usize,
) -> Result<(usize, Vec<usize>, Vec<f64>, Vec<f64>), String> {
    let m_obs = cl_raw.len();
    if m_obs < 2 {
        return Err(format!("too few active observations ({m_obs}) for clustered SCPC"));
    }
    let mut cluster_id = Vec::with_capacity(m_obs);
    let mut max_id = 0usize;
    for &val in cl_raw {
        if !val.is_finite() || val < 1.0 || (val - val.round()).abs() > 0.0 {
            return Err(format!("invalid cluster index '{val}' in plugin varlist"));
        }
        let g = (val as usize) - 1;
        max_id = max_id.max(g);
        cluster_id.push(g);
    }
    let n_clust = max_id + 1;
    if n_clust < 2 {
        return Err(format!("too few active clusters ({n_clust}) for SCPC"));
    }

    let mut seen = vec![false; n_clust];
    let mut coords_cl = vec![0.0; n_clust * d];
    let mut scores_cl = vec![0.0; n_clust * p];
    for (i, &g) in cluster_id.iter().enumerate() {
        let obs_c = &coords_obs[i * d..(i + 1) * d];
        let cl_c = &mut coords_cl[g * d..(g + 1) * d];
        if !seen[g] {
            seen[g] = true;
            cl_c.copy_from_slice(obs_c);
            scores_cl[g * p..(g + 1) * p].copy_from_slice(&scores_obs[i * p..(i + 1) * p]);
        } else {
            for k in 0..d {
                if (obs_c[k] - cl_c[k]).abs() > 1e-8 {
                    return Err("coordinates are not constant within clusters".to_string());
                }
            }
        }
    }
    if seen.iter().any(|&s| !s) {
        return Err("cluster index vector contains empty clusters".to_string());
    }
    Ok((n_clust, cluster_id, coords_cl, scores_cl))
}

struct ExactScpcArgs<'a> {
    d: usize,
    p: usize,
    k_use: usize,
    has_cons: bool,
    avc: f64,
    latlong: bool,
    uncond: bool,
    cvs: bool,
    nthreads: Option<usize>,
    b_mat_name: &'a str,
    bread_mat_name: &'a str,
    stats_mat_name: &'a str,
    cvs_mat_name: &'a str,
    c0_scal_name: &'a str,
    cv_scal_name: &'a str,
    q_scal_name: &'a str,
    tied_scal_name: Option<&'a str>,
    has_cluster: bool,
    n_endog: usize,
    n_z_excl: usize,
}

impl<'a> ExactScpcArgs<'a> {
    fn parse(args: &'a [&'a str]) -> Result<Self, String> {
        if !(16..=18).contains(&args.len()) && args.len() != 20 {
            return Err(format!(
                "exact_scpc expects 16, 17, 18, or 20 arguments after subcommand, got {}",
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
        let tied_scal_name = args.get(16).copied();
        let has_cluster = match args.get(17) {
            Some(flag) => parse_bool_flag(flag, "has_cluster")?,
            None => false,
        };
        let n_endog = match args.get(18) {
            Some(raw) => parse_usize(raw, "n_endog")?,
            None => 0,
        };
        let n_z_excl = match args.get(19) {
            Some(raw) => parse_usize(raw, "n_z_excl")?,
            None => 0,
        };

        if d == 0 || p == 0 || k_use == 0 || k_use > p {
            return Err(format!(
                "invalid dimensions: d={d}, p={p}, k_use={k_use}"
            ));
        }
        let nx = p - usize::from(has_cons);
        if n_endog > nx {
            return Err(format!("invalid n_endog={n_endog} > nx={nx}"));
        }
        if latlong && d != 2 {
            return Err(format!("latlong requires d=2 coordinates, got d={d}"));
        }
        if !avc.is_finite() || avc <= 0.001 || avc >= 0.99 {
            return Err(format!("avc must lie strictly in (0.001, 0.99), got {avc}"));
        }
        if b_mat_name.is_empty()
            || bread_mat_name.is_empty()
            || stats_mat_name.is_empty()
            || cvs_mat_name.is_empty()
            || c0_scal_name.is_empty()
            || cv_scal_name.is_empty()
            || q_scal_name.is_empty()
            || tied_scal_name.is_some_and(str::is_empty)
        {
            return Err("matrix and scalar names passed to exact_scpc must be non-empty".to_string());
        }

        Ok(Self {
            d,
            p,
            k_use,
            has_cons,
            avc,
            latlong,
            uncond,
            cvs,
            nthreads,
            b_mat_name,
            bread_mat_name,
            stats_mat_name,
            cvs_mat_name,
            c0_scal_name,
            cv_scal_name,
            q_scal_name,
            tied_scal_name,
            has_cluster,
            n_endog,
            n_z_excl,
        })
    }
}

struct LoadedStataScpc {
    n_sites: usize,
    metric: Metric,
    coords: Vec<f64>,
    b: Vec<f64>,
    cluster_id: Option<Vec<usize>>,
    inputs: PreparedInputs,
}

fn load_stata_scpc_inputs(parsed: &ExactScpcArgs<'_>) -> Result<LoadedStataScpc, String> {
    let nx = parsed.p - usize::from(parsed.has_cons);
    let cl_vars = usize::from(parsed.has_cluster);
    let iv_vars = if !parsed.uncond && parsed.n_endog > 0 {
        parsed.n_z_excl
    } else {
        0
    };
    let expected_vars =
        parsed.d + parsed.p + cl_vars + if parsed.uncond { 0 } else { nx + iv_vars };
    let actual_vars = plugin_var_count();
    if actual_vars != expected_vars {
        return Err(format!(
            "scpc plugin expected {expected_vars} variables in plugin varlist, got {actual_vars}"
        ));
    }

    let obs = active_obs();
    let m_obs = obs.len();
    if m_obs < 2 {
        return Err(format!("too few active observations ({m_obs}) for SCPC"));
    }

    // Read all Stata dataset variables and matrices on the calling thread before Rayon.
    let coords_obs = read_var_matrix(&obs, 1, parsed.d, "coordinates (s_*)")?;
    let scores_obs = read_var_matrix(&obs, 1 + parsed.d, parsed.p, "scores (scpc_score*)")?;
    let cl_raw = if parsed.has_cluster {
        read_var_matrix(&obs, 1 + parsed.d + parsed.p, 1, "cluster index (cl_idx)")?
    } else {
        Vec::new()
    };
    let x_raw = if !parsed.uncond && nx > 0 {
        read_var_matrix(
            &obs,
            1 + parsed.d + parsed.p + cl_vars,
            nx,
            "conditional regressors",
        )?
    } else {
        Vec::new()
    };
    let z_excl = if !parsed.uncond && parsed.n_endog > 0 && parsed.n_z_excl > 0 {
        read_var_matrix(
            &obs,
            1 + parsed.d + parsed.p + cl_vars + nx,
            parsed.n_z_excl,
            "excluded IV instruments",
        )?
    } else {
        Vec::new()
    };

    let (n_sites, cluster_id_vec, coords, scores) = if parsed.has_cluster {
        let (nc, cid, cc, sc) =
            aggregate_cluster_inputs(&cl_raw, &coords_obs, &scores_obs, parsed.d, parsed.p)?;
        (nc, Some(cid), cc, sc)
    } else {
        (m_obs, None, coords_obs, scores_obs)
    };

    let b = read_matrix(parsed.b_mat_name, 1, parsed.p)?;
    let bread = read_matrix(parsed.bread_mat_name, parsed.p, parsed.p)?;

    let inputs = if parsed.n_endog > 0 {
        prepare_iv_inputs(
            n_sites,
            m_obs,
            parsed.p,
            parsed.k_use,
            parsed.has_cons,
            parsed.uncond,
            parsed.n_endog,
            parsed.n_z_excl,
            &b,
            &bread,
            &scores,
            &x_raw,
            &z_excl,
            cluster_id_vec.as_deref(),
        )?
    } else {
        prepare_ols_inputs(
            n_sites,
            m_obs,
            parsed.p,
            parsed.k_use,
            parsed.has_cons,
            parsed.uncond,
            &b,
            &bread,
            &scores,
            &x_raw,
            cluster_id_vec.as_deref(),
        )?
    };
    let metric = if parsed.latlong {
        Metric::GreatCircleFractionLatLon
    } else {
        Metric::Euclidean
    };

    Ok(LoadedStataScpc {
        n_sites,
        metric,
        coords,
        b,
        cluster_id: cluster_id_vec,
        inputs,
    })
}

fn write_stata_scpc_outputs(parsed: &ExactScpcArgs<'_>, out: &ExactScpcOutput) -> Result<(), String> {
    write_matrix(parsed.stats_mat_name, parsed.k_use, 6, &out.scpcstats)?;
    if parsed.cvs {
        let cvs_mat = out
            .scpccvs
            .as_deref()
            .ok_or_else(|| "missing scpccvs output when cvs=1".to_string())?;
        write_matrix(parsed.cvs_mat_name, parsed.k_use, 4, cvs_mat)?;
    }
    save_scalar(parsed.c0_scal_name, out.c0)?;
    save_scalar(parsed.cv_scal_name, out.cvfin)?;
    save_scalar(parsed.q_scal_name, out.q as f64)?;
    if let Some(tied_scal_name) = parsed.tied_scal_name {
        save_scalar(tied_scal_name, if out.tied { 1.0 } else { 0.0 })?;
    }
    Ok(())
}

fn handle_exact_scpc(args: &[&str]) -> Result<(), String> {
    let parsed = ExactScpcArgs::parse(args)?;
    let loaded = load_stata_scpc_inputs(&parsed)?;
    let out = run_exact_scpc(&ExactScpcRequest {
        coords: &loaded.coords,
        n: loaded.n_sites,
        d: parsed.d,
        metric: loaded.metric,
        avc: parsed.avc,
        uncond: parsed.uncond,
        cvs: parsed.cvs,
        coef: &loaded.b[..parsed.k_use],
        wj: &loaded.inputs.wj,
        cluster_id: loaded.cluster_id.as_deref(),
        fe_id: None,
        xj: loaded.inputs.xj.as_deref(),
        xjs: loaded.inputs.xjs.as_deref(),
        proj_u: loaded.inputs.proj_u.as_deref(),
        proj_v: loaded.inputs.proj_v.as_deref(),
        proj_rank: loaded.inputs.proj_rank,
        cv_tolerance: DEFAULT_CV_TOL,
        nthreads: parsed.nthreads,
    })
    .map_err(format_core_error)?;
    write_stata_scpc_outputs(&parsed, &out)
}

const DEFAULT_NYSTROM_LANDMARKS: usize = 1000;
const DEFAULT_NYSTROM_SAMPLES: usize = 20;
const DEFAULT_NYSTROM_PAIR_COUNT: usize = 1_000_000;

struct PermutedApproxInputs {
    coords: Vec<f64>,
    wj: Vec<f64>,
    cluster_id: Option<Vec<usize>>,
    xj: Option<Vec<f64>>,
    xjs: Option<Vec<f64>>,
    proj_u: Option<Vec<f64>>,
    proj_v: Option<Vec<f64>>,
}

#[allow(clippy::too_many_arguments)]
fn prepare_permuted_approx_inputs(
    coords: &[f64],
    n_sites: usize,
    d: usize,
    metric: Metric,
    k_use: usize,
    uncond: bool,
    cluster_id: Option<&[usize]>,
    inputs: &PreparedInputs,
) -> Result<PermutedApproxInputs, String> {
    let (coords_norm, perm) =
        normalize_coords(coords, n_sites, d, metric).map_err(format_core_error)?;
    let wj = permute_rows(&inputs.wj, &perm, k_use);

    if uncond {
        return Ok(PermutedApproxInputs {
            coords: coords_norm,
            wj,
            cluster_id: None,
            xj: None,
            xjs: None,
            proj_u: None,
            proj_v: None,
        });
    }

    if let Some(cid) = cluster_id {
        let mut inv_perm = vec![0usize; n_sites];
        for (new_idx, &old_idx) in perm.iter().enumerate() {
            inv_perm[old_idx] = new_idx;
        }
        let remapped = cid.iter().map(|&g| inv_perm[g]).collect();
        return Ok(PermutedApproxInputs {
            coords: coords_norm,
            wj,
            cluster_id: Some(remapped),
            xj: inputs.xj.clone(),
            xjs: inputs.xjs.clone(),
            proj_u: inputs.proj_u.clone(),
            proj_v: inputs.proj_v.clone(),
        });
    }

    let xj = inputs.xj.as_deref().map(|v| permute_rows(v, &perm, k_use));
    let xjs = inputs.xjs.as_deref().map(|v| permute_rows(v, &perm, k_use));
    let proj_u = if inputs.proj_rank > 0 {
        inputs
            .proj_u
            .as_deref()
            .map(|v| permute_rows(v, &perm, inputs.proj_rank))
    } else {
        None
    };
    let proj_v = if inputs.proj_rank > 0 {
        inputs
            .proj_v
            .as_deref()
            .map(|v| permute_rows(v, &perm, inputs.proj_rank))
    } else {
        None
    };
    Ok(PermutedApproxInputs {
        coords: coords_norm,
        wj,
        cluster_id: None,
        xj,
        xjs,
        proj_u,
        proj_v,
    })
}

fn handle_approx_scpc(args: &[&str]) -> Result<(), String> {
    if !(20..=21).contains(&args.len()) {
        return Err(format!(
            "approx_scpc expects 20 or 21 arguments after subcommand, got {}",
            args.len()
        ));
    }
    let parsed = ExactScpcArgs::parse(&args[..20])?;
    let seed = match args.get(20) {
        Some(raw) => parse_seed(raw)?,
        None => 1,
    };
    let loaded = load_stata_scpc_inputs(&parsed)?;
    let perm_inputs = prepare_permuted_approx_inputs(
        &loaded.coords,
        loaded.n_sites,
        parsed.d,
        loaded.metric,
        parsed.k_use,
        parsed.uncond,
        loaded.cluster_id.as_deref(),
        &loaded.inputs,
    )?;

    match run_approx_scpc(&ApproxScpcRequest {
        coords: &perm_inputs.coords,
        n: loaded.n_sites,
        d: parsed.d,
        metric: loaded.metric,
        avc: parsed.avc,
        uncond: parsed.uncond,
        cvs: parsed.cvs,
        coef: &loaded.b[..parsed.k_use],
        wj: &perm_inputs.wj,
        cluster_id: perm_inputs.cluster_id.as_deref(),
        fe_id: None,
        xj: perm_inputs.xj.as_deref(),
        xjs: perm_inputs.xjs.as_deref(),
        proj_u: perm_inputs.proj_u.as_deref(),
        proj_v: perm_inputs.proj_v.as_deref(),
        proj_rank: loaded.inputs.proj_rank,
        cv_tolerance: DEFAULT_CV_TOL,
        nthreads: parsed.nthreads,
        seed,
        landmarks: DEFAULT_NYSTROM_LANDMARKS,
        samples: DEFAULT_NYSTROM_SAMPLES,
        pair_count: DEFAULT_NYSTROM_PAIR_COUNT,
    }) {
        Ok(out) => write_stata_scpc_outputs(&parsed, &out),
        Err(err) if err.kind() == ErrorKind::TiedEigenspace => {
            if let Some(tied_scal_name) = parsed.tied_scal_name {
                save_scalar(tied_scal_name, 1.0)?;
                Ok(())
            } else {
                Err(format_core_error(err))
            }
        }
        Err(err) => Err(format_core_error(err)),
    }
}

fn dispatch(args: &[&str]) -> Result<(), String> {
    let subcmd = *args
        .first()
        .ok_or_else(|| "missing plugin subcommand".to_string())?;
    match subcmd {
        "ping" | "ping_v2" | "ping_v3" | "ping_v4" | "ping_v5" | "ping_spur_v1" | "ping_spur_v2" => Ok(()),
        "exact_scpc" => handle_exact_scpc(&args[1..]),
        "approx_scpc" => handle_approx_scpc(&args[1..]),
        "spur_i0" => spur_stata::handle_spur_test(&args[1..], spur_stata::SpurTestFlavor::I0),
        "spur_i1" => spur_stata::handle_spur_test(&args[1..], spur_stata::SpurTestFlavor::I1),
        "spur_i0_resid" => {
            spur_stata::handle_spur_test(&args[1..], spur_stata::SpurTestFlavor::I0Resid)
        }
        "spur_i1_resid" => {
            spur_stata::handle_spur_test(&args[1..], spur_stata::SpurTestFlavor::I1Resid)
        }
        "spur_halflife" => spur_stata::handle_spur_halflife(&args[1..]),
        "spur_transform" => spur_stata::handle_spur_transform(&args[1..]),
        "spur_pipeline" => spur_stata::handle_spur_pipeline(&args[1..]),
        other => Err(format!("unknown spatial-stata subcommand '{other}'")),
    }
}

/// Entry point invoked from `c/stata_bridge.c`.
///
/// # Safety
/// `argv` must point to at least `argc` valid NUL-terminated C string pointers when `argc > 0`.
#[no_mangle]
pub unsafe extern "C" fn rs_stata_call(argc: c_int, argv: *const *const c_char) -> c_int {
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

        let inputs = prepare_ols_inputs(n, n, p, p, true, false, &b, &bread, &scores, &x_raw, None).unwrap();
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
            cluster_id: None,
            fe_id: None,
            xj: inputs.xj.as_deref(),
            xjs: inputs.xjs.as_deref(),
            proj_u: inputs.proj_u.as_deref(),
            proj_v: inputs.proj_v.as_deref(),
            proj_rank: inputs.proj_rank,
            cv_tolerance: DEFAULT_CV_TOL,
            nthreads: Some(2),
        })
        .unwrap();
        assert_eq!(out_euclid.q, 19);
        let dist_euclid =
            spatial_core::distance::distance_matrix(&coords_lon_lat, n, 2, Metric::Euclidean)
                .unwrap();
        let mut dvec_euclid = Vec::with_capacity(n * (n - 1) / 2);
        for i in 1..n {
            dvec_euclid.extend_from_slice(&dist_euclid[i * n..i * n + i]);
        }
        assert!(
            (spatial_core::calibration::avg_correlation(out_euclid.c0, &dvec_euclid) - 0.03).abs()
                < 1e-12
        );
        assert!((out_euclid.c0 - 3.500409719601247).abs() < 1e-3);
        assert!((out_euclid.cvfin - 2.1008365544083087).abs() < 1e-3);
        assert!(!out_euclid.tied);
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
            cluster_id: None,
            fe_id: None,
            xj: inputs.xj.as_deref(),
            xjs: inputs.xjs.as_deref(),
            proj_u: inputs.proj_u.as_deref(),
            proj_v: inputs.proj_v.as_deref(),
            proj_rank: inputs.proj_rank,
            cv_tolerance: DEFAULT_CV_TOL,
            nthreads: Some(2),
        })
        .unwrap();
        assert_eq!(out_geo.q, 18);
        let dist_geo = spatial_core::distance::distance_matrix(
            &coords_lat_lon,
            n,
            2,
            Metric::GreatCircleFractionLatLon,
        )
        .unwrap();
        let mut dvec_geo = Vec::with_capacity(n * (n - 1) / 2);
        for i in 1..n {
            dvec_geo.extend_from_slice(&dist_geo[i * n..i * n + i]);
        }
        assert!(
            (spatial_core::calibration::avg_correlation(out_geo.c0, &dvec_geo) - 0.03).abs()
                < 1e-12
        );
        assert!(
            (out_geo.c0 - 1396.4193247376622).abs() < 1.0,
            "geo c0={}",
            out_geo.c0
        );
        assert!(
            (out_geo.cvfin - 2.1009220402410382).abs() < 5e-3,
            "geo cvfin={}",
            out_geo.cvfin
        );
        assert!(!out_geo.tied);

        // First 20 rows of bench_200 exhibit tied eigenspace for automatic fallback.
        let n_tie = 20;
        let mut xtx_20 = vec![0.0; p * p];
        let mut xty_20 = vec![0.0; p];
        for i in 0..n_tie {
            for r in 0..p {
                xty_20[r] += x_design[i * p + r] * y[i];
                for c in 0..p {
                    xtx_20[r * p + c] += x_design[i * p + r] * x_design[i * p + c];
                }
            }
        }
        let bread_20 = invert_symmetric(&xtx_20, p);
        let mut b_20 = vec![0.0; p];
        for r in 0..p {
            b_20[r] = (0..p).map(|c| bread_20[r * p + c] * xty_20[c]).sum();
        }
        let mut scores_20 = vec![0.0; n_tie * p];
        for i in 0..n_tie {
            let yhat: f64 = (0..p).map(|c| x_design[i * p + c] * b_20[c]).sum();
            let resid = y[i] - yhat;
            for c in 0..p {
                scores_20[i * p + c] = x_design[i * p + c] * resid;
            }
        }
        let inputs_20 = prepare_ols_inputs(
            n_tie,
            n_tie,
            p,
            p,
            true,
            false,
            &b_20,
            &bread_20,
            &scores_20,
            &x_raw[..n_tie * nx],
            None,
        )
        .unwrap();
        let out_20 = run_exact_scpc(&ExactScpcRequest {
            coords: &coords_lon_lat[..n_tie * 2],
            n: n_tie,
            d: 2,
            metric: Metric::Euclidean,
            avc: 0.03,
            uncond: false,
            cvs: false,
            coef: &b_20,
            wj: &inputs_20.wj,
            cluster_id: None,
            fe_id: None,
            xj: inputs_20.xj.as_deref(),
            xjs: inputs_20.xjs.as_deref(),
            proj_u: inputs_20.proj_u.as_deref(),
            proj_v: inputs_20.proj_v.as_deref(),
            proj_rank: inputs_20.proj_rank,
            cv_tolerance: DEFAULT_CV_TOL,
            nthreads: Some(1),
        })
        .unwrap();
        assert_eq!(out_20.q, 7);
        assert!(out_20.tied);
    }

    #[test]
    fn orthonormal_basis_is_scale_invariant_and_drops_collinear_columns() {
        let n = 5;
        let base = [1.0, -2.0, 0.5, 3.0, -1.5];
        // Small-unit regressor 1e-12 * x and extreme tiny 1e-200 * x stay full rank with _cons.
        for scale in [1e-12, 1e-200, 1e12] {
            let x_raw: Vec<f64> = base.iter().map(|&v| v * scale).collect();
            let (q, rank) = orthonormal_basis_mgs(&x_raw, n, 1, true);
            assert_eq!(rank, 2, "scale={scale}");
            for r1 in 0..2 {
                for r2 in 0..2 {
                    let dot: f64 = (0..n).map(|i| q[i * 2 + r1] * q[i * 2 + r2]).sum();
                    let exp = if r1 == r2 { 1.0 } else { 0.0 };
                    assert!((dot - exp).abs() < 1e-14, "scale={scale}, dot={dot}");
                }
            }
        }

        // Duplicate/collinear and zero columns drop cleanly while scaled regressor remains.
        let mut x_raw = Vec::with_capacity(n * 3);
        for &v in &base {
            x_raw.extend_from_slice(&[1e-12 * v, 2e-12 * v, 0.0]);
        }
        let (q, rank) = orthonormal_basis_mgs(&x_raw, n, 3, true);
        assert_eq!(rank, 2);
        for r1 in 0..2 {
            for r2 in 0..2 {
                let dot: f64 = (0..n).map(|i| q[i * 2 + r1] * q[i * 2 + r2]).sum();
                let exp = if r1 == r2 { 1.0 } else { 0.0 };
                assert!((dot - exp).abs() < 1e-14);
            }
        }
    }

    #[test]
    fn dispatch_and_exact_scpc_args_support_v3_v4_and_cluster_flag() {
        assert!(dispatch(&["ping"]).is_ok());
        assert!(dispatch(&["ping_v2"]).is_ok());
        assert!(dispatch(&["ping_v3"]).is_ok());
        assert!(dispatch(&["ping_v4"]).is_ok());
        assert!(dispatch(&[]).is_err());
        assert!(dispatch(&["unknown"]).is_err());

        let base16 = [
            "2", "3", "3", "1", "0.03", "0", "0", "1", "-1", "b", "bread", "stats", "cvs", "c0",
            "cv", "q",
        ];
        let parsed16 = ExactScpcArgs::parse(&base16).unwrap();
        assert_eq!(parsed16.tied_scal_name, None);
        assert!(!parsed16.has_cluster);

        let mut args17 = base16.to_vec();
        args17.push("tied");
        let parsed17 = ExactScpcArgs::parse(&args17).unwrap();
        assert_eq!(parsed17.tied_scal_name, Some("tied"));
        assert!(!parsed17.has_cluster);

        let mut args18_ok = args17.clone();
        args18_ok.push("1");
        let parsed18 = ExactScpcArgs::parse(&args18_ok).unwrap();
        assert!(parsed18.has_cluster);

        // 17 args enters handle_exact_scpc and reaches plugin_var_count() check (0 under test stubs).
        let err = handle_exact_scpc(&args17).unwrap_err();
        assert!(err.contains("variables in plugin varlist"), "err={err}");

        // Wrong argument counts or invalid cluster flags fail cleanly.
        assert!(handle_exact_scpc(&base16[..15]).is_err());
        let mut args18_bad = args17.clone();
        args18_bad.push("extra");
        assert!(handle_exact_scpc(&args18_bad).is_err());
        let mut args19 = args18_ok.clone();
        args19.push("extra");
        assert!(handle_exact_scpc(&args19).is_err());

        // Empty tied scalar name or invalid avc/dims fail cleanly.
        let mut bad_tied = base16.to_vec();
        bad_tied.push("");
        assert!(ExactScpcArgs::parse(&bad_tied).is_err());

        let mut bad_latlong = args17.clone();
        bad_latlong[0] = "3";
        bad_latlong[5] = "1";
        assert!(ExactScpcArgs::parse(&bad_latlong).is_err());

        let mut bad_avc = args17;
        bad_avc[4] = "0.99";
        assert!(ExactScpcArgs::parse(&bad_avc).is_err());
    }

    #[test]
    fn error_mapping_and_c_entrypoint_behave_safely() {
        assert_eq!(
            format_core_error(CoreError::calibration("bad target")),
            "calibration failed: bad target"
        );
        assert_eq!(
            format_core_error(CoreError::degenerate("zero var")),
            "degenerate variance: zero var"
        );
        assert_eq!(
            format_core_error(CoreError::tied("tie")),
            "tied eigenspace: tie"
        );
        assert_eq!(
            format_core_error(CoreError::unsupported("unsupported")),
            "unsupported case: unsupported"
        );
        assert_eq!(
            format_core_error(CoreError::invalid_input("invalid")),
            "invalid input: invalid"
        );
        assert_eq!(
            format_core_error(CoreError::internal("boom")),
            "internal failure: boom"
        );

        let ping_cstr = std::ffi::CString::new("ping_v3").unwrap();
        let argv_ok = [ping_cstr.as_ptr()];
        assert_eq!(unsafe { rs_stata_call(1, argv_ok.as_ptr()) }, STATA_OK);
        assert_eq!(unsafe { rs_stata_call(0, std::ptr::null()) }, STATA_ERR);
        let argv_null = [std::ptr::null()];
        assert_eq!(unsafe { rs_stata_call(1, argv_null.as_ptr()) }, STATA_ERR);
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

    #[test]
    fn iv_preparation_matches_direct_2sls_projection() {
        let m_obs = 8;
        let p = 3; // [x_endog, w_exog, _cons]
        let nx = 2;
        let mut x_raw = Vec::with_capacity(m_obs * nx);
        let mut z_excl = Vec::with_capacity(m_obs);
        let mut y = Vec::with_capacity(m_obs);
        for i in 0..m_obs {
            let fi = i as f64;
            let wi = (fi * 0.7).sin();
            let zi = (fi * 1.1).cos() + 0.2 * fi;
            let xi = 0.8 * zi + 0.3 * wi + 0.25 * (fi * 1.9).sin();
            let yi = 1.0 + 1.4 * xi + 0.5 * wi + 0.3 * (fi * 1.9).sin();
            x_raw.extend_from_slice(&[xi, wi]);
            z_excl.push(zi);
            y.push(yi);
        }
        let b = [1.4, 0.5, 1.0];
        let mut bread = vec![0.0; p * p];
        for i in 0..p {
            bread[i * p + i] = 1.0;
        }
        let mut scores = vec![0.0; m_obs * p];
        for i in 0..m_obs {
            scores[i * p] = z_excl[i] * 0.1;
            scores[i * p + 1] = x_raw[i * nx + 1] * 0.1;
            scores[i * p + 2] = 0.1;
        }
        let inputs = prepare_iv_inputs(
            m_obs, m_obs, p, p, true, false, 1, 1, &b, &bread, &scores, &x_raw, &z_excl, None,
        )
        .unwrap();
        let u = inputs.proj_u.as_ref().unwrap();
        let v = inputs.proj_v.as_ref().unwrap();
        let xj = inputs.xj.as_ref().unwrap();
        let xjs = inputs.xjs.as_ref().unwrap();
        assert_eq!(inputs.proj_rank, p);
        // Verify factor orientation corresponds to M' = I - U V' where U = P_Z X A^{-1}, V = X:
        // V' * U == I_p and each xjs[i] == sign(xj[i]).
        for r in 0..p {
            for c in 0..p {
                let dot: f64 = (0..m_obs).map(|i| v[i * p + r] * u[i * p + c]).sum();
                let exp = if r == c { 1.0 } else { 0.0 };
                assert!((dot - exp).abs() < 1e-11, "V'U[{r},{c}]={dot}");
            }
        }
        // Check M' projects onto P_Z X's nullspace via column orthogonality: for any A, X' (A - U V' A) == 0.
        for i in 0..m_obs {
            for j in 0..p {
                let expected_sign = if xj[i * p + j] > 0.0 {
                    1.0
                } else if xj[i * p + j] < 0.0 {
                    -1.0
                } else {
                    0.0
                };
                assert_eq!(xjs[i * p + j], expected_sign);
            }
        }
        // V columns store raw X (intercept last in Stata), whereas U stores P_Z X A^{-1}.
        for i in 0..m_obs {
            assert_eq!(v[i * p], x_raw[i * nx]);
            assert_eq!(v[i * p + 1], x_raw[i * nx + 1]);
            assert_eq!(v[i * p + 2], 1.0);
        }
        assert!(dispatch(&["ping_v5"]).is_ok());
    }

    #[test]
    fn approx_scpc_args_normalization_permutation_and_cluster_remapping() {
        assert!(parse_seed("1").unwrap() == 1);
        assert!(parse_seed("4294967295").unwrap() == u32::MAX);
        assert!(parse_seed("-1").is_err());
        assert!(parse_seed("4294967296").is_err());

        let base20 = [
            "2", "3", "3", "1", "0.05", "0", "0", "1", "-1", "b", "bread", "stats", "cvs", "c0",
            "cv", "q", "tied", "1", "0", "0",
        ];
        assert!(handle_approx_scpc(&base20[..19]).is_err());
        let err = handle_approx_scpc(&base20).unwrap_err();
        assert!(err.contains("variables in plugin varlist"), "err={err}");

        let mut args21 = base20.to_vec();
        args21.push("7");
        let err21 = handle_approx_scpc(&args21).unwrap_err();
        assert!(err21.contains("variables in plugin varlist"), "err21={err21}");

        // Synthetic clustered OLS through prepare_ols_inputs + prepare_permuted_approx_inputs
        let n_sites = 12usize;
        let per_cl = 3usize;
        let m_obs = n_sites * per_cl;
        let p = 2usize;
        let mut coords_cl = Vec::with_capacity(n_sites * 2);
        for g in 0..n_sites {
            // Reversed order so normalize_coords produces a non-identity permutation.
            let gf = (n_sites - 1 - g) as f64;
            coords_cl.push(gf * 0.25);
            coords_cl.push((gf * 1.7).sin() * 0.3);
        }
        let mut cid = Vec::with_capacity(m_obs);
        let mut x_raw = Vec::with_capacity(m_obs);
        let mut scores_cl = vec![0.0; n_sites * p];
        for g in 0..n_sites {
            for r in 0..per_cl {
                let obs_i = g * per_cl + r;
                let xv = (obs_i as f64 * 0.41).cos() + 0.2 * (g as f64);
                let u = ((obs_i + 1) as f64 * 0.73).sin() * 0.5;
                cid.push(g);
                x_raw.push(xv);
                scores_cl[g * p] += u * xv;
                scores_cl[g * p + 1] += u;
            }
        }
        let b = [1.2, 0.4];
        let bread = [1.0 / m_obs as f64, 0.0, 0.0, 1.0 / m_obs as f64];
        let inputs = prepare_ols_inputs(
            n_sites, m_obs, p, p, true, false, &b, &bread, &scores_cl, &x_raw, Some(&cid),
        )
        .unwrap();
        let perm_cl = prepare_permuted_approx_inputs(
            &coords_cl,
            n_sites,
            2,
            Metric::GreatCircleFractionLatLon,
            p,
            false,
            Some(&cid),
            &inputs,
        )
        .unwrap();
        let out_1 = run_approx_scpc(&ApproxScpcRequest {
            coords: &perm_cl.coords,
            n: n_sites,
            d: 2,
            metric: Metric::GreatCircleFractionLatLon,
            avc: 0.1,
            uncond: false,
            cvs: true,
            coef: &b,
            wj: &perm_cl.wj,
            cluster_id: perm_cl.cluster_id.as_deref(),
            fe_id: None,
            xj: perm_cl.xj.as_deref(),
            xjs: perm_cl.xjs.as_deref(),
            proj_u: perm_cl.proj_u.as_deref(),
            proj_v: perm_cl.proj_v.as_deref(),
            proj_rank: inputs.proj_rank,
            cv_tolerance: DEFAULT_CV_TOL,
            nthreads: Some(1),
            seed: 5,
            landmarks: DEFAULT_NYSTROM_LANDMARKS,
            samples: DEFAULT_NYSTROM_SAMPLES,
            pair_count: DEFAULT_NYSTROM_PAIR_COUNT,
        })
        .unwrap();
        let out_2 = run_approx_scpc(&ApproxScpcRequest {
            coords: &perm_cl.coords,
            n: n_sites,
            d: 2,
            metric: Metric::GreatCircleFractionLatLon,
            avc: 0.1,
            uncond: false,
            cvs: true,
            coef: &b,
            wj: &perm_cl.wj,
            cluster_id: perm_cl.cluster_id.as_deref(),
            fe_id: None,
            xj: perm_cl.xj.as_deref(),
            xjs: perm_cl.xjs.as_deref(),
            proj_u: perm_cl.proj_u.as_deref(),
            proj_v: perm_cl.proj_v.as_deref(),
            proj_rank: inputs.proj_rank,
            cv_tolerance: DEFAULT_CV_TOL,
            nthreads: Some(2),
            seed: 5,
            landmarks: DEFAULT_NYSTROM_LANDMARKS,
            samples: DEFAULT_NYSTROM_SAMPLES,
            pair_count: DEFAULT_NYSTROM_PAIR_COUNT,
        })
        .unwrap();
        assert_eq!(out_1.q, out_2.q);
        for (&a, &b) in out_1.scpcstats.iter().zip(out_2.scpcstats.iter()) {
            assert!((a - b).abs() < 1e-12);
        }
    }
}
