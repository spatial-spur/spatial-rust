//! Stata SPI v3 handlers for SPUR diagnostics, half-life CI, spatial transforms, and pipeline.

use crate::stata_spi::{
    active_obs, plugin_var_count, read_f64_bin_file, read_var_matrix, save_scalar,
    save_scalar_or_missing, write_matrix, write_matrix_with_infinity, write_var_matrix,
};
use crate::{format_core_error, parse_bool_flag, parse_nthreads, parse_usize};
use spatial_core::{
    run_spur_halflife, run_spur_i0, run_spur_i0_resid, run_spur_i1, run_spur_i1_resid,
    run_spur_pipeline, run_spur_transform, Metric, SpurHalfLifeRequest, SpurPipelineRequest,
    SpurTestRequest, SpurTransformKind, SpurTransformRequest,
};
use std::collections::HashMap;

/// SPUR host convention for converting great-circle fractions to meters.
use spatial_core::EARTH_CIRCUMFERENCE_M;
const HALF_LIFE_COMPONENT_ROWS: usize = 65;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum SpurTestFlavor {
    I0,
    I1,
    I0Resid,
    I1Resid,
}

impl SpurTestFlavor {
    fn is_residual(self) -> bool {
        matches!(self, Self::I0Resid | Self::I1Resid)
    }

    fn label(self) -> &'static str {
        match self {
            Self::I0 => "spur_i0",
            Self::I1 => "spur_i1",
            Self::I0Resid => "spur_i0_resid",
            Self::I1Resid => "spur_i1_resid",
        }
    }
}

fn parse_spur_metric(d: usize, latlong: bool) -> Result<Metric, String> {
    if d == 0 {
        return Err("coordinate dimension d must be at least 1".to_string());
    }
    if latlong {
        if d != 2 {
            return Err(format!("latlong requires d=2 coordinates, got d={d}"));
        }
        Ok(Metric::GreatCircleFractionLatLon)
    } else {
        Ok(Metric::Euclidean)
    }
}

fn require_non_empty(val: &str, name: &str, context: &str) -> Result<(), String> {
    if val.is_empty() {
        return Err(format!("{name} passed to {context} must be non-empty"));
    }
    Ok(())
}

fn read_emat_checked(path: &str, q: usize, nrep: usize) -> Result<Vec<f64>, String> {
    let expected = q
        .checked_mul(nrep)
        .ok_or_else(|| "q * nrep size overflow".to_string())?;
    if expected == 0 {
        return Err("q and nrep must be at least 1".to_string());
    }
    read_f64_bin_file(path, expected, "emat")
}

/// Prepend an intercept column of `1.0` when `include_intercept` is true, matching Stata `[1, X_raw]`.
fn build_x_in(
    x_raw: &[f64],
    n: usize,
    k_x_raw: usize,
    include_intercept: bool,
) -> Result<(Vec<f64>, usize), String> {
    let k_x = k_x_raw + usize::from(include_intercept);
    if k_x == 0 {
        return Err("residual SPUR test requires at least one regressor or intercept".to_string());
    }
    let expected_raw = n
        .checked_mul(k_x_raw)
        .ok_or_else(|| "x_raw size overflow".to_string())?;
    if x_raw.len() != expected_raw {
        return Err(format!(
            "x_raw length {} does not match {n}x{k_x_raw}",
            x_raw.len()
        ));
    }
    let out_len = n
        .checked_mul(k_x)
        .ok_or_else(|| "X_in size overflow".to_string())?;
    let mut out = Vec::with_capacity(out_len);
    for i in 0..n {
        if include_intercept {
            out.push(1.0);
        }
        if k_x_raw > 0 {
            out.extend_from_slice(&x_raw[i * k_x_raw..(i + 1) * k_x_raw]);
        }
    }
    Ok((out, k_x))
}

/// Map arbitrary finite cluster values (including gapped group IDs from `egen = group()` with `if`)
/// into contiguous 0-based cluster indices `0..n_clusters`.
fn normalize_cluster_ids(raw: &[f64]) -> Result<Vec<usize>, String> {
    let mut map: HashMap<u64, usize> = HashMap::new();
    let mut ids = Vec::with_capacity(raw.len());
    for &val in raw {
        if !val.is_finite() {
            return Err("missing or non-finite cluster value".to_string());
        }
        let key = if val == 0.0 {
            0.0_f64.to_bits()
        } else {
            val.to_bits()
        };
        let next_id = map.len();
        let id = *map.entry(key).or_insert(next_id);
        ids.push(id);
    }
    Ok(ids)
}

struct SpurTestSummary {
    lr: f64,
    pvalue: f64,
    cvalue: Vec<f64>,
    ha_param: f64,
    info: [f64; 3],
}

fn diagnostic_info(tied: bool, all: bool, best: Option<f64>) -> [f64; 3] {
    [
        f64::from(tied),
        f64::from(all),
        best.unwrap_or(f64::INFINITY),
    ]
}

pub(crate) fn handle_spur_test(args: &[&str], flavor: SpurTestFlavor) -> Result<(), String> {
    let ctx = flavor.label();
    if args.len() != 13 {
        return Err(format!(
            "{ctx} expects 13 arguments after subcommand, got {}",
            args.len()
        ));
    }

    let d = parse_usize(args[0], "d")?;
    let k_x_raw = parse_usize(args[1], "k_x_raw")?;
    let include_intercept = parse_bool_flag(args[2], "include_intercept")?;
    let latlong = parse_bool_flag(args[3], "latlong")?;
    let q = parse_usize(args[4], "q")?;
    let nrep = parse_usize(args[5], "nrep")?;
    let nthreads = parse_nthreads(args[6])?;
    let emat_path = args[7];
    let cv_mat_name = args[8];
    let stat_scal_name = args[9];
    let p_scal_name = args[10];
    let ha_scal_name = args[11];
    let info_mat_name = args[12];

    let metric = parse_spur_metric(d, latlong)?;
    if q == 0 || nrep == 0 {
        return Err(format!("invalid simulation dimensions: q={q}, nrep={nrep}"));
    }
    if !flavor.is_residual() && (k_x_raw != 0 || include_intercept) {
        return Err(format!("{ctx} requires k_x_raw=0 and include_intercept=0"));
    }
    if flavor.is_residual() && k_x_raw + usize::from(include_intercept) == 0 {
        return Err(format!(
            "{ctx} requires at least one regressor or intercept"
        ));
    }
    require_non_empty(emat_path, "emat_path", ctx)?;
    require_non_empty(cv_mat_name, "cv_mat_name", ctx)?;
    require_non_empty(stat_scal_name, "stat_scal_name", ctx)?;
    require_non_empty(p_scal_name, "p_scal_name", ctx)?;
    require_non_empty(ha_scal_name, "ha_scal_name", ctx)?;
    require_non_empty(info_mat_name, "info_mat_name", ctx)?;

    let expected_vars = d + 1 + k_x_raw;
    let actual_vars = plugin_var_count();
    if actual_vars != expected_vars {
        return Err(format!(
            "{ctx} expected {expected_vars} variables in plugin varlist, got {actual_vars}"
        ));
    }

    let obs = active_obs();
    let n = obs.len();
    if n < 3 {
        return Err(format!("too few active observations ({n}) for {ctx}"));
    }

    let coords = read_var_matrix(&obs, 1, d, "coordinates (s_*)")?;
    let y = read_var_matrix(&obs, 1 + d, 1, "response (y)")?;
    let x_raw = if k_x_raw > 0 {
        read_var_matrix(&obs, 2 + d, k_x_raw, "regressors (X)")?
    } else {
        Vec::new()
    };
    let (x_in_buf, k_x) = if flavor.is_residual() {
        let (buf, kx) = build_x_in(&x_raw, n, k_x_raw, include_intercept)?;
        (Some(buf), kx)
    } else {
        (None, 0)
    };
    let emat = read_emat_checked(emat_path, q.checked_mul(2).ok_or("q overflow")?, nrep)?;

    let req = SpurTestRequest {
        coords: &coords,
        n,
        d,
        metric,
        y: &y,
        n_y: 1,
        x_in: x_in_buf.as_deref(),
        k_x,
        emat: &emat,
        q,
        nrep,
        nthreads,
    };

    let summary = match flavor {
        SpurTestFlavor::I0 => {
            let out = run_spur_i0(&req).map_err(format_core_error)?;
            SpurTestSummary {
                lr: out.lr[0],
                pvalue: out.pvalue[0],
                cvalue: out.cvalue,
                ha_param: out.ha_param,
                info: diagnostic_info(out.tied, out.all_directions, out.best_tuning_power),
            }
        }
        SpurTestFlavor::I1 => {
            let out = run_spur_i1(&req).map_err(format_core_error)?;
            SpurTestSummary {
                lr: out.lr[0],
                pvalue: out.pvalue[0],
                cvalue: out.cvalue,
                ha_param: out.ha_param,
                info: diagnostic_info(out.tied, out.all_directions, out.best_tuning_power),
            }
        }
        SpurTestFlavor::I0Resid => {
            let out = run_spur_i0_resid(&req).map_err(format_core_error)?;
            SpurTestSummary {
                lr: out.lr[0],
                pvalue: out.pvalue[0],
                cvalue: out.cvalue,
                ha_param: out.ha_param,
                info: diagnostic_info(out.tied, out.all_directions, out.best_tuning_power),
            }
        }
        SpurTestFlavor::I1Resid => {
            let out = run_spur_i1_resid(&req).map_err(format_core_error)?;
            SpurTestSummary {
                lr: out.lr[0],
                pvalue: out.pvalue[0],
                cvalue: out.cvalue,
                ha_param: out.ha_param,
                info: diagnostic_info(out.tied, out.all_directions, out.best_tuning_power),
            }
        }
    };

    write_matrix_with_infinity(info_mat_name, 1, 3, &summary.info)?;
    write_matrix_with_infinity(cv_mat_name, 3, 1, &summary.cvalue)?;
    save_scalar(stat_scal_name, summary.lr)?;
    save_scalar(p_scal_name, summary.pvalue)?;
    save_scalar(ha_scal_name, summary.ha_param)?;
    Ok(())
}

pub(crate) fn handle_spur_halflife(args: &[&str]) -> Result<(), String> {
    const CTX: &str = "spur_halflife";
    if args.len() != 13 {
        return Err(format!(
            "{CTX} expects 13 arguments after subcommand, got {}",
            args.len()
        ));
    }

    let d = parse_usize(args[0], "d")?;
    let latlong = parse_bool_flag(args[1], "latlong")?;
    let q = parse_usize(args[2], "q")?;
    if q > 100 {
        return Err(format!(
            "q={q} exceeds the validated half-life quadrature limit of 100"
        ));
    }
    let nrep = parse_usize(args[3], "nrep")?;
    let level_pct = args[4]
        .parse::<f64>()
        .map_err(|_| format!("invalid level '{}'", args[4]))?;
    let nthreads = parse_nthreads(args[5])?;
    let emat_path = args[6];
    let ci_l_scal_name = args[7];
    let ci_u_scal_name = args[8];
    let max_dist_scal_name = args[9];
    let info_mat_name = args[10];
    let components_name = args[11];
    let n_components_name = args[12];

    let metric = parse_spur_metric(d, latlong)?;
    if q == 0 || nrep == 0 {
        return Err(format!("invalid simulation dimensions: q={q}, nrep={nrep}"));
    }
    if !level_pct.is_finite() || level_pct <= 0.0 || level_pct >= 100.0 {
        return Err(format!(
            "level must lie strictly in (0, 100), got {level_pct}"
        ));
    }
    require_non_empty(emat_path, "emat_path", CTX)?;
    require_non_empty(ci_l_scal_name, "ci_l_scal_name", CTX)?;
    require_non_empty(ci_u_scal_name, "ci_u_scal_name", CTX)?;
    require_non_empty(max_dist_scal_name, "max_dist_scal_name", CTX)?;
    require_non_empty(info_mat_name, "info_mat_name", CTX)?;

    let expected_vars = d + 1;
    let actual_vars = plugin_var_count();
    if actual_vars != expected_vars {
        return Err(format!(
            "{CTX} expected {expected_vars} variables in plugin varlist, got {actual_vars}"
        ));
    }

    let obs = active_obs();
    let n = obs.len();
    if n < 3 {
        return Err(format!("too few active observations ({n}) for {CTX}"));
    }

    let coords = read_var_matrix(&obs, 1, d, "coordinates (s_*)")?;
    let y = read_var_matrix(&obs, 1 + d, 1, "response (y)")?;
    let emat = read_emat_checked(emat_path, q, nrep)?;

    let out = run_spur_halflife(&SpurHalfLifeRequest {
        coords: &coords,
        n,
        d,
        metric,
        y: &y,
        emat: &emat,
        q,
        nrep,
        level: level_pct / 100.0,
        nthreads,
    })
    .map_err(format_core_error)?;

    if out.components.len() > HALF_LIFE_COMPONENT_ROWS {
        return Err(format!(
            "detected {} half-life components, exceeding Stata's output capacity of {HALF_LIFE_COMPONENT_ROWS}",
            out.components.len()
        ));
    }
    write_matrix_with_infinity(
        info_mat_name,
        1,
        3,
        &diagnostic_info(out.tied, out.all_directions, None),
    )?;
    let max_dist_out = if latlong {
        out.max_dist * EARTH_CIRCUMFERENCE_M
    } else {
        out.max_dist
    };
    save_scalar(max_dist_scal_name, max_dist_out)?;
    save_scalar_or_missing(
        ci_l_scal_name,
        out.ci_lower.is_finite().then_some(out.ci_lower),
    )?;
    save_scalar_or_missing(
        ci_u_scal_name,
        out.ci_upper.is_finite().then_some(out.ci_upper),
    )?;
    save_scalar(n_components_name, out.components.len() as f64)?;
    let mut data = vec![f64::INFINITY; HALF_LIFE_COMPONENT_ROWS * 2];
    for (row, bounds) in data.chunks_exact_mut(2).zip(&out.components) {
        row.copy_from_slice(bounds);
    }
    write_matrix_with_infinity(components_name, HALF_LIFE_COMPONENT_ROWS, 2, &data)
}

pub(crate) fn handle_spur_transform(args: &[&str]) -> Result<(), String> {
    const CTX: &str = "spur_transform";
    if args.len() != 7 {
        return Err(format!(
            "{CTX} expects 7 arguments after subcommand, got {}",
            args.len()
        ));
    }

    let d = parse_usize(args[0], "d")?;
    let kind_str = args[1];
    let latlong = parse_bool_flag(args[2], "latlong")?;
    let radius = args[3]
        .parse::<f64>()
        .map_err(|_| format!("invalid radius '{}'", args[3]))?;
    let has_cluster = parse_bool_flag(args[4], "has_cluster")?;
    let n_y = parse_usize(args[5], "n_y")?;
    let nthreads = parse_nthreads(args[6])?;

    if n_y == 0 {
        return Err("n_y must be at least 1".to_string());
    }

    let kind = match kind_str {
        "nn" => {
            if has_cluster {
                return Err("has_cluster=1 is only allowed for cluster transform".to_string());
            }
            SpurTransformKind::Nn
        }
        "iso" => {
            if has_cluster {
                return Err("has_cluster=1 is only allowed for cluster transform".to_string());
            }
            if !radius.is_finite() || radius <= 0.0 {
                return Err(format!(
                    "iso radius must be positive and finite, got {radius}"
                ));
            }
            let eff_radius = if latlong {
                radius / EARTH_CIRCUMFERENCE_M
            } else {
                radius
            };
            SpurTransformKind::Iso { radius: eff_radius }
        }
        "cluster" => {
            if !has_cluster {
                return Err("cluster transform requires has_cluster=1".to_string());
            }
            SpurTransformKind::Cluster
        }
        "lbmgls" => {
            if has_cluster {
                return Err("has_cluster=1 is only allowed for cluster transform".to_string());
            }
            SpurTransformKind::LbmGls
        }
        other => {
            return Err(format!(
                "invalid transformation '{other}' (expected nn, iso, cluster, or lbmgls)"
            ))
        }
    };

    let metric = if has_cluster {
        Metric::Euclidean
    } else {
        parse_spur_metric(d, latlong)?
    };
    let cl_vars = usize::from(has_cluster);
    let expected_vars = d + cl_vars + 2 * n_y;
    let actual_vars = plugin_var_count();
    if actual_vars != expected_vars {
        return Err(format!(
            "{CTX} expected {expected_vars} variables in plugin varlist, got {actual_vars}"
        ));
    }

    let obs = active_obs();
    let n = obs.len();

    let coords = read_var_matrix(&obs, 1, d, "coordinates (s_*)")?;
    let cluster_ids = if has_cluster {
        let cl_raw = read_var_matrix(&obs, 1 + d, 1, "cluster index")?;
        Some(normalize_cluster_ids(&cl_raw)?)
    } else {
        None
    };
    let y_in = read_var_matrix(&obs, 1 + d + cl_vars, n_y, "input variables (Y)")?;

    let transformed = run_spur_transform(&SpurTransformRequest {
        coords: (d > 0).then_some(&coords),
        n,
        d,
        metric,
        kind,
        cluster_id: cluster_ids.as_deref(),
        y: &y_in,
        n_y,
        nthreads,
    })
    .map_err(format_core_error)?;

    write_var_matrix(
        &obs,
        1 + d + cl_vars + n_y,
        n_y,
        &transformed,
        "transformed variables",
    )
}

pub(crate) fn handle_spur_pipeline(args: &[&str]) -> Result<(), String> {
    const CTX: &str = "spur_pipeline";
    if args.len() != 11 {
        return Err(format!(
            "{CTX} expects 11 arguments after subcommand, got {}",
            args.len()
        ));
    }

    let d = parse_usize(args[0], "d")?;
    let k_x_raw = parse_usize(args[1], "k_x_raw")?;
    let include_intercept = parse_bool_flag(args[2], "include_intercept")?;
    let latlong = parse_bool_flag(args[3], "latlong")?;
    let q = parse_usize(args[4], "q")?;
    let nrep = parse_usize(args[5], "nrep")?;
    let nthreads = parse_nthreads(args[6])?;
    let emat_path = args[7];
    let diag_mat_name = args[8];
    let cvs_mat_name = args[9];
    let info_mat_name = args[10];

    let metric = parse_spur_metric(d, latlong)?;
    if k_x_raw == 0 || q == 0 || nrep == 0 {
        return Err(format!(
            "invalid pipeline dimensions: k_x_raw={k_x_raw}, q={q}, nrep={nrep}"
        ));
    }
    require_non_empty(emat_path, "emat_path", CTX)?;
    require_non_empty(diag_mat_name, "diag_mat_name", CTX)?;
    require_non_empty(cvs_mat_name, "cvs_mat_name", CTX)?;
    require_non_empty(info_mat_name, "info_mat_name", CTX)?;

    let n_vars_tr = 1 + k_x_raw;
    let expected_vars = d + 2 * n_vars_tr;
    let actual_vars = plugin_var_count();
    if actual_vars != expected_vars {
        return Err(format!(
            "{CTX} expected {expected_vars} variables in plugin varlist, got {actual_vars}"
        ));
    }

    let obs = active_obs();
    let n = obs.len();
    if n < 3 {
        return Err(format!("too few active observations ({n}) for {CTX}"));
    }

    let coords = read_var_matrix(&obs, 1, d, "coordinates (s_*)")?;
    let y_and_x = read_var_matrix(&obs, 1 + d, n_vars_tr, "spur variables (depvar indepvars)")?;
    let emat = read_emat_checked(emat_path, q.checked_mul(2).ok_or("q overflow")?, nrep)?;

    let k_x = k_x_raw + usize::from(include_intercept);
    let mut y = Vec::with_capacity(n);
    let mut x_in = Vec::with_capacity(n * k_x);
    for i in 0..n {
        let row = &y_and_x[i * n_vars_tr..(i + 1) * n_vars_tr];
        y.push(row[0]);
        if include_intercept {
            x_in.push(1.0);
        }
        x_in.extend_from_slice(&row[1..]);
    }

    let out = run_spur_pipeline(&SpurPipelineRequest {
        coords: &coords,
        n,
        d,
        metric,
        y: &y,
        x_in: &x_in,
        k_x,
        emat: &emat,
        q,
        nrep,
        transform_vars: &y_and_x,
        n_transform_vars: n_vars_tr,
        nthreads,
    })
    .map_err(format_core_error)?;

    let info = [
        diagnostic_info(out.i0.tied, out.i0.all_directions, out.i0.best_tuning_power),
        diagnostic_info(out.i1.tied, out.i1.all_directions, out.i1.best_tuning_power),
        diagnostic_info(
            out.i0resid.tied,
            out.i0resid.all_directions,
            out.i0resid.best_tuning_power,
        ),
        diagnostic_info(
            out.i1resid.tied,
            out.i1resid.all_directions,
            out.i1resid.best_tuning_power,
        ),
    ];
    write_matrix_with_infinity(info_mat_name, 4, 3, info.as_flattened())?;

    let diag = [
        out.i0.lr[0],
        out.i0.pvalue[0],
        out.i0.ha_param,
        out.i1.lr[0],
        out.i1.pvalue[0],
        out.i1.ha_param,
        out.i0resid.lr[0],
        out.i0resid.pvalue[0],
        out.i0resid.ha_param,
        out.i1resid.lr[0],
        out.i1resid.pvalue[0],
        out.i1resid.ha_param,
    ];
    let mut cvs = Vec::with_capacity(12);
    cvs.extend_from_slice(&out.i0.cvalue);
    cvs.extend_from_slice(&out.i1.cvalue);
    cvs.extend_from_slice(&out.i0resid.cvalue);
    cvs.extend_from_slice(&out.i1resid.cvalue);

    write_matrix(diag_mat_name, 4, 3, &diag)?;
    write_matrix_with_infinity(cvs_mat_name, 4, 3, &cvs)?;
    write_var_matrix(
        &obs,
        1 + d + n_vars_tr,
        n_vars_tr,
        &out.transformed_vars,
        "transformed spur variables",
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dispatch;

    #[test]
    fn ping_and_spur_arg_validation_work() {
        assert!(dispatch(&["ping_unified_v1"]).is_ok());
        for retired in ["ping_spur_v1", "ping_spur_v2", "ping_spur_v3"] {
            assert!(dispatch(&[retired]).is_err());
        }

        // Wrong arity errors cleanly
        assert!(dispatch(&["spur_i0", "2"]).is_err());
        assert!(dispatch(&["spur_i1", "2"]).is_err());
        assert!(dispatch(&["spur_i0_resid", "2"]).is_err());
        assert!(dispatch(&["spur_i1_resid", "2"]).is_err());
        assert!(dispatch(&["spur_halflife", "2"]).is_err());
        assert!(dispatch(&["spur_transform", "2"]).is_err());
        assert!(dispatch(&["spur_pipeline", "2"]).is_err());

        // Non-residual test rejects k_x_raw > 0 or include_intercept = 1
        assert!(dispatch(&[
            "spur_i0", "2", "1", "0", "1", "15", "100", "-1", "e.bin", "cv", "st", "p", "ha", "t"
        ])
        .is_err());
        assert!(dispatch(&[
            "spur_i1", "2", "0", "1", "1", "15", "100", "-1", "e.bin", "cv", "st", "p", "ha", "t"
        ])
        .is_err());

        // Residual test with empty regressor + no intercept errors cleanly
        assert!(dispatch(&[
            "spur_i0_resid",
            "2",
            "0",
            "0",
            "1",
            "15",
            "100",
            "-1",
            "e.bin",
            "cv",
            "st",
            "p",
            "ha",
            "t"
        ])
        .is_err());

        // Invalid halflife level
        assert!(dispatch(&[
            "spur_halflife",
            "2",
            "1",
            "15",
            "100",
            "100",
            "-1",
            "e.bin",
            "cil",
            "ciu",
            "md",
            "info",
            "comp",
            "cnt"
        ])
        .unwrap_err()
        .contains("level must lie strictly in (0, 100)"));

        // Invalid transform combinations
        assert!(dispatch(&["spur_transform", "2", "iso", "1", "0", "0", "1", "-1"]).is_err());
        assert!(dispatch(&["spur_transform", "2", "cluster", "0", "-1", "0", "1", "-1"]).is_err());
        assert!(dispatch(&["spur_transform", "2", "bogus", "0", "-1", "0", "1", "-1"]).is_err());

        // Invalid pipeline k_x_raw=0
        assert!(dispatch(&[
            "spur_pipeline",
            "2",
            "0",
            "1",
            "1",
            "15",
            "100",
            "-1",
            "e.bin",
            "d",
            "c",
            "t"
        ])
        .is_err());
    }

    #[test]
    fn build_x_in_and_cluster_normalization_match_stata_conventions() {
        let x_raw = [2.0, 3.0, 4.0, 5.0];
        let (x_with_int, k1) = build_x_in(&x_raw, 2, 2, true).unwrap();
        assert_eq!(k1, 3);
        assert_eq!(x_with_int, vec![1.0, 2.0, 3.0, 1.0, 4.0, 5.0]);

        let (x_only_int, k0) = build_x_in(&[], 3, 0, true).unwrap();
        assert_eq!(k0, 1);
        assert_eq!(x_only_int, vec![1.0, 1.0, 1.0]);

        let (x_no_int, k2) = build_x_in(&x_raw, 2, 2, false).unwrap();
        assert_eq!(k2, 2);
        assert_eq!(x_no_int, vec![2.0, 3.0, 4.0, 5.0]);

        // Gapped / negative cluster labels map densely while preserving equality groups.
        let cl = normalize_cluster_ids(&[10.0, -5.0, 10.0, 0.0, -0.0, -5.0]).unwrap();
        assert_eq!(cl, vec![0, 1, 0, 2, 2, 1]);
        assert!(normalize_cluster_ids(&[1.0, f64::NAN]).is_err());
    }
}
