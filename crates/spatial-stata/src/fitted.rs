//! Stata postestimation adapter. Stata supplies the active estimation sample,
//! materialized model columns and structural residuals.

use std::fs::File;
use std::io::{BufWriter, Write};

use spatial_core::distance::canonical_coordinates;
use spatial_core::preparation::regression_scores;
use spatial_core::{
    run_fitted_scpc, Metric, ScpcGeometry, ScpcInference, ScpcOptions, ScpcPreparation,
};

use crate::stata_spi::{active_obs, plugin_var_count, read_matrix, read_var_matrix, write_matrix};
use crate::{
    cluster_coordinates, format_core_error, parse_bool_flag, parse_nthreads, parse_seed,
    parse_usize, StataDesign,
};

pub fn run_scpc(args: &[&str]) -> Result<(), String> {
    let [d, regressors, report, has_cons, endog, excluded, cluster, geographic, avc, uncond, cvs, nthreads, method, seed, landmarks, samples, tolerance, estimates, table, cutoffs, effective, meta, weights_file] =
        *args
    else {
        return Err("invalid fitted SCPC plugin argument count".into());
    };
    let d = parse_usize(d, "coordinate width")?;
    let p = parse_usize(regressors, "regressor count")?;
    let k = parse_usize(report, "coefficient count")?;
    let has_cons = parse_bool_flag(has_cons, "constant")?;
    let endog = parse_usize(endog, "endogenous columns")?;
    let excluded = parse_usize(excluded, "excluded instruments")?;
    let cluster = parse_bool_flag(cluster, "cluster")?;
    let geographic = parse_bool_flag(geographic, "latlong")?;
    let avc = avc.parse().map_err(|_| "invalid avc")?;
    let uncond = parse_bool_flag(uncond, "unconditional")?;
    let cvs = parse_bool_flag(cvs, "cvs")?;
    let nthreads = parse_nthreads(nthreads)?;
    let geometry = match method {
        "exact" => ScpcGeometry::Exact,
        "nystrom" => ScpcGeometry::Nystrom {
            seed: parse_seed(seed)?,
            landmarks: parse_usize(landmarks, "landmarks")?,
            samples: parse_usize(samples, "landmark samples")?,
        },
        _ => return Err("method must be exact or nystrom".into()),
    };
    let cv_tolerance = tolerance
        .parse()
        .map_err(|_| "invalid critical-value tolerance")?;
    let inference = if uncond {
        ScpcInference::Unconditional
    } else if endog > 0 {
        ScpcInference::ConditionalTwoStageLeastSquares
    } else {
        ScpcInference::ConditionalOls
    };
    let metric = if geographic {
        Metric::GreatCircleFractionLatLon
    } else {
        Metric::Euclidean
    };

    let rows = active_obs();
    if d == 0 || k == 0 || k > p || p > rows.len() {
        return Err("invalid fitted dimensions".into());
    }
    let n_obs = rows.len();
    let nx = p - usize::from(has_cons);
    let cluster_cols = usize::from(cluster);
    if plugin_var_count() != d + 1 + cluster_cols + nx + excluded {
        return Err("fitted input variables and regression columns disagree".into());
    }
    let coords = read_var_matrix(&rows, 1, d, "spatial coordinates")?;
    let residuals = read_var_matrix(&rows, d + 1, 1, "fitted residuals")?;
    let (n_sites, membership, coords) = if cluster {
        let ids = read_var_matrix(&rows, d + 2, 1, "spatial cluster index")?;
        let coords = if geographic {
            canonical_coordinates(&coords, n_obs, d, metric).map_err(format_core_error)?
        } else {
            coords
        };
        let (n, ids, coords) = cluster_coordinates(&ids, &coords, d)?;
        (n, Some(ids), coords)
    } else {
        (n_obs, None, coords)
    };
    let x = read_var_matrix(&rows, d + 2 + cluster_cols, nx, "structural regressors")?;
    let z_excl = read_var_matrix(
        &rows,
        d + 2 + cluster_cols + nx,
        excluded,
        "excluded instruments",
    )?;
    let design = StataDesign::from_columns(n_obs, p, has_cons, endog, excluded, &x, &z_excl)?;
    let (scores, bread) = regression_scores(
        &design.x,
        &residuals,
        p,
        design.z.as_deref().map(|z| (z, design.z_cols)),
    )
    .map_err(format_core_error)?;
    let estimates = read_matrix(estimates, 1, p)?;
    let output = run_fitted_scpc(
        ScpcPreparation {
            n_sites,
            n_obs,
            n_scores: p,
            coefficients: &estimates[..k],
            bread: &bread[..k * p],
            scores: &scores,
            design: Some(&design.x),
            design_cols: p,
            instruments: design.z.as_deref(),
            instrument_cols: design.z_cols,
            direction_design: None,
            cluster_ids: membership.as_deref(),
            fixed_effects: None,
        },
        &coords,
        d,
        metric,
        ScpcOptions {
            inference,
            geometry,
            avc,
            cvs,
            cv_tolerance,
            nthreads,
        },
        true,
        None,
    )
    .map_err(format_core_error)?;
    write_matrix(table, k, 6, &output.scpcstats)?;
    if let Some(cvs) = &output.scpccvs {
        write_matrix(cutoffs, k, 4, cvs)?;
    }
    write_matrix(effective, k, 1, &output.critical_95)?;
    write_matrix(
        meta,
        1,
        8,
        &[
            output.c0,
            output.cmax,
            output.cvfin,
            output.q as f64,
            if output.tied { 1.0 } else { 0.0 },
            n_obs as f64,
            n_sites as f64,
            output.distance_scale,
        ],
    )?;
    if weights_file != "_" {
        write_f64(weights_file, &output.wfin)
            .map_err(|error| format!("cannot write spatial weights: {error}"))?;
    }
    Ok(())
}

fn write_f64(path: &str, values: &[f64]) -> std::io::Result<()> {
    let mut writer = BufWriter::new(File::create(path)?);
    for value in values {
        writer.write_all(&value.to_le_bytes())?;
    }
    writer.flush()
}
