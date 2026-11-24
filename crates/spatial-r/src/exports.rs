use std::panic::{catch_unwind, AssertUnwindSafe};

use spatial_core::{
    run_fitted_scpc, run_spur_halflife, run_spur_i0, run_spur_i0_resid, run_spur_i1,
    run_spur_i1_resid, run_spur_pipeline, run_spur_transform, Metric, ScpcGeometry, ScpcInference,
    ScpcOptions, ScpcOutput, ScpcPreparation, SpurHalfLifeOutput, SpurHalfLifeRequest,
    SpurI0Output, SpurI1Output, SpurPipelineRequest, SpurTestRequest, SpurTransformKind,
    SpurTransformRequest,
};

use crate::r_ffi::{
    as_bool, as_f64, as_f64_matrix, as_f64_slice, as_i32, as_opt_f64, as_opt_f64_matrix,
    as_opt_function, as_opt_i32, as_opt_i32_slice, as_str, call_residualizer, make_error_obj,
    make_named_list, make_real_matrix_col_major_from_row_major, make_real_vec, make_scalar_bool,
    make_scalar_int, make_scalar_real, make_scalar_str, r_nil, ProtectGuard, SEXP,
};
use crate::{col_major_to_row_major, parse_r_index_slice, parse_r_nthreads, ScpcError};

unsafe fn wrap_call(body: impl FnOnce() -> Result<SEXP, ScpcError>) -> SEXP {
    match catch_unwind(AssertUnwindSafe(body)) {
        Ok(Ok(res)) => res,
        Ok(Err(err)) => make_error_obj(&err.message, err.class),
        Err(_) => make_error_obj("native spatial Rust call panicked", None),
    }
}

unsafe fn build_r_scpc_list(out: ScpcOutput, n: usize, k_use: usize) -> Result<SEXP, String> {
    let mut guard = ProtectGuard::new();
    let scpcstats = guard.protect(make_real_matrix_col_major_from_row_major(
        &out.scpcstats,
        k_use,
        6,
    )?);
    let scpccvs = match out.scpccvs {
        Some(ref cvs_rm) => guard.protect(make_real_matrix_col_major_from_row_major(
            cvs_rm, k_use, 4,
        )?),
        None => r_nil(),
    };
    let wfin = guard.protect(make_real_matrix_col_major_from_row_major(
        &out.wfin,
        n,
        out.q + 1,
    )?);
    let q_i32 =
        i32::try_from(out.q).map_err(|_| "selected q exceeds 32-bit R integer".to_string())?;
    let c0 = guard.protect(make_scalar_real(out.c0));
    let cmax = guard.protect(make_scalar_real(out.cmax));
    let cvfin = guard.protect(make_scalar_real(out.cvfin));
    let critical_95 = guard.protect(make_real_vec(&out.critical_95)?);
    let distance_scale = guard.protect(make_scalar_real(out.distance_scale));
    let q_sexp = guard.protect(make_scalar_int(q_i32));
    let tied_sexp = guard.protect(make_scalar_bool(out.tied));

    make_named_list(&[
        ("scpcstats", scpcstats),
        ("scpccvs", scpccvs),
        ("wfin", wfin),
        ("c0", c0),
        ("cmax", cmax),
        ("cvfin", cvfin),
        ("critical_95", critical_95),
        ("distance_scale", distance_scale),
        ("q", q_sexp),
        ("tied", tied_sexp),
    ])
}

unsafe fn build_r_spur_i0_list(out: SpurI0Output, n_y: usize) -> Result<SEXP, String> {
    let mut guard = ProtectGuard::new();
    let n_rho = out.rho_grid.len();
    let lr = guard.protect(make_real_vec(&out.lr)?);
    let pvalue = guard.protect(make_real_vec(&out.pvalue)?);
    let cvalue = guard.protect(make_real_vec(&out.cvalue)?);
    let ha_parm = guard.protect(make_scalar_real(out.ha_param));
    let best_tuning_power = guard.protect(make_scalar_real(
        out.best_tuning_power.unwrap_or(f64::NAN),
    ));
    let cvalue_mat = guard.protect(make_real_matrix_col_major_from_row_major(
        &out.cvalue_mat,
        n_rho,
        3,
    )?);
    let pvalue_mat = guard.protect(make_real_matrix_col_major_from_row_major(
        &out.pvalue_mat,
        n_rho,
        n_y,
    )?);
    let rho_grid = guard.protect(make_real_vec(&out.rho_grid)?);
    let tied = guard.protect(make_scalar_bool(out.tied));
    let all_directions = guard.protect(make_scalar_bool(out.all_directions));

    make_named_list(&[
        ("LR", lr),
        ("pvalue", pvalue),
        ("cvalue", cvalue),
        ("ha_parm", ha_parm),
        ("best_tuning_power", best_tuning_power),
        ("cvalue_mat", cvalue_mat),
        ("pvalue_mat", pvalue_mat),
        ("rho_grid", rho_grid),
        ("tied", tied),
        ("all_directions", all_directions),
    ])
}

unsafe fn build_r_spur_i1_list(out: SpurI1Output) -> Result<SEXP, String> {
    let mut guard = ProtectGuard::new();
    let lr = guard.protect(make_real_vec(&out.lr)?);
    let pvalue = guard.protect(make_real_vec(&out.pvalue)?);
    let cv_vec = guard.protect(make_real_vec(&out.cvalue)?);
    let ha_parm = guard.protect(make_scalar_real(out.ha_param));
    let best_tuning_power = guard.protect(make_scalar_real(
        out.best_tuning_power.unwrap_or(f64::NAN),
    ));
    let tied = guard.protect(make_scalar_bool(out.tied));
    let all_directions = guard.protect(make_scalar_bool(out.all_directions));

    make_named_list(&[
        ("LR", lr),
        ("pvalue", pvalue),
        ("cv_vec", cv_vec),
        ("ha_parm", ha_parm),
        ("best_tuning_power", best_tuning_power),
        ("tied", tied),
        ("all_directions", all_directions),
    ])
}

unsafe fn build_r_spur_hl_list(out: SpurHalfLifeOutput) -> Result<SEXP, String> {
    let mut guard = ProtectGuard::new();
    let n_components = out.components.len();
    let bounds = out.components.into_iter().flatten().collect::<Vec<_>>();
    let hl_ci_lower = guard.protect(make_scalar_real(out.ci_lower));
    let hl_ci_upper = guard.protect(make_scalar_real(out.ci_upper));
    let max_dist = guard.protect(make_scalar_real(out.max_dist));
    let pv_mat = guard.protect(make_real_vec(&out.pv_mat)?);
    let hl_grid_ho = guard.protect(make_real_vec(&out.hl_grid_ho)?);
    let components = guard.protect(make_real_matrix_col_major_from_row_major(
        &bounds,
        n_components,
        2,
    )?);
    let tied = guard.protect(make_scalar_bool(out.tied));
    let all_directions = guard.protect(make_scalar_bool(out.all_directions));

    make_named_list(&[
        ("hl_ci_lower", hl_ci_lower),
        ("hl_ci_upper", hl_ci_upper),
        ("max_dist", max_dist),
        ("pv_mat", pv_mat),
        ("hl_grid_ho", hl_grid_ho),
        ("components", components),
        ("tied", tied),
        ("all_directions", all_directions),
    ])
}

unsafe fn run_spur_test_r_impl(
    coords: SEXP,
    latlong: SEXP,
    y: SEXP,
    x_in: Option<SEXP>,
    emat: SEXP,
    nthreads: SEXP,
    is_i0: bool,
) -> Result<SEXP, ScpcError> {
    let coords = as_f64_matrix(coords, "coords")?;
    let latlong = as_bool(latlong, "latlong")?;
    let y = as_f64_matrix(y, "y")?;
    let emat = as_f64_matrix(emat, "emat")?;
    let nthreads_usize = parse_r_nthreads(as_opt_i32(nthreads, "nthreads")?)?;

    let n = coords.nrows;
    let d = coords.ncols;
    if y.nrows != n || y.ncols == 0 {
        return Err(format!(
            "Y dimensions ({} x {}) incompatible with coords rows n = {n}",
            y.nrows, y.ncols
        )
        .into());
    }
    let n_y = y.ncols;
    let draw_rows = emat.nrows;
    if draw_rows % 2 != 0 {
        return Err("SPUR diagnostic emat must have 2*q rows".into());
    }
    let q = draw_rows / 2;
    let nrep = emat.ncols;
    let metric = if latlong {
        Metric::GreatCircleFractionLonLat
    } else {
        Metric::Euclidean
    };

    let coords_rm = col_major_to_row_major(coords.data, n, d)?;
    let y_rm = col_major_to_row_major(y.data, n, n_y)?;
    let emat_rm = col_major_to_row_major(emat.data, draw_rows, nrep)?;
    let (x_rm, k_x) = match x_in {
        None => (None, 0usize),
        Some(xs) => {
            let xm = as_f64_matrix(xs, "x_in")?;
            if xm.nrows != n || xm.ncols == 0 {
                return Err(format!(
                    "X_in dimensions ({} x {}) incompatible with n = {n}",
                    xm.nrows, xm.ncols
                )
                .into());
            }
            let kx = xm.ncols;
            (Some(col_major_to_row_major(xm.data, n, kx)?), kx)
        }
    };

    let req = SpurTestRequest {
        coords: &coords_rm,
        n,
        d,
        metric,
        y: &y_rm,
        n_y,
        x_in: x_rm.as_deref(),
        k_x,
        emat: &emat_rm,
        q,
        nrep,
        nthreads: nthreads_usize,
    };

    if is_i0 {
        let out = if req.x_in.is_some() {
            run_spur_i0_resid(&req)
        } else {
            run_spur_i0(&req)
        }?;
        Ok(build_r_spur_i0_list(out, n_y)?)
    } else {
        let out = if req.x_in.is_some() {
            run_spur_i1_resid(&req)
        } else {
            run_spur_i1(&req)
        }?;
        Ok(build_r_spur_i1_list(out)?)
    }
}

#[no_mangle]
pub unsafe extern "C" fn wrap__rust_build_info() -> SEXP {
    wrap_call(|| {
        let mut guard = ProtectGuard::new();
        let version = guard.protect(make_scalar_str(spatial_core::VERSION)?);
        let package_version = guard.protect(make_scalar_str(env!("CARGO_PKG_VERSION"))?);
        let spec_version = guard.protect(make_scalar_str(spatial_core::SPEC_VERSION)?);
        let earth_circumference_m =
            guard.protect(make_scalar_real(spatial_core::EARTH_CIRCUMFERENCE_M));
        let os = guard.protect(make_scalar_str(std::env::consts::OS)?);
        let arch = guard.protect(make_scalar_str(std::env::consts::ARCH)?);
        let linalg_provider = guard.protect(make_scalar_str(spatial_core::linalg_provider())?);
        let accel = guard.protect(make_scalar_bool(
            spatial_core::threads::inherited_accelerate_state(),
        ));
        let scpc_api_version =
            guard.protect(make_scalar_int(spatial_core::SCPC_API_VERSION as i32));
        let spur_api_version =
            guard.protect(make_scalar_int(spatial_core::SPUR_API_VERSION as i32));

        Ok(make_named_list(&[
            ("version", version),
            ("package_version", package_version),
            ("spec_version", spec_version),
            ("earth_circumference_m", earth_circumference_m),
            ("os", os),
            ("arch", arch),
            ("linalg_provider", linalg_provider),
            ("inherited_accelerate_state", accel),
            ("scpc_api_version", scpc_api_version),
            ("spur_api_version", spur_api_version),
        ])?)
    })
}

#[no_mangle]
pub unsafe extern "C" fn wrap__rust_run_spur_i0(
    coords: SEXP,
    latlong: SEXP,
    y: SEXP,
    emat: SEXP,
    nthreads: SEXP,
) -> SEXP {
    wrap_call(|| run_spur_test_r_impl(coords, latlong, y, None, emat, nthreads, true))
}

#[no_mangle]
pub unsafe extern "C" fn wrap__rust_run_spur_i1(
    coords: SEXP,
    latlong: SEXP,
    y: SEXP,
    emat: SEXP,
    nthreads: SEXP,
) -> SEXP {
    wrap_call(|| run_spur_test_r_impl(coords, latlong, y, None, emat, nthreads, false))
}

#[no_mangle]
pub unsafe extern "C" fn wrap__rust_run_spur_i0_resid(
    coords: SEXP,
    latlong: SEXP,
    y: SEXP,
    x_in: SEXP,
    emat: SEXP,
    nthreads: SEXP,
) -> SEXP {
    wrap_call(|| run_spur_test_r_impl(coords, latlong, y, Some(x_in), emat, nthreads, true))
}

#[no_mangle]
pub unsafe extern "C" fn wrap__rust_run_spur_i1_resid(
    coords: SEXP,
    latlong: SEXP,
    y: SEXP,
    x_in: SEXP,
    emat: SEXP,
    nthreads: SEXP,
) -> SEXP {
    wrap_call(|| run_spur_test_r_impl(coords, latlong, y, Some(x_in), emat, nthreads, false))
}

#[no_mangle]
pub unsafe extern "C" fn wrap__rust_run_spur_halflife(
    coords: SEXP,
    latlong: SEXP,
    y: SEXP,
    emat: SEXP,
    level: SEXP,
    nthreads: SEXP,
) -> SEXP {
    wrap_call(|| {
        let coords = as_f64_matrix(coords, "coords")?;
        let latlong = as_bool(latlong, "latlong")?;
        let y_slice = as_f64_slice(y, "y")?;
        let emat = as_f64_matrix(emat, "emat")?;
        let level = as_f64(level, "level")?;
        let nthreads_usize = parse_r_nthreads(as_opt_i32(nthreads, "nthreads")?)?;

        let n = coords.nrows;
        let d = coords.ncols;
        if y_slice.len() != n {
            return Err(format!(
                "y length {} does not match coords rows n = {n}",
                y_slice.len()
            )
            .into());
        }
        let q = emat.nrows;
        let nrep = emat.ncols;
        let metric = if latlong {
            Metric::GreatCircleFractionLonLat
        } else {
            Metric::Euclidean
        };
        let coords_rm = col_major_to_row_major(coords.data, n, d)?;
        let emat_rm = col_major_to_row_major(emat.data, q, nrep)?;
        let out = run_spur_halflife(&SpurHalfLifeRequest {
            coords: &coords_rm,
            n,
            d,
            metric,
            y: y_slice,
            emat: &emat_rm,
            q,
            nrep,
            level,
            nthreads: nthreads_usize,
        })?;
        Ok(build_r_spur_hl_list(out)?)
    })
}

#[no_mangle]
pub unsafe extern "C" fn wrap__rust_run_spur_transform(
    kind: SEXP,
    n: SEXP,
    y: SEXP,
    coords: SEXP,
    latlong: SEXP,
    radius: SEXP,
    cluster_id: SEXP,
    nthreads: SEXP,
) -> SEXP {
    wrap_call(|| {
        let kind = as_str(kind, "kind")?;
        let n = as_i32(n, "n")?;
        let n_usize = usize::try_from(n)
            .map_err(|_| "transform requires a nonnegative sample size".to_string())?;
        let y = as_f64_matrix(y, "y")?;
        let coords_opt = as_opt_f64_matrix(coords, "coords")?;
        let latlong = as_bool(latlong, "latlong")?;
        let radius_opt = as_opt_f64(radius, "radius")?;
        let cluster_slice = as_opt_i32_slice(cluster_id, "cluster_id")?;
        let nthreads_usize = parse_r_nthreads(as_opt_i32(nthreads, "nthreads")?)?;

        let tr_kind = match kind {
            "nn" => SpurTransformKind::Nn,
            "iso" => {
                let r = radius_opt.ok_or_else(|| "iso transform requires radius".to_string())?;
                SpurTransformKind::Iso { radius: r }
            }
            "cluster" => SpurTransformKind::Cluster,
            "lbmgls" => SpurTransformKind::LbmGls,
            other => {
                return Err(format!(
                    "unknown transform kind '{other}' (expected nn, iso, cluster, or lbmgls)"
                )
                .into());
            }
        };

        let (coords_rm, d, metric) = match coords_opt {
            None => (None, 0usize, Metric::Euclidean),
            Some(cm) => {
                if cm.nrows != n_usize {
                    return Err(format!(
                        "coords rows {} does not match expected n = {n_usize}",
                        cm.nrows
                    )
                    .into());
                }
                let dc = cm.ncols;
                let met = if latlong {
                    Metric::GreatCircleFractionLonLat
                } else {
                    Metric::Euclidean
                };
                (Some(col_major_to_row_major(cm.data, n_usize, dc)?), dc, met)
            }
        };

        let cluster_vec = parse_r_index_slice(cluster_slice, n_usize, "cluster_id", "observations")?;
        let n_y = y.ncols;
        if y.nrows != n_usize || n_y == 0 {
            return Err(format!(
                "Y dimensions ({} x {n_y}) incompatible with n = {n_usize}",
                y.nrows
            )
            .into());
        }
        let y_rm = col_major_to_row_major(y.data, n_usize, n_y)?;

        let out = run_spur_transform(&SpurTransformRequest {
            coords: coords_rm.as_deref(),
            n: n_usize,
            d,
            metric,
            kind: tr_kind,
            cluster_id: cluster_vec.as_deref(),
            y: &y_rm,
            n_y,
            nthreads: nthreads_usize,
        })?;

        let mut guard = ProtectGuard::new();
        let transformed =
            guard.protect(make_real_matrix_col_major_from_row_major(&out, n_usize, n_y)?);
        Ok(make_named_list(&[("transformed", transformed)])?)
    })
}

#[no_mangle]
pub unsafe extern "C" fn wrap__rust_run_spur_pipeline(
    coords: SEXP,
    latlong: SEXP,
    y: SEXP,
    x_in: SEXP,
    emat: SEXP,
    transform_vars: SEXP,
    nthreads: SEXP,
) -> SEXP {
    wrap_call(|| {
        let coords = as_f64_matrix(coords, "coords")?;
        let latlong = as_bool(latlong, "latlong")?;
        let y_slice = as_f64_slice(y, "y")?;
        let x_in = as_f64_matrix(x_in, "x_in")?;
        let emat = as_f64_matrix(emat, "emat")?;
        let transform_vars = as_f64_matrix(transform_vars, "transform_vars")?;
        let nthreads_usize = parse_r_nthreads(as_opt_i32(nthreads, "nthreads")?)?;

        let n = coords.nrows;
        let d = coords.ncols;
        if y_slice.len() != n {
            return Err(format!(
                "y length {} does not match coords rows n = {n}",
                y_slice.len()
            )
            .into());
        }
        if x_in.nrows != n || x_in.ncols == 0 {
            return Err(format!(
                "X_in dimensions ({} x {}) incompatible with n = {n}",
                x_in.nrows, x_in.ncols
            )
            .into());
        }
        if transform_vars.nrows != n || transform_vars.ncols == 0 {
            return Err(format!(
                "transform_vars dimensions ({} x {}) incompatible with n = {n}",
                transform_vars.nrows, transform_vars.ncols
            )
            .into());
        }
        let k_x = x_in.ncols;
        let n_tr = transform_vars.ncols;
        let draw_rows = emat.nrows;
        if draw_rows % 2 != 0 {
            return Err("SPUR diagnostic emat must have 2*q rows".into());
        }
        let q = draw_rows / 2;
        let nrep = emat.ncols;
        let metric = if latlong {
            Metric::GreatCircleFractionLonLat
        } else {
            Metric::Euclidean
        };

        let coords_rm = col_major_to_row_major(coords.data, n, d)?;
        let x_rm = col_major_to_row_major(x_in.data, n, k_x)?;
        let emat_rm = col_major_to_row_major(emat.data, draw_rows, nrep)?;
        let tr_rm = col_major_to_row_major(transform_vars.data, n, n_tr)?;

        let out = run_spur_pipeline(&SpurPipelineRequest {
            coords: &coords_rm,
            n,
            d,
            metric,
            y: y_slice,
            x_in: &x_rm,
            k_x,
            emat: &emat_rm,
            q,
            nrep,
            transform_vars: &tr_rm,
            n_transform_vars: n_tr,
            nthreads: nthreads_usize,
        })?;

        let mut guard = ProtectGuard::new();
        let i0 = guard.protect(build_r_spur_i0_list(out.i0, 1)?);
        let i1 = guard.protect(build_r_spur_i1_list(out.i1)?);
        let i0resid = guard.protect(build_r_spur_i0_list(out.i0resid, 1)?);
        let i1resid = guard.protect(build_r_spur_i1_list(out.i1resid)?);
        let tr_mat = guard.protect(make_real_matrix_col_major_from_row_major(
            &out.transformed_vars,
            n,
            n_tr,
        )?);
        Ok(make_named_list(&[
            ("i0", i0),
            ("i1", i1),
            ("i0resid", i0resid),
            ("i1resid", i1resid),
            ("transformed_vars", tr_mat),
        ])?)
    })
}

#[no_mangle]
#[allow(clippy::too_many_arguments)]
pub unsafe extern "C" fn wrap__rust_run_model_scpc(
    coords: SEXP,
    latlong: SEXP,
    avc: SEXP,
    cvs: SEXP,
    coef: SEXP,
    scores: SEXP,
    bread: SEXP,
    design: SEXP,
    instruments: SEXP,
    cluster_id: SEXP,
    fe_id: SEXP,
    nthreads: SEXP,
    inference: SEXP,
    method: SEXP,
    seed: SEXP,
    landmarks: SEXP,
    samples: SEXP,
    cv_tolerance: SEXP,
    residualizer: SEXP,
) -> SEXP {
    wrap_call(|| {
        let coords = as_f64_matrix(coords, "coords")?;
        let latlong = as_bool(latlong, "latlong")?;
        let avc = as_f64(avc, "avc")?;
        let cvs = as_bool(cvs, "cvs")?;
        let coefficients = as_f64_slice(coef, "coef")?;
        let scores = as_f64_matrix(scores, "scores")?;
        let bread = as_f64_matrix(bread, "bread")?;
        let design = as_opt_f64_matrix(design, "design")?;
        let instruments = as_opt_f64_matrix(instruments, "instruments")?;
        let cluster_slice = as_opt_i32_slice(cluster_id, "cluster_id")?;
        let fe_slice = as_opt_i32_slice(fe_id, "fe_id")?;
        let nthreads = as_opt_i32(nthreads, "nthreads")?;
        let inference = as_str(inference, "inference")?;
        let method = as_str(method, "method")?;
        let seed = as_opt_f64(seed, "seed")?;
        let landmarks = as_i32(landmarks, "landmarks")?;
        let samples = as_i32(samples, "samples")?;
        let cv_tolerance = as_f64(cv_tolerance, "cv_tolerance")?;
        let residualizer_fn = as_opt_function(residualizer, "residualizer")?;

        let n = coords.nrows;
        let d = coords.ncols;
        let p = scores.ncols;
        let k = coefficients.len();
        if bread.nrows != k || bread.ncols != p {
            return Err("score or bread dimensions do not match the model".into());
        }
        let cluster = parse_r_index_slice(cluster_slice, n, "cluster_id", "clusters")?;
        let n_obs = cluster.as_ref().map_or(n, Vec::len);
        if scores.nrows != n_obs {
            return Err("scores must have one row per observation".into());
        }
        let fe = parse_r_index_slice(fe_slice, n_obs, "fe_id", "observations")?;
        let read_design =
            |matrix: Option<crate::r_ffi::RealMatrixView<'_>>| -> Result<(Option<Vec<f64>>, usize), String> {
                match matrix {
                    None => Ok((None, 0)),
                    Some(m) => {
                        if m.nrows != n_obs {
                            return Err(
                                "a model matrix has the wrong number of observations".to_string()
                            );
                        }
                        Ok((Some(col_major_to_row_major(m.data, n_obs, m.ncols)?), m.ncols))
                    }
                }
            };
        let (x, design_cols) = read_design(design)?;
        let (z, instrument_cols) = read_design(instruments)?;
        if !cv_tolerance.is_finite() || cv_tolerance <= 0.0 {
            return Err("cv_tolerance must be finite and strictly positive".into());
        }
        let seed = match seed {
            None => None,
            Some(s) if (0.0..4_294_967_296.0).contains(&s) && s.fract() == 0.0 => Some(s as u32),
            _ => return Err("seed must be an integer in [0, 2^32)".into()),
        };
        let inference = match inference {
            "conditional_ols" => ScpcInference::ConditionalOls,
            "conditional_2sls" => ScpcInference::ConditionalTwoStageLeastSquares,
            "unconditional" => ScpcInference::Unconditional,
            _ => return Err("unknown SCPC inference procedure".into()),
        };
        let geometry = match (method, seed) {
            ("exact", None) => ScpcGeometry::Exact,
            ("nystrom", Some(seed)) if landmarks >= 2 && samples >= 1 => ScpcGeometry::Nystrom {
                seed,
                landmarks: landmarks as usize,
                samples: samples as usize,
            },
            _ => return Err("invalid SCPC method or approximation settings".into()),
        };
        let coord_data = col_major_to_row_major(coords.data, n, d)?;
        let obs_scores = col_major_to_row_major(scores.data, n_obs, p)?;
        let bread_rows = col_major_to_row_major(bread.data, k, p)?;
        let model = ScpcPreparation {
            n_sites: n,
            n_obs,
            n_scores: p,
            coefficients,
            bread: &bread_rows,
            scores: &obs_scores,
            design: x.as_deref(),
            design_cols,
            instruments: z.as_deref(),
            instrument_cols,
            direction_design: None,
            cluster_ids: cluster.as_deref(),
            fixed_effects: fe.as_deref(),
        };
        let project = |values: &mut [f64], width| {
            let func = residualizer_fn.expect("validated callback");
            let projected = call_residualizer(func, values, n_obs, width).map_err(|error| {
                spatial_core::CoreError::invalid_input("fixed-effect projection failed")
                    .with_message(format!("Fixed-effect projection failed: {error}"))
            })?;
            values.copy_from_slice(&projected);
            Ok(())
        };
        let out = run_fitted_scpc(
            model,
            &coord_data,
            d,
            if latlong {
                Metric::GreatCircleFractionLonLat
            } else {
                Metric::Euclidean
            },
            ScpcOptions {
                avc,
                cvs,
                nthreads: parse_r_nthreads(nthreads)?,
                geometry,
                inference,
                cv_tolerance,
            },
            true,
            residualizer_fn.map(|_| &project as &spatial_core::Residualize<'_>),
        )?;
        Ok(build_r_scpc_list(out, n, k)?)
    })
}
