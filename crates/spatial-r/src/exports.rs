use crate as spatial_r;
use extendr_api::prelude::*;
use spatial_r::{
    run_spur_halflife_r, run_spur_i0_r, run_spur_i0_resid_r, run_spur_i1_r, run_spur_i1_resid_r,
    run_spur_pipeline_r, run_spur_transform_r,
};

#[extendr]
fn rust_run_spur_i0(
    coords: RMatrix<f64>,
    latlong: bool,
    y: RMatrix<f64>,
    emat: RMatrix<f64>,
    nthreads: Nullable<i32>,
) -> extendr_api::Result<List> {
    run_spur_i0_r(coords, latlong, y, emat, nthreads)
}

#[extendr]
fn rust_run_spur_i1(
    coords: RMatrix<f64>,
    latlong: bool,
    y: RMatrix<f64>,
    emat: RMatrix<f64>,
    nthreads: Nullable<i32>,
) -> extendr_api::Result<List> {
    run_spur_i1_r(coords, latlong, y, emat, nthreads)
}

#[extendr]
fn rust_run_spur_i0_resid(
    coords: RMatrix<f64>,
    latlong: bool,
    y: RMatrix<f64>,
    x_in: RMatrix<f64>,
    emat: RMatrix<f64>,
    nthreads: Nullable<i32>,
) -> extendr_api::Result<List> {
    run_spur_i0_resid_r(coords, latlong, y, x_in, emat, nthreads)
}

#[extendr]
fn rust_run_spur_i1_resid(
    coords: RMatrix<f64>,
    latlong: bool,
    y: RMatrix<f64>,
    x_in: RMatrix<f64>,
    emat: RMatrix<f64>,
    nthreads: Nullable<i32>,
) -> extendr_api::Result<List> {
    run_spur_i1_resid_r(coords, latlong, y, x_in, emat, nthreads)
}

#[extendr]
fn rust_run_spur_halflife(
    coords: RMatrix<f64>,
    latlong: bool,
    y: Doubles,
    emat: RMatrix<f64>,
    level: f64,
    nthreads: Nullable<i32>,
) -> extendr_api::Result<List> {
    run_spur_halflife_r(coords, latlong, y, emat, level, nthreads)
}

#[extendr]
#[allow(clippy::too_many_arguments)]
fn rust_run_spur_transform(
    kind: &str,
    n: i32,
    y: RMatrix<f64>,
    coords: Nullable<RMatrix<f64>>,
    latlong: bool,
    radius: Nullable<f64>,
    cluster_id: Nullable<Integers>,
    nthreads: Nullable<i32>,
) -> extendr_api::Result<List> {
    run_spur_transform_r(kind, n, y, coords, latlong, radius, cluster_id, nthreads)
}

#[extendr]
fn rust_run_spur_pipeline(
    coords: RMatrix<f64>,
    latlong: bool,
    y: Doubles,
    x_in: RMatrix<f64>,
    emat: RMatrix<f64>,
    transform_vars: RMatrix<f64>,
    nthreads: Nullable<i32>,
) -> extendr_api::Result<List> {
    run_spur_pipeline_r(coords, latlong, y, x_in, emat, transform_vars, nthreads)
}

#[extendr]
fn rust_build_info() -> List {
    spatial_r::build_info_r()
}

#[extendr]
#[allow(clippy::too_many_arguments)]
fn rust_run_model_scpc(
    coords: RMatrix<f64>,
    latlong: bool,
    avc: f64,
    cvs: bool,
    coef: Doubles,
    scores: RMatrix<f64>,
    bread: RMatrix<f64>,
    design: Nullable<RMatrix<f64>>,
    instruments: Nullable<RMatrix<f64>>,
    cluster_id: Nullable<Integers>,
    fe_id: Nullable<Integers>,
    nthreads: Nullable<i32>,
    inference: &str,
    method: &str,
    seed: Nullable<f64>,
    landmarks: i32,
    samples: i32,
    cv_tolerance: f64,
    residualizer: Nullable<Function>,
) -> std::result::Result<List, spatial_r::ScpcError> {
    spatial_r::rust_run_model_scpc(
        coords,
        latlong,
        avc,
        cvs,
        coef,
        scores,
        bread,
        design,
        instruments,
        cluster_id,
        fe_id,
        nthreads,
        inference,
        method,
        seed,
        landmarks,
        samples,
        cv_tolerance,
        residualizer,
    )
}
