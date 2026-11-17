use numpy::ndarray::{Array1, Array2};
use numpy::{
    IntoPyArray, PyArray1, PyArray2, PyArrayMethods, PyReadonlyArray1, PyReadonlyArray2,
    PyUntypedArrayMethods,
};
use pyo3::exceptions::{PyNotImplementedError, PyRuntimeError, PyTypeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyBool, PyDict, PyTuple};
use spatial_core::{
    run_approx_scpc as core_run_approx_scpc, run_exact_scpc as core_run_exact_scpc,
    run_spur_halflife as core_run_spur_halflife, run_spur_i0 as core_run_spur_i0,
    run_spur_i0_resid as core_run_spur_i0_resid, run_spur_i1 as core_run_spur_i1,
    run_spur_i1_resid as core_run_spur_i1_resid, run_spur_pipeline as core_run_spur_pipeline,
    run_spur_transform as core_run_spur_transform, ApproxScpcRequest, CoreError, ErrorKind,
    ExactScpcOutput, ExactScpcRequest, Metric, SpurHalfLifeOutput, SpurHalfLifeRequest,
    SpurI0Output, SpurI1Output, SpurPipelineRequest, SpurTestRequest, SpurTransformKind,
    SpurTransformRequest,
};
use std::mem::align_of;

fn map_core_err(err: CoreError) -> PyErr {
    match err.kind() {
        ErrorKind::InvalidInput
        | ErrorKind::CalibrationFailed
        | ErrorKind::DegenerateVariance => PyValueError::new_err(err.message()),
        ErrorKind::UnsupportedCase => PyNotImplementedError::new_err(err.message()),
        ErrorKind::TiedEigenspace => {
            PyRuntimeError::new_err(format!("tied eigenspace: {}", err.message()))
        }
        ErrorKind::InternalFailure => PyRuntimeError::new_err(err.message()),
    }
}

fn require_c_aligned_f64<'py>(
    require_fn: &Bound<'py, PyAny>,
    kwargs: &Bound<'py, PyDict>,
    obj: &Bound<'py, PyAny>,
    name: &str,
) -> PyResult<Bound<'py, PyAny>> {
    if obj.is_none() {
        return Err(PyTypeError::new_err(format!("{name} must not be None")));
    }
    require_fn.call((obj,), Some(kwargs))
}

fn extract_1d_f64<'py>(
    require_fn: &Bound<'py, PyAny>,
    kwargs: &Bound<'py, PyDict>,
    obj: &Bound<'py, PyAny>,
    name: &str,
) -> PyResult<Vec<f64>> {
    let aligned = require_c_aligned_f64(require_fn, kwargs, obj, name)?;
    let arr: PyReadonlyArray1<'py, f64> = aligned.extract()?;
    let data = if arr.is_empty() {
        Vec::new()
    } else {
        if !arr.is_c_contiguous() || (arr.data() as usize) % align_of::<f64>() != 0 {
            return Err(PyValueError::new_err(format!(
                "{name} must be aligned and C-contiguous"
            )));
        }
        arr.as_slice()
            .map_err(|_| PyValueError::new_err(format!("{name} must be C-contiguous")))?
            .to_vec()
    };
    if data.iter().any(|v| !v.is_finite()) {
        return Err(PyValueError::new_err(format!(
            "{name} must contain only finite values"
        )));
    }
    Ok(data)
}

fn extract_2d_row_major_f64<'py>(
    require_fn: &Bound<'py, PyAny>,
    kwargs: &Bound<'py, PyDict>,
    obj: &Bound<'py, PyAny>,
    name: &str,
) -> PyResult<(Vec<f64>, usize, usize)> {
    let aligned = require_c_aligned_f64(require_fn, kwargs, obj, name)?;
    let arr: PyReadonlyArray2<'py, f64> = aligned.extract()?;
    let shape = arr.shape();
    let (rows, cols) = (shape[0], shape[1]);
    let data = if arr.is_empty() {
        Vec::new()
    } else {
        if !arr.is_c_contiguous() || (arr.data() as usize) % align_of::<f64>() != 0 {
            return Err(PyValueError::new_err(format!(
                "{name} must be aligned and C-contiguous"
            )));
        }
        arr.as_slice()
            .map_err(|_| PyValueError::new_err(format!("{name} must be C-contiguous")))?
            .to_vec()
    };
    if data.iter().any(|v| !v.is_finite()) {
        return Err(PyValueError::new_err(format!(
            "{name} must contain only finite values"
        )));
    }
    Ok((data, rows, cols))
}

fn vec_to_pyarray1<'py>(py: Python<'py>, data: Vec<f64>) -> Bound<'py, PyArray1<f64>> {
    Array1::from_vec(data).into_pyarray_bound(py)
}

fn vec_to_pyarray2<'py>(
    py: Python<'py>,
    data: Vec<f64>,
    rows: usize,
    cols: usize,
) -> PyResult<Bound<'py, PyArray2<f64>>> {
    let array = Array2::from_shape_vec((rows, cols), data)
        .map_err(|err| PyValueError::new_err(format!("invalid output shape: {err}")))?;
    Ok(array.into_pyarray_bound(py))
}

fn parse_nthreads(nthreads: Option<Bound<'_, PyAny>>) -> PyResult<Option<usize>> {
    let Some(obj) = nthreads else {
        return Ok(None);
    };
    if obj.is_instance_of::<PyBool>() {
        return Err(PyTypeError::new_err(
            "nthreads must be an integer >= 1 or None",
        ));
    }
    let val: i64 = obj.extract()?;
    if val < 1 {
        return Err(PyValueError::new_err("nthreads must be >= 1"));
    }
    Ok(Some(val as usize))
}

fn extract_1d_index_id<'py>(
    py: Python<'py>,
    require_fn: &Bound<'py, PyAny>,
    obj: &Bound<'py, PyAny>,
    upper_bound: usize,
    name: &str,
    bound_label: &str,
) -> PyResult<Vec<usize>> {
    if obj.is_none() {
        return Err(PyTypeError::new_err(format!("{name} must not be None")));
    }
    let np_mod = py.import_bound("numpy")?;
    let raw_arr = np_mod.getattr("asarray")?.call1((obj,))?;
    let dtype = raw_arr.getattr("dtype")?;
    let kind: String = dtype.getattr("kind")?.extract()?;
    if kind != "i" && kind != "u" {
        return Err(PyTypeError::new_err(format!(
            "{name} must have an integer dtype (int/uint)"
        )));
    }
    let int64_dtype = np_mod.getattr("int64")?;
    let reqs = PyTuple::new_bound(py, ["C", "A"]);
    let kwargs = PyDict::new_bound(py);
    kwargs.set_item("dtype", int64_dtype)?;
    kwargs.set_item("requirements", reqs)?;
    let aligned = require_fn.call((raw_arr,), Some(&kwargs))?;
    let arr: PyReadonlyArray1<'py, i64> = aligned.extract()?;
    if !arr.is_empty() && (!arr.is_c_contiguous() || (arr.data() as usize) % align_of::<i64>() != 0)
    {
        return Err(PyValueError::new_err(format!(
            "{name} must be aligned and C-contiguous"
        )));
    }
    let slice = arr
        .as_slice()
        .map_err(|_| PyValueError::new_err(format!("{name} must be C-contiguous")))?;
    let mut out = Vec::with_capacity(slice.len());
    for &val in slice {
        if val < 0 || (val as usize) >= upper_bound {
            return Err(PyValueError::new_err(format!(
                "{name} entry {val} out of bounds for {upper_bound} {bound_label}"
            )));
        }
        out.push(val as usize);
    }
    Ok(out)
}

struct ParsedScpcInputs {
    coords: Vec<f64>,
    n: usize,
    d: usize,
    metric: Metric,
    coef: Vec<f64>,
    wj: Vec<f64>,
    cluster_id: Option<Vec<usize>>,
    fe_id: Option<Vec<usize>>,
    xj: Option<Vec<f64>>,
    xjs: Option<Vec<f64>>,
    proj_u: Option<Vec<f64>>,
    proj_v: Option<Vec<f64>>,
    proj_rank: usize,
    nthreads: Option<usize>,
}

#[allow(clippy::too_many_arguments)]
fn parse_scpc_inputs<'py>(
    py: Python<'py>,
    coords: &Bound<'py, PyAny>,
    coef: &Bound<'py, PyAny>,
    wj: &Bound<'py, PyAny>,
    avc: f64,
    latlong: bool,
    uncond: bool,
    xj: Option<Bound<'py, PyAny>>,
    xjs: Option<Bound<'py, PyAny>>,
    proj_u: Option<Bound<'py, PyAny>>,
    proj_v: Option<Bound<'py, PyAny>>,
    nthreads: Option<Bound<'py, PyAny>>,
    cluster_id: Option<Bound<'py, PyAny>>,
    fe_id: Option<Bound<'py, PyAny>>,
) -> PyResult<ParsedScpcInputs> {
    let nthreads_usize = parse_nthreads(nthreads)?;

    if !avc.is_finite() || avc <= 0.001 || avc >= 0.99 {
        return Err(PyValueError::new_err(
            "avc must lie strictly in (0.001, 0.99)",
        ));
    }

    let np_mod = py.import_bound("numpy")?;
    let require_fn = np_mod.getattr("require")?;
    let float64_dtype = np_mod.getattr("float64")?;
    let reqs = PyTuple::new_bound(py, ["C", "A"]);
    let kwargs = PyDict::new_bound(py);
    kwargs.set_item("dtype", float64_dtype)?;
    kwargs.set_item("requirements", reqs)?;

    let (coords_vec, n, d) = extract_2d_row_major_f64(&require_fn, &kwargs, coords, "coords")?;
    if n < 2 || d == 0 {
        return Err(PyValueError::new_err(
            "coords must have at least 2 rows and 1 column",
        ));
    }
    if latlong && d != 2 {
        return Err(PyValueError::new_err(
            "geodesic coordinates require exactly 2 columns (lon, lat)",
        ));
    }

    let coef_vec = extract_1d_f64(&require_fn, &kwargs, coef, "coef")?;
    let k_use = coef_vec.len();
    if k_use == 0 {
        return Err(PyValueError::new_err("coef must be non-empty"));
    }

    let (wj_vec, wj_rows, wj_cols) = extract_2d_row_major_f64(&require_fn, &kwargs, wj, "wj")?;
    if wj_rows != n || wj_cols != k_use {
        return Err(PyValueError::new_err(format!(
            "wj shape ({wj_rows}, {wj_cols}) does not match expected ({n}, {k_use})"
        )));
    }

    let cluster_vec: Option<Vec<usize>> = match cluster_id.as_ref() {
        Some(obj) => Some(extract_1d_index_id(
            py,
            &require_fn,
            obj,
            n,
            "cluster_id",
            "clusters",
        )?),
        None => None,
    };
    let m_obs = cluster_vec.as_ref().map_or(n, Vec::len);
    let fe_vec: Option<Vec<usize>> = match fe_id.as_ref() {
        Some(obj) => Some(extract_1d_index_id(
            py,
            &require_fn,
            obj,
            m_obs,
            "fe_id",
            "observations",
        )?),
        None => None,
    };

    let mut xj_vec: Option<Vec<f64>> = None;
    let mut xjs_vec: Option<Vec<f64>> = None;
    let mut proj_u_vec: Option<Vec<f64>> = None;
    let mut proj_v_vec: Option<Vec<f64>> = None;
    let mut proj_rank = 0usize;

    if !uncond {
        let xj_arr = xj.as_ref().ok_or_else(|| {
            PyValueError::new_err("conditional SCPC requires xj direction matrix")
        })?;
        let xjs_arr = xjs.as_ref().ok_or_else(|| {
            PyValueError::new_err("conditional SCPC requires xjs scaling matrix")
        })?;

        let (xj_data, xj_rows, xj_cols) =
            extract_2d_row_major_f64(&require_fn, &kwargs, xj_arr, "xj")?;
        if xj_rows != m_obs || xj_cols != k_use {
            return Err(PyValueError::new_err(format!(
                "xj shape ({xj_rows}, {xj_cols}) does not match expected ({m_obs}, {k_use})"
            )));
        }
        let (xjs_data, xjs_rows, xjs_cols) =
            extract_2d_row_major_f64(&require_fn, &kwargs, xjs_arr, "xjs")?;
        if xjs_rows != m_obs || xjs_cols != k_use {
            return Err(PyValueError::new_err(format!(
                "xjs shape ({xjs_rows}, {xjs_cols}) does not match expected ({m_obs}, {k_use})"
            )));
        }
        xj_vec = Some(xj_data);
        xjs_vec = Some(xjs_data);

        match (proj_u.as_ref(), proj_v.as_ref()) {
            (None, None) => {}
            (Some(u_arr), Some(v_arr)) => {
                let (u_data, u_rows, u_cols) =
                    extract_2d_row_major_f64(&require_fn, &kwargs, u_arr, "proj_u")?;
                let (v_data, v_rows, v_cols) =
                    extract_2d_row_major_f64(&require_fn, &kwargs, v_arr, "proj_v")?;
                if u_rows != m_obs || v_rows != m_obs || u_cols != v_cols {
                    return Err(PyValueError::new_err(format!(
                        "proj_u shape ({u_rows}, {u_cols}) and proj_v shape ({v_rows}, {v_cols}) must both be ({m_obs}, r)"
                    )));
                }
                if u_cols > 0 {
                    proj_rank = u_cols;
                    proj_u_vec = Some(u_data);
                    proj_v_vec = Some(v_data);
                }
            }
            _ => {
                return Err(PyValueError::new_err(
                    "proj_u and proj_v must both be provided or both be None",
                ));
            }
        }
    }

    let metric = if latlong {
        Metric::GreatCircleFractionLonLat
    } else {
        Metric::Euclidean
    };

    Ok(ParsedScpcInputs {
        coords: coords_vec,
        n,
        d,
        metric,
        coef: coef_vec,
        wj: wj_vec,
        cluster_id: cluster_vec,
        fe_id: fe_vec,
        xj: xj_vec,
        xjs: xjs_vec,
        proj_u: proj_u_vec,
        proj_v: proj_v_vec,
        proj_rank,
        nthreads: nthreads_usize,
    })
}

fn build_scpc_result_dict<'py>(
    py: Python<'py>,
    output: ExactScpcOutput,
    n: usize,
    k_use: usize,
) -> PyResult<Bound<'py, PyDict>> {
    let result = PyDict::new_bound(py);
    result.set_item("scpcstats", vec_to_pyarray2(py, output.scpcstats, k_use, 6)?)?;
    match output.scpccvs {
        Some(cvs_data) => {
            result.set_item("scpccvs", vec_to_pyarray2(py, cvs_data, k_use, 4)?)?;
        }
        None => {
            result.set_item("scpccvs", py.None())?;
        }
    }
    result.set_item("wfin", vec_to_pyarray2(py, output.wfin, n, output.q + 1)?)?;
    result.set_item("c0", output.c0)?;
    result.set_item("cmax", output.cmax)?;
    result.set_item("cv", output.cvfin)?;
    result.set_item("q", output.q)?;
    result.set_item("tied", output.tied)?;
    Ok(result)
}

#[pyfunction]
#[pyo3(signature = (
    coords,
    coef,
    wj,
    avc,
    latlong = false,
    uncond = false,
    cvs = false,
    xj = None,
    xjs = None,
    proj_u = None,
    proj_v = None,
    nthreads = None,
    cluster_id = None,
    fe_id = None
))]
#[allow(clippy::too_many_arguments)]
fn run_exact_scpc<'py>(
    py: Python<'py>,
    coords: &Bound<'py, PyAny>,
    coef: &Bound<'py, PyAny>,
    wj: &Bound<'py, PyAny>,
    avc: f64,
    latlong: bool,
    uncond: bool,
    cvs: bool,
    xj: Option<Bound<'py, PyAny>>,
    xjs: Option<Bound<'py, PyAny>>,
    proj_u: Option<Bound<'py, PyAny>>,
    proj_v: Option<Bound<'py, PyAny>>,
    nthreads: Option<Bound<'py, PyAny>>,
    cluster_id: Option<Bound<'py, PyAny>>,
    fe_id: Option<Bound<'py, PyAny>>,
) -> PyResult<Bound<'py, PyDict>> {
    let inputs = parse_scpc_inputs(
        py, coords, coef, wj, avc, latlong, uncond, xj, xjs, proj_u, proj_v, nthreads,
        cluster_id, fe_id,
    )?;
    let n = inputs.n;
    let k_use = inputs.coef.len();
    let output = py
        .allow_threads(|| {
            let req = ExactScpcRequest {
                coords: &inputs.coords,
                n,
                d: inputs.d,
                metric: inputs.metric,
                avc,
                uncond,
                cvs,
                coef: &inputs.coef,
                wj: &inputs.wj,
                cluster_id: inputs.cluster_id.as_deref(),
                fe_id: inputs.fe_id.as_deref(),
                xj: inputs.xj.as_deref(),
                xjs: inputs.xjs.as_deref(),
                proj_u: inputs.proj_u.as_deref(),
                proj_v: inputs.proj_v.as_deref(),
                proj_rank: inputs.proj_rank,
                cv_tolerance: 0.001,
                nthreads: inputs.nthreads,
            };
            core_run_exact_scpc(&req)
        })
        .map_err(map_core_err)?;
    build_scpc_result_dict(py, output, n, k_use)
}

#[pyfunction]
#[pyo3(signature = (
    coords,
    coef,
    wj,
    avc,
    latlong = false,
    uncond = false,
    cvs = false,
    xj = None,
    xjs = None,
    proj_u = None,
    proj_v = None,
    nthreads = None,
    cluster_id = None,
    fe_id = None,
    seed = 1,
    landmarks = 1000,
    samples = 20,
    pair_count = 1_000_000
))]
#[allow(clippy::too_many_arguments)]
fn run_approx_scpc<'py>(
    py: Python<'py>,
    coords: &Bound<'py, PyAny>,
    coef: &Bound<'py, PyAny>,
    wj: &Bound<'py, PyAny>,
    avc: f64,
    latlong: bool,
    uncond: bool,
    cvs: bool,
    xj: Option<Bound<'py, PyAny>>,
    xjs: Option<Bound<'py, PyAny>>,
    proj_u: Option<Bound<'py, PyAny>>,
    proj_v: Option<Bound<'py, PyAny>>,
    nthreads: Option<Bound<'py, PyAny>>,
    cluster_id: Option<Bound<'py, PyAny>>,
    fe_id: Option<Bound<'py, PyAny>>,
    seed: u32,
    landmarks: usize,
    samples: usize,
    pair_count: usize,
) -> PyResult<Bound<'py, PyDict>> {
    let inputs = parse_scpc_inputs(
        py, coords, coef, wj, avc, latlong, uncond, xj, xjs, proj_u, proj_v, nthreads,
        cluster_id, fe_id,
    )?;
    let n = inputs.n;
    let k_use = inputs.coef.len();
    let output = py
        .allow_threads(|| {
            let req = ApproxScpcRequest {
                coords: &inputs.coords,
                n,
                d: inputs.d,
                metric: inputs.metric,
                avc,
                uncond,
                cvs,
                coef: &inputs.coef,
                wj: &inputs.wj,
                cluster_id: inputs.cluster_id.as_deref(),
                fe_id: inputs.fe_id.as_deref(),
                xj: inputs.xj.as_deref(),
                xjs: inputs.xjs.as_deref(),
                proj_u: inputs.proj_u.as_deref(),
                proj_v: inputs.proj_v.as_deref(),
                proj_rank: inputs.proj_rank,
                cv_tolerance: 0.001,
                nthreads: inputs.nthreads,
                seed,
                landmarks,
                samples,
                pair_count,
            };
            core_run_approx_scpc(&req)
        })
        .map_err(map_core_err)?;
    build_scpc_result_dict(py, output, n, k_use)
}

fn make_float64_kwargs<'py>(
    py: Python<'py>,
) -> PyResult<(Bound<'py, PyAny>, Bound<'py, PyDict>)> {
    let np_mod = py.import_bound("numpy")?;
    let require_fn = np_mod.getattr("require")?;
    let float64_dtype = np_mod.getattr("float64")?;
    let reqs = PyTuple::new_bound(py, ["C", "A"]);
    let kwargs = PyDict::new_bound(py);
    kwargs.set_item("dtype", float64_dtype)?;
    kwargs.set_item("requirements", reqs)?;
    Ok((require_fn, kwargs))
}

fn extract_y_matrix_f64<'py>(
    require_fn: &Bound<'py, PyAny>,
    kwargs: &Bound<'py, PyDict>,
    obj: &Bound<'py, PyAny>,
    expected_n: Option<usize>,
    name: &str,
) -> PyResult<(Vec<f64>, usize, usize, bool)> {
    let aligned = require_c_aligned_f64(require_fn, kwargs, obj, name)?;
    let ndim: usize = aligned.getattr("ndim")?.extract()?;
    match ndim {
        1 => {
            let v = extract_1d_f64(require_fn, kwargs, &aligned, name)?;
            let n = v.len();
            if let Some(exp) = expected_n {
                if n != exp {
                    return Err(PyValueError::new_err(format!(
                        "{name} length {n} does not match expected {exp}"
                    )));
                }
            }
            Ok((v, n, 1, true))
        }
        2 => {
            let (v, rows, cols) = extract_2d_row_major_f64(require_fn, kwargs, &aligned, name)?;
            if let Some(exp) = expected_n {
                if rows != exp {
                    return Err(PyValueError::new_err(format!(
                        "{name} shape ({rows}, {cols}) does not match expected rows {exp}"
                    )));
                }
            }
            Ok((v, rows, cols, false))
        }
        _ => Err(PyValueError::new_err(format!(
            "{name} must be a 1D or 2D float array"
        ))),
    }
}

fn parse_spur_coords<'py>(
    require_fn: &Bound<'py, PyAny>,
    kwargs: &Bound<'py, PyDict>,
    coords: &Bound<'py, PyAny>,
    latlong: bool,
) -> PyResult<(Vec<f64>, usize, usize, Metric)> {
    let (coords_vec, n, d) = extract_2d_row_major_f64(require_fn, kwargs, coords, "coords")?;
    if n < 2 || d == 0 {
        return Err(PyValueError::new_err(
            "coords must have at least 2 rows and 1 column",
        ));
    }
    if latlong && d != 2 {
        return Err(PyValueError::new_err(
            "geodesic coordinates require exactly 2 columns (lon, lat)",
        ));
    }
    let metric = if latlong {
        Metric::GreatCircleFractionLonLat
    } else {
        Metric::Euclidean
    };
    Ok((coords_vec, n, d, metric))
}

fn build_spur_i0_dict<'py>(
    py: Python<'py>,
    out: SpurI0Output,
    n_y: usize,
) -> PyResult<Bound<'py, PyDict>> {
    let result = PyDict::new_bound(py);
    result.set_item("lr", vec_to_pyarray1(py, out.lr))?;
    result.set_item("pvalue", vec_to_pyarray1(py, out.pvalue))?;
    result.set_item("cvalue", vec_to_pyarray1(py, out.cvalue))?;
    result.set_item("ha_param", out.ha_param)?;
    let n_rho = out.rho_grid.len();
    result.set_item("pvalue_mat", vec_to_pyarray2(py, out.pvalue_mat, n_rho, n_y)?)?;
    result.set_item("cvalue_mat", vec_to_pyarray2(py, out.cvalue_mat, n_rho, 3)?)?;
    result.set_item("rho_grid", vec_to_pyarray1(py, out.rho_grid))?;
    result.set_item("tied", out.tied)?;
    Ok(result)
}

fn build_spur_i1_dict<'py>(py: Python<'py>, out: SpurI1Output) -> PyResult<Bound<'py, PyDict>> {
    let result = PyDict::new_bound(py);
    result.set_item("lr", vec_to_pyarray1(py, out.lr))?;
    result.set_item("pvalue", vec_to_pyarray1(py, out.pvalue))?;
    result.set_item("cvalue", vec_to_pyarray1(py, out.cvalue))?;
    result.set_item("ha_param", out.ha_param)?;
    result.set_item("tied", out.tied)?;
    Ok(result)
}

fn build_spur_halflife_dict<'py>(
    py: Python<'py>,
    out: SpurHalfLifeOutput,
) -> PyResult<Bound<'py, PyDict>> {
    let result = PyDict::new_bound(py);
    result.set_item("ci_lower", out.ci_lower)?;
    result.set_item("ci_upper", out.ci_upper)?;
    result.set_item("max_dist", out.max_dist)?;
    result.set_item("pv_mat", vec_to_pyarray1(py, out.pv_mat))?;
    result.set_item("hl_grid_ho", vec_to_pyarray1(py, out.hl_grid_ho))?;
    result.set_item("tied", out.tied)?;
    Ok(result)
}

fn run_spur_test_impl<'py>(
    py: Python<'py>,
    coords: &Bound<'py, PyAny>,
    y: &Bound<'py, PyAny>,
    emat: &Bound<'py, PyAny>,
    q: usize,
    nrep: usize,
    latlong: bool,
    nthreads: Option<Bound<'py, PyAny>>,
    x_in: Option<&Bound<'py, PyAny>>,
    is_i0: bool,
) -> PyResult<Bound<'py, PyDict>> {
    let nthreads_usize = parse_nthreads(nthreads)?;
    let (require_fn, kwargs) = make_float64_kwargs(py)?;
    let (coords_vec, n, d, metric) = parse_spur_coords(&require_fn, &kwargs, coords, latlong)?;
    let (y_vec, _, n_y, _) = extract_y_matrix_f64(&require_fn, &kwargs, y, Some(n), "y")?;
    let (emat_vec, emat_rows, emat_cols) =
        extract_2d_row_major_f64(&require_fn, &kwargs, emat, "emat")?;
    if emat_rows != q || emat_cols != nrep {
        return Err(PyValueError::new_err(format!(
            "emat shape ({emat_rows}, {emat_cols}) does not match expected ({q}, {nrep})"
        )));
    }
    let (x_vec, k_x) = match x_in {
        Some(x_obj) => {
            let (xv, xr, xc) = extract_2d_row_major_f64(&require_fn, &kwargs, x_obj, "x_in")?;
            if xr != n {
                return Err(PyValueError::new_err(format!(
                    "x_in rows {xr} does not match coords rows {n}"
                )));
            }
            (Some(xv), xc)
        }
        None => (None, 0),
    };

    if is_i0 {
        let output = py
            .allow_threads(|| {
                let req = SpurTestRequest {
                    coords: &coords_vec,
                    n,
                    d,
                    metric,
                    y: &y_vec,
                    n_y,
                    x_in: x_vec.as_deref(),
                    k_x,
                    emat: &emat_vec,
                    q,
                    nrep,
                    nthreads: nthreads_usize,
                };
                if req.x_in.is_some() {
                    core_run_spur_i0_resid(&req)
                } else {
                    core_run_spur_i0(&req)
                }
            })
            .map_err(map_core_err)?;
        build_spur_i0_dict(py, output, n_y)
    } else {
        let output = py
            .allow_threads(|| {
                let req = SpurTestRequest {
                    coords: &coords_vec,
                    n,
                    d,
                    metric,
                    y: &y_vec,
                    n_y,
                    x_in: x_vec.as_deref(),
                    k_x,
                    emat: &emat_vec,
                    q,
                    nrep,
                    nthreads: nthreads_usize,
                };
                if req.x_in.is_some() {
                    core_run_spur_i1_resid(&req)
                } else {
                    core_run_spur_i1(&req)
                }
            })
            .map_err(map_core_err)?;
        build_spur_i1_dict(py, output)
    }
}

#[pyfunction]
#[pyo3(signature = (coords, y, emat, q, nrep, latlong = false, nthreads = None))]
#[allow(clippy::too_many_arguments)]
fn run_spur_i0<'py>(
    py: Python<'py>,
    coords: &Bound<'py, PyAny>,
    y: &Bound<'py, PyAny>,
    emat: &Bound<'py, PyAny>,
    q: usize,
    nrep: usize,
    latlong: bool,
    nthreads: Option<Bound<'py, PyAny>>,
) -> PyResult<Bound<'py, PyDict>> {
    run_spur_test_impl(py, coords, y, emat, q, nrep, latlong, nthreads, None, true)
}

#[pyfunction]
#[pyo3(signature = (coords, y, emat, q, nrep, latlong = false, nthreads = None))]
#[allow(clippy::too_many_arguments)]
fn run_spur_i1<'py>(
    py: Python<'py>,
    coords: &Bound<'py, PyAny>,
    y: &Bound<'py, PyAny>,
    emat: &Bound<'py, PyAny>,
    q: usize,
    nrep: usize,
    latlong: bool,
    nthreads: Option<Bound<'py, PyAny>>,
) -> PyResult<Bound<'py, PyDict>> {
    run_spur_test_impl(py, coords, y, emat, q, nrep, latlong, nthreads, None, false)
}

#[pyfunction]
#[pyo3(signature = (coords, y, x_in, emat, q, nrep, latlong = false, nthreads = None))]
#[allow(clippy::too_many_arguments)]
fn run_spur_i0_resid<'py>(
    py: Python<'py>,
    coords: &Bound<'py, PyAny>,
    y: &Bound<'py, PyAny>,
    x_in: &Bound<'py, PyAny>,
    emat: &Bound<'py, PyAny>,
    q: usize,
    nrep: usize,
    latlong: bool,
    nthreads: Option<Bound<'py, PyAny>>,
) -> PyResult<Bound<'py, PyDict>> {
    run_spur_test_impl(
        py,
        coords,
        y,
        emat,
        q,
        nrep,
        latlong,
        nthreads,
        Some(x_in),
        true,
    )
}

#[pyfunction]
#[pyo3(signature = (coords, y, x_in, emat, q, nrep, latlong = false, nthreads = None))]
#[allow(clippy::too_many_arguments)]
fn run_spur_i1_resid<'py>(
    py: Python<'py>,
    coords: &Bound<'py, PyAny>,
    y: &Bound<'py, PyAny>,
    x_in: &Bound<'py, PyAny>,
    emat: &Bound<'py, PyAny>,
    q: usize,
    nrep: usize,
    latlong: bool,
    nthreads: Option<Bound<'py, PyAny>>,
) -> PyResult<Bound<'py, PyDict>> {
    run_spur_test_impl(
        py,
        coords,
        y,
        emat,
        q,
        nrep,
        latlong,
        nthreads,
        Some(x_in),
        false,
    )
}

#[pyfunction]
#[pyo3(signature = (coords, y, emat, q, nrep, level = 0.95, latlong = false, nthreads = None))]
#[allow(clippy::too_many_arguments)]
fn run_spur_halflife<'py>(
    py: Python<'py>,
    coords: &Bound<'py, PyAny>,
    y: &Bound<'py, PyAny>,
    emat: &Bound<'py, PyAny>,
    q: usize,
    nrep: usize,
    level: f64,
    latlong: bool,
    nthreads: Option<Bound<'py, PyAny>>,
) -> PyResult<Bound<'py, PyDict>> {
    let nthreads_usize = parse_nthreads(nthreads)?;
    let (require_fn, kwargs) = make_float64_kwargs(py)?;
    let (coords_vec, n, d, metric) = parse_spur_coords(&require_fn, &kwargs, coords, latlong)?;
    let y_vec = extract_1d_f64(&require_fn, &kwargs, y, "y")?;
    if y_vec.len() != n {
        return Err(PyValueError::new_err(format!(
            "y length {} does not match coords rows {n}",
            y_vec.len()
        )));
    }
    let (emat_vec, emat_rows, emat_cols) =
        extract_2d_row_major_f64(&require_fn, &kwargs, emat, "emat")?;
    if emat_rows != q || emat_cols != nrep {
        return Err(PyValueError::new_err(format!(
            "emat shape ({emat_rows}, {emat_cols}) does not match expected ({q}, {nrep})"
        )));
    }
    let output = py
        .allow_threads(|| {
            core_run_spur_halflife(&SpurHalfLifeRequest {
                coords: &coords_vec,
                n,
                d,
                metric,
                y: &y_vec,
                emat: &emat_vec,
                q,
                nrep,
                level,
                nthreads: nthreads_usize,
            })
        })
        .map_err(map_core_err)?;
    build_spur_halflife_dict(py, output)
}

#[pyfunction]
#[pyo3(signature = (
    kind,
    n,
    y = None,
    coords = None,
    latlong = false,
    radius = None,
    strict_lt = true,
    exclude_zero_distances = false,
    cluster_id = None,
    return_matrix = false,
    demean_lbmgls = false,
    nthreads = None
))]
#[allow(clippy::too_many_arguments)]
fn run_spur_transform<'py>(
    py: Python<'py>,
    kind: &str,
    n: usize,
    y: Option<Bound<'py, PyAny>>,
    coords: Option<Bound<'py, PyAny>>,
    latlong: bool,
    radius: Option<f64>,
    strict_lt: bool,
    exclude_zero_distances: bool,
    cluster_id: Option<Bound<'py, PyAny>>,
    return_matrix: bool,
    demean_lbmgls: bool,
    nthreads: Option<Bound<'py, PyAny>>,
) -> PyResult<Bound<'py, PyDict>> {
    let nthreads_usize = parse_nthreads(nthreads)?;
    let (require_fn, kwargs) = make_float64_kwargs(py)?;
    let transform_kind = match kind {
        "nn" => SpurTransformKind::Nn {
            exclude_zero_distances,
        },
        "iso" => {
            let r = radius.ok_or_else(|| PyValueError::new_err("iso requires radius"))?;
            SpurTransformKind::Iso {
                radius: r,
                strict_lt,
            }
        }
        "cluster" => SpurTransformKind::Cluster,
        "lbmgls" => SpurTransformKind::LbmGls,
        _ => {
            return Err(PyValueError::new_err(format!(
                "unknown transform kind '{kind}' (expected nn, iso, cluster, or lbmgls)"
            )));
        }
    };

    let (coords_opt, d, metric) = match coords.as_ref() {
        Some(c_obj) => {
            let (cv, cn, cd, cm) = parse_spur_coords(&require_fn, &kwargs, c_obj, latlong)?;
            if cn != n {
                return Err(PyValueError::new_err(format!(
                    "coords rows {cn} does not match expected n={n}"
                )));
            }
            (Some(cv), cd, cm)
        }
        None => (None, 0, Metric::Euclidean),
    };

    let cluster_opt = match cluster_id.as_ref() {
        Some(cl_obj) => Some(extract_1d_index_id(
            py,
            &require_fn,
            cl_obj,
            n,
            "cluster_id",
            "observations",
        )?),
        None => None,
    };

    let (y_opt, n_y, y_was_1d) = match y.as_ref() {
        Some(y_obj) => {
            let (yv, _, ny, is_1d) =
                extract_y_matrix_f64(&require_fn, &kwargs, y_obj, Some(n), "y")?;
            (Some(yv), ny, is_1d)
        }
        None => (None, 0, false),
    };

    let output = py
        .allow_threads(|| {
            core_run_spur_transform(&SpurTransformRequest {
                coords: coords_opt.as_deref(),
                n,
                d,
                metric,
                kind: transform_kind,
                cluster_id: cluster_opt.as_deref(),
                y: y_opt.as_deref(),
                n_y,
                return_matrix,
                demean_lbmgls,
                nthreads: nthreads_usize,
            })
        })
        .map_err(map_core_err)?;

    let result = PyDict::new_bound(py);
    match output.matrix {
        Some(mat) => {
            result.set_item("matrix", vec_to_pyarray2(py, mat, n, n)?)?;
        }
        None => {
            result.set_item("matrix", py.None())?;
        }
    }
    match output.transformed {
        Some(tr) => {
            if y_was_1d {
                result.set_item("transformed", vec_to_pyarray1(py, tr))?;
            } else {
                result.set_item("transformed", vec_to_pyarray2(py, tr, n, n_y)?)?;
            }
        }
        None => {
            result.set_item("transformed", py.None())?;
        }
    }
    Ok(result)
}

#[pyfunction]
#[pyo3(signature = (
    coords,
    y,
    x_in,
    emat,
    q,
    nrep,
    transform_vars,
    latlong = false,
    demean_lbmgls = false,
    nthreads = None
))]
#[allow(clippy::too_many_arguments)]
fn run_spur_pipeline<'py>(
    py: Python<'py>,
    coords: &Bound<'py, PyAny>,
    y: &Bound<'py, PyAny>,
    x_in: &Bound<'py, PyAny>,
    emat: &Bound<'py, PyAny>,
    q: usize,
    nrep: usize,
    transform_vars: &Bound<'py, PyAny>,
    latlong: bool,
    demean_lbmgls: bool,
    nthreads: Option<Bound<'py, PyAny>>,
) -> PyResult<Bound<'py, PyDict>> {
    let nthreads_usize = parse_nthreads(nthreads)?;
    let (require_fn, kwargs) = make_float64_kwargs(py)?;
    let (coords_vec, n, d, metric) = parse_spur_coords(&require_fn, &kwargs, coords, latlong)?;
    let y_vec = extract_1d_f64(&require_fn, &kwargs, y, "y")?;
    if y_vec.len() != n {
        return Err(PyValueError::new_err(format!(
            "y length {} does not match coords rows {n}",
            y_vec.len()
        )));
    }
    let (x_vec, xr, k_x) = extract_2d_row_major_f64(&require_fn, &kwargs, x_in, "x_in")?;
    if xr != n {
        return Err(PyValueError::new_err(format!(
            "x_in rows {xr} does not match coords rows {n}"
        )));
    }
    let (emat_vec, emat_rows, emat_cols) =
        extract_2d_row_major_f64(&require_fn, &kwargs, emat, "emat")?;
    if emat_rows != q || emat_cols != nrep {
        return Err(PyValueError::new_err(format!(
            "emat shape ({emat_rows}, {emat_cols}) does not match expected ({q}, {nrep})"
        )));
    }
    let (tr_vec, tr_rows, n_tr) =
        extract_2d_row_major_f64(&require_fn, &kwargs, transform_vars, "transform_vars")?;
    if tr_rows != n {
        return Err(PyValueError::new_err(format!(
            "transform_vars rows {tr_rows} does not match coords rows {n}"
        )));
    }

    let output = py
        .allow_threads(|| {
            core_run_spur_pipeline(&SpurPipelineRequest {
                coords: &coords_vec,
                n,
                d,
                metric,
                y: &y_vec,
                x_in: &x_vec,
                k_x,
                emat: &emat_vec,
                q,
                nrep,
                transform_vars: &tr_vec,
                n_transform_vars: n_tr,
                demean_lbmgls,
                nthreads: nthreads_usize,
            })
        })
        .map_err(map_core_err)?;

    let result = PyDict::new_bound(py);
    result.set_item("i0", build_spur_i0_dict(py, output.i0, 1)?)?;
    result.set_item("i1", build_spur_i1_dict(py, output.i1)?)?;
    result.set_item("i0resid", build_spur_i0_dict(py, output.i0resid, 1)?)?;
    result.set_item("i1resid", build_spur_i1_dict(py, output.i1resid)?)?;
    result.set_item(
        "transformed_vars",
        vec_to_pyarray2(py, output.transformed_vars, n, n_tr)?,
    )?;
    Ok(result)
}

#[pymodule]
fn spatial_spur_rust(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("API_VERSION", 5usize)?;
    m.add("SPUR_API_VERSION", 2usize)?;
    m.add_function(wrap_pyfunction!(run_exact_scpc, m)?)?;
    m.add_function(wrap_pyfunction!(run_approx_scpc, m)?)?;
    m.add_function(wrap_pyfunction!(run_spur_i0, m)?)?;
    m.add_function(wrap_pyfunction!(run_spur_i1, m)?)?;
    m.add_function(wrap_pyfunction!(run_spur_i0_resid, m)?)?;
    m.add_function(wrap_pyfunction!(run_spur_i1_resid, m)?)?;
    m.add_function(wrap_pyfunction!(run_spur_halflife, m)?)?;
    m.add_function(wrap_pyfunction!(run_spur_transform, m)?)?;
    m.add_function(wrap_pyfunction!(run_spur_pipeline, m)?)?;
    Ok(())
}
