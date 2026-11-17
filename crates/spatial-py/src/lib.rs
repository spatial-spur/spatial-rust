use numpy::ndarray::Array2;
use numpy::{
    IntoPyArray, PyArray2, PyArrayMethods, PyReadonlyArray1, PyReadonlyArray2,
    PyUntypedArrayMethods,
};
use pyo3::exceptions::{PyNotImplementedError, PyRuntimeError, PyTypeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyBool, PyDict, PyTuple};
use spatial_core::{
    run_approx_scpc as core_run_approx_scpc, run_exact_scpc as core_run_exact_scpc,
    ApproxScpcRequest, CoreError, ErrorKind, ExactScpcOutput, ExactScpcRequest, Metric,
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

#[pymodule]
fn spatial_spur_rust(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("API_VERSION", 5usize)?;
    m.add_function(wrap_pyfunction!(run_exact_scpc, m)?)?;
    m.add_function(wrap_pyfunction!(run_approx_scpc, m)?)?;
    Ok(())
}
