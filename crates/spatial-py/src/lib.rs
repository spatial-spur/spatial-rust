use numpy::ndarray::{Array1, Array2};
use numpy::{
    IntoPyArray, PyArray1, PyArray2, PyArrayMethods, PyReadonlyArray1, PyReadonlyArray2,
    PyUntypedArrayMethods,
};
use pyo3::exceptions::{PyNotImplementedError, PyRuntimeError, PyTypeError, PyValueError};
use pyo3::marker::Ungil;
use pyo3::prelude::*;
use pyo3::types::{PyBool, PyDict, PyTuple};
use spatial_core::{
    run_spur_halflife as core_run_spur_halflife, run_spur_i0 as core_run_spur_i0,
    run_spur_i0_resid as core_run_spur_i0_resid, run_spur_i1 as core_run_spur_i1,
    run_spur_i1_resid as core_run_spur_i1_resid, run_spur_pipeline as core_run_spur_pipeline,
    run_spur_transform as core_run_spur_transform, CoreError, ErrorKind, Metric, ScpcOutput,
    SpurHalfLifeOutput, SpurHalfLifeRequest, SpurI0Output, SpurI1Output, SpurPipelineRequest,
    SpurTestRequest, SpurTransformKind, SpurTransformRequest,
};
use std::mem::align_of;

// PyO3 0.22's exception macro checks its own `gil-refs` feature in the caller.
#[allow(unexpected_cfgs)]
mod exceptions {
    use pyo3::exceptions::PyValueError;
    pyo3::create_exception!(
        spatial_spur_rust,
        InvalidSampledCovarianceError,
        PyValueError
    );
    pyo3::create_exception!(spatial_spur_rust, TiedEigenspaceError, PyValueError);
}
use exceptions::{InvalidSampledCovarianceError, TiedEigenspaceError};

fn map_core_err(err: CoreError) -> PyErr {
    match err.kind() {
        ErrorKind::InvalidInput | ErrorKind::CalibrationFailed | ErrorKind::DegenerateVariance => {
            PyValueError::new_err(err.message())
        }
        ErrorKind::UnsupportedCase => PyNotImplementedError::new_err(err.message()),
        ErrorKind::TiedEigenspace => TiedEigenspaceError::new_err(err.message()),
        ErrorKind::InvalidSampledCovariance => {
            InvalidSampledCovarianceError::new_err(err.message())
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

#[cfg(unix)]
static NATIVE_PID: std::sync::atomic::AtomicU32 = std::sync::atomic::AtomicU32::new(0);

#[cfg(unix)]
fn fork_error() -> PyErr {
    PyRuntimeError::new_err(
        "Spatial Rust was used before this process forked. Use the multiprocessing 'spawn' or 'forkserver' start method for native calls.",
    )
}

fn check_not_forked() -> PyResult<()> {
    #[cfg(unix)]
    {
        use std::sync::atomic::Ordering;
        let parent = NATIVE_PID.load(Ordering::Acquire);
        if parent != 0 && parent != std::process::id() {
            return Err(fork_error());
        }
    }
    Ok(())
}

fn mark_native_used() -> PyResult<()> {
    #[cfg(unix)]
    {
        use std::sync::atomic::Ordering;
        let pid = std::process::id();
        if let Err(parent) =
            NATIVE_PID.compare_exchange(0, pid, Ordering::AcqRel, Ordering::Acquire)
        {
            if parent != pid {
                return Err(fork_error());
            }
        }
    }
    Ok(())
}

fn allow_threads_checked<T, F>(py: Python<'_>, f: F) -> PyResult<T>
where
    F: Ungil + FnOnce() -> T,
    T: Ungil,
{
    mark_native_used()?;
    Ok(py.allow_threads(f))
}

fn extract_integer<'py, T: FromPyObject<'py>>(value: &Bound<'py, PyAny>) -> PyResult<T> {
    if value.is_instance_of::<PyBool>() {
        return Err(PyTypeError::new_err("expected an integer, received bool"));
    }
    value.extract()
}

fn parse_nthreads(nthreads: Option<Bound<'_, PyAny>>) -> PyResult<Option<usize>> {
    check_not_forked()?;
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
    let count =
        Some(usize::try_from(val).map_err(|_| PyValueError::new_err("nthreads is too large"))?);
    spatial_core::threads::validate_nthreads(count).map_err(map_core_err)?;
    Ok(count)
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

fn build_scpc_result_dict<'py>(
    py: Python<'py>,
    output: ScpcOutput,
    n: usize,
    k_use: usize,
) -> PyResult<Bound<'py, PyDict>> {
    let result = PyDict::new_bound(py);
    result.set_item(
        "scpcstats",
        vec_to_pyarray2(py, output.scpcstats, k_use, 6)?,
    )?;
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

/// Run SCPC on fitted-model inputs. All numerical preparation and row permutations
/// take place inside the same native call as inference.
#[pyfunction]
#[pyo3(signature = (coords, coef, bread, scores, avc, *, design=None, instruments=None,
    direction_design=None, cluster_id=None, fe_id=None, latlong=false, cvs=false,
    nthreads=None, seed=None))]
#[allow(clippy::too_many_arguments)]
fn run_model_scpc<'py>(
    py: Python<'py>,
    coords: &Bound<'py, PyAny>,
    coef: &Bound<'py, PyAny>,
    bread: &Bound<'py, PyAny>,
    scores: &Bound<'py, PyAny>,
    avc: f64,
    design: Option<Bound<'py, PyAny>>,
    instruments: Option<Bound<'py, PyAny>>,
    direction_design: Option<Bound<'py, PyAny>>,
    cluster_id: Option<Bound<'py, PyAny>>,
    fe_id: Option<Bound<'py, PyAny>>,
    latlong: bool,
    cvs: bool,
    nthreads: Option<Bound<'py, PyAny>>,
    #[pyo3(from_py_with = "extract_integer")] seed: Option<u32>,
) -> PyResult<Bound<'py, PyDict>> {
    use spatial_core::{run_model_scpc as run, ScpcPreparation};

    let nthreads = parse_nthreads(nthreads)?;
    if !avc.is_finite() || avc <= 0.001 || avc >= 0.99 {
        return Err(PyValueError::new_err(
            "avc must lie strictly in (0.001, 0.99)",
        ));
    }
    let (require, kwargs) = make_float64_kwargs(py)?;
    let (coords, n_sites, d, metric) = parse_coords(&require, &kwargs, coords, latlong)?;
    let coefficients = extract_1d_f64(&require, &kwargs, coef, "coef")?;
    let k = coefficients.len();
    let (scores, n_obs, n_scores) = extract_2d_row_major_f64(&require, &kwargs, scores, "scores")?;
    let (bread, rows, cols) = extract_2d_row_major_f64(&require, &kwargs, bread, "bread")?;
    if k == 0 || n_scores < k || (rows, cols) != (k, n_scores) {
        return Err(PyValueError::new_err(
            "scores, bread, coefficients and coordinates have incompatible dimensions",
        ));
    }
    let cluster_ids = cluster_id
        .as_ref()
        .map(|x| extract_1d_index_id(py, &require, x, n_sites, "cluster_id", "sites"))
        .transpose()?;
    if cluster_ids.as_ref().map_or(n_sites, Vec::len) != n_obs {
        return Err(PyValueError::new_err(
            "scores must have one row per model observation",
        ));
    }
    let fixed_effects = fe_id
        .as_ref()
        .map(|x| extract_1d_index_id(py, &require, x, n_obs, "fe_id", "observations"))
        .transpose()?;
    if fixed_effects.as_ref().is_some_and(|fe| fe.len() != n_obs) {
        return Err(PyValueError::new_err(
            "fe_id must have one entry per observation",
        ));
    }
    let read_design =
        |value: Option<Bound<'py, PyAny>>, name: &str| -> PyResult<(Option<Vec<f64>>, usize)> {
            let Some(value) = value else {
                return Ok((None, 0));
            };
            let (data, rows, cols) = extract_2d_row_major_f64(&require, &kwargs, &value, name)?;
            if rows != n_obs || cols == 0 {
                return Err(PyValueError::new_err(format!(
                    "{name} must have one row per observation and at least one column",
                )));
            }
            Ok((Some(data), cols))
        };
    let (design, design_cols) = read_design(design, "design")?;
    let (instruments, instrument_cols) = read_design(instruments, "instruments")?;
    let (direction_design, direction_cols) = read_design(direction_design, "direction_design")?;
    if (direction_design.is_some() && direction_cols != n_scores)
        || (design.is_none()
            && (instruments.is_some() || direction_design.is_some() || fixed_effects.is_some()))
    {
        return Err(PyValueError::new_err(
            "conditional model matrices have incompatible dimensions",
        ));
    }
    let output = allow_threads_checked(py, || {
        run(
            ScpcPreparation {
                n_sites,
                n_obs,
                n_scores,
                coefficients: &coefficients,
                bread: &bread,
                scores: &scores,
                design: design.as_deref(),
                design_cols,
                instruments: instruments.as_deref(),
                instrument_cols,
                direction_design: direction_design.as_deref(),
                cluster_ids: cluster_ids.as_deref(),
                fixed_effects: fixed_effects.as_deref(),
            },
            &coords,
            d,
            metric,
            avc,
            cvs,
            nthreads,
            seed,
        )
    })?
    .map_err(map_core_err)?;
    build_scpc_result_dict(py, output, n_sites, k)
}

fn make_float64_kwargs<'py>(py: Python<'py>) -> PyResult<(Bound<'py, PyAny>, Bound<'py, PyDict>)> {
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

fn parse_coords<'py>(
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
    result.set_item(
        "pvalue_mat",
        vec_to_pyarray2(py, out.pvalue_mat, n_rho, n_y)?,
    )?;
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
    let (coords_vec, n, d, metric) = parse_coords(&require_fn, &kwargs, coords, latlong)?;
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
        let output = allow_threads_checked(py, || {
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
        })?
        .map_err(map_core_err)?;
        build_spur_i0_dict(py, output, n_y)
    } else {
        let output = allow_threads_checked(py, || {
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
        })?
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
    #[pyo3(from_py_with = "extract_integer")] q: usize,
    #[pyo3(from_py_with = "extract_integer")] nrep: usize,
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
    #[pyo3(from_py_with = "extract_integer")] q: usize,
    #[pyo3(from_py_with = "extract_integer")] nrep: usize,
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
    #[pyo3(from_py_with = "extract_integer")] q: usize,
    #[pyo3(from_py_with = "extract_integer")] nrep: usize,
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
    #[pyo3(from_py_with = "extract_integer")] q: usize,
    #[pyo3(from_py_with = "extract_integer")] nrep: usize,
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
    #[pyo3(from_py_with = "extract_integer")] q: usize,
    #[pyo3(from_py_with = "extract_integer")] nrep: usize,
    level: f64,
    latlong: bool,
    nthreads: Option<Bound<'py, PyAny>>,
) -> PyResult<Bound<'py, PyDict>> {
    let nthreads_usize = parse_nthreads(nthreads)?;
    let (require_fn, kwargs) = make_float64_kwargs(py)?;
    let (coords_vec, n, d, metric) = parse_coords(&require_fn, &kwargs, coords, latlong)?;
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
    let output = allow_threads_checked(py, || {
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
    })?
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
    #[pyo3(from_py_with = "extract_integer")] n: usize,
    y: Option<Bound<'py, PyAny>>,
    coords: Option<Bound<'py, PyAny>>,
    latlong: bool,
    radius: Option<Bound<'py, PyAny>>,
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
            let r_obj = radius.ok_or_else(|| PyValueError::new_err("iso requires radius"))?;
            if r_obj.is_instance_of::<PyBool>() {
                return Err(PyTypeError::new_err("radius must be a float, not bool"));
            }
            let r: f64 = r_obj.extract()?;
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
            let (cv, cn, cd, cm) = parse_coords(&require_fn, &kwargs, c_obj, latlong)?;
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

    let output = allow_threads_checked(py, || {
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
    })?
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
    #[pyo3(from_py_with = "extract_integer")] q: usize,
    #[pyo3(from_py_with = "extract_integer")] nrep: usize,
    transform_vars: &Bound<'py, PyAny>,
    latlong: bool,
    demean_lbmgls: bool,
    nthreads: Option<Bound<'py, PyAny>>,
) -> PyResult<Bound<'py, PyDict>> {
    let nthreads_usize = parse_nthreads(nthreads)?;
    let (require_fn, kwargs) = make_float64_kwargs(py)?;
    let (coords_vec, n, d, metric) = parse_coords(&require_fn, &kwargs, coords, latlong)?;
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

    let output = allow_threads_checked(py, || {
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
    })?
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

#[pyfunction]
fn build_info<'py>(py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
    let info = PyDict::new_bound(py);
    info.set_item("version", spatial_core::VERSION)?;
    info.set_item("spec_version", spatial_core::SPEC_VERSION)?;
    info.set_item("os", std::env::consts::OS)?;
    info.set_item("arch", std::env::consts::ARCH)?;
    info.set_item("linalg_provider", spatial_core::linalg_provider())?;
    info.set_item("api_version", spatial_core::SCPC_API_VERSION)?;
    info.set_item("spur_api_version", spatial_core::SPUR_API_VERSION)?;
    Ok(info)
}

#[pymodule]
fn spatial_spur_rust(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("API_VERSION", spatial_core::SCPC_API_VERSION)?;
    m.add("SPUR_API_VERSION", spatial_core::SPUR_API_VERSION)?;
    m.add("LINALG_PROVIDER", spatial_core::linalg_provider())?;
    m.add(
        "InvalidSampledCovarianceError",
        m.py().get_type_bound::<InvalidSampledCovarianceError>(),
    )?;
    m.add(
        "TiedEigenspaceError",
        m.py().get_type_bound::<TiedEigenspaceError>(),
    )?;
    m.add_function(wrap_pyfunction!(build_info, m)?)?;
    m.add_function(wrap_pyfunction!(run_model_scpc, m)?)?;
    m.add_function(wrap_pyfunction!(run_spur_i0, m)?)?;
    m.add_function(wrap_pyfunction!(run_spur_i1, m)?)?;
    m.add_function(wrap_pyfunction!(run_spur_i0_resid, m)?)?;
    m.add_function(wrap_pyfunction!(run_spur_i1_resid, m)?)?;
    m.add_function(wrap_pyfunction!(run_spur_halflife, m)?)?;
    m.add_function(wrap_pyfunction!(run_spur_transform, m)?)?;
    m.add_function(wrap_pyfunction!(run_spur_pipeline, m)?)?;
    Ok(())
}
