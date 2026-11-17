use numpy::ndarray::Array2;
use numpy::{IntoPyArray, PyArray2, PyReadonlyArray1, PyReadonlyArray2};
use pyo3::exceptions::{PyTypeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyBool, PyDict};
use spatial_core::{run_exact_scpc as core_run_exact_scpc, ExactScpcRequest, Metric};

fn extract_1d_f64(arr: &PyReadonlyArray1<'_, f64>, name: &str) -> PyResult<Vec<f64>> {
    let view = arr.as_array();
    let data = if let Some(slice) = view.as_slice() {
        slice.to_vec()
    } else {
        let std_layout = view.as_standard_layout();
        std_layout.iter().copied().collect()
    };
    if data.iter().any(|v| !v.is_finite()) {
        return Err(PyValueError::new_err(format!(
            "{name} must contain only finite values"
        )));
    }
    Ok(data)
}

fn extract_2d_row_major_f64(
    arr: &PyReadonlyArray2<'_, f64>,
    name: &str,
) -> PyResult<(Vec<f64>, usize, usize)> {
    let view = arr.as_array();
    let (rows, cols) = view.dim();
    let data = if let Some(slice) = view.as_slice() {
        slice.to_vec()
    } else {
        let owned = view.as_standard_layout().into_owned();
        owned.iter().copied().collect()
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
    nthreads = None
))]
#[allow(clippy::too_many_arguments)]
fn run_exact_scpc<'py>(
    py: Python<'py>,
    coords: PyReadonlyArray2<'py, f64>,
    coef: PyReadonlyArray1<'py, f64>,
    wj: PyReadonlyArray2<'py, f64>,
    avc: f64,
    latlong: bool,
    uncond: bool,
    cvs: bool,
    xj: Option<PyReadonlyArray2<'py, f64>>,
    xjs: Option<PyReadonlyArray2<'py, f64>>,
    proj_u: Option<PyReadonlyArray2<'py, f64>>,
    proj_v: Option<PyReadonlyArray2<'py, f64>>,
    nthreads: Option<Bound<'py, PyAny>>,
) -> PyResult<Bound<'py, PyDict>> {
    let nthreads_usize = parse_nthreads(nthreads)?;

    if !avc.is_finite() || avc <= 0.001 || avc >= 0.99 {
        return Err(PyValueError::new_err(
            "avc must lie strictly in (0.001, 0.99)",
        ));
    }

    let (coords_vec, n, d) = extract_2d_row_major_f64(&coords, "coords")?;
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

    let coef_vec = extract_1d_f64(&coef, "coef")?;
    let k_use = coef_vec.len();
    if k_use == 0 {
        return Err(PyValueError::new_err("coef must be non-empty"));
    }

    let (wj_vec, wj_rows, wj_cols) = extract_2d_row_major_f64(&wj, "wj")?;
    if wj_rows != n || wj_cols != k_use {
        return Err(PyValueError::new_err(format!(
            "wj shape ({wj_rows}, {wj_cols}) does not match expected ({n}, {k_use})"
        )));
    }

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

        let (xj_data, xj_rows, xj_cols) = extract_2d_row_major_f64(xj_arr, "xj")?;
        if xj_rows != n || xj_cols != k_use {
            return Err(PyValueError::new_err(format!(
                "xj shape ({xj_rows}, {xj_cols}) does not match expected ({n}, {k_use})"
            )));
        }
        let (xjs_data, xjs_rows, xjs_cols) = extract_2d_row_major_f64(xjs_arr, "xjs")?;
        if xjs_rows != n || xjs_cols != k_use {
            return Err(PyValueError::new_err(format!(
                "xjs shape ({xjs_rows}, {xjs_cols}) does not match expected ({n}, {k_use})"
            )));
        }
        xj_vec = Some(xj_data);
        xjs_vec = Some(xjs_data);

        match (proj_u.as_ref(), proj_v.as_ref()) {
            (None, None) => {}
            (Some(u_arr), Some(v_arr)) => {
                let (u_data, u_rows, u_cols) = extract_2d_row_major_f64(u_arr, "proj_u")?;
                let (v_data, v_rows, v_cols) = extract_2d_row_major_f64(v_arr, "proj_v")?;
                if u_rows != n || v_rows != n || u_cols != v_cols {
                    return Err(PyValueError::new_err(format!(
                        "proj_u shape ({u_rows}, {u_cols}) and proj_v shape ({v_rows}, {v_cols}) must both be ({n}, r)"
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

    let output = py
        .allow_threads(|| {
            let req = ExactScpcRequest {
                coords: &coords_vec,
                n,
                d,
                metric,
                avc,
                uncond,
                cvs,
                coef: &coef_vec,
                wj: &wj_vec,
                xj: xj_vec.as_deref(),
                xjs: xjs_vec.as_deref(),
                proj_u: proj_u_vec.as_deref(),
                proj_v: proj_v_vec.as_deref(),
                proj_rank,
                cv_tolerance: 0.001,
                nthreads: nthreads_usize,
            };
            core_run_exact_scpc(&req)
        })
        .map_err(|err| PyValueError::new_err(err.0))?;

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

#[pymodule]
fn spatial_spur_rust(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("API_VERSION", 2usize)?;
    m.add_function(wrap_pyfunction!(run_exact_scpc, m)?)?;
    Ok(())
}
