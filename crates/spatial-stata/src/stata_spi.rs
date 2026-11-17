//! Safe Rust wrappers over Stata SPI v3 C shims (`c/stata_bridge.c`).

use std::ffi::{CStr, CString};
use std::os::raw::{c_char, c_double, c_int, c_uchar};

pub const STATA_OK: c_int = 0;
pub const STATA_ERR: c_int = 198;

#[cfg(not(test))]
#[allow(dead_code)]
extern "C" {
    fn st_error(message: *const c_char) -> c_int;
    fn st_display(message: *const c_char) -> c_int;
    fn st_in1() -> c_int;
    fn st_in2() -> c_int;
    fn st_nobs() -> c_int;
    fn st_nvar() -> c_int;
    fn st_nvars() -> c_int;
    fn st_ifobs(obs: c_int) -> c_uchar;
    fn st_is_missing(z: c_double) -> c_uchar;
    fn st_vdata(var: c_int, obs: c_int, z: *mut c_double) -> c_int;
    fn st_row(mat: *const c_char) -> c_int;
    fn st_col(mat: *const c_char) -> c_int;
    fn st_mat_el(mat: *const c_char, row: c_int, col: c_int, z: *mut c_double) -> c_int;
    fn st_mat_store(mat: *const c_char, row: c_int, col: c_int, z: c_double) -> c_int;
    fn st_scal_save(scal: *const c_char, z: c_double) -> c_int;
}

#[cfg(test)]
#[allow(dead_code)]
mod test_stubs {
    use super::*;

    pub(super) unsafe fn st_error(_message: *const c_char) -> c_int {
        0
    }
    pub(super) unsafe fn st_display(_message: *const c_char) -> c_int {
        0
    }
    pub(super) unsafe fn st_in1() -> c_int {
        1
    }
    pub(super) unsafe fn st_in2() -> c_int {
        0
    }
    pub(super) unsafe fn st_nobs() -> c_int {
        0
    }
    pub(super) unsafe fn st_nvar() -> c_int {
        0
    }
    pub(super) unsafe fn st_nvars() -> c_int {
        0
    }
    pub(super) unsafe fn st_ifobs(_obs: c_int) -> c_uchar {
        0
    }
    pub(super) unsafe fn st_is_missing(z: c_double) -> c_uchar {
        u8::from(!z.is_finite())
    }
    pub(super) unsafe fn st_vdata(_var: c_int, _obs: c_int, _z: *mut c_double) -> c_int {
        STATA_ERR
    }
    pub(super) unsafe fn st_row(_mat: *const c_char) -> c_int {
        0
    }
    pub(super) unsafe fn st_col(_mat: *const c_char) -> c_int {
        0
    }
    pub(super) unsafe fn st_mat_el(
        _mat: *const c_char,
        _row: c_int,
        _col: c_int,
        _z: *mut c_double,
    ) -> c_int {
        STATA_ERR
    }
    pub(super) unsafe fn st_mat_store(
        _mat: *const c_char,
        _row: c_int,
        _col: c_int,
        _z: c_double,
    ) -> c_int {
        STATA_ERR
    }
    pub(super) unsafe fn st_scal_save(_scal: *const c_char, _z: c_double) -> c_int {
        STATA_ERR
    }
}

#[cfg(test)]
use test_stubs::*;

fn to_cstring(text: &str, label: &str) -> Result<CString, String> {
    CString::new(text).map_err(|_| format!("{label} contains interior NUL byte"))
}

fn ensure_newline_cstring(message: &str) -> CString {
    let sanitized = message.replace('\0', " ");
    let with_newline = if sanitized.ends_with('\n') {
        sanitized
    } else {
        format!("{sanitized}\n")
    };
    CString::new(with_newline).expect("sanitized message contains no interior NUL")
}

/// Print an informational message to the Stata console.
#[allow(dead_code)]
pub fn display(message: &str) {
    let cstr = ensure_newline_cstring(message);
    unsafe {
        let _ = st_display(cstr.as_ptr());
    }
}

/// Print an error message to the Stata console.
pub fn error(message: &str) {
    let cstr = ensure_newline_cstring(message);
    unsafe {
        let _ = st_error(cstr.as_ptr());
    }
}

/// Return total dataset observation count (`SF_nobs`).
#[allow(dead_code)]
pub fn total_obs() -> usize {
    usize::try_from(unsafe { st_nobs() }).unwrap_or(0)
}

/// Return total dataset variable count (`SF_nvar`).
#[allow(dead_code)]
pub fn dataset_var_count() -> usize {
    usize::try_from(unsafe { st_nvar() }).unwrap_or(0)
}

/// Return number of variables in the `plugin call` varlist (`SF_nvars`).
pub fn plugin_var_count() -> usize {
    usize::try_from(unsafe { st_nvars() }).unwrap_or(0)
}

/// Collect 1-based Stata observation indices in `SF_in1()..=SF_in2()` passing `SF_ifobs`.
pub fn active_obs() -> Vec<c_int> {
    let in1 = unsafe { st_in1() };
    let in2 = unsafe { st_in2() };
    if in1 < 1 || in2 < in1 {
        return Vec::new();
    }
    let cap = usize::try_from(i64::from(in2) - i64::from(in1) + 1).unwrap_or(0);
    let mut obs = Vec::with_capacity(cap);
    for row in in1..=in2 {
        if unsafe { st_ifobs(row) } != 0 {
            obs.push(row);
        }
    }
    obs
}

/// Read `ncols` contiguous variables from the `plugin call` varlist starting at
/// 1-based `var_1based_start` across `obs`, returning row-major `obs.len() x ncols`.
pub fn read_var_matrix(
    obs: &[c_int],
    var_1based_start: usize,
    ncols: usize,
    label: &str,
) -> Result<Vec<f64>, String> {
    if ncols == 0 {
        return Ok(Vec::new());
    }
    let nvars = plugin_var_count();
    let last_var = var_1based_start
        .checked_add(ncols)
        .and_then(|end| end.checked_sub(1))
        .ok_or_else(|| format!("{label} variable index overflow"))?;
    if var_1based_start == 0 || last_var > nvars {
        return Err(format!(
            "{label} requires plugin varlist indices {var_1based_start}..={last_var}, but only {nvars} variables were passed"
        ));
    }
    if c_int::try_from(last_var).is_err() {
        return Err(format!("{label} variable index {last_var} exceeds c_int range"));
    }

    let n = obs.len();
    let total = n
        .checked_mul(ncols)
        .ok_or_else(|| format!("{label} buffer size overflow"))?;
    let mut out = vec![0.0; total];

    for (row_idx, &st_obs) in obs.iter().enumerate() {
        for col_idx in 0..ncols {
            let var_idx = c_int::try_from(var_1based_start + col_idx)
                .map_err(|_| format!("{label} variable index exceeds c_int range"))?;
            let mut val: c_double = 0.0;
            let rc = unsafe { st_vdata(var_idx, st_obs, &mut val) };
            if rc != STATA_OK {
                return Err(format!(
                    "SF_vdata failed for {label} var {var_idx} obs {st_obs} (rc={rc})"
                ));
            }
            if unsafe { st_is_missing(val) } != 0 || !val.is_finite() {
                return Err(format!(
                    "missing or non-finite value in {label} at obs {st_obs}, col {}",
                    col_idx + 1
                ));
            }
            out[row_idx * ncols + col_idx] = val;
        }
    }
    Ok(out)
}

fn matrix_dims_cstr(name: &str, c_name: &CStr) -> Result<(c_int, c_int), String> {
    let rows = unsafe { st_row(c_name.as_ptr()) };
    let cols = unsafe { st_col(c_name.as_ptr()) };
    if rows <= 0 || cols <= 0 {
        return Err(format!(
            "Stata matrix '{name}' is missing or empty ({rows}x{cols})"
        ));
    }
    Ok((rows, cols))
}

/// Query Stata matrix dimensions `(rows, cols)`.
#[allow(dead_code)]
pub fn matrix_dims(name: &str) -> Result<(usize, usize), String> {
    let c_name = to_cstring(name, "matrix name")?;
    let (rows, cols) = matrix_dims_cstr(name, &c_name)?;
    let rows_u = usize::try_from(rows).map_err(|_| format!("invalid row count {rows} for '{name}'"))?;
    let cols_u = usize::try_from(cols).map_err(|_| format!("invalid col count {cols} for '{name}'"))?;
    Ok((rows_u, cols_u))
}

/// Read a Stata matrix into a row-major `expected_rows x expected_cols` `Vec<f64>`.
///
/// Stata SPI v3 `SF_mat_el` uses 1-based row/col indices (`1..=rows`, `1..=cols`).
pub fn read_matrix(
    name: &str,
    expected_rows: usize,
    expected_cols: usize,
) -> Result<Vec<f64>, String> {
    let expected_rows_i32 = c_int::try_from(expected_rows)
        .map_err(|_| format!("Stata matrix '{name}' expected rows {expected_rows} exceeds i32::MAX"))?;
    let expected_cols_i32 = c_int::try_from(expected_cols)
        .map_err(|_| format!("Stata matrix '{name}' expected cols {expected_cols} exceeds i32::MAX"))?;
    let c_name = to_cstring(name, "matrix name")?;
    let (mat_rows, mat_cols) = matrix_dims_cstr(name, &c_name)?;
    if mat_rows != expected_rows_i32 || mat_cols != expected_cols_i32 {
        return Err(format!(
            "Stata matrix '{name}' has dimensions {mat_rows}x{mat_cols}, expected {expected_rows_i32}x{expected_cols_i32}"
        ));
    }
    let total = expected_rows
        .checked_mul(expected_cols)
        .ok_or_else(|| format!("Stata matrix '{name}' size overflow"))?;
    let mut out = vec![0.0; total];
    for r in 0..expected_rows {
        let r_i32 = c_int::try_from(r + 1)
            .map_err(|_| format!("Stata matrix '{name}' row index overflow"))?;
        for c in 0..expected_cols {
            let c_i32 = c_int::try_from(c + 1)
                .map_err(|_| format!("Stata matrix '{name}' col index overflow"))?;
            let mut val: c_double = 0.0;
            let rc = unsafe { st_mat_el(c_name.as_ptr(), r_i32, c_i32, &mut val) };
            if rc != STATA_OK {
                return Err(format!(
                    "SF_mat_el failed for '{name}'[{r_i32},{c_i32}] (rc={rc})"
                ));
            }
            if unsafe { st_is_missing(val) } != 0 || !val.is_finite() {
                return Err(format!(
                    "missing or non-finite entry in matrix '{name}'[{r_i32},{c_i32}]"
                ));
            }
            out[r * expected_cols + c] = val;
        }
    }
    Ok(out)
}

/// Store a row-major `rows x cols` slice into an existing Stata matrix.
///
/// Stata SPI v3 `SF_mat_store` uses 1-based row/col indices (`1..=rows`, `1..=cols`).
pub fn write_matrix(name: &str, rows: usize, cols: usize, data: &[f64]) -> Result<(), String> {
    if rows == 0 || cols == 0 {
        return Err(format!(
            "target Stata matrix '{name}' dimensions must be positive (got {rows}x{cols})"
        ));
    }
    let rows_i32 = c_int::try_from(rows)
        .map_err(|_| format!("Stata matrix '{name}' row count {rows} exceeds i32::MAX"))?;
    let cols_i32 = c_int::try_from(cols)
        .map_err(|_| format!("Stata matrix '{name}' col count {cols} exceeds i32::MAX"))?;
    let expected_len = rows
        .checked_mul(cols)
        .ok_or_else(|| format!("Stata matrix '{name}' size overflow"))?;
    if data.len() != expected_len {
        return Err(format!(
            "output buffer length {} does not match {rows}x{cols} for matrix '{name}'",
            data.len()
        ));
    }
    for (idx, &val) in data.iter().enumerate() {
        if !val.is_finite() {
            return Err(format!(
                "cannot store non-finite value into matrix '{name}'[{},{}]",
                idx / cols + 1,
                idx % cols + 1
            ));
        }
    }
    let c_name = to_cstring(name, "matrix name")?;
    let mat_rows = unsafe { st_row(c_name.as_ptr()) };
    let mat_cols = unsafe { st_col(c_name.as_ptr()) };
    if mat_rows <= 0 || mat_cols <= 0 {
        return Err(format!(
            "Stata matrix '{name}' is missing or empty ({mat_rows}x{mat_cols})"
        ));
    }
    if mat_rows != rows_i32 || mat_cols != cols_i32 {
        return Err(format!(
            "target Stata matrix '{name}' has dimensions {mat_rows}x{mat_cols}, expected {rows_i32}x{cols_i32}"
        ));
    }
    for r in 0..rows {
        let r_i32 = c_int::try_from(r + 1)
            .map_err(|_| format!("Stata matrix '{name}' row index overflow"))?;
        for c in 0..cols {
            let c_i32 = c_int::try_from(c + 1)
                .map_err(|_| format!("Stata matrix '{name}' col index overflow"))?;
            let val = data[r * cols + c];
            let rc = unsafe { st_mat_store(c_name.as_ptr(), r_i32, c_i32, val) };
            if rc != STATA_OK {
                return Err(format!(
                    "SF_mat_store failed for '{name}'[{r_i32},{c_i32}] (rc={rc})"
                ));
            }
        }
    }
    Ok(())
}

/// Save a finite scalar value into a Stata scalar.
pub fn save_scalar(name: &str, value: f64) -> Result<(), String> {
    if !value.is_finite() {
        return Err(format!("cannot store non-finite scalar '{name}'"));
    }
    let c_name = to_cstring(name, "scalar name")?;
    let rc = unsafe { st_scal_save(c_name.as_ptr(), value) };
    if rc != STATA_OK {
        return Err(format!("SF_scal_save failed for '{name}' (rc={rc})"));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn write_matrix_validates_dims_and_inputs_before_store() {
        assert!(write_matrix("m", 0, 2, &[]).is_err());
        assert!(write_matrix("m", usize::MAX, 1, &[]).is_err());
        assert!(write_matrix("m", 2, 2, &[1.0, 2.0, 3.0]).is_err());
        assert!(write_matrix("m", 2, 2, &[1.0, f64::NAN, 3.0, 4.0]).is_err());
        // test_stubs::st_row/st_col return 0, so dimension check rejects missing matrix
        let err = write_matrix("m", 2, 2, &[1.0, 2.0, 3.0, 4.0]).unwrap_err();
        assert!(err.contains("missing or empty"));
        let err_read = read_matrix("m", usize::MAX, 1).unwrap_err();
        assert!(err_read.contains("i32::MAX"));
    }
}
