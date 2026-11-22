//! Safe Rust wrappers over Stata SPI v3 C shims (`c/stata_bridge.c`).

use std::ffi::{CStr, CString};
use std::fs::File;
use std::io::{BufReader, Read};
use std::os::raw::{c_char, c_double, c_int, c_uchar};

pub const STATA_OK: c_int = 0;
pub const STATA_ERR: c_int = 198;

fn storable(value: f64) -> bool {
    value.is_finite() && unsafe { st_is_missing(value) } == 0
}

#[cfg(not(test))]
extern "C" {
    fn st_error(message: *const c_char) -> c_int;
    fn st_display(message: *const c_char) -> c_int;
    fn st_in1() -> c_int;
    fn st_in2() -> c_int;
    fn st_nvars() -> c_int;
    fn st_ifobs(obs: c_int) -> c_uchar;
    fn st_is_missing(z: c_double) -> c_uchar;
    fn st_missval() -> c_double;
    fn st_vdata(var: c_int, obs: c_int, z: *mut c_double) -> c_int;
    fn st_vstore(var: c_int, obs: c_int, z: c_double) -> c_int;
    fn st_row(mat: *const c_char) -> c_int;
    fn st_col(mat: *const c_char) -> c_int;
    fn st_mat_el(mat: *const c_char, row: c_int, col: c_int, z: *mut c_double) -> c_int;
    fn st_mat_store(mat: *const c_char, row: c_int, col: c_int, z: c_double) -> c_int;
    fn st_scal_save(scal: *const c_char, z: c_double) -> c_int;
}

#[cfg(test)]
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
    pub(super) unsafe fn st_nvars() -> c_int {
        0
    }
    pub(super) unsafe fn st_ifobs(_obs: c_int) -> c_uchar {
        0
    }
    pub(super) unsafe fn st_is_missing(z: c_double) -> c_uchar {
        u8::from(!z.is_finite())
    }
    pub(super) unsafe fn st_missval() -> c_double {
        f64::NAN
    }
    pub(super) unsafe fn st_vdata(_var: c_int, _obs: c_int, _z: *mut c_double) -> c_int {
        STATA_ERR
    }
    pub(super) unsafe fn st_vstore(_var: c_int, _obs: c_int, _z: c_double) -> c_int {
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
        return Err(format!(
            "{label} variable index {last_var} exceeds c_int range"
        ));
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

/// Read a Stata matrix into a row-major `expected_rows x expected_cols` `Vec<f64>`.
///
/// Stata SPI v3 `SF_mat_el` uses 1-based row/col indices (`1..=rows`, `1..=cols`).
pub fn read_matrix(
    name: &str,
    expected_rows: usize,
    expected_cols: usize,
) -> Result<Vec<f64>, String> {
    let expected_rows_i32 = c_int::try_from(expected_rows).map_err(|_| {
        format!("Stata matrix '{name}' expected rows {expected_rows} exceeds i32::MAX")
    })?;
    let expected_cols_i32 = c_int::try_from(expected_cols).map_err(|_| {
        format!("Stata matrix '{name}' expected cols {expected_cols} exceeds i32::MAX")
    })?;
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

/// Write a row-major `obs.len() x ncols` `f64` buffer into `ncols` contiguous
/// dataset variables in the `plugin call` varlist starting at `var_1based_start`.
pub fn write_var_matrix(
    obs: &[c_int],
    var_1based_start: usize,
    ncols: usize,
    data: &[f64],
    label: &str,
) -> Result<(), String> {
    if ncols == 0 {
        return Ok(());
    }
    let nvars = plugin_var_count();
    let last_var = var_1based_start
        .checked_add(ncols)
        .and_then(|end| end.checked_sub(1))
        .ok_or_else(|| format!("{label} output variable index overflow"))?;
    if var_1based_start == 0 || last_var > nvars {
        return Err(format!(
            "{label} requires plugin output varlist indices {var_1based_start}..={last_var}, but only {nvars} variables were passed"
        ));
    }
    let expected = obs
        .len()
        .checked_mul(ncols)
        .ok_or_else(|| format!("{label} output size overflow"))?;
    if data.len() != expected {
        return Err(format!(
            "{label} output length {} does not match {}x{ncols}",
            data.len(),
            obs.len()
        ));
    }
    c_int::try_from(last_var)
        .map_err(|_| format!("{label} output variable index exceeds c_int range"))?;
    if obs.iter().any(|&row| row <= 0) {
        return Err(format!(
            "{label} output observation indices must be positive"
        ));
    }
    if let Some(index) = data.iter().position(|&value| !storable(value)) {
        return Err(format!(
            "invalid Stata output value in {label} at obs {}, col {}",
            obs[index / ncols],
            index % ncols + 1
        ));
    }
    for (row_idx, &st_obs) in obs.iter().enumerate() {
        for col_idx in 0..ncols {
            let var_idx = (var_1based_start + col_idx) as c_int;
            let val = data[row_idx * ncols + col_idx];
            let rc = unsafe { st_vstore(var_idx, st_obs, val) };
            if rc != STATA_OK {
                return Err(format!(
                    "SF_vstore failed for {label} var {var_idx} obs {st_obs} (rc={rc})"
                ));
            }
        }
    }
    Ok(())
}

/// Read `expected_len` little-endian IEEE-754 `f64` values written by Mata `fbufput(..., "%8z", ...)`.
pub fn read_f64_bin_file(path: &str, expected_len: usize, label: &str) -> Result<Vec<f64>, String> {
    let expected_bytes = expected_len
        .checked_mul(8)
        .ok_or_else(|| format!("{label} byte length overflow"))?;
    let file =
        File::open(path).map_err(|e| format!("failed to open {label} file '{path}': {e}"))?;
    let metadata = file
        .metadata()
        .map_err(|e| format!("failed to inspect {label} file '{path}': {e}"))?;
    if !metadata.is_file() || metadata.len() != expected_bytes as u64 {
        return Err(format!(
            "{label} binary file length {} bytes does not match expected {expected_bytes} bytes",
            metadata.len()
        ));
    }
    let mut reader = BufReader::new(file);
    let mut out = Vec::new();
    out.try_reserve_exact(expected_len)
        .map_err(|e| format!("failed to allocate {label} inputs: {e}"))?;
    let mut chunk = [0; 8];
    for idx in 0..expected_len {
        reader
            .read_exact(&mut chunk)
            .map_err(|e| format!("failed to read {label} file '{path}': {e}"))?;
        let val = f64::from_le_bytes(chunk);
        if !storable(val) {
            return Err(format!(
                "invalid Stata numeric entry in {label} file at index {idx}"
            ));
        }
        out.push(val);
    }
    if reader
        .read(&mut chunk[..1])
        .map_err(|e| format!("failed to read {label} file '{path}': {e}"))?
        != 0
    {
        return Err(format!(
            "unexpected trailing bytes in {label} file '{path}'"
        ));
    }
    Ok(out)
}

/// Store a row-major `rows x cols` slice into an existing Stata matrix.
///
/// Stata SPI v3 `SF_mat_store` uses 1-based row/col indices (`1..=rows`, `1..=cols`).
pub fn write_matrix(name: &str, rows: usize, cols: usize, data: &[f64]) -> Result<(), String> {
    write_matrix_values(name, rows, cols, data, false)
}

/// Stata represents an unbounded endpoint or unreachable cutoff with a missing value.
pub fn write_matrix_with_infinity(
    name: &str,
    rows: usize,
    cols: usize,
    data: &[f64],
) -> Result<(), String> {
    write_matrix_values(name, rows, cols, data, true)
}

fn write_matrix_values(
    name: &str,
    rows: usize,
    cols: usize,
    data: &[f64],
    allow_infinity: bool,
) -> Result<(), String> {
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
        if !storable(val) && !(allow_infinity && val == f64::INFINITY) {
            return Err(format!(
                "cannot store non-finite or missing value into matrix '{name}'[{},{}]",
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
            let val = if val == f64::INFINITY {
                unsafe { st_missval() }
            } else {
                val
            };
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
    if !storable(value) {
        return Err(format!(
            "cannot store non-finite or missing scalar '{name}'"
        ));
    }
    let c_name = to_cstring(name, "scalar name")?;
    let rc = unsafe { st_scal_save(c_name.as_ptr(), value) };
    if rc != STATA_OK {
        return Err(format!("SF_scal_save failed for '{name}' (rc={rc})"));
    }
    Ok(())
}

/// Save a finite scalar value, or Stata's missing value (`SV_missval`) when `None`.
pub fn save_scalar_or_missing(name: &str, value: Option<f64>) -> Result<(), String> {
    let raw = match value {
        Some(v) => {
            if !storable(v) {
                return Err(format!(
                    "cannot store non-finite or missing scalar '{name}'"
                ));
            }
            v
        }
        None => unsafe { st_missval() },
    };
    let c_name = to_cstring(name, "scalar name")?;
    let rc = unsafe { st_scal_save(c_name.as_ptr(), raw) };
    if rc != STATA_OK {
        return Err(format!("SF_scal_save failed for '{name}' (rc={rc})"));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

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
        assert!(save_scalar_or_missing("s", Some(f64::NAN)).is_err());
    }

    #[test]
    fn write_var_matrix_validates_bounds_and_inputs() {
        let obs = [1, 2];
        assert!(write_var_matrix(&obs, 1, 0, &[], "out").is_ok());
        let err_zero = write_var_matrix(&obs, 0, 1, &[1.0, 2.0], "out").unwrap_err();
        assert!(err_zero.contains("indices"));
        let err_var = write_var_matrix(&obs, 1, 1, &[1.0, 2.0], "out").unwrap_err();
        assert!(err_var.contains("varlist indices"));
        let err_len = write_var_matrix(&obs, 1, 1, &[1.0], "out").unwrap_err();
        assert!(err_len.contains("varlist indices") || err_len.contains("length"));
    }

    #[test]
    fn read_f64_bin_file_round_trips_and_validates() {
        let dir = std::env::temp_dir();
        let path = dir.join(format!("spatial_stata_spi_test_{}.bin", std::process::id()));
        let path_str = path.to_str().unwrap();

        let vals = [1.5_f64, -2.25_f64, 0.0_f64, 42.125_f64];
        {
            let mut f = File::create(&path).unwrap();
            for v in vals {
                f.write_all(&v.to_le_bytes()).unwrap();
            }
        }
        let parsed = read_f64_bin_file(path_str, vals.len(), "emat").unwrap();
        assert_eq!(parsed, vals);
        assert!(read_f64_bin_file(path_str, vals.len() + 1, "emat").is_err());

        {
            let mut f = File::create(&path).unwrap();
            f.write_all(&f64::NAN.to_le_bytes()).unwrap();
        }
        assert!(read_f64_bin_file(path_str, 1, "emat").is_err());
        let _ = std::fs::remove_file(&path);
        assert!(read_f64_bin_file(path_str, 1, "emat").is_err());
    }
}
