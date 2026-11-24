//! Minimal direct R C-API bindings using only official public `<Rinternals.h>` symbols.

use std::ffi::{c_char, c_int, c_uint, CStr};
use std::slice;

/// Opaque pointer to an R `SEXPREC`.
#[repr(C)]
pub struct SEXPREC {
    _private: [u8; 0],
}

pub type SEXP = *mut SEXPREC;
type Rboolean = c_int;
#[allow(non_camel_case_types)]
type R_xlen_t = isize;
type SEXPTYPE = c_uint;

const NILSXP: c_int = 0;
const LGLSXP: SEXPTYPE = 10;
const INTSXP: SEXPTYPE = 13;
const REALSXP: SEXPTYPE = 14;
const STRSXP: SEXPTYPE = 16;
const VECSXP: SEXPTYPE = 19;
const CE_UTF8: c_int = 1;

extern "C" {
    fn TYPEOF(x: SEXP) -> c_int;
    fn Rf_xlength(x: SEXP) -> R_xlen_t;
    fn Rf_nrows(x: SEXP) -> c_int;
    fn Rf_ncols(x: SEXP) -> c_int;
    fn Rf_isMatrix(x: SEXP) -> Rboolean;
    fn Rf_isFunction(x: SEXP) -> Rboolean;
    fn REAL_RO(x: SEXP) -> *const f64;
    fn REAL(x: SEXP) -> *mut f64;
    fn INTEGER_RO(x: SEXP) -> *const i32;
    fn LOGICAL_RO(x: SEXP) -> *const i32;
    fn STRING_ELT(x: SEXP, i: R_xlen_t) -> SEXP;
    fn R_CHAR(x: SEXP) -> *const c_char;
    fn Rf_allocVector(tp: SEXPTYPE, n: R_xlen_t) -> SEXP;
    fn Rf_allocMatrix(tp: SEXPTYPE, nrow: c_int, ncol: c_int) -> SEXP;
    fn Rf_ScalarReal(v: f64) -> SEXP;
    fn Rf_ScalarInteger(v: i32) -> SEXP;
    fn Rf_ScalarLogical(v: c_int) -> SEXP;
    fn Rf_ScalarString(s: SEXP) -> SEXP;
    fn Rf_mkCharLenCE(s: *const c_char, len: c_int, enc: c_int) -> SEXP;
    fn SET_STRING_ELT(x: SEXP, i: R_xlen_t, v: SEXP);
    fn SET_VECTOR_ELT(x: SEXP, i: R_xlen_t, v: SEXP) -> SEXP;
    fn Rf_setAttrib(vec: SEXP, name: SEXP, val: SEXP) -> SEXP;
    fn Rf_protect(x: SEXP) -> SEXP;
    fn Rf_unprotect(count: c_int);
    fn Rf_lang2(fun: SEXP, arg1: SEXP) -> SEXP;
    fn R_tryEvalSilent(expr: SEXP, env: SEXP, error_occurred: *mut c_int) -> SEXP;

    fn sp_r_nil() -> SEXP;
    fn sp_r_names_sym() -> SEXP;
    fn sp_r_class_sym() -> SEXP;
    fn sp_r_global_env() -> SEXP;
}

/// Stack-based RAII protection guard ensuring balanced `Rf_unprotect` even if an early return or panic happens.
pub struct ProtectGuard {
    count: c_int,
}

impl ProtectGuard {
    pub fn new() -> Self {
        Self { count: 0 }
    }

    pub unsafe fn protect(&mut self, sexp: SEXP) -> SEXP {
        Rf_protect(sexp);
        self.count = self.count.saturating_add(1);
        sexp
    }
}

impl Drop for ProtectGuard {
    fn drop(&mut self) {
        if self.count > 0 {
            unsafe { Rf_unprotect(self.count) };
        }
    }
}

pub unsafe fn r_nil() -> SEXP {
    sp_r_nil()
}

pub unsafe fn is_null(sexp: SEXP) -> bool {
    sexp.is_null() || TYPEOF(sexp) == NILSXP
}

pub struct RealMatrixView<'a> {
    pub data: &'a [f64],
    pub nrows: usize,
    pub ncols: usize,
}

pub unsafe fn as_f64_matrix<'a>(sexp: SEXP, name: &str) -> Result<RealMatrixView<'a>, String> {
    if TYPEOF(sexp) != REALSXP as c_int || Rf_isMatrix(sexp) == 0 {
        return Err(format!("`{name}` must be a double-precision numeric matrix"));
    }
    let nrows = Rf_nrows(sexp);
    let ncols = Rf_ncols(sexp);
    if nrows < 0 || ncols < 0 {
        return Err(format!("`{name}` has invalid matrix dimensions"));
    }
    let len = Rf_xlength(sexp);
    let nrows_usize = nrows as usize;
    let ncols_usize = ncols as usize;
    let expected = nrows_usize
        .checked_mul(ncols_usize)
        .ok_or_else(|| format!("`{name}` matrix dimensions overflow"))?;
    if len != expected as R_xlen_t {
        return Err(format!("`{name}` length does not match matrix dimensions"));
    }
    let ptr = REAL_RO(sexp);
    let data = if expected == 0 {
        &[]
    } else {
        if ptr.is_null() {
            return Err(format!("`{name}` data pointer is null"));
        }
        slice::from_raw_parts(ptr, expected)
    };
    Ok(RealMatrixView {
        data,
        nrows: nrows_usize,
        ncols: ncols_usize,
    })
}

pub unsafe fn as_opt_f64_matrix<'a>(
    sexp: SEXP,
    name: &str,
) -> Result<Option<RealMatrixView<'a>>, String> {
    if is_null(sexp) {
        Ok(None)
    } else {
        Ok(Some(as_f64_matrix(sexp, name)?))
    }
}

pub unsafe fn as_f64_slice<'a>(sexp: SEXP, name: &str) -> Result<&'a [f64], String> {
    if TYPEOF(sexp) != REALSXP as c_int {
        return Err(format!("`{name}` must be a double-precision vector"));
    }
    let len = Rf_xlength(sexp);
    if len < 0 {
        return Err(format!("`{name}` has negative length"));
    }
    let len_usize = len as usize;
    if len_usize == 0 {
        return Ok(&[]);
    }
    let ptr = REAL_RO(sexp);
    if ptr.is_null() {
        return Err(format!("`{name}` data pointer is null"));
    }
    Ok(slice::from_raw_parts(ptr, len_usize))
}

pub unsafe fn as_opt_i32_slice<'a>(sexp: SEXP, name: &str) -> Result<Option<&'a [i32]>, String> {
    if is_null(sexp) {
        return Ok(None);
    }
    if TYPEOF(sexp) != INTSXP as c_int {
        return Err(format!("`{name}` must be an integer vector"));
    }
    let len = Rf_xlength(sexp);
    if len < 0 {
        return Err(format!("`{name}` has negative length"));
    }
    let len_usize = len as usize;
    if len_usize == 0 {
        return Ok(Some(&[]));
    }
    let ptr = INTEGER_RO(sexp);
    if ptr.is_null() {
        return Err(format!("`{name}` data pointer is null"));
    }
    Ok(Some(slice::from_raw_parts(ptr, len_usize)))
}

pub unsafe fn as_bool(sexp: SEXP, name: &str) -> Result<bool, String> {
    if TYPEOF(sexp) != LGLSXP as c_int || Rf_xlength(sexp) != 1 {
        return Err(format!("`{name}` must be a scalar logical"));
    }
    let val = *LOGICAL_RO(sexp);
    if val == i32::MIN {
        return Err(format!("`{name}` must not be NA"));
    }
    Ok(val != 0)
}

pub unsafe fn as_f64(sexp: SEXP, name: &str) -> Result<f64, String> {
    if TYPEOF(sexp) != REALSXP as c_int || Rf_xlength(sexp) != 1 {
        return Err(format!("`{name}` must be a scalar numeric value"));
    }
    Ok(*REAL_RO(sexp))
}

pub unsafe fn as_opt_f64(sexp: SEXP, name: &str) -> Result<Option<f64>, String> {
    if is_null(sexp) {
        Ok(None)
    } else {
        Ok(Some(as_f64(sexp, name)?))
    }
}

pub unsafe fn as_i32(sexp: SEXP, name: &str) -> Result<i32, String> {
    if TYPEOF(sexp) != INTSXP as c_int || Rf_xlength(sexp) != 1 {
        return Err(format!("`{name}` must be a scalar integer"));
    }
    let val = *INTEGER_RO(sexp);
    if val == i32::MIN {
        return Err(format!("`{name}` must not be NA"));
    }
    Ok(val)
}

pub unsafe fn as_opt_i32(sexp: SEXP, name: &str) -> Result<Option<i32>, String> {
    if is_null(sexp) {
        Ok(None)
    } else {
        Ok(Some(as_i32(sexp, name)?))
    }
}

pub unsafe fn as_str<'a>(sexp: SEXP, name: &str) -> Result<&'a str, String> {
    if TYPEOF(sexp) != STRSXP as c_int || Rf_xlength(sexp) != 1 {
        return Err(format!("`{name}` must be a single character string"));
    }
    let elt = STRING_ELT(sexp, 0);
    let cptr = R_CHAR(elt);
    if cptr.is_null() {
        return Err(format!("`{name}` string pointer is null"));
    }
    CStr::from_ptr(cptr)
        .to_str()
        .map_err(|_| format!("`{name}` is not valid UTF-8"))
}

pub unsafe fn as_opt_function(sexp: SEXP, name: &str) -> Result<Option<SEXP>, String> {
    if is_null(sexp) {
        return Ok(None);
    }
    if Rf_isFunction(sexp) == 0 {
        return Err(format!("`{name}` must be a function"));
    }
    Ok(Some(sexp))
}

pub unsafe fn call_residualizer(
    func: SEXP,
    row_major: &[f64],
    nrows: usize,
    ncols: usize,
) -> Result<Vec<f64>, String> {
    let mut guard = ProtectGuard::new();
    let mat = guard.protect(make_real_matrix_col_major_from_row_major(
        row_major, nrows, ncols,
    )?);
    let call = guard.protect(Rf_lang2(func, mat));
    let mut err_flag: c_int = 0;
    let res = R_tryEvalSilent(call, sp_r_global_env(), &mut err_flag);
    if err_flag != 0 {
        return Err("host residualizer evaluation failed".to_string());
    }
    let res = guard.protect(res);
    let view = as_f64_matrix(res, "within directions")?;
    if view.nrows != nrows || view.ncols != ncols {
        return Err("host residualizer changed the dimensions of its input".to_string());
    }
    crate::col_major_to_row_major(view.data, nrows, ncols)
}

unsafe fn mk_charsxp(s: &str) -> Result<SEXP, String> {
    let len = c_int::try_from(s.len()).map_err(|_| "string length exceeds R limit".to_string())?;
    Ok(Rf_mkCharLenCE(s.as_ptr().cast::<c_char>(), len, CE_UTF8))
}

pub unsafe fn make_scalar_str(s: &str) -> Result<SEXP, String> {
    let mut guard = ProtectGuard::new();
    let ch = guard.protect(mk_charsxp(s)?);
    Ok(Rf_ScalarString(ch))
}

pub unsafe fn make_scalar_real(v: f64) -> SEXP {
    Rf_ScalarReal(v)
}

pub unsafe fn make_scalar_int(v: i32) -> SEXP {
    Rf_ScalarInteger(v)
}

pub unsafe fn make_scalar_bool(v: bool) -> SEXP {
    Rf_ScalarLogical(c_int::from(v))
}

pub unsafe fn make_real_vec(values: &[f64]) -> Result<SEXP, String> {
    let len =
        R_xlen_t::try_from(values.len()).map_err(|_| "vector length exceeds R limit".to_string())?;
    let mut guard = ProtectGuard::new();
    let vec = guard.protect(Rf_allocVector(REALSXP, len));
    if !values.is_empty() {
        slice::from_raw_parts_mut(REAL(vec), values.len()).copy_from_slice(values);
    }
    Ok(vec)
}

pub unsafe fn make_real_matrix_col_major_from_row_major(
    row_major: &[f64],
    nrows: usize,
    ncols: usize,
) -> Result<SEXP, String> {
    let expected = nrows
        .checked_mul(ncols)
        .ok_or_else(|| "matrix dimensions overflow usize".to_string())?;
    if row_major.len() != expected {
        return Err(format!(
            "row-major buffer length {} does not match {} x {}",
            row_major.len(),
            nrows,
            ncols
        ));
    }
    let nr = c_int::try_from(nrows)
        .map_err(|_| "matrix row count exceeds 32-bit R integer".to_string())?;
    let nc = c_int::try_from(ncols)
        .map_err(|_| "matrix column count exceeds 32-bit R integer".to_string())?;
    let mut guard = ProtectGuard::new();
    let mat = guard.protect(Rf_allocMatrix(REALSXP, nr, nc));
    if expected > 0 {
        let dst = slice::from_raw_parts_mut(REAL(mat), expected);
        for col in 0..ncols {
            let dst_col = &mut dst[col * nrows..(col + 1) * nrows];
            for row in 0..nrows {
                dst_col[row] = row_major[row * ncols + col];
            }
        }
    }
    Ok(mat)
}

pub unsafe fn make_str_vec(items: &[&str]) -> Result<SEXP, String> {
    let len =
        R_xlen_t::try_from(items.len()).map_err(|_| "vector length exceeds R limit".to_string())?;
    let mut guard = ProtectGuard::new();
    let vec = guard.protect(Rf_allocVector(STRSXP, len));
    for (idx, item) in items.iter().enumerate() {
        let ch = guard.protect(mk_charsxp(item)?);
        SET_STRING_ELT(vec, idx as R_xlen_t, ch);
    }
    Ok(vec)
}

pub unsafe fn set_class(obj: SEXP, classes: &[&str]) -> Result<(), String> {
    let mut guard = ProtectGuard::new();
    let cls_vec = guard.protect(make_str_vec(classes)?);
    Rf_setAttrib(obj, sp_r_class_sym(), cls_vec);
    Ok(())
}

/// Build an R named list. Caller is responsible for protecting elements passed inside `entries` if
/// constructed prior to this call, or better yet build entries inside a single `ProtectGuard` before assembling.
pub unsafe fn make_named_list(entries: &[(&str, SEXP)]) -> Result<SEXP, String> {
    let len = R_xlen_t::try_from(entries.len())
        .map_err(|_| "list length exceeds R limit".to_string())?;
    let mut guard = ProtectGuard::new();
    let list = guard.protect(Rf_allocVector(VECSXP, len));
    let names = guard.protect(Rf_allocVector(STRSXP, len));
    for (idx, &(key, val)) in entries.iter().enumerate() {
        SET_VECTOR_ELT(list, idx as R_xlen_t, val);
        let ch = guard.protect(mk_charsxp(key)?);
        SET_STRING_ELT(names, idx as R_xlen_t, ch);
    }
    Rf_setAttrib(list, sp_r_names_sym(), names);
    Ok(list)
}

pub unsafe fn make_error_obj(message: &str, class: Option<&str>) -> SEXP {
    let mut guard = ProtectGuard::new();
    let msg_sexp = match make_scalar_str(message) {
        Ok(s) => guard.protect(s),
        Err(_) => r_nil(),
    };
    let cls_sexp = match class {
        Some(c) => match make_scalar_str(c) {
            Ok(s) => guard.protect(s),
            Err(_) => r_nil(),
        },
        None => r_nil(),
    };
    let Ok(obj) = make_named_list(&[("message", msg_sexp), ("class", cls_sexp)]) else {
        return r_nil();
    };
    let obj = guard.protect(obj);
    let _ = set_class(obj, &["spatial_native_error"]);
    obj
}
