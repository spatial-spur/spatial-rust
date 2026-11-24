//! Native bindings for the spatialSpur R distribution using direct official R C-API entrypoints.

#[cfg(not(test))]
mod r_ffi;
#[cfg(not(test))]
mod exports;

use spatial_core::{CoreError, ErrorKind};

fn core_err(err: CoreError) -> String {
    match err.kind() {
        ErrorKind::InvalidInput | ErrorKind::TiedEigenspace => err.message().to_string(),
        ErrorKind::UnsupportedCase => format!("unsupported case: {}", err.message()),
        ErrorKind::CalibrationFailed => format!("calibration failed: {}", err.message()),
        ErrorKind::DegenerateVariance => format!("degenerate variance: {}", err.message()),
        ErrorKind::InternalFailure => format!("internal failure: {}", err.message()),
    }
}

/// Native error payload converted into a `spatial_native_error` R object at the FFI boundary.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ScpcError {
    pub message: String,
    pub class: Option<&'static str>,
}

impl From<String> for ScpcError {
    fn from(message: String) -> Self {
        Self {
            message,
            class: None,
        }
    }
}

impl From<&str> for ScpcError {
    fn from(message: &str) -> Self {
        Self {
            message: message.to_string(),
            class: None,
        }
    }
}

impl From<CoreError> for ScpcError {
    fn from(error: CoreError) -> Self {
        let class = match error.kind() {
            ErrorKind::TiedEigenspace => Some("spatial_tied_eigenspace"),
            _ => None,
        };
        Self {
            message: core_err(error),
            class,
        }
    }
}

/// Convert an R column-major matrix buffer (`nrows x ncols`) into a contiguous row-major `Vec<f64>`.
pub fn col_major_to_row_major(data: &[f64], nrows: usize, ncols: usize) -> Result<Vec<f64>, String> {
    let expected = nrows
        .checked_mul(ncols)
        .ok_or_else(|| "matrix dimensions overflow usize".to_string())?;
    if data.len() != expected {
        return Err(format!(
            "column-major buffer length {} does not match {} x {}",
            data.len(),
            nrows,
            ncols
        ));
    }
    let mut row_major = vec![0.0; expected];
    for col in 0..ncols {
        let src_col = &data[col * nrows..(col + 1) * nrows];
        for row in 0..nrows {
            row_major[row * ncols + col] = src_col[row];
        }
    }
    Ok(row_major)
}

fn parse_r_nthreads(nthreads: Option<i32>) -> Result<Option<usize>, String> {
    match nthreads {
        None => Ok(None),
        Some(t) => {
            if t == i32::MIN || t <= 0 {
                return Err("`nthreads` must be a positive integer.".to_string());
            }
            let count = Some(t as usize);
            spatial_core::threads::validate_nthreads(count).map_err(core_err)?;
            Ok(count)
        }
    }
}

fn parse_r_index_slice(
    ids: Option<&[i32]>,
    upper_bound: usize,
    name: &str,
    label: &str,
) -> Result<Option<Vec<usize>>, String> {
    let Some(slice) = ids else {
        return Ok(None);
    };
    let mut out = Vec::with_capacity(slice.len());
    for &val in slice {
        if val == i32::MIN || val < 0 || (val as usize) >= upper_bound {
            return Err(format!(
                "`{name}` entries must be non-NA integers in 0..{} for {label} (got {val})",
                upper_bound
            ));
        }
        out.push(val as usize);
    }
    Ok(Some(out))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn col_major_to_row_major_transposes_correctly() {
        let col_major = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        let row_major = col_major_to_row_major(&col_major, 2, 3).unwrap();
        assert_eq!(row_major, vec![1.0, 3.0, 5.0, 2.0, 4.0, 6.0]);
    }

    #[test]
    fn col_major_to_row_major_rejects_mismatched_length() {
        assert!(col_major_to_row_major(&[1.0, 2.0], 2, 2).is_err());
    }

    #[test]
    fn core_err_maps_all_error_kinds() {
        let cases = [
            (CoreError::invalid_input("bad arg"), "bad arg"),
            (CoreError::unsupported("approx"), "unsupported case: approx"),
            (CoreError::tied("ambiguous"), "ambiguous"),
            (
                CoreError::calibration("no root"),
                "calibration failed: no root",
            ),
            (
                CoreError::degenerate("zero var"),
                "degenerate variance: zero var",
            ),
            (CoreError::internal("lapack"), "internal failure: lapack"),
        ];
        for (err, expected) in cases {
            assert_eq!(core_err(err), expected);
        }
    }

    #[test]
    fn approx_scpc_error_sets_expected_classes_and_messages() {
        let tied = ScpcError::from(CoreError::tied("tied basis"));
        assert_eq!(tied.class, Some("spatial_tied_eigenspace"));
        assert_eq!(tied.message, "tied basis");

        let other_core = ScpcError::from(CoreError::degenerate("zero var"));
        assert_eq!(other_core.class, None);
        assert_eq!(other_core.message, "degenerate variance: zero var");

        let direct = ScpcError::from("direct msg");
        assert_eq!(direct.class, None);
        assert_eq!(direct.message, "direct msg");
    }

    #[test]
    fn parse_r_index_slice_rejects_na_and_out_of_bounds() {
        assert!(parse_r_index_slice(Some(&[i32::MIN]), 5, "cluster_id", "clusters").is_err());
        assert!(parse_r_index_slice(Some(&[-1]), 5, "cluster_id", "clusters").is_err());
        assert!(parse_r_index_slice(Some(&[5]), 5, "cluster_id", "clusters").is_err());
        assert_eq!(
            parse_r_index_slice(Some(&[0, 4]), 5, "cluster_id", "clusters").unwrap(),
            Some(vec![0, 4])
        );
        assert!(parse_r_nthreads(Some(0)).is_err());
        assert!(parse_r_nthreads(Some(i32::MIN)).is_err());
        assert_eq!(parse_r_nthreads(None).unwrap(), None);
        assert_eq!(parse_r_nthreads(Some(1)).unwrap(), Some(1));
    }
}
