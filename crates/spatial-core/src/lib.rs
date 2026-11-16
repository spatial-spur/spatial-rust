//! Shared numerical core for spatial SCPC and SPUR.
//!
//! Matrices at this boundary are contiguous and row-major.

pub mod basis;
pub mod calibration;
mod covariance;
pub mod distance;
mod linalg;
mod quadrature;
pub mod scpc_op;
mod student_t;
pub mod threads;

pub use distance::Metric;
pub use scpc_op::{run_exact_scpc, ExactScpcOutput, ExactScpcRequest};

use std::fmt;

/// A numerical input or calculation error, translated by each host adapter.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct CoreError(pub &'static str);

impl fmt::Display for CoreError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.0.fmt(f)
    }
}

impl std::error::Error for CoreError {}

pub(crate) fn checked_len(dimensions: &[usize]) -> Result<usize, CoreError> {
    dimensions
        .iter()
        .try_fold(1usize, |size, dimension| size.checked_mul(*dimension))
        .filter(|&size| size <= isize::MAX as usize / size_of::<f64>())
        .ok_or(CoreError("matrix dimensions are too large"))
}
