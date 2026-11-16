//! Released SCPC decay calibration (`getc0fromavc`).

use crate::{linalg, CoreError};
use rayon::prelude::*;

const CHUNK_SIZE: usize = 32_768;

/// Evaluate average pairwise correlation `mean(exp(-c * dist))`.
pub fn avg_correlation(c: f64, distances: &[f64]) -> f64 {
    if distances.is_empty() {
        return f64::NAN;
    }
    let sum = if distances.len() >= 4 * CHUNK_SIZE {
        distances
            .par_chunks(CHUNK_SIZE)
            .map_init(Vec::new, |scratch, chunk| {
                scratch.resize(chunk.len(), 0.0);
                linalg::exp_scaled(chunk, -c, scratch);
                scratch.iter().sum::<f64>()
            })
            .sum::<f64>()
    } else {
        let mut scratch = vec![0.0; distances.len()];
        linalg::exp_scaled(distances, -c, &mut scratch);
        scratch.iter().sum::<f64>()
    };
    sum / distances.len() as f64
}

/// Solve for exponential decay parameter `c` matching `target_avc` using
/// the exact geometric-midpoint bracketing and bisection of released `main`.
pub fn solve_c0_from_avc(distances: &[f64], target_avc: f64) -> Result<f64, CoreError> {
    if distances.is_empty() {
        return Err(CoreError("at least two locations are required"));
    }
    if !target_avc.is_finite() || target_avc <= 0.0 || target_avc >= 1.0 {
        return Err(CoreError(
            "target average correlation must lie strictly in (0, 1)",
        ));
    }
    if distances.iter().any(|&d| !d.is_finite() || d < 0.0) {
        return Err(CoreError("distances must be finite and nonnegative"));
    }

    let mut c0 = 10.0_f64;
    let mut c1 = 10.0_f64;
    while avg_correlation(c0, distances) < target_avc {
        c1 = c0;
        c0 *= 0.5;
        if c0 < 1e-300 {
            return Err(CoreError(
                "decay calibration lower bracket underflowed; distances may be degenerate",
            ));
        }
    }
    while avg_correlation(c1, distances) > target_avc && c1 < 5000.0 {
        c0 = c1;
        c1 *= 2.0;
    }

    let mut iter = 0usize;
    loop {
        let c = (c0 * c1).sqrt();
        if c == c0 || c == c1 {
            return Ok(c);
        }
        if avg_correlation(c, distances) > target_avc {
            c0 = c;
        } else {
            c1 = c;
        }
        iter += 1;
        if c1 - c0 <= 0.001 || iter >= 1024 {
            return Ok(c);
        }
    }
}
