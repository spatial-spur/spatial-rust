//! Landmark Nyström basis construction for large-n SCPC.

use crate::distance::{Distances, Metric};
use crate::{basis, calibration, checked_len, linalg, CoreError};
use rayon::prelude::*;

const ROW_TILE: usize = 128;

#[derive(Clone, Copy)]
struct Lcg(u32);

impl Lcg {
    #[inline]
    fn next_u(&mut self) -> f64 {
        self.0 = self.0.wrapping_mul(64389).wrapping_add(1);
        f64::from(self.0) / 4_294_967_296.0
    }

    #[inline]
    fn index(&mut self, count: usize) -> usize {
        if count <= 1 {
            self.0 = self.0.wrapping_mul(64389).wrapping_add(1);
            return 0;
        }
        let raw = (self.next_u() * count as f64) as usize;
        raw.min(count - 1)
    }
}

/// Put sites in canonical order before seeded landmark sampling. For Euclidean
/// inputs, a common translation and scale condition the arithmetic without
/// changing distance ratios. The return value also records the scale divisor.
pub fn normalize_coords(
    coords: &[f64],
    n: usize,
    d: usize,
    metric: Metric,
) -> Result<(Vec<f64>, Vec<usize>, f64), CoreError> {
    let distance = Distances::new(coords, n, d, metric)?;
    let (mut values, scale) = if let Some(canonical) = distance.canonical_coordinates() {
        (canonical, 1.0)
    } else {
        let mut bounds = vec![(f64::INFINITY, f64::NEG_INFINITY); d];
        for row in coords.chunks_exact(d) {
            for (&value, (min, max)) in row.iter().zip(&mut bounds) {
                *min = min.min(value);
                *max = max.max(value);
            }
        }
        let center: Vec<_> = bounds.iter().map(|(a, b)| a.midpoint(*b)).collect();
        let scaled: Vec<_> = coords
            .chunks_exact(d)
            .flat_map(|row| row.iter().zip(&center).map(|(x, c)| x - c))
            .collect();
        let divisor = scaled.iter().copied().map(f64::abs).fold(0.0, f64::max);
        (scaled, if divisor == 0.0 { 1.0 } else { divisor })
    };
    for value in &mut values {
        *value = if *value == 0.0 { 0.0 } else { *value / scale };
    }
    let order = distance.coordinate_order();
    let sorted = order
        .iter()
        .flat_map(|&i| values[i * d..(i + 1) * d].iter().copied())
        .collect();
    Ok((sorted, order, scale))
}

fn landmarks(n: usize, size: usize, samples: usize, seed: u32) -> Vec<Vec<usize>> {
    let mut rng = Lcg(seed);
    let mut order: Vec<usize> = (0..n).collect();
    (0..samples)
        .map(|_| {
            for i in 0..size {
                let j = i + rng.index(n - i);
                order.swap(i, j);
            }
            order[..size].to_vec()
        })
        .collect()
}

/// Lower and upper SCPC decay scales with co-location floor adjustment.
pub fn decay_limits(distances: &[f64], avc: f64, min_avc: f64) -> Result<[f64; 2], CoreError> {
    let floor = calibration::correlation_floor(distances)?;
    if avc <= floor || !avc.is_finite() || avc >= 1.0 {
        return Err(CoreError::calibration("unattainable sampled SCPC decay").with_message(format!(
            "Target correlation avc={avc} must exceed the co-located-pair fraction {floor:.6}."
        )));
    }
    let min_avc_eff = floor + ((1.0 - floor) * min_avc).min(0.1 * (avc - floor));
    let c0 = calibration::solve_c0_from_avc(distances, avc)?;
    let cmax = calibration::solve_c0_from_avc(distances, min_avc_eff)?;
    Ok([c0, cmax])
}

fn streamed_decay_limits(
    distance: &Distances<'_>,
    n: usize,
    avc: f64,
    min_avc: f64,
) -> Result<[f64; 2], CoreError> {
    let summarize = |i: usize| {
        let mut zeros = 0usize;
        let mut sum_pos = 0.0_f64;
        let mut max_d = 0.0_f64;
        for j in 0..i {
            let d = distance.pair(i, j);
            if !d.is_finite() || d < 0.0 {
                return Err(CoreError::invalid_input(
                    "distances must be finite and nonnegative",
                ));
            }
            if d == 0.0 {
                zeros += 1;
            } else {
                sum_pos += d;
                max_d = max_d.max(d);
            }
        }
        Ok((zeros, sum_pos, max_d))
    };
    let row_stats: Vec<(usize, f64, f64)> = if !crate::threads::is_serial() && n >= 64 {
        (1..n)
            .into_par_iter()
            .map(summarize)
            .collect::<Result<Vec<_>, _>>()?
    } else {
        (1..n).map(summarize).collect::<Result<Vec<_>, _>>()?
    };
    let (zeros, sum_pos, max_d) =
        row_stats
            .into_iter()
            .fold((0usize, 0.0_f64, 0.0_f64), |acc, item| {
                (acc.0 + item.0, acc.1 + item.1, acc.2.max(item.2))
            });
    let pairs = n * (n - 1) / 2;
    let floor = zeros as f64 / pairs as f64;
    if !avc.is_finite() || avc <= floor || avc >= 1.0 {
        return Err(CoreError::calibration("unattainable SCPC decay").with_message(format!(
            "Target correlation avc={avc} must exceed the co-located-pair fraction {floor:.6}."
        )));
    }
    let min_avc_eff = floor + ((1.0 - floor) * min_avc).min(0.1 * (avc - floor));
    let mean_pos = sum_pos / (pairs - zeros) as f64;
    let solve_target = |target: f64| {
        let near_one = target > 0.9999;
        calibration::solve_decay(target, floor, max_d, mean_pos, |rate, _| {
            let eval_row = |i: usize| {
                let row = (0..i).map(|j| distance.pair(i, j)).collect::<Vec<_>>();
                calibration::eval_residual_and_slope_chunk(&row, rate, near_one, &mut Vec::new())
            };
            let sums: Vec<_> = if !crate::threads::is_serial() && n >= 64 {
                (1..n).into_par_iter().map(eval_row).collect()
            } else {
                (1..n).map(eval_row).collect()
            };
            let (mass, derivative) = sums
                .into_iter()
                .fold((0.0_f64, 0.0_f64), |acc, (a, b)| (acc.0 + a, acc.1 + b));
            let count = pairs as f64;
            let residual = if near_one {
                (1.0 - target) + mass / count
            } else {
                mass / count - target
            };
            (residual, -derivative / count)
        })
    };
    let c0 = solve_target(avc)?;
    let cmax = solve_target(min_avc_eff)?;
    Ok([c0, cmax])
}

/// Controls the landmark Nyström approximation.
pub struct LandmarkRequest {
    pub avc: f64,
    pub min_avc: f64,
    pub components: usize,
    pub landmarks: usize,
    pub samples: usize,
    pub seed: u32,
}

/// Row-major Nyström basis (`n x (q + 1)`), and decay limits `[c0, cmax]`.
pub struct LandmarkBasis {
    pub weights: Vec<f64>,
    pub decay: [f64; 2],
    pub candidates: Vec<usize>,
    pub tied: bool,
    pub cutoff_tied: bool,
    pub more_positive: bool,
}

/// Embed complete local eigenspaces, then take the leading subspace of the
/// pooled projections. Orthogonal rotations of a retained local block leave
/// the pooled covariance invariant.
pub fn landmark_basis(
    coords: &[f64],
    n: usize,
    width: usize,
    metric: Metric,
    request: &LandmarkRequest,
) -> Result<LandmarkBasis, CoreError> {
    landmark_basis_with_shortcut(coords, n, width, metric, request, true)
}

fn landmark_basis_with_shortcut(
    coords: &[f64],
    n: usize,
    width: usize,
    metric: Metric,
    request: &LandmarkRequest,
    allow_exact_shortcut: bool,
) -> Result<LandmarkBasis, CoreError> {
    let distance = Distances::new(coords, n, width, metric)?;
    let q = request.components;
    let m = request.landmarks.min(n);
    if q == 0 || q >= m || request.samples == 0 || n > i32::MAX as usize {
        return Err(CoreError::invalid_input("invalid SCPC landmark dimensions"));
    }
    let pairs = checked_len(&[m, m - 1])? / 2;
    let subsets = if m == n {
        vec![(0..n).collect()]
    } else {
        landmarks(n, m, request.samples, request.seed)
    };
    let samples = subsets.len();

    let mut distances = vec![0.0; checked_len(&[samples, pairs])?];
    if crate::threads::is_serial() || samples == 1 {
        for (block, subset) in distances.chunks_exact_mut(pairs).zip(&subsets) {
            let mut k = 0;
            for j in 0..m - 1 {
                for i in j + 1..m {
                    block[k] = distance.pair(subset[i], subset[j]);
                    k += 1;
                }
            }
        }
    } else {
        distances
            .par_chunks_exact_mut(pairs)
            .zip(&subsets)
            .for_each(|(block, subset)| {
                let mut k = 0;
                for j in 0..m - 1 {
                    for i in j + 1..m {
                        block[k] = distance.pair(subset[i], subset[j]);
                        k += 1;
                    }
                }
            });
    }

    let decay = streamed_decay_limits(&distance, n, request.avc, request.min_avc)?;

    if allow_exact_shortcut && m == n {
        let dense = crate::distance::distance_matrix(coords, n, width, metric)?;
        let exact = basis::scpc_basis(&dense, n, decay[0], q)?;
        return Ok(LandmarkBasis {
            weights: exact.weights,
            decay,
            candidates: exact.candidates,
            tied: exact.tied,
            cutoff_tied: exact.cutoff_tied,
            more_positive: exact.more_positive,
        });
    }

    let extend_one =
        |(triangle, subset): (&[f64], &Vec<usize>)| -> Result<Option<(Vec<f64>, usize, bool)>, CoreError> {
            let mut dense = crate::zeroed(&[m, m], "landmark distances")?;
            let mut k = 0;
            for j in 0..m - 1 {
                for i in j + 1..m {
                    dense[i * m + j] = triangle[k];
                    dense[j * m + i] = triangle[k];
                    k += 1;
                }
            }
            let mut requested = q;
            let (local, retained) = loop {
                let local = match basis::scpc_basis(&dense, m, decay[0], requested) {
                    Ok(basis) => basis,
                    Err(err) if err.kind() == crate::ErrorKind::DegenerateVariance => {
                        return Ok(None);
                    }
                    Err(err) => return Err(err),
                };
                if local.cutoff_tied && requested < m - 1 {
                    requested = (requested + requested.div_ceil(2)).min(m - 1);
                    continue;
                }
                let Some(&retained) = local.candidates.last() else {
                    if local.cutoff_tied {
                        return Err(CoreError::degenerate(
                            "landmark sample splits the leading eigenspace without an interior cutoff",
                        ));
                    }
                    return Ok(None);
                };
                break (local, retained);
            };
            let mut inverse_roots = Vec::with_capacity(retained);
            let mut start = 0;
            for &end in &local.candidates {
                // A common eigenvalue scale leaves a numerically tied block
                // invariant to the provider's orientation of its eigenvectors.
                let mean = local.eigenvalues[start..end].iter().sum::<f64>() / (end - start) as f64;
                inverse_roots.resize(end, mean.sqrt().recip());
                start = end;
            }
            let weights = local
                .weights
                .chunks_exact(requested + 1)
                .flat_map(|row| {
                    row[1..=retained]
                        .iter()
                        .zip(&inverse_roots)
                        .map(|(v, scale)| v * scale)
                })
                .collect::<Vec<_>>();
            let projected = project_landmarks(&distance, n, subset, &weights, retained, decay[0])?;
            Ok(Some((projected, retained, local.more_positive)))
        };

    let raw_projected: Vec<Option<(Vec<f64>, usize, bool)>> =
        if crate::threads::is_serial() || samples == 1 {
            distances
                .chunks_exact(pairs)
                .zip(&subsets)
                .map(extend_one)
                .collect::<Result<Vec<_>, _>>()?
        } else {
            distances
                .par_chunks_exact(pairs)
                .zip(&subsets)
                .map(extend_one)
                .collect::<Result<Vec<_>, _>>()?
        };
    drop(distances);

    let projected: Vec<_> = raw_projected.into_iter().flatten().collect();
    if projected.is_empty() {
        return Err(CoreError::degenerate(
            "landmark samples have no positive spatial directions",
        ));
    }
    let any_local_more_positive = projected.iter().any(|(_, _, more)| *more);
    let total_columns = projected.iter().map(|(_, width, _)| width).sum();
    let mut wall = crate::zeroed(&[n, total_columns], "embedded landmark directions")?;
    let mut offset = 0;
    for (projection, width, _) in projected {
        for (destination, source) in wall
            .chunks_exact_mut(total_columns)
            .zip(projection.chunks_exact(width))
        {
            destination[offset..offset + width].copy_from_slice(source);
        }
        offset += width;
    }
    let (weights, blocks) = combine_landmarks(&wall, n, q, total_columns, any_local_more_positive)?;
    Ok(LandmarkBasis {
        weights,
        decay,
        candidates: blocks.ends,
        tied: blocks.tied,
        cutoff_tied: blocks.cutoff_tied,
        more_positive: blocks.more_positive,
    })
}

fn project_landmarks(
    distance: &Distances<'_>,
    n: usize,
    subset: &[usize],
    weights: &[f64],
    q: usize,
    decay: f64,
) -> Result<Vec<f64>, CoreError> {
    let m = subset.len();
    let mut projected = vec![0.0; n * q];
    let mut d = vec![0.0; ROW_TILE.min(n) * m];
    let mut kernel = d.clone();
    for (block_index, output) in projected.chunks_mut(ROW_TILE * q).enumerate() {
        let rows = output.len() / q;
        for (row, distances) in d[..rows * m].chunks_exact_mut(m).enumerate() {
            let i = block_index * ROW_TILE + row;
            for (value, &j) in distances.iter_mut().zip(subset) {
                let dij = distance.pair(i, j);
                if !dij.is_finite() {
                    return Err(CoreError::invalid_input(
                        "non-finite pairwise distance in Nystrom landmark projection",
                    ));
                }
                *value = dij;
            }
        }
        linalg::exp_scaled(&d[..rows * m], -decay, &mut kernel[..rows * m]);
        linalg::multiply(&kernel, weights, output, rows, m, q, false);
    }
    // Center and use one Frobenius scale per landmark sample. Normalizing
    // columns individually would assign different weights after a block rotation.
    let mut scale_squared = 0.0;
    for j in 0..q {
        let mean = projected.chunks_exact(q).map(|row| row[j]).sum::<f64>() / n as f64;
        for row in projected.chunks_exact_mut(q) {
            row[j] -= mean;
            scale_squared += row[j] * row[j];
        }
    }
    if scale_squared > 0.0 && scale_squared.is_finite() {
        let scale = scale_squared.sqrt().recip();
        for value in &mut projected {
            *value *= scale;
        }
    }
    Ok(projected)
}

fn combine_landmarks(
    wall: &[f64],
    n: usize,
    q: usize,
    columns: usize,
    any_local_more_positive: bool,
) -> Result<(Vec<f64>, basis::PositiveBlocks), CoreError> {
    let mut gram = crate::zeroed(&[columns, columns], "landmark Gram matrix")?;
    linalg::multiply(wall, wall, &mut gram, columns, n, columns, true);
    let retained = q.min(columns);
    let eigen = basis::symmetric_eigenpairs(&gram, columns, retained)?;
    let mut blocks = basis::positive_blocks(&eigen.values, eigen.next_value, 1.6e-11);
    let rank = blocks.ends.last().copied().unwrap_or(0);
    let mut weights = vec![(n as f64).sqrt().recip(); checked_len(&[n, q + 1])?];
    for row in weights.chunks_exact_mut(q + 1) {
        row[1..].fill(0.0);
    }
    if rank > 0 {
        // The working matrix is n × columns. The right singular vectors
        // arrive one vector per row and must become columns before W V.
        let mut right = crate::zeroed(&[columns, rank], "right landmark directions")?;
        for (j, vector) in eigen.vectors[..rank * columns]
            .chunks_exact(columns)
            .enumerate()
        {
            let inverse = eigen.values[j].sqrt().recip();
            for (row, value) in right.chunks_exact_mut(rank).zip(vector) {
                row[j] = value * inverse;
            }
        }
        let mut directions = crate::zeroed(&[n, rank], "pooled landmark basis")?;
        linalg::multiply(wall, &right, &mut directions, n, columns, rank, false);

        // Forming W'W squares the condition number; re-center and run two-pass
        // Modified Gram-Schmidt to restore centering and orthogonality.
        let mut ortho_cols: Vec<Vec<f64>> = (0..rank)
            .map(|c| {
                let mut col: Vec<f64> = directions.chunks_exact(rank).map(|r| r[c]).collect();
                let mean = col.iter().sum::<f64>() / n as f64;
                for v in &mut col {
                    *v -= mean;
                }
                col
            })
            .collect();
        let mut kept = 0usize;
        for c in 0..rank {
            let mut col = std::mem::take(&mut ortho_cols[c]);
            for _ in 0..2 {
                let mean = col.iter().sum::<f64>() / n as f64;
                for v in &mut col {
                    *v -= mean;
                }
                for prev in &ortho_cols[..kept] {
                    let dot: f64 = prev.iter().zip(&col).map(|(a, b)| a * b).sum();
                    for (v, &p) in col.iter_mut().zip(prev) {
                        *v -= dot * p;
                    }
                }
            }
            let norm = col.iter().map(|v| v * v).sum::<f64>().sqrt();
            if !norm.is_finite() || norm <= 1e-12 {
                break;
            }
            for v in &mut col {
                *v /= norm;
            }
            ortho_cols[kept] = col;
            kept += 1;
        }
        let mut final_rank = blocks.ends.last().copied().unwrap_or(0);
        if kept < rank {
            blocks.ends.retain(|&end| end <= kept);
            blocks.cutoff_tied = blocks.ends.is_empty();
            blocks.more_positive = false;
            final_rank = blocks.ends.last().copied().unwrap_or(0);
        } else if final_rank == q && any_local_more_positive {
            blocks.more_positive = true;
        }
        for (i, row) in weights.chunks_exact_mut(q + 1).enumerate() {
            for (c, col) in ortho_cols[..final_rank].iter().enumerate() {
                row[1 + c] = col[i];
            }
        }
    }
    Ok((weights, blocks))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::distance::distance_matrix;

    #[test]
    fn landmark_sampling_is_reproducible_and_without_replacement() {
        let subsets = landmarks(16, 8, 3, 42);
        assert_eq!(subsets, landmarks(16, 8, 3, 42));
        assert_ne!(subsets, landmarks(16, 8, 3, 43));
        for subset in subsets {
            assert_eq!(
                subset
                    .iter()
                    .copied()
                    .collect::<std::collections::HashSet<_>>()
                    .len(),
                8
            );
        }
    }

    #[test]
    fn pooled_landmarks_recover_left_singular_directions() {
        // Centered, orthogonal embedded columns with singular values sqrt(2) and sqrt(18).
        let wall = [1., 0., -1., 0., 0., 3., 0., -3., 0., 0.];
        let (out, blocks) = combine_landmarks(&wall, 5, 2, 2, false).unwrap();
        assert_eq!(blocks.ends, [1, 2]);
        let large = [0., 0., 1. / 2.0_f64.sqrt(), -1. / 2.0_f64.sqrt(), 0.];
        let small = [1. / 2.0_f64.sqrt(), -1. / 2.0_f64.sqrt(), 0., 0., 0.];
        for (column, expected) in [(1, large), (2, small)] {
            let dot: f64 = out
                .chunks_exact(3)
                .zip(expected)
                .map(|(r, e)| r[column] * e)
                .sum();
            assert!((dot.abs() - 1.0).abs() < 1e-14);
        }
    }

    #[test]
    fn full_geographic_landmark_extension_spans_exact_subspace() {
        let n = 36;
        let mut coords = Vec::with_capacity(n * 2);
        for i in 0..n {
            let fi = i as f64;
            coords.push(-170.0 + (fi * 37.3).rem_euclid(340.0));
            coords.push(-65.0 + (fi * 19.7).rem_euclid(130.0));
        }
        let lb = landmark_basis_with_shortcut(
            &coords,
            n,
            2,
            Metric::GreatCircleFractionLonLat,
            &LandmarkRequest {
                avc: 0.1,
                min_avc: 1e-4,
                components: 3,
                landmarks: n,
                samples: 1,
                seed: 7,
            },
            false,
        )
        .unwrap();
        let dm = distance_matrix(&coords, n, 2, Metric::GreatCircleFractionLonLat).unwrap();
        let exact = basis::scpc_basis(&dm, n, lb.decay[0], 3).unwrap();
        // Check column-space projector P_approx == P_exact on columns 1..=3
        let mut p_approx = vec![0.0; n * n];
        let mut p_exact = vec![0.0; n * n];
        let a_sub: Vec<f64> = lb
            .weights
            .chunks_exact(4)
            .flat_map(|r| r[1..].iter().copied())
            .collect();
        let e_sub: Vec<f64> = exact
            .weights
            .chunks_exact(4)
            .flat_map(|r| r[1..].iter().copied())
            .collect();
        let mut a_sub_t = vec![0.0; 3 * n];
        let mut e_sub_t = vec![0.0; 3 * n];
        for i in 0..n {
            for c in 0..3 {
                a_sub_t[c * n + i] = a_sub[i * 3 + c];
                e_sub_t[c * n + i] = e_sub[i * 3 + c];
            }
        }
        linalg::multiply(&a_sub, &a_sub_t, &mut p_approx, n, 3, n, false);
        linalg::multiply(&e_sub, &e_sub_t, &mut p_exact, n, 3, n, false);
        for (a, b) in p_approx.iter().zip(&p_exact) {
            assert!((a - b).abs() < 1e-7, "projector diff a={a}, b={b}");
        }
    }

    #[test]
    fn normalize_coords_is_scale_invariant_at_extreme_magnitudes() {
        let base = vec![-1.5, 0.25, 0.0, -0.75, 2.0, 1.5, -0.5, 2.25, 1.25, -1.0];
        let (norm_ref, perm_ref, _) = normalize_coords(&base, 5, 2, Metric::Euclidean).unwrap();
        for factor in [1e-200_f64, 1e200_f64] {
            let scaled: Vec<f64> = base.iter().map(|&v| v * factor).collect();
            let (norm, perm, _) = normalize_coords(&scaled, 5, 2, Metric::Euclidean).unwrap();
            assert_eq!(perm, perm_ref);
            for (a, b) in norm.iter().zip(&norm_ref) {
                assert!(a.is_finite() && (a - b).abs() < 1e-12, "a={a}, b={b}");
            }
        }
        let shifted: Vec<_> = base.iter().map(|&v| 1e308 + v * 1e306).collect();
        let (norm, perm, _) = normalize_coords(&shifted, 5, 2, Metric::Euclidean).unwrap();
        assert_eq!(perm, perm_ref);
        for (a, b) in norm.iter().zip(&norm_ref) {
            assert!((a - b).abs() < 1e-12, "a={a}, b={b}");
        }
    }

    #[test]
    fn colocated_landmark_sample_is_skipped_without_injecting_roundoff_noise() {
        // 8 co-located sites out of 9 (floor = 28/36 < 0.9), landmarks = 7.
        // All-colocated landmark draws must be skipped cleanly; informative samples give rank 1.
        let mut coords = vec![0.0_f64; 16];
        coords.extend_from_slice(&[1.0, 0.0]);
        let lb = landmark_basis(
            &coords,
            9,
            2,
            Metric::Euclidean,
            &LandmarkRequest {
                avc: 0.9,
                min_avc: 1e-4,
                components: 2,
                landmarks: 7,
                samples: 200,
                seed: 1,
            },
        )
        .unwrap();
        assert_eq!(lb.candidates, vec![1]);
        assert!(!lb.more_positive);

        // When all sites coincide globally, calibration rejects avc <= 1.0 cleanly.
        let all_same = vec![0.0_f64; 18];
        let err = landmark_basis(
            &all_same,
            9,
            2,
            Metric::Euclidean,
            &LandmarkRequest {
                avc: 0.9,
                min_avc: 1e-4,
                components: 2,
                landmarks: 7,
                samples: 10,
                seed: 1,
            },
        )
        .err()
        .expect("all-colocated global geometry must reject avc <= floor");
        assert_eq!(err.kind(), crate::ErrorKind::CalibrationFailed);
    }

    #[test]
    fn single_sample_nystrom_preserves_more_positive_and_matches_full_decay() {
        let n = 48;
        let coords: Vec<f64> = (0..n)
            .flat_map(|i| {
                let f = i as f64;
                [-105.0 + 0.31 * f + (f * 0.7).sin(), 30.0 + 0.19 * f]
            })
            .collect();
        let lb = landmark_basis(
            &coords,
            n,
            2,
            Metric::GreatCircleFractionLonLat,
            &LandmarkRequest {
                avc: 0.03,
                min_avc: 1e-5,
                components: 10,
                landmarks: 27,
                samples: 1,
                seed: 7,
            },
        )
        .unwrap();
        assert!(lb.more_positive);
        let dm = distance_matrix(&coords, n, 2, Metric::GreatCircleFractionLonLat).unwrap();
        let full_c0 = calibration::solve_c0_from_avc(&crate::distance::lower_triangle(&dm, n), 0.03).unwrap();
        assert!((lb.decay[0] - full_c0).abs() <= 1e-10 * full_c0);
    }

    #[test]
    fn geographic_nystrom_lonlat_and_latlon_orderings_agree_exactly() {
        let n = 36;
        let lonlat: Vec<f64> = (0..n)
            .flat_map(|i| {
                let f = i as f64;
                [-120.0 + (f * 23.0) % 45.0, 25.0 + (f * 13.0) % 30.0]
            })
            .collect();
        let (norm_lonlat, perm_lonlat, _) =
            normalize_coords(&lonlat, n, 2, Metric::GreatCircleFractionLonLat).unwrap();
        let latlon: Vec<f64> = lonlat
            .chunks_exact(2)
            .flat_map(|pair| [pair[1], pair[0]])
            .collect();
        let (norm_latlon, perm_latlon, _) =
            normalize_coords(&latlon, n, 2, Metric::GreatCircleFractionLatLon).unwrap();
        assert_eq!(perm_lonlat, perm_latlon);
        let req = LandmarkRequest {
            avc: 0.05,
            min_avc: 1e-5,
            components: 8,
            landmarks: 20,
            samples: 4,
            seed: 42,
        };
        let b1 = landmark_basis(&norm_lonlat, n, 2, Metric::GreatCircleFractionLonLat, &req).unwrap();
        let b2 = landmark_basis(&norm_latlon, n, 2, Metric::GreatCircleFractionLatLon, &req).unwrap();
        assert_eq!(b1.candidates, b2.candidates);
        assert_eq!(b1.decay, b2.decay);
        for (a, b) in b1.weights.iter().zip(&b2.weights) {
            assert!((a - b).abs() < 1e-12, "a={a}, b={b}");
        }
    }
}
