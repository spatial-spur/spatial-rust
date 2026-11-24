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
    let mut order: Vec<_> = (0..n).collect();
    order.sort_by(|&i, &j| {
        values[i * d..(i + 1) * d]
            .iter()
            .zip(&values[j * d..(j + 1) * d])
            .map(|(a, b)| a.total_cmp(b))
            .find(|cmp| !cmp.is_eq())
            .unwrap_or(std::cmp::Ordering::Equal)
    });
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
    if avc <= floor {
        return Err(CoreError::calibration("unattainable sampled SCPC decay").with_message(format!(
            "Target correlation avc={avc} must exceed the pooled landmark co-located-pair fraction {floor:.6}."
        )));
    }
    let c0 = calibration::solve_c0_from_avc(distances, avc)?;
    let cmax = calibration::solve_c0_from_avc(distances, floor + (1.0 - floor) * min_avc)?;
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

    let decay = decay_limits(&distances, request.avc, request.min_avc)?;

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
        |(triangle, subset): (&[f64], &Vec<usize>)| -> Result<(Vec<f64>, usize), CoreError> {
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
                let local = basis::scpc_basis(&dense, m, decay[0], requested)?;
                if local.cutoff_tied {
                    requested = (requested + requested.div_ceil(2)).min(m - 1);
                    continue;
                }
                let retained = *local.candidates.last().ok_or_else(|| {
                    CoreError::degenerate("landmark sample has no positive spatial directions")
                })?;
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
            Ok((projected, retained))
        };

    let projected: Vec<(Vec<f64>, usize)> = if crate::threads::is_serial() || samples == 1 {
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

    let total_columns = projected.iter().map(|(_, width)| width).sum();
    let mut wall = crate::zeroed(&[n, total_columns], "embedded landmark directions")?;
    let mut offset = 0;
    for (projection, width) in projected {
        for (destination, source) in wall
            .chunks_exact_mut(total_columns)
            .zip(projection.chunks_exact(width))
        {
            destination[offset..offset + width].copy_from_slice(source);
        }
        offset += width;
    }
    let (weights, blocks) = combine_landmarks(&wall, n, q, total_columns)?;
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
        if kept < rank {
            blocks.ends.retain(|&end| end <= kept);
            blocks.cutoff_tied = blocks.ends.is_empty();
            blocks.more_positive = false;
        }
        let final_rank = blocks.ends.last().copied().unwrap_or(0);
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
        let (out, blocks) = combine_landmarks(&wall, 5, 2, 2).unwrap();
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
    fn colocated_landmark_sample_is_rejected_without_injecting_roundoff_noise() {
        // 8 co-located sites out of 9 (floor = 28/36 < 0.9), landmarks = 7.
        // With 200 samples, an all-co-located subset is drawn with overwhelming probability.
        let mut coords = vec![0.0_f64; 16];
        coords.extend_from_slice(&[1.0, 0.0]);
        let err = landmark_basis(
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
        .err()
        .expect("an all-colocated landmark subset must not accept roundoff directions");
        assert!(
            err.message().contains("no positive spatial directions"),
            "got: {err:?}"
        );
    }
}
