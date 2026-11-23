use crate::{checked_len, linalg, CoreError};
use rayon::prelude::*;

const TILE: usize = 256;

/// The exponential family at infinite decay. Combine rows at identical
/// coordinates before forming the Gram matrix, including different declared sites.
pub(crate) fn independence_limit(
    locations: &crate::distance::Locations,
    basis: &[f64],
    width: usize,
) -> Result<Vec<f64>, CoreError> {
    let n = locations.representatives.len();
    let mut sums = crate::zeroed(&[n, width], "co-location projection")?;
    for (&g, row) in locations.ids.iter().zip(basis.chunks_exact(width)) {
        for (total, &value) in sums[g * width..(g + 1) * width].iter_mut().zip(row) {
            *total += value;
        }
    }
    let mut grid = crate::zeroed(&[width, width], "infinite-decay covariance")?;
    linalg::multiply(&sums, &sums, &mut grid, width, n, width, true);
    for i in 0..width {
        for j in 0..i {
            let value = grid[i * width + j].midpoint(grid[j * width + i]);
            grid[i * width + j] = value;
            grid[j * width + i] = value;
        }
    }
    Ok(grid)
}

/// The supported centered covariance families. Infinite decay is handled by
/// independence_limit so equivalent sites contribute once via their row sums.
#[derive(Clone, Copy)]
pub(crate) enum Kernel {
    Exponential(f64),
    CenteredExponential(f64),
    HalfLife(f64),
    Brownian,
}

impl Kernel {
    fn fill(self, distance: &[f64], out: &mut [f64]) {
        match self {
            Self::Exponential(rate) => linalg::exp_scaled(distance, -rate, out),
            Self::CenteredExponential(rate) => {
                for (value, &d) in out.iter_mut().zip(distance) {
                    *value = (-rate * d).exp_m1();
                }
            }
            Self::HalfLife(half_life) => {
                let rate = std::f64::consts::LN_2 / half_life;
                for (value, &d) in out.iter_mut().zip(distance) {
                    *value = if rate.is_infinite() {
                        (-std::f64::consts::LN_2 * (d / half_life)).exp()
                    } else {
                        0.5 * (-rate * d).exp_m1() / rate
                    };
                }
            }
            Self::Brownian => {
                for (value, &d) in out.iter_mut().zip(distance) {
                    *value = -0.5 * d;
                }
            }
        }
    }
}

/// Project a symmetric covariance kernel for several bases and decay scales.
pub(crate) fn project(
    n: usize,
    bases: &[&[f64]],
    widths: &[usize],
    scales: &[f64],
    fill_tile: impl Fn(usize, usize, usize, usize, &mut [f64]) + Sync,
) -> Result<Vec<Vec<f64>>, CoreError> {
    if scales
        .iter()
        .any(|&scale| !scale.is_finite() || scale < 0.0)
    {
        return Err(CoreError::invalid_input(
            "covariance scales must be finite and nonnegative",
        ));
    }
    let kernels = scales
        .iter()
        .copied()
        .map(Kernel::Exponential)
        .collect::<Vec<_>>();
    project_kernels(n, bases, widths, &kernels, fill_tile)
}

pub(crate) fn project_kernels(
    n: usize,
    bases: &[&[f64]],
    widths: &[usize],
    kernels: &[Kernel],
    fill_tile: impl Fn(usize, usize, usize, usize, &mut [f64]) + Sync,
) -> Result<Vec<Vec<f64>>, CoreError> {
    if n == 0 || bases.is_empty() || bases.len() != widths.len() || kernels.is_empty() {
        return Err(CoreError::invalid_input(
            "expected positive covariance-grid dimensions",
        ));
    }
    let mut offsets = Vec::with_capacity(widths.len());
    let mut scale_width = 0usize;
    for (&basis, &width) in bases.iter().zip(widths) {
        if width == 0 || basis.len() != checked_len(&[n, width])? {
            return Err(CoreError::invalid_input(
                "covariance-grid basis dimensions do not match",
            ));
        }
        if basis.iter().any(|value| !value.is_finite()) {
            return Err(CoreError::invalid_input(
                "covariance-grid basis must be finite",
            ));
        }
        offsets.push(scale_width);
        scale_width = scale_width
            .checked_add(checked_len(&[width, width])?)
            .ok_or(CoreError::invalid_input("covariance grid is too large"))?;
    }
    let size = checked_len(&[kernels.len(), scale_width])?;
    let blocks = n.div_ceil(TILE);
    let mut tiles = Vec::new();
    for i in 0..blocks {
        for j in i..blocks {
            tiles.push((i * TILE, j * TILE));
        }
    }
    let grid = Projection {
        n,
        bases,
        widths,
        offsets: &offsets,
        scale_width,
        fill_tile: &fill_tile,
    };
    let mut result = crate::zeroed(&[size], "projected covariance grid")?;
    // Keep shard boundaries independent of the execution budget so summation
    // does not change with the number of workers.
    if tiles.len() > 8 && size <= 2_000_000 {
        let stride = tiles.len().div_ceil(8);
        let eval = |chunk: &[(usize, usize)]| {
            let mut shard = vec![0.0; size];
            grid.evaluate(chunk, kernels, &mut shard);
            shard
        };
        let shards: Vec<_> = if crate::threads::is_serial() {
            tiles.chunks(stride).map(eval).collect()
        } else {
            tiles.par_chunks(stride).map(eval).collect()
        };
        for shard in shards {
            for (total, value) in result.iter_mut().zip(shard) {
                *total += value;
            }
        }
    } else if crate::threads::is_serial() {
        grid.evaluate(&tiles, kernels, &mut result);
    } else {
        result
            .par_chunks_mut(scale_width)
            .zip(kernels.par_iter())
            .for_each(|(out, kernel)| {
                grid.evaluate(&tiles, std::slice::from_ref(kernel), out);
            });
    }

    Ok(widths
        .iter()
        .zip(offsets)
        .map(|(&width, start)| {
            let matrix_size = width * width;
            result
                .chunks_exact(scale_width)
                .flat_map(|scale| scale[start..start + matrix_size].iter().copied())
                .collect()
        })
        .collect())
}

struct Projection<'a, F> {
    n: usize,
    bases: &'a [&'a [f64]],
    widths: &'a [usize],
    offsets: &'a [usize],
    scale_width: usize,
    fill_tile: &'a F,
}

impl<F: Fn(usize, usize, usize, usize, &mut [f64])> Projection<'_, F> {
    fn evaluate(&self, tiles: &[(usize, usize)], kernels: &[Kernel], output: &mut [f64]) {
        let width_max = *self.widths.iter().max().expect("at least one basis");
        let block_size = self.n.min(TILE);
        let mut distance = vec![0.0; block_size * block_size];
        let mut kernel = vec![0.0; block_size * block_size];
        let mut product = vec![0.0; block_size * width_max];
        let mut contribution = vec![0.0; width_max * width_max];

        for &(i, j) in tiles {
            let ni = (self.n - i).min(TILE);
            let nj = (self.n - j).min(TILE);
            let distance = &mut distance[..ni * nj];
            let kernel = &mut kernel[..ni * nj];
            (self.fill_tile)(i, j, ni, nj, distance);
            for (&scale, out) in kernels
                .iter()
                .zip(output.chunks_exact_mut(self.scale_width))
            {
                scale.fill(distance, kernel);
                for ((&basis, &q), &offset) in self.bases.iter().zip(self.widths).zip(self.offsets)
                {
                    linalg::multiply(
                        kernel,
                        &basis[j * q..(j + nj) * q],
                        &mut product,
                        ni,
                        nj,
                        q,
                        false,
                    );
                    linalg::multiply(
                        &basis[i * q..(i + ni) * q],
                        &product,
                        &mut contribution,
                        q,
                        ni,
                        q,
                        true,
                    );
                    for a in 0..q {
                        for b in a..q {
                            let ab = a * q + b;
                            let ba = b * q + a;
                            let value = if i == j {
                                0.5 * (contribution[ab] + contribution[ba])
                            } else {
                                contribution[ab] + contribution[ba]
                            };
                            out[offset + ab] += value;
                            if ab != ba {
                                out[offset + ba] += value;
                            }
                        }
                    }
                }
            }
        }
    }
}
