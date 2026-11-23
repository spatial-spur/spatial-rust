//! Geographic and Euclidean distances in the host's requested units.

use crate::{checked_len, CoreError};
use rayon::prelude::*;
use std::collections::HashMap;

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Metric {
    Euclidean,
    /// Great-circle distance as fraction of 2 * pi * R (equivalently asin(...) / pi).
    /// Input columns are `(lat, lon)` in degrees.
    GreatCircleFractionLatLon,
    /// Same as `GreatCircleFractionLatLon`, with input columns `(lon, lat)` in degrees.
    GreatCircleFractionLonLat,
    /// Great-circle distance in the units of `radius`. Input columns are `(lat, lon)` in degrees.
    GreatCircle {
        radius: f64,
    },
}

#[derive(Clone, Copy)]
struct GeographicPoint {
    lat: f64,
    lon: f64,
    cos_lat: f64,
}

impl GeographicPoint {
    fn new(lat: f64, lon: f64) -> Result<Self, CoreError> {
        if lat.abs() > 90.0 {
            return Err(CoreError::invalid_input(
                "latitude must be between -90 and 90 degrees",
            ));
        }
        let (lat, lon, pole) = (
            positive_zero(lat),
            canonical_longitude(lon),
            lat.abs() == 90.0,
        );
        Ok(Self {
            lat,
            lon: if pole { 0.0 } else { lon },
            cos_lat: if pole { 0.0 } else { lat.to_radians().cos() },
        })
    }

    fn half_angle(self, other: Self) -> f64 {
        // Difference in degrees first. Converting positions separately breaks
        // symmetry on grids, and subtracting 360 near the seam loses precision.
        let lo = self.lon.min(other.lon);
        let hi = self.lon.max(other.lon);
        let delta = if hi - lo > 180.0 {
            (180.0 - hi) + (lo + 180.0)
        } else {
            hi - lo
        };
        let lat = (0.5 * (self.lat - other.lat).abs()).to_radians().sin();
        let lon = (0.5 * delta).to_radians().sin();
        lat.hypot((self.cos_lat * other.cos_lat).sqrt() * lon)
            .min(1.0)
            .asin()
    }
}

fn positive_zero(value: f64) -> f64 {
    if value == 0.0 {
        0.0
    } else {
        value
    }
}

fn canonical_longitude(lon: f64) -> f64 {
    let lon = lon % 360.0;
    positive_zero(if lon >= 180.0 {
        lon - 360.0
    } else if lon < -180.0 {
        lon + 360.0
    } else {
        lon
    })
}

/// The distinction between repeated observations at a location and separately
/// declared regression sites belongs to the caller. This map is for geometry.
pub(crate) struct Locations {
    pub ids: Vec<usize>,
    pub representatives: Vec<usize>,
}

impl Locations {
    pub fn pair_fraction(&self) -> f64 {
        let mut counts = vec![0.0; self.representatives.len()];
        let mut matching = 0.0;
        for &g in &self.ids {
            matching += counts[g];
            counts[g] += 1.0;
        }
        let n = self.ids.len() as f64;
        matching / (n * (n - 1.0) * 0.5)
    }
}

pub(crate) struct Distances<'a> {
    coordinates: &'a [f64],
    width: usize,
    geographic: Option<Vec<GeographicPoint>>,
    multiplier: f64,
    metric: Metric,
}

impl<'a> Distances<'a> {
    pub(crate) fn new(
        coordinates: &'a [f64],
        n: usize,
        width: usize,
        metric: Metric,
    ) -> Result<Self, CoreError> {
        if width == 0 || coordinates.len() != checked_len(&[n, width])? {
            return Err(CoreError::invalid_input(
                "coordinate dimensions do not match",
            ));
        }
        if coordinates.iter().any(|v| !v.is_finite()) {
            return Err(CoreError::invalid_input("coordinates must be finite"));
        }
        let geographic = match metric {
            Metric::Euclidean => None,
            _ => {
                if width != 2 {
                    return Err(CoreError::invalid_input(
                        "geographic coordinates require two columns",
                    ));
                }
                let lat = usize::from(metric == Metric::GreatCircleFractionLonLat);
                Some(
                    coordinates
                        .chunks_exact(2)
                        .map(|row| GeographicPoint::new(row[lat], row[1 - lat]))
                        .collect::<Result<_, _>>()?,
                )
            }
        };
        let multiplier = match metric {
            Metric::Euclidean => 1.0,
            Metric::GreatCircleFractionLatLon | Metric::GreatCircleFractionLonLat => {
                std::f64::consts::FRAC_1_PI
            }
            Metric::GreatCircle { radius } if radius.is_finite() && radius > 0.0 => 2.0 * radius,
            Metric::GreatCircle { .. } => {
                return Err(CoreError::invalid_input(
                    "Earth radius must be finite and positive",
                ))
            }
        };
        Ok(Self {
            coordinates,
            width,
            geographic,
            multiplier,
            metric,
        })
    }

    pub fn canonical_coordinates(&self) -> Option<Vec<f64>> {
        self.geographic.as_ref().map(|points| {
            let lat_first = self.metric != Metric::GreatCircleFractionLonLat;
            points
                .iter()
                .flat_map(|p| {
                    if lat_first {
                        [p.lat, p.lon]
                    } else {
                        [p.lon, p.lat]
                    }
                })
                .collect()
        })
    }

    pub fn coordinate_order(&self) -> Vec<usize> {
        let mut rows: Vec<_> = (0..self.coordinates.len() / self.width).collect();
        rows.sort_by(|&i, &j| {
            if let Some(points) = &self.geographic {
                return points[i]
                    .lat
                    .total_cmp(&points[j].lat)
                    .then(points[i].lon.total_cmp(&points[j].lon));
            }
            self.coordinates[i * self.width..(i + 1) * self.width]
                .iter()
                .zip(&self.coordinates[j * self.width..(j + 1) * self.width])
                .map(|(&a, &b)| positive_zero(a).total_cmp(&positive_zero(b)))
                .find(|ord| !ord.is_eq())
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        rows
    }

    pub fn locations(&self) -> Locations {
        let mut ids = Vec::with_capacity(self.coordinates.len() / self.width);
        let mut representatives = Vec::new();
        let mut groups = HashMap::new();
        for (i, row) in self.coordinates.chunks_exact(self.width).enumerate() {
            let point;
            let coords = if let Some(geographic) = &self.geographic {
                point = [geographic[i].lat, geographic[i].lon];
                &point[..]
            } else {
                row
            };
            let key: Vec<_> = coords.iter().map(|&v| positive_zero(v).to_bits()).collect();
            let id = groups.entry(key).or_insert_with(|| {
                representatives.push(i);
                representatives.len() - 1
            });
            ids.push(*id);
        }
        Locations {
            ids,
            representatives,
        }
    }

    pub fn rounding_resolution(&self, i: usize) -> f64 {
        if self.geographic.is_some() {
            return 0.0;
        }
        let row = &self.coordinates[i * self.width..(i + 1) * self.width];
        row.iter().fold(0.0_f64, |norm, &v| {
            let bits = v.abs().to_bits();
            let next = f64::from_bits(if v.abs() == f64::MAX {
                bits - 1
            } else {
                bits + 1
            });
            norm.hypot((next - v.abs()).abs())
        })
    }

    pub(crate) fn pair(&self, i: usize, j: usize) -> f64 {
        if let Some(points) = &self.geographic {
            self.multiplier * points[i].half_angle(points[j])
        } else {
            let a = &self.coordinates[i * self.width..(i + 1) * self.width];
            let b = &self.coordinates[j * self.width..(j + 1) * self.width];
            let squares = a
                .iter()
                .zip(b)
                .map(|(left, right)| (left - right).powi(2))
                .sum::<f64>();
            if squares.is_finite() && squares >= f64::MIN_POSITIVE {
                squares.sqrt()
            } else {
                // Norms remain representable even when their squares underflow
                // or overflow. Keep the fast path for ordinary coordinates.
                a.iter()
                    .zip(b)
                    .fold(0.0_f64, |norm, (left, right)| norm.hypot(left - right))
            }
        }
    }

    pub(crate) fn fill_tile(&self, i: usize, j: usize, nj: usize, scale: f64, tile: &mut [f64]) {
        for (row, values) in tile.chunks_exact_mut(nj).enumerate() {
            for (col, value) in values.iter_mut().enumerate() {
                *value = self.pair(i + row, j + col) / scale;
            }
        }
    }
}

/// Canonical coordinate encoding used for exact co-location checks.
pub fn canonical_coordinates(
    coordinates: &[f64],
    n: usize,
    width: usize,
    metric: Metric,
) -> Result<Vec<f64>, CoreError> {
    let distances = Distances::new(coordinates, n, width, metric)?;
    Ok(distances
        .canonical_coordinates()
        .unwrap_or_else(|| coordinates.to_vec()))
}

/// The full symmetric distance matrix in row-major order.
pub fn distance_matrix(
    coordinates: &[f64],
    n: usize,
    width: usize,
    metric: Metric,
) -> Result<Vec<f64>, CoreError> {
    let distances = Distances::new(coordinates, n, width, metric)?;
    let mut matrix = crate::zeroed(&[n, n], "dense distance matrix")?;
    if n == 0 {
        return Ok(matrix);
    }
    if crate::threads::is_serial() || n < 64 {
        for (i, row) in matrix.chunks_exact_mut(n).enumerate() {
            for (j, value) in row[..i].iter_mut().enumerate() {
                *value = distances.pair(i, j);
            }
        }
    } else {
        matrix.par_chunks_mut(n).enumerate().for_each(|(i, row)| {
            for (j, value) in row[..i].iter_mut().enumerate() {
                *value = distances.pair(i, j);
            }
        });
    }
    for i in 0..n {
        for j in 0..i {
            matrix[j * n + i] = matrix[i * n + j];
        }
    }
    if matrix.iter().any(|d| !d.is_finite()) {
        return Err(CoreError::invalid_input(
            "pairwise distance exceeds the floating-point range",
        ));
    }
    Ok(matrix)
}

pub(crate) fn distinct_location_count(distances: &[f64], n: usize) -> usize {
    distances
        .chunks_exact(n)
        .enumerate()
        .filter(|(i, row)| !row[..*i].contains(&0.0))
        .count()
}

/// Extract strict lower triangle from a symmetric `n x n` matrix.
pub(crate) fn lower_triangle(distances: &[f64], n: usize) -> Vec<f64> {
    let pairs = n.saturating_mul(n.saturating_sub(1)) / 2;
    let mut out = Vec::with_capacity(pairs);
    for (i, row) in distances.chunks_exact(n).enumerate() {
        out.extend_from_slice(&row[..i]);
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn euclidean_distance_handles_extreme_units() {
        for scale in [1e-200, 1e-163, 1.0, 1e200] {
            let coords = [0.0, 0.0, 3.0 * scale, 4.0 * scale];
            let distances = distance_matrix(&coords, 2, 2, Metric::Euclidean).unwrap();
            assert!((distances[1] / scale - 5.0).abs() < 2e-15);
            assert_eq!(distances[1], distances[2]);
        }
        assert!(distance_matrix(&[f64::MAX, -f64::MAX], 2, 1, Metric::Euclidean).is_err());
    }
}
