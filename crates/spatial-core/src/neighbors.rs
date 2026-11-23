//! Exact metric searches with conservative pruning near floating-point ties.

use crate::distance::{Distances, Locations};
use crate::CoreError;

struct Node {
    point: usize,
    split: f64,
    near: Option<Box<Node>>,
    far: Option<Box<Node>>,
}

pub(crate) struct NeighborIndex<'a> {
    distances: Distances<'a>,
    locations: Locations,
    observations: Vec<Vec<usize>>,
    resolution: Vec<f64>,
    max_resolution: f64,
    root: Option<Box<Node>>,
}

fn checked_distance(distances: &Distances<'_>, i: usize, j: usize) -> Result<f64, CoreError> {
    let distance = distances.pair(i, j);
    if distance.is_finite() {
        Ok(distance)
    } else {
        Err(CoreError::invalid_input(
            "pairwise distance exceeds the floating-point range",
        ))
    }
}

impl Node {
    fn build(
        points: &mut [usize],
        distances: &Distances<'_>,
    ) -> Result<Option<Box<Self>>, CoreError> {
        let Some((point, children)) = points.split_last_mut() else {
            return Ok(None);
        };
        let point = *point;
        let mut node = Self {
            point,
            split: 0.0,
            near: None,
            far: None,
        };
        if !children.is_empty() {
            let mut ranked = children
                .iter()
                .map(|&j| checked_distance(distances, point, j).map(|distance| (j, distance)))
                .collect::<Result<Vec<_>, _>>()?;
            let mid = ranked.len() / 2;
            ranked.select_nth_unstable_by(mid, |a, b| a.1.total_cmp(&b.1).then(a.0.cmp(&b.0)));
            node.split = ranked[mid].1;
            for (dest, (j, _)) in children.iter_mut().zip(ranked) {
                *dest = j;
            }
            let (near, far) = children.split_at_mut(mid);
            node.near = Self::build(near, distances)?;
            node.far = Self::build(far, distances)?;
        }
        Ok(Some(Box::new(node)))
    }

    fn search(
        &self,
        distances: &Distances<'_>,
        query: usize,
        limit: &mut f64,
        visit: &mut impl FnMut(usize, f64, &mut f64),
    ) -> Result<(), CoreError> {
        let distance = checked_distance(distances, query, self.point)?;
        visit(self.point, distance, limit);
        let (first, second) = if distance <= self.split {
            (&self.near, &self.far)
        } else {
            (&self.far, &self.near)
        };
        if let Some(first) = first {
            first.search(distances, query, limit, visit)?;
        }
        // The trigonometric distance is ill-conditioned near antipodes. Use
        // a conservative pruning allowance, then compare actual distances.
        let margin = 4.0 * f64::EPSILON.sqrt() * distance.max(self.split);
        if (distance - self.split).abs() <= *limit + margin {
            if let Some(second) = second {
                second.search(distances, query, limit, visit)?;
            }
        }
        Ok(())
    }
}

impl<'a> NeighborIndex<'a> {
    pub fn new(distances: Distances<'a>) -> Result<Self, CoreError> {
        let locations = distances.locations();
        let mut points = locations.representatives.clone();
        let root = Node::build(&mut points, &distances)?;
        let mut observations = vec![Vec::new(); points.len()];
        for (i, &site) in locations.ids.iter().enumerate() {
            observations[site].push(i);
        }
        let resolution: Vec<_> = locations
            .representatives
            .iter()
            .map(|&i| distances.rounding_resolution(i))
            .collect();
        let max_resolution = resolution.iter().copied().fold(0.0, f64::max);
        Ok(Self {
            distances,
            locations,
            observations,
            resolution,
            max_resolution,
            root,
        })
    }

    pub fn nearest(&self, i: usize, result: &mut Vec<usize>) -> Result<(), CoreError> {
        if self.observations.len() < 2 {
            return Err(CoreError::invalid_input(
                "nearest-neighbor transform requires at least two distinct locations",
            ));
        }
        let query = self.locations.ids[i];
        let i = self.locations.representatives[query];
        let mut best = (f64::INFINITY, query);
        let mut limit = f64::INFINITY;
        self.root.as_ref().expect("nonempty index").search(
            &self.distances,
            i,
            &mut limit,
            &mut |j, d, limit| {
                let site = self.locations.ids[j];
                if site != query {
                    if d < best.0
                        || (d == best.0 && self.resolution[site] > self.resolution[best.1])
                    {
                        best = (d, site);
                        *limit = d
                            + 16.0 * f64::EPSILON * d
                            + self.resolution[query]
                            + self.max_resolution;
                    }
                    if d <= *limit {
                        result.push(site);
                    }
                }
            },
        )?;
        let tolerance = 16.0 * f64::EPSILON * best.0 + self.resolution[query];
        result.retain(|&site| {
            let d = self.distances.pair(i, self.locations.representatives[site]);
            d - best.0 <= tolerance + self.resolution[site].max(self.resolution[best.1])
        });
        let sites = std::mem::take(result);
        result.extend(
            sites
                .into_iter()
                .flat_map(|site| self.observations[site].iter().copied()),
        );
        result.sort_unstable();
        Ok(())
    }

    pub fn within(&self, i: usize, radius: f64, result: &mut Vec<usize>) -> Result<(), CoreError> {
        let site = self.locations.ids[i];
        let query = self.locations.representatives[site];
        let mut limit = radius;
        self.root.as_ref().expect("nonempty index").search(
            &self.distances,
            query,
            &mut limit,
            &mut |j, d, _| {
                if d < radius {
                    result.extend(
                        self.observations[self.locations.ids[j]]
                            .iter()
                            .copied()
                            .filter(|&k| k != i),
                    );
                }
            },
        )?;
        result.sort_unstable();
        Ok(())
    }
}
