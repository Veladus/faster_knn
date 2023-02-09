use good_lp::{constraint, default_solver, variables, Solution, SolverModel};

use super::{
    quadratic::InvertedQuadraticFunction,
};
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Serialize, Deserialize)]
struct ApproximationSegment {
    // end of segment (exclusive)
    end_segment: usize,
    inversion: InvertedQuadraticFunction,
    max_distance: f64,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
/// An approximation for a set of points by a piecewise quadratic function
pub struct QuadraticApproximation {
    approximation_segments: Vec<ApproximationSegment>,
}

impl QuadraticApproximation {
    fn first_overestimation_index(cur_index: usize, threshold: f64) -> usize {
        let how_many_smaller_are_allowed = (cur_index as f64 + 1.0) * (1.0 + threshold);
        if ((how_many_smaller_are_allowed + 0.5).fract() - 0.5)
            .abs()
            < 1e-8
        {
            // "exactly" at a threshold
            // so return how_many_smaller_allowed
            // due to rounding we must add some stuff
            (how_many_smaller_are_allowed + 0.5) as usize
        } else {
            how_many_smaller_are_allowed as usize
        }
    }

    fn approximate_start(
        // i-th, self_dist, threshold dist
        points: &[(usize, f64, Option<f64>)],
    ) -> ([f64; 3], usize)
    {
        let empty_dist: f64 = points.last().unwrap().1 * 2.0 + 1.0;
        let find_start_approx = |number_of_points: usize| {
            variables! {
                vars:
                    a;
                    b;
                    c;
                    0 <= eps_strict <= 1000.0;
            }
            let mut problem = vars.maximise(eps_strict).using(default_solver);
            for index in 0..number_of_points {
                let (ith, self_dist, threshold_opt) = points[index];
                let threshold_dist = threshold_opt.unwrap_or(empty_dist);
                let ithf64 = ith as f64;
                problem.add_constraint(
                    constraint! { a * (self_dist * self_dist) + b * self_dist + c + 0.1 <= ithf64},
                );
                problem.add_constraint(
                    constraint! { a * (threshold_dist * threshold_dist) + b * threshold_dist + c - eps_strict >= ithf64},
                );
            }

            match problem.solve() {
                Ok(solution) => {
                    let eps_val = solution.eval(eps_strict);
                    if eps_val < 1e-4 {
                        None
                    } else {
                        let a_val = solution.eval(a);
                        let b_val = solution.eval(b);
                        let c_val = solution.eval(c);

                        Some([
                            a_val,
                            b_val,
                            c_val,
                        ])
                    }
                }
                Err(_) => None,
            }
        };
        let mut separation = find_start_approx(1).unwrap_or_else(||
            [0.0, 1.0, points[0].1 - points[0].0 as f64]
        );
        let mut separation_len = 1;

        let mut separation_bound = loop {
            let testing_range = points.len().min(2 * separation_len);
            match find_start_approx(testing_range) {
                Some(plane) => {
                    separation = plane;
                    separation_len = testing_range;

                    if testing_range == points.len() {
                        return (separation, points.len());
                    }
                }
                None => {
                    break testing_range;
                }
            }
        };

        while separation_len + 1 < separation_bound {
            let middle = (separation_len + separation_bound) / 2;

            match find_start_approx(middle) {
                Some(plane) => {
                    separation = plane;
                    separation_len = middle;
                }
                None => {
                    separation_bound = middle;
                }
            }
        }

        return (separation, separation_len);
    }

    fn from_distances_inner(
        distances: &[f64],
        threshold: f64,
        num_segments: usize,
    ) -> Option<Self> {
        let enriched_distances =
            distances
                .iter()
                .enumerate()
                .fold(Some(Vec::new()), |vec_opt, (index, dist)| {
                    vec_opt.and_then(|mut vec| {
                        let first_overestimation =
                            Self::first_overestimation_index(index, threshold);

                        if let Some(one_past_overestimation) = distances.get(first_overestimation) {
                            if *one_past_overestimation <= dist + 0.2 {
                                // such an approximation can't exist
                                return None;
                            }
                        }

                        let distance_threshold_opt = distances.get(first_overestimation).copied();
                        let entry = (index + 1, *dist, distance_threshold_opt);
                        vec.push(entry);

                        Some(vec)
                    })
                })?;

        let mut approximation_segments = Vec::new();
        let mut remaining_points = enriched_distances.len();

        while remaining_points > 0 && approximation_segments.len() < num_segments {
            let start_index = enriched_distances.len() - remaining_points;
            let slice = &enriched_distances[start_index..];
            let (separation, points) = Self::approximate_start(slice);

            let start_segment = enriched_distances.len() - remaining_points;
            let end_segment = start_segment + points;
            let max_distance = distances[end_segment - 1];

            remaining_points -= points;

            let [a, b, c] = separation;
            // println!("{} * x^2 + {} * x + {}", a, b, c);
            let inversion =
                if let Some((left_arc, right_arc, _)) = InvertedQuadraticFunction::with(a, b, c) {
                    // only one of the sides can be sensible
                    if let (ith, _, Some(upper_bound)) = enriched_distances[start_index] {
                        let right_approx_opt = right_arc.x_value_for(ith as f64);
                        if right_approx_opt.is_none() {
                            println!("{} * x^2 + {} * x + {} at {}", a, b, c, ith);
                        }
                        let right_approx = right_approx_opt.unwrap();

                        if upper_bound > right_approx {
                            right_arc
                        } else {
                            left_arc
                        }
                    } else {
                        right_arc
                    }
                } else {
                    InvertedQuadraticFunction::constant_output(max_distance)
                };

            approximation_segments.push(ApproximationSegment {
                end_segment,
                inversion,
                max_distance,
            });
        }

        if approximation_segments.last().unwrap().end_segment == distances.len() {
            Some(Self {
                approximation_segments,
            })
        } else {
            None
        }
    }

    /// Creates a [QuadraticApproximation] from the given distances.
    /// For each query point the number of distances within the approximation is at most (1+threshold)
    /// times the number of actual distances at most that of the query point.
    pub fn from_distances_with_overestimation_threshold(
        distances: &[f64],
        threshold: f64,
    ) -> Option<Self> {
        Self::from_distances_inner(distances, threshold, distances.len())
    }

    /// Creates a [QuadraticApproximation] from the given distances.
    /// The approximation will be done with at most `num_segments` quadratic functions.
    pub fn from_distances_with_segment_threshold(distances: &[f64], num_segments: usize) -> Self {
        let mut left_factor = 0.0;
        let mut right_factor = distances.len() as f64;

        const PRECISION: f64 = 1e-5;
        let num_iterations = (((distances.len() as f64 + 1.0) / PRECISION).log2()).ceil() as usize;
        for _ in 0..num_iterations {
            let middle_factor = (left_factor + right_factor) / 2.0;
            let opt =
                Self::from_distances_inner(distances, middle_factor, num_segments);

            if opt.is_some() {
                right_factor = middle_factor;
            } else {
                left_factor = middle_factor;
            }
        }

        Self::from_distances_inner(distances, right_factor, num_segments).unwrap()
    }

    /// Gives the approximate distance to the k-th nearest neighbor.
    /// This approximation always gives an upper bound.
    pub fn approximate_knn_distance(&self, k: usize) -> Option<f64> {
        let matching_index = self
            .approximation_segments
            .partition_point(|segment| segment.end_segment <= k);
        self.approximation_segments
            .get(matching_index)
            .map(|segment| {
                segment
                    .inversion
                    .x_value_for((k + 1) as f64)
                    // .map(|inverted_value| inverted_value.min(segment.max_distance))
                    .unwrap_or(segment.max_distance)
            })
            .map(|approx| approx.min(self.approximation_segments.last().unwrap().max_distance))
    }
}

#[cfg(test)]
mod tests {
    use proptest::prelude::*;

    use super::*;

    macro_rules! check_overestimation_threshold {
        ($assert:ident!, $distances:expr, $threshold:expr, $approx:expr) => {
            for (i, real_dist) in $distances.iter().enumerate() {
                // an approximation should be given for all supplied values
                let current_approx = $approx.approximate_knn_distance(i).unwrap();
                $assert!(
                    *real_dist <= current_approx + 1.0,
                    "Approximation underestimated\nreal value: {}\napproximation: {}",
                    real_dist,
                    current_approx,
                );

                let first_above_approximation =
                    $distances.partition_point(|dist| *dist <= current_approx - 1.0);
                    let last_allowed = QuadraticApproximation::first_overestimation_index(i, $threshold);
                $assert!(
                    first_above_approximation
                        <= last_allowed,
                    "Approximation overestimated too much\nlast allowed overestimation index: {} with dist {:?}\nreal first overestimation index: {} with approximation {}",
                    last_allowed,
                    $distances.get(last_allowed),
                    first_above_approximation,
                    current_approx,
                );
            }
        };
    }

    macro_rules! build_and_check_overestimation_threshold {
        (prop_assume!, $assert:ident!, $distances:expr, $threshold:expr) => {
            {
                let mut tests = $distances.clone();
                tests.sort_by(|a, b| a.partial_cmp(&b).unwrap());
                for i in (0..$distances.len()) {
                    let first_overestimation_index = QuadraticApproximation::first_overestimation_index(i, $threshold);
                    if first_overestimation_index < tests.len() {
                        prop_assume!(tests[first_overestimation_index] > tests[i] + 1.0);
                    }
                }

                build_and_check_overestimation_threshold!($assert!, $distances, $threshold);
            }
        };
        ($assert:ident!, $distances:expr, $threshold:expr) => {
            {
                let approx = QuadraticApproximation::from_distances_with_overestimation_threshold(
                    &$distances, $threshold,
                )
                .unwrap();

                println!("{:#?}", approx);
                check_overestimation_threshold!(assert!, $distances, $threshold, approx);

                approx
            }
        };
    }

    #[test]
    fn points_on_quadratic_function_are_approximable_by_one_quadratic() {
        const LENGTH: usize = 50;
        let factor: f64 = 100.0 / (((LENGTH + 1) as f64).sqrt() - (LENGTH as f64).sqrt());
        let distances: Vec<_> = (0..LENGTH)
            .map(|ith| factor * ((ith + 1) as f64).sqrt())
            .collect();

        let approx = build_and_check_overestimation_threshold!(assert!, distances, 0.0);

        assert_eq!(approx.approximation_segments.len(), 1);
    }

    #[test]
    fn very_simple_example() {
        let distances = vec![10.0, 20.0, 30.0, 40.0];
        let threshold = 0.0;

        build_and_check_overestimation_threshold!(assert!, distances, threshold);
    }

    #[test]
    fn regression_take_correct_arc_right() {
        let distances = vec![
            10.0,
            83.61978122286034,
            88.75568242934371,
            91.39771593241738,
            97.1991578286616,
        ];
        let threshold = 0.0;

        build_and_check_overestimation_threshold!(assert!, distances, threshold);
    }

    #[test]
    fn to_debug() {
        let distances = vec![
            487.3301383201039,
            551.6888748096712,
            794.0422431188146,
            857.403583646506,
        ];
        let threshold = 1.8864806403957868;

        build_and_check_overestimation_threshold!(assert!, distances, threshold);
    }

    #[test]
    fn regression_take_correct_arc_left() {
        let distances = vec![
            10.0,
            125.96804771721295,
            140.43388941986376,
            147.53283066538626,
            148.6883783075379,
            166.82191213358652,
            284.5548124371531,
            510.7414546749714,
        ];
        let threshold = 0.5457547195489343;

        build_and_check_overestimation_threshold!(assert!, distances, threshold);
    }

    #[test]
    fn regression_dont_overestimate_into_next_segment() {
        let distances = vec![
            10.0,
            125.96804771721295,
            140.43388941986376,
            147.53283066538626,
            148.6883783075379,
            166.82191213358652,
            284.5548124371531,
            510.7414546749714,
        ];
        let threshold = 0.0;

        build_and_check_overestimation_threshold!(assert!, distances, threshold);
    }

    #[test]
    fn regression_wrong_number_of_approximated_points_after_4_points() {
        let distances = vec![
            10.0,
            88.08168025079267,
            297.713902626059,
            310.3028595110715,
            321.7086686223136,
            544.0320839017126,
        ];
        let threshold = 0.0;

        build_and_check_overestimation_threshold!(assert!, distances, threshold);
    }

    proptest! {
        #[test]
        fn with_threshold_overestimation_respects_bounds(
            mut distances in prop::collection::vec(10f64..1000f64, 1..64),
            threshold in 0.0..10.0,
        ) {
            distances.sort_by(|a, b| a.partial_cmp(b).unwrap());

            for i in 0..distances.len() {
                let overestimation_index = QuadraticApproximation::first_overestimation_index(i, threshold);
                if let Some(overestimation_dist) = distances.get(overestimation_index) {
                    // otherwise there is no such approximation
                    prop_assume!(distances[i] < *overestimation_dist);
                }
            }

            build_and_check_overestimation_threshold!(prop_assume!, prop_assert!, distances, threshold);
        }
    }
}
