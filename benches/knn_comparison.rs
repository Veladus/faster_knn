
use std::{
    collections::hash_map::DefaultHasher,
    hash::{Hash, Hasher},
    io::Write,
    panic,
    path::PathBuf,
    sync::atomic::{AtomicUsize, Ordering},
    time::Duration,
};


use approximate_knn::QuadraticApproximation;
use kd_tree::{KdMap, KdTree2};
use rand_xorshift::XorShiftRng;
use rstar::{primitives::PointWithData, PointDistance, RTree};

use criterion::{black_box, criterion_group, criterion_main, Criterion};
use rand::{Rng, SeedableRng};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::fs;
use std::fs::File;
use std::io::{BufReader, BufWriter};
use std::path::Path;


const RNG_SEED: [u8; 16] = *b"0123456789abcdef";

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
struct Point([i64; 2]);

impl Point {
    pub fn distance_2(&self, other: &Self) -> i64 {
        let dx = self.0[0] - other.0[0];
        let dy = self.0[1] - other.0[1];
        dx * dx + dy * dy
    }
}

type RStarPoint = PointWithData<(), [i64; 2]>;

impl From<RStarPoint> for Point {
    fn from(point: RStarPoint) -> Self {
        Self([point.position()[0], point.position()[1]])
    }
}

impl From<Point> for RStarPoint {
    fn from(seg_tree_point: Point) -> Self {
        RStarPoint::new((), [seg_tree_point.0[0], seg_tree_point.0[1]])
    }
}

type KdTreePoint = [i64; 2];

impl From<Point> for KdTreePoint {
    fn from(point: Point) -> Self {
        point.0
    }
}

impl From<KdTreePoint> for Point {
    fn from(point: KdTreePoint) -> Self {
        Self(point)
    }
}

fn default_num_trees() -> usize {
    1
}

#[derive(Debug, Clone, Deserialize)]
enum ApproximationQuality {
    Segments(usize),
    Ratio(f64),
}

impl ApproximationQuality {
    fn to_signature(&self) -> String {
        match self {
            ApproximationQuality::Segments(num_segments) => { format!("S{}", num_segments) }
            ApproximationQuality::Ratio(ratio) => { format!("R{:.4}", ratio) }
        }
    }

    fn approximate(&self, distances: &[f64]) -> QuadraticApproximation {
        match self {
            &ApproximationQuality::Segments(num_segments) => {
                QuadraticApproximation::from_distances_with_segment_threshold(
                    distances,
                    num_segments,
                )
            }
            &ApproximationQuality::Ratio(quality_ratio) => {
                let ratio = quality_ratio - 1.0; // QuadraticApproximation expects overestimation
                let non_adjusted_approximation = QuadraticApproximation::from_distances_with_overestimation_threshold(distances, ratio);
                non_adjusted_approximation.unwrap_or_else(|| {
                    let adjusted_distances: Vec<_> = distances.iter().scan(-1.0, |prev_dist, cur_dist| {
                        *prev_dist = cur_dist.max(*prev_dist + 0.3);
                        Some(*prev_dist)
                    }).collect();
                    QuadraticApproximation::from_distances_with_overestimation_threshold(&adjusted_distances, ratio)
                        .expect("Did not find an approximation with adjusted distances")
                })
            }
        }
    }
}

#[derive(Debug, Clone, Deserialize)]
struct Config {
    approximation_upto: usize,
    approximation_quality: ApproximationQuality,
    #[serde(default = "default_num_trees")]
    num_trees: usize,
    data_points_path: String,
    approximation_points_path: String,
}

impl Config {
    fn get_approximation_tree(
        &self,
        query_positions: &[Point],
        all_points: &RTree<RStarPoint>,
    ) -> KdMap<KdTreePoint, QuadraticApproximation> {
        let query_hash = {
            let mut hashable: Vec<_> = query_positions.iter().collect();
            hashable.sort_unstable();

            let mut hasher = DefaultHasher::new();
            hashable.hash(&mut hasher);
            hasher.finish()
        };

        let points_hash = {
            let mut hashable: Vec<_> = all_points.iter().collect();
            hashable.sort_unstable();

            let mut hasher = DefaultHasher::new();
            hashable.hash(&mut hasher);
            hasher.finish()
        };

        let mut path = PathBuf::from("benches/data");
        fs::create_dir_all(&path).unwrap();

        let file_name = format!(
            "{}-{}-{}-{}.dat",
            query_hash, points_hash, self.approximation_upto, self.approximation_quality.to_signature()
        );
        path.push(file_name);
        eprintln!("looking for approximations in {}", path.to_str().unwrap());

        if path.exists() {
            eprintln!("Found approximations");
            let points: Vec<(Point, QuadraticApproximation)> = serde_json::de::from_reader(
                BufReader::new(fs::File::open(path).expect("Could not open file")),
            )
                .expect("Error deserializing approximtations");

            KdMap::build(
                points
                    .into_iter()
                    .map(|(point, data)| (KdTreePoint::from(point), data))
                    .collect::<Vec<_>>(),
            )
        } else {
            eprintln!("Calculating");
            let cnt = AtomicUsize::new(0);
            let approximation_points: Vec<_> = query_positions
                .par_iter()
                .cloned()
                .filter_map(|point| {
                    match panic::catch_unwind(|| {
                        let distances: Vec<_> = all_points
                            .nearest_neighbor_iter(&point.0)
                            .map(|other_point| {
                                (point.distance_2(&Point(other_point.position().clone())) as f64 + 0.1)
                                    .sqrt()
                            })
                            .take(self.approximation_upto)
                            .collect();

                        let result = (
                            point,
                            self.approximation_quality.approximate(&distances),
                        );
                        let done = cnt.fetch_add(1, Ordering::SeqCst) + 1;
                        if done > 0 && done % 100 == 0 {
                            eprintln!("Completed {} approximations", done);
                        }
                        result
                    }) {
                        Ok(out) => Some(out),
                        Err(_) => {
                            eprintln!("Error during approximation calculation");
                            None
                        }
                    }
                })
                .collect();

            let mut file = fs::File::create(path).unwrap();
            serde_json::ser::to_writer(BufWriter::new(&mut file), &approximation_points).unwrap();
            file.flush().unwrap();

            KdMap::build(
                approximation_points
                    .into_iter()
                    .map(|(point, data)| (KdTreePoint::from(point), data))
                    .collect::<Vec<_>>(),
            )
        }
    }
}

#[derive(Debug, Clone, Deserialize, Serialize)]
struct CsvDataPoint {
    x: f64,
    y: f64,
    data: Option<String>,
}

fn load_point_set_from_csv(path: impl AsRef<Path>) -> Vec<Point> {
    eprintln!("loading csv data from {}", path.as_ref().to_string_lossy());
    let mut reader = csv::Reader::from_path(path).expect("Could not open CSV reader to file");

    reader
        .deserialize::<CsvDataPoint>()
        .map(|csv_point_res| {
            let csv_point = csv_point_res.expect("Error reading point");
            Point([csv_point.x.round() as i64, csv_point.y.round() as i64])
        })
        .collect()
}

fn get_point_sets<R: Rng>(config: &Config, _rng: &mut R) -> Vec<Vec<Point>> {
    vec![load_point_set_from_csv(&config.data_points_path); config.num_trees]
}

fn get_approximation_point_sets<R: Rng>(config: &Config, _rng: &mut R) -> Vec<Vec<Point>> {
    vec![load_point_set_from_csv(&config.approximation_points_path); config.num_trees]
}

fn get_rtrees(all_points: &[Vec<Point>]) -> Vec<RTree<RStarPoint>> {
    all_points
        .iter()
        .map(|cur_point_set| {
            let points = cur_point_set
                .clone()
                .into_iter()
                .map(|point| point.into())
                .collect();
            RTree::bulk_load(points)
        })
        .collect()
}

fn get_kd_trees(all_points: &[Vec<Point>]) -> Vec<KdTree2<KdTreePoint>> {
    all_points
        .iter()
        .map(|cur_point_set| {
            let points: Vec<KdTreePoint> = cur_point_set
                .clone()
                .into_iter()
                .map(|point| point.into())
                .collect();
            KdTree2::build(points)
        })
        .collect()
}

fn get_approximation_trees(
    config: &Config,
    query_points: &[Vec<Point>],
    rtrees: &[RTree<RStarPoint>],
) -> Vec<KdMap<KdTreePoint, QuadraticApproximation>> {
    query_points
        .iter()
        .enumerate()
        .map(|(index, _)| config.get_approximation_tree(&query_points[index], &rtrees[index]))
        .collect()
}

fn get_random_query_with_count<R: Rng>(
    neighbors: usize,
    possible_positions: &[Point],
    rng: &mut R,
) -> (Point, usize) {
    (
        possible_positions[rng.gen_range(0..possible_positions.len())],
        neighbors,
    )
}

fn get_random_query_with_distance_exact_2<R: Rng>(
    neighbours: usize,
    possible_positions: &[Point],
    tree: &RTree<RStarPoint>,
    rng: &mut R,
) -> (Point, i64) {
    let point = possible_positions[rng.gen_range(0..possible_positions.len())];
    let rstar_point: RStarPoint = point.into();
    (
        point,
        tree.nearest_neighbor_iter_with_distance_2(rstar_point.position())
            .take(neighbours)
            .last()
            .unwrap()
            .1,
    )
}

fn enrich_with_distance_approximation(
    point: Point,
    size: usize,
    approx_tree: &KdMap<KdTreePoint, QuadraticApproximation>,
) -> (Point, f64) {
    let approx_point_opt = approx_tree.nearest(&point.0);
    if let Some(approx_point) = approx_point_opt {
        let approx_opt = &approx_point.item.1;
        (
            point,
            approx_opt
                .approximate_knn_distance(size)
                .unwrap_or_else(|| {
                    eprintln!("broken no approx");
                    0.0
                })
                + (point.distance_2(&Point(approx_point.item.0.clone())) as f64 + 1.0).sqrt(),
        )
    } else {
        eprintln!("broken no point");
        (point, 1.0)
    }
}

fn rstar_query_bench_nearest(query: (Point, usize), container: &RTree<RStarPoint>) -> usize {
    let (point, count) = query;
    let point: RStarPoint = point.into();
    black_box(
        container
            .nearest_neighbor_iter(point.position())
            .take(count),
    )
        .count()
}

fn rstar_query_bench_within_distance(query: (Point, f64), container: &RTree<RStarPoint>) -> usize {
    let (point, distance) = query;
    let point: RStarPoint = point.into();
    black_box(
        container.locate_within_distance(*point.position(), (distance.ceil().powi(2) + 1.0) as i64),
    )
        .count()
}

fn rstar_query_bench_within_distance_exact(
    query: (Point, f64),
    count: usize,
    container: &RTree<RStarPoint>,
) -> usize {
    let (point, distance) = query;
    let point: RStarPoint = point.into();
    let mut all: Vec<_> = container
        .locate_within_distance(*point.position(), (distance.ceil().powi(2) + 1.0) as i64)
        .collect();
    let partition_point = count.min(all.len()).saturating_sub(1);
    if partition_point > 0 {
        all.select_nth_unstable_by_key(partition_point, |other| {
            point.position().distance_2(other.position())
        });
    }
    black_box(all.len())
}

fn kd_tree_query_bench_nearest(query: (Point, usize), container: &KdTree2<KdTreePoint>) -> usize {
    let (point, count) = query;
    let point: KdTreePoint = point.into();
    let all = container.nearests(&point, count);
    all.len()
}

fn kd_tree_query_bench_within_distance(
    query: (Point, f64),
    container: &KdTree2<KdTreePoint>,
) -> usize {
    let (point, distance) = query;
    let point: KdTreePoint = point.into();
    let all = container.within_radius(&point, (distance.ceil() + 1.0) as i64);
    black_box(all.len())
}

fn kd_tree_query_bench_within_distance_exact(
    query: (Point, f64),
    count: usize,
    container: &KdTree2<KdTreePoint>,
) -> usize {
    let (point, distance) = query;
    let kd_point: KdTreePoint = point.into();
    let mut all: Vec<_> = container.within_radius(&kd_point, (distance.ceil() + 1.0) as i64);
    let partition_point = count.min(all.len()).saturating_sub(1);
    if partition_point > 0 {
        all.select_nth_unstable_by_key(partition_point, |other| {
            point.distance_2(&(**other).into())
        });
    }
    black_box(all.len())
}

fn size_tree_group<C>(
    _config: &Config,
    c: &mut Criterion,
    rng: &mut impl Rng,
    size: usize,
    approximations: &Vec<KdMap<KdTreePoint, QuadraticApproximation>>,
    all_points: &Vec<Vec<Point>>,
    r_trees: &Vec<RTree<RStarPoint>>,
    container_name: &'static str,
    get_container: fn(&[Vec<Point>]) -> Vec<C>,
    nearest_query_opt: Option<fn((Point, usize), &C) -> usize>,
    distance_query: fn((Point, f64), &C) -> usize,
    distance_post_selection_query: fn((Point, f64), usize, &C) -> usize,
) {
    let containers = get_container(&all_points);

    let mut group = c.benchmark_group(format!("{}/{}", size, container_name));
    let running_time = std::env::var("RUNNING_TIME").ok().and_then(|s| s.parse().ok()).unwrap_or(10);
    group.measurement_time(Duration::new(running_time, 0));

    group.bench_function("exact_distance", |b| {
        b.iter_with_setup(
            || {
                let index = rng.gen_range(0..containers.len());
                let (point, dist) = get_random_query_with_distance_exact_2(
                    size,
                    &all_points[index],
                    &r_trees[index],
                    rng,
                );
                ((point, (dist as f64).sqrt() + 1.0), index)
            },
            |(query, index)| assert!(distance_query(query, &containers[index]) >= size),
        );
    });

    group.bench_function("approximate_distance", |b| {
        b.iter_with_setup(
            || {
                let index = rng.gen_range(0..containers.len());
                let point = all_points[index][rng.gen_range(0..all_points[index].len())];
                (point, index)
            },
            |(point, index)| {
                let query = enrich_with_distance_approximation(point, size, &approximations[index]);
                let found_size = distance_query(query, &containers[index]);
                if found_size < size {
                    eprintln!(
                        "broken approx dist {} < {} {} {:?}",
                        found_size, size, index, query
                    );
                }
            },
        );
    });

    group.bench_function("exact_distance_post_selection", |b| {
        b.iter_with_setup(
            || {
                let index = rng.gen_range(0..containers.len());
                let (point, dist) = get_random_query_with_distance_exact_2(
                    size,
                    &all_points[index],
                    &r_trees[index],
                    rng,
                );
                ((point, (dist as f64).sqrt() + 1.0), index)
            },
            |(query, index)| {
                assert!(distance_post_selection_query(query, size, &containers[index]) >= size)
            },
        );
    });

    group.bench_function("approximate_distance_post_selection", |b| {
        b.iter_with_setup(
            || {
                let index = rng.gen_range(0..containers.len());
                let point = all_points[index][rng.gen_range(0..all_points[index].len())];
                (point, index)
            },
            |(point, index)| {
                let query = enrich_with_distance_approximation(point, size, &approximations[index]);
                let found_size = distance_post_selection_query(query, size, &containers[index]);
                if found_size < size {
                    eprintln!(
                        "broken approx dist select {} < {} {} {:?}",
                        found_size, size, index, query
                    );
                }
            },
        );
    });

    if let Some(nearest_query) = nearest_query_opt {
        group.bench_function("plain-knn", |b| {
            b.iter_with_setup(
                || {
                    let index = rng.gen_range(0..containers.len());
                    (
                        get_random_query_with_count(size, &all_points[index], rng),
                        index,
                    )
                },
                |(query, idx)| assert!(nearest_query(query, &containers[idx]) >= size),
            );
        });
    }

    group.finish();
}

fn size_groups(c: &mut Criterion, size: usize) {
    let config_path = std::env::var("CONFIG_FILE")
        .expect("You need to provide a CONFIG_FILE as environment variable");
    let config: Config =
        serde_json::from_reader(File::open(config_path).expect("Could not open CONFIG_FILE"))
            .expect("CONFIG_FILE could not be read");

    let rng = XorShiftRng::from_seed(RNG_SEED);
    let all_points = get_point_sets(&config, &mut rng.clone());
    let r_trees = get_rtrees(&all_points);

    let all_approximation_points = get_approximation_point_sets(&config, &mut rng.clone());
    let approximations = get_approximation_trees(&config, &all_approximation_points, &r_trees);

    if std::env::var("COMPUTE_ONLY").map_or(false, |value| value.to_lowercase() == "true") {
        std::process::exit(0);
    }

    size_tree_group(
        &config,
        c,
        &mut rng.clone(),
        size,
        &approximations,
        &all_points,
        &r_trees,
        "rstar",
        get_rtrees,
        Some(rstar_query_bench_nearest),
        rstar_query_bench_within_distance,
        rstar_query_bench_within_distance_exact,
    );

    size_tree_group(
        &config,
        c,
        &mut rng.clone(),
        size,
        &approximations,
        &all_points,
        &r_trees,
        "kd",
        get_kd_trees,
        Some(kd_tree_query_bench_nearest),
        kd_tree_query_bench_within_distance,
        kd_tree_query_bench_within_distance_exact,
    );
}

fn all_groups(c: &mut Criterion) {
    for size in &[10, 30, 100, 300, 1000, 3000] {
        size_groups(c, *size);
    }
}

criterion_group!(all_groups_wrapper, all_groups);
criterion_main!(all_groups_wrapper);
