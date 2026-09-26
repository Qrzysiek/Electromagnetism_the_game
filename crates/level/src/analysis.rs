//! Difficulty analysis of a level (SPEC §7.6).
//!
//! A good puzzle has a large configuration space, a small solution set (so guessing or
//! brute force is hopeless), and a "distance to solution" that changes fairly smoothly
//! under small changes (so a player can learn from attempts). This module measures:
//!
//! - **configuration space**: number of distinct player placements within the limits;
//! - **random solution rate**: fraction of uniformly random placements that solve the
//!   level (pure guessing);
//! - **heuristic search**: how often, and how fast, a local search that only sees the
//!   distance objective (like a player watching the preview) finds a solution;
//! - **smoothness**: rank correlation of the distance between a placement and a
//!   neighbouring one (one small move).

use rayon::prelude::*;

use crate::solve::{Rng, is_verified_solution, objective, single_element_options};
use crate::{Element, ElementKind, Level};

#[derive(Clone, Debug)]
pub struct Analysis {
    /// log10 of the number of distinct placements.
    pub config_space_log10: f64,
    pub random_samples: usize,
    pub random_solutions: usize,
    pub search_runs: usize,
    pub search_successes: usize,
    /// Mean number of objective evaluations of the successful searches.
    pub search_mean_evaluations: f64,
    /// Evaluation budget per search run.
    pub search_budget: u32,
    /// Spearman rank correlation of the objective across one-move neighbours.
    pub smoothness: f64,
}

impl Analysis {
    /// Upper estimate of the random solution rate (rule of three when none were found).
    pub fn random_rate(&self) -> f64 {
        #[allow(clippy::cast_precision_loss)]
        let (s, n) = (self.random_solutions as f64, self.random_samples as f64);
        if self.random_solutions == 0 {
            3.0 / n
        } else {
            s / n
        }
    }

    pub fn search_rate(&self) -> f64 {
        #[allow(clippy::cast_precision_loss)]
        let r = self.search_successes as f64 / self.search_runs.max(1) as f64;
        r
    }

    /// Expected objective evaluations for the heuristic search to find a solution
    /// (restarting after each failed run); infinite if it never succeeded.
    pub fn expected_search_effort(&self) -> f64 {
        let p = self.search_rate();
        if p == 0.0 {
            return f64::INFINITY;
        }
        // Failed runs cost the full budget; the successful one its mean.
        f64::from(self.search_budget) * (1.0 - p) / p + self.search_mean_evaluations
    }

    /// Expected evaluations for random guessing to find a solution.
    pub fn expected_guess_effort(&self) -> f64 {
        1.0 / self.random_rate()
    }
}

fn ln_choose(n: u64, k: u64) -> f64 {
    #[allow(clippy::cast_precision_loss)]
    (0..k)
        .map(|i| ((n - i) as f64).ln() - ((i + 1) as f64).ln())
        .sum()
}

/// log10 of the number of placements: subsets of at most `max` options of each kind.
fn config_space_log10(level: &Level, options: &[Element]) -> f64 {
    let count = |kind: ElementKind| options.iter().filter(|e| e.kind == kind).count() as u64;
    let subsets = |n: u64, max: u64| -> f64 {
        // ln Σ_{j=0..max} C(n, j), via log-sum-exp.
        let terms: Vec<f64> = (0..=max.min(n)).map(|j| ln_choose(n, j)).collect();
        let m = terms.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        m + terms.iter().map(|t| (t - m).exp()).sum::<f64>().ln()
    };
    let ln = subsets(
        count(ElementKind::Charge),
        u64::from(level.limits.max_charges),
    ) + subsets(
        count(ElementKind::Magnet),
        u64::from(level.limits.max_magnets),
    ) + subsets(
        count(ElementKind::Antenna),
        u64::from(level.limits.max_antennas),
    );
    ln / std::f64::consts::LN_10
}

/// A random valid placement: a random number of elements of each kind (at least one in
/// total), each a random allowed option.
fn random_placement(level: &Level, options: &[Element], rng: &mut Rng) -> Vec<Element> {
    let by_kind = |k: ElementKind| -> Vec<Element> {
        options.iter().copied().filter(|e| e.kind == k).collect()
    };
    let (charges, magnets, antennas) = (
        by_kind(ElementKind::Charge),
        by_kind(ElementKind::Magnet),
        by_kind(ElementKind::Antenna),
    );
    loop {
        let nc = rng.below(level.limits.max_charges as usize + 1);
        let nm = rng.below(level.limits.max_magnets as usize + 1);
        let na = rng.below(level.limits.max_antennas as usize + 1);
        if nc + nm + na == 0 {
            continue;
        }
        let mut p: Vec<Element> = Vec::new();
        for (pool, n) in [(&charges, nc), (&magnets, nm), (&antennas, na)] {
            let mut tries = 0;
            while p.iter().filter(|e| pool.contains(e)).count() < n && tries < 100 {
                tries += 1;
                let e = pool[rng.below(pool.len())];
                let mut trial = p.clone();
                trial.push(e);
                if level.check_placement(&trial).is_ok() {
                    p = trial;
                }
            }
        }
        if !p.is_empty() {
            return p;
        }
    }
}

/// One small change: move an element by one or two nodes, change its value to a
/// neighbouring allowed one (or rotate an antenna by one step), flip its sign, or
/// add/remove an element.
fn neighbour(level: &Level, p: &[Element], options: &[Element], rng: &mut Rng) -> Vec<Element> {
    for _ in 0..50 {
        let mut t = p.to_vec();
        match rng.below(5) {
            0 | 1 if !t.is_empty() => {
                let i = rng.below(t.len());
                let d = |r: &mut Rng| i64::try_from(r.below(5)).unwrap() - 2;
                let n = t[i].node;
                t[i].node = [n[0] + d(rng), n[1] + d(rng), n[2]];
            }
            2 if t.iter().any(|e| e.kind == ElementKind::Antenna) && rng.below(2) == 0 => {
                // Rotate an antenna by one allowed step.
                let idx: Vec<usize> = (0..t.len())
                    .filter(|&j| t[j].kind == ElementKind::Antenna)
                    .collect();
                let i = idx[rng.below(idx.len())];
                let angles = crate::ANTENNA_ANGLES;
                let k = angles
                    .iter()
                    .position(|a| a.to_bits() == t[i].angle_deg.to_bits())
                    .unwrap_or(0);
                let k2 = if rng.below(2) == 0 {
                    (k + angles.len() - 1) % angles.len()
                } else {
                    (k + 1) % angles.len()
                };
                t[i].angle_deg = angles[k2];
            }
            2 if !t.is_empty() => {
                let i = rng.below(t.len());
                let list = match t[i].kind {
                    ElementKind::Charge => &level.limits.magnitudes,
                    ElementKind::Magnet => &level.limits.magnet_strengths,
                    ElementKind::Antenna => &level.limits.antenna_amplitudes,
                };
                let k = list
                    .iter()
                    .position(|m| m.to_bits() == t[i].value.abs().to_bits())
                    .unwrap_or(0);
                let k2 = if rng.below(2) == 0 {
                    k.saturating_sub(1)
                } else {
                    (k + 1).min(list.len() - 1)
                };
                t[i].value = t[i].value.signum() * list[k2];
            }
            3 if !t.is_empty() => {
                let i = rng.below(t.len());
                t[i].value = -t[i].value;
            }
            _ => {
                if !t.is_empty() && rng.below(2) == 0 {
                    let i = rng.below(t.len());
                    t.remove(i);
                } else {
                    t.push(options[rng.below(options.len())]);
                }
            }
        }
        if !t.is_empty() && t != p && level.check_placement(&t).is_ok() {
            return t;
        }
    }
    p.to_vec()
}

/// Heuristic search that only sees the distance objective: greedy local search with
/// occasional acceptance of worse moves. Returns the evaluations used if it found a
/// verified solution within `budget`.
fn search(level: &Level, options: &[Element], budget: u32, seed: u64) -> Option<u32> {
    let mut rng = Rng::new(seed);
    let mut current = random_placement(level, options, &mut rng);
    let (mut score, mut outcome) = objective(level, &current);
    for eval in 1..=budget {
        if outcome == physics::trajectory::Outcome::Arrived && is_verified_solution(level, &current)
        {
            return Some(eval);
        }
        let trial = neighbour(level, &current, options, &mut rng);
        let (s, o) = objective(level, &trial);
        let temperature = 0.5 * (1.0 - f64::from(eval) / f64::from(budget));
        let accept = s <= score || rng.unit() < (-(s - score) / temperature.max(1e-9)).exp();
        if accept {
            current = trial;
            score = s;
            outcome = o;
        }
    }
    None
}

fn ranks(v: &[f64]) -> Vec<f64> {
    let mut idx: Vec<usize> = (0..v.len()).collect();
    idx.sort_by(|&a, &b| v[a].total_cmp(&v[b]));
    let mut r = vec![0.0; v.len()];
    for (rank, &i) in idx.iter().enumerate() {
        #[allow(clippy::cast_precision_loss)]
        {
            r[i] = rank as f64;
        }
    }
    r
}

fn spearman(a: &[f64], b: &[f64]) -> f64 {
    let (ra, rb) = (ranks(a), ranks(b));
    #[allow(clippy::cast_precision_loss)]
    let n = ra.len() as f64;
    let (ma, mb) = (ra.iter().sum::<f64>() / n, rb.iter().sum::<f64>() / n);
    let cov: f64 = ra.iter().zip(&rb).map(|(x, y)| (x - ma) * (y - mb)).sum();
    let va: f64 = ra.iter().map(|x| (x - ma).powi(2)).sum();
    let vb: f64 = rb.iter().map(|y| (y - mb).powi(2)).sum();
    cov / (va * vb).sqrt()
}

pub fn analyze(level: &Level, samples: usize, runs: usize, budget: u32, seed: u64) -> Analysis {
    let options = single_element_options(level);
    let config_space_log10 = config_space_log10(level, &options);

    // Random guessing.
    let random_solutions = (0..samples)
        .into_par_iter()
        .filter(|&i| {
            let mut rng = Rng::new(seed ^ (i as u64).wrapping_mul(0x9E37_79B9));
            let p = random_placement(level, &options, &mut rng);
            objective(level, &p).1 == physics::trajectory::Outcome::Arrived
                && is_verified_solution(level, &p)
        })
        .count();

    // Heuristic search.
    let results: Vec<Option<u32>> = (0..runs)
        .into_par_iter()
        .map(|r| search(level, &options, budget, seed.wrapping_add(1000 + r as u64)))
        .collect();
    let successes: Vec<u32> = results.into_iter().flatten().collect();
    #[allow(clippy::cast_precision_loss)]
    let search_mean_evaluations = if successes.is_empty() {
        f64::NAN
    } else {
        successes.iter().map(|&e| f64::from(e)).sum::<f64>() / successes.len() as f64
    };

    // Smoothness.
    let pairs: Vec<(f64, f64)> = (0..samples.min(400))
        .into_par_iter()
        .map(|i| {
            let mut rng = Rng::new(seed ^ 0xABCD ^ (i as u64).wrapping_mul(31));
            let p = random_placement(level, &options, &mut rng);
            let q = neighbour(level, &p, &options, &mut rng);
            (objective(level, &p).0, objective(level, &q).0)
        })
        .collect();
    let (a, b): (Vec<f64>, Vec<f64>) = pairs.into_iter().unzip();

    Analysis {
        config_space_log10,
        random_samples: samples,
        random_solutions,
        search_runs: runs,
        search_successes: successes.len(),
        search_mean_evaluations,
        search_budget: budget,
        smoothness: spearman(&a, &b),
    }
}
