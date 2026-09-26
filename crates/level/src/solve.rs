//! Searching for player placements that solve a level (SPEC §7.2).

use crate::{Element, ElementKind, Level, Node};
use physics::trajectory::{Outcome, RunSettings, run};
use rayon::prelude::*;

/// Deterministic pseudo-random numbers (SplitMix64).
pub struct Rng(u64);

impl Rng {
    pub fn new(seed: u64) -> Self {
        Self(seed)
    }

    pub fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    /// Uniform in [0, 1) with 53 random bits.
    #[allow(clippy::cast_precision_loss)] // 53-bit integers are exact in f64
    pub fn unit(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }

    /// Uniform integer in `0..n`.
    pub fn below(&mut self, n: usize) -> usize {
        usize::try_from(self.next_u64() % n as u64).expect("fits")
    }
}

/// Every single element the player may place (all free nodes × allowed signed values of
/// every allowed kind).
pub fn single_element_options(level: &Level) -> Vec<Element> {
    let m = level.grid.max_node();
    let z_range = if level.grid.is_2d() { 0..=0 } else { 0..=m[2] };
    let mut kinds: Vec<(ElementKind, f64, f64, Option<f64>)> = Vec::new();
    for &q in &level.limits.magnitudes {
        if level.limits.allow_positive {
            kinds.push((ElementKind::Charge, q, 0.0, None));
        }
        if level.limits.allow_negative {
            kinds.push((ElementKind::Charge, -q, 0.0, None));
        }
    }
    if level.limits.max_magnets > 0 {
        for &mu in &level.limits.magnet_strengths {
            kinds.push((ElementKind::Magnet, mu, 0.0, None));
            kinds.push((ElementKind::Magnet, -mu, 0.0, None));
        }
    }
    if level.limits.max_antennas > 0 {
        let omegas: Vec<Option<f64>> = if level.limits.antenna_omegas.is_empty() {
            vec![None]
        } else {
            level
                .limits
                .antenna_omegas
                .iter()
                .copied()
                .map(Some)
                .collect()
        };
        for &p in &level.limits.antenna_amplitudes {
            for a in crate::ANTENNA_ANGLES {
                for &w in &omegas {
                    kinds.push((ElementKind::Antenna, p, a, w));
                    kinds.push((ElementKind::Antenna, -p, a, w));
                }
            }
        }
    }
    let mut out = Vec::new();
    for z in z_range {
        for y in 0..=m[1] {
            for x in 0..=m[0] {
                for &(kind, value, angle_deg, omega) in &kinds {
                    let e = Element {
                        node: [x, y, z],
                        kind,
                        value,
                        angle_deg,
                        omega,
                    };
                    if level.check_placement(&[e]).is_ok() {
                        out.push(e);
                    }
                }
            }
        }
    }
    out
}

/// Search objective at the preview tolerance, summed over all shots: 0 for a shot that
/// arrives, otherwise its closest approach to its detector before the flight ended (a
/// continuous measure of how near the placement is to a solution). The outcome is
/// `Arrived` only if every shot arrives; otherwise that of the first shot that does not.
pub fn objective(level: &Level, placement: &[Element]) -> (f64, Outcome) {
    let mut rs = RunSettings::with_tolerance(level.physics.tolerances.preview);
    rs.record = false;
    let mut score = 0.0;
    let mut outcome = Outcome::Arrived;
    for scn in level.scenarios(placement) {
        let tr = run(&scn, &rs);
        if tr.outcome != Outcome::Arrived {
            let d = tr
                .margins
                .as_ref()
                .and_then(|m| m.detector)
                .unwrap_or(f64::INFINITY);
            score += d.max(0.0);
            // In the detector but outside its acceptance: how far outside (radians or
            // relative energy), so that the search can improve on it.
            if tr.outcome == Outcome::Rejected {
                let a = tr
                    .margins
                    .as_ref()
                    .and_then(|m| m.acceptance)
                    .unwrap_or(1.0);
                score += (-a).max(0.0) + 1e-3;
            }
            if outcome == Outcome::Arrived {
                outcome = tr.outcome;
            }
        }
    }
    (score, outcome)
}

/// Whether a placement is a verified solution: every shot arrives, verified.
pub fn is_verified_solution(level: &Level, placement: &[Element]) -> bool {
    level.check_placement(placement).is_ok()
        && level
            .verify_flights(placement)
            .iter()
            .all(|v| v.outcome() == Outcome::Arrived && v.status.is_verified())
}

/// All verified single-element solutions (exhaustive).
pub fn single_element_solutions(level: &Level) -> Vec<Element> {
    single_element_options(level)
        .into_par_iter()
        .filter(|c| objective(level, &[*c]).1 == Outcome::Arrived)
        .filter(|c| is_verified_solution(level, &[*c]))
        .collect()
}

fn canonical(mut p: Vec<Element>) -> Vec<Element> {
    p.sort_by(|a, b| a.node.cmp(&b.node).then(a.value.total_cmp(&b.value)));
    p
}

/// Random-restart simulated annealing for placements of exactly `k` elements. Returns the
/// distinct verified solutions found.
pub fn anneal(
    level: &Level,
    k: usize,
    restarts: u64,
    iterations: u32,
    seed: u64,
) -> Vec<Vec<Element>> {
    let options = single_element_options(level);
    let mut found: Vec<Vec<Element>> = (0..restarts)
        .into_par_iter()
        .filter_map(|r| anneal_once(level, &options, k, iterations, seed.wrapping_add(r)))
        .map(canonical)
        .collect();
    found.sort_by(|a, b| {
        let key = |p: &Vec<Element>| {
            p.iter()
                .map(|c| (c.node, c.value.to_bits()))
                .collect::<Vec<_>>()
        };
        key(a).cmp(&key(b))
    });
    found.dedup();
    found
}

fn anneal_once(
    level: &Level,
    options: &[Element],
    k: usize,
    iterations: u32,
    seed: u64,
) -> Option<Vec<Element>> {
    let mut rng = Rng::new(seed);
    let mut current: Vec<Element> = Vec::with_capacity(k);
    while current.len() < k {
        let c = options[rng.below(options.len())];
        let mut trial = current.clone();
        trial.push(c);
        if level.check_placement(&trial).is_ok() {
            current = trial;
        }
    }
    let (mut score, _) = objective(level, &current);
    // Allowed (value, orientation, frequency) combinations per kind.
    let values = |kind: ElementKind| -> Vec<(f64, f64, Option<f64>)> {
        let mut v: Vec<(f64, f64, Option<f64>)> = options
            .iter()
            .filter(|c| c.kind == kind)
            .map(|c| (c.value, c.angle_deg, c.omega))
            .collect();
        let key =
            |x: &(f64, f64, Option<f64>)| (x.0.to_bits(), x.1.to_bits(), x.2.map(f64::to_bits));
        v.sort_by(|a, b| {
            a.0.total_cmp(&b.0)
                .then(a.1.total_cmp(&b.1))
                .then(key(a).2.cmp(&key(b).2))
        });
        v.dedup_by(|a, b| key(a) == key(b));
        v
    };
    let charge_values = values(ElementKind::Charge);
    let magnet_values = values(ElementKind::Magnet);
    let antenna_values = values(ElementKind::Antenna);
    for it in 0..iterations {
        if score == 0.0 && is_verified_solution(level, &current) {
            return Some(current);
        }
        let temperature = 2.0 * (1.0 - f64::from(it) / f64::from(iterations));
        let mut trial = current.clone();
        let i = rng.below(k);
        match rng.below(3) {
            0 => {
                // Move by up to 3 nodes in x and y.
                let dx = i64::try_from(rng.below(7)).unwrap() - 3;
                let dy = i64::try_from(rng.below(7)).unwrap() - 3;
                let n: Node = trial[i].node;
                trial[i].node = [n[0] + dx, n[1] + dy, n[2]];
            }
            1 => {
                let v = match trial[i].kind {
                    ElementKind::Charge => &charge_values,
                    ElementKind::Magnet => &magnet_values,
                    ElementKind::Antenna => &antenna_values,
                };
                (trial[i].value, trial[i].angle_deg, trial[i].omega) = v[rng.below(v.len())];
            }
            _ => trial[i] = options[rng.below(options.len())],
        }
        if level.check_placement(&trial).is_err() {
            continue;
        }
        let (s, _) = objective(level, &trial);
        let accept = s <= score || {
            let u = rng.unit();
            temperature > 0.0 && u < (-(s - score) / temperature.max(1e-9)).exp()
        };
        if accept {
            current = trial;
            score = s;
        }
    }
    (score == 0.0 && is_verified_solution(level, &current)).then_some(current)
}
