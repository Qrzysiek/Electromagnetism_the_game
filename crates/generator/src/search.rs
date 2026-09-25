//! Searching for player placements that solve a level (SPEC §7.2).

use level::{Charge, Level, Node};
use physics::trajectory::{Outcome, RunSettings, run};
use physics::verify::verify;
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

/// Every single charge the player may place (all free nodes × allowed signed magnitudes).
pub fn single_charge_options(level: &Level) -> Vec<Charge> {
    let m = level.grid.max_node();
    let z_range = if level.grid.is_2d() { 0..=0 } else { 0..=m[2] };
    let mut signed = Vec::new();
    for &q in &level.limits.magnitudes {
        if level.limits.allow_positive {
            signed.push(q);
        }
        if level.limits.allow_negative {
            signed.push(-q);
        }
    }
    let mut out = Vec::new();
    for z in z_range {
        for y in 0..=m[1] {
            for x in 0..=m[0] {
                for &q in &signed {
                    let c = Charge {
                        node: [x, y, z],
                        charge: q,
                    };
                    if level.check_placement(&[c]).is_ok() {
                        out.push(c);
                    }
                }
            }
        }
    }
    out
}

/// Search objective at the preview tolerance: 0 when the particle arrives, otherwise the
/// closest approach to the detector before the flight ended (a continuous measure of how
/// near the placement is to a solution).
pub fn objective(level: &Level, placement: &[Charge]) -> (f64, Outcome) {
    let scn = level.scenario(placement);
    let mut rs = RunSettings::with_tolerance(level.physics.tolerances.preview);
    rs.record = false;
    let tr = run(&scn, &rs);
    let d = tr
        .margins
        .as_ref()
        .and_then(|m| m.detector)
        .unwrap_or(f64::INFINITY);
    let score = if tr.outcome == Outcome::Arrived {
        0.0
    } else {
        d.max(0.0)
    };
    (score, tr.outcome)
}

/// Whether a placement is a verified solution.
pub fn is_verified_solution(level: &Level, placement: &[Charge]) -> bool {
    level.check_placement(placement).is_ok() && {
        let v = verify(&level.scenario(placement), level.tolerances());
        v.outcome() == Outcome::Arrived && v.status.is_verified()
    }
}

/// All verified single-charge solutions (exhaustive).
pub fn single_charge_solutions(level: &Level) -> Vec<Charge> {
    single_charge_options(level)
        .into_par_iter()
        .filter(|c| objective(level, &[*c]).1 == Outcome::Arrived)
        .filter(|c| is_verified_solution(level, &[*c]))
        .collect()
}

fn canonical(mut p: Vec<Charge>) -> Vec<Charge> {
    p.sort_by(|a, b| a.node.cmp(&b.node).then(a.charge.total_cmp(&b.charge)));
    p
}

/// Random-restart simulated annealing for placements of exactly `k` charges. Returns the
/// distinct verified solutions found.
pub fn anneal(
    level: &Level,
    k: usize,
    restarts: u64,
    iterations: u32,
    seed: u64,
) -> Vec<Vec<Charge>> {
    let options = single_charge_options(level);
    let mut found: Vec<Vec<Charge>> = (0..restarts)
        .into_par_iter()
        .filter_map(|r| anneal_once(level, &options, k, iterations, seed.wrapping_add(r)))
        .map(canonical)
        .collect();
    found.sort_by(|a, b| {
        let key = |p: &Vec<Charge>| {
            p.iter()
                .map(|c| (c.node, c.charge.to_bits()))
                .collect::<Vec<_>>()
        };
        key(a).cmp(&key(b))
    });
    found.dedup();
    found
}

fn anneal_once(
    level: &Level,
    options: &[Charge],
    k: usize,
    iterations: u32,
    seed: u64,
) -> Option<Vec<Charge>> {
    let mut rng = Rng::new(seed);
    let mut current: Vec<Charge> = Vec::with_capacity(k);
    while current.len() < k {
        let c = options[rng.below(options.len())];
        let mut trial = current.clone();
        trial.push(c);
        if level.check_placement(&trial).is_ok() {
            current = trial;
        }
    }
    let (mut score, _) = objective(level, &current);
    let signed: Vec<f64> = {
        let mut v: Vec<f64> = options.iter().map(|c| c.charge).collect();
        v.sort_by(f64::total_cmp);
        v.dedup();
        v
    };
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
            1 => trial[i].charge = signed[rng.below(signed.len())],
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
