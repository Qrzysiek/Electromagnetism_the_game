//! A floating contact's noise against equipartition (prototypes/noise/EQUIPARTITION.md):
//! copy into crates/physics/tests/ and run with --ignored; RESERVOIR=1 (Poisson), 2 (a regular
//! train), 3 (Poisson with each release's Boltzmann factor, the potential evaluated) or 4 (the
//! same to first order, the face's field over δ) needs reservoir_contacts.patch applied to
//! crates/physics/src/tube/emission.rs (`git apply`).
//! Parameters by environment variable: L, W, CL, CPA, DT, T0, T1, SEEDS, LAG, SEG, SOFT, MARGIN,
//! FIXED (a fixed right contact: explore_nyquist).
#![allow(clippy::all, clippy::pedantic, clippy::disallowed_methods)]

use physics::DVec3;
use physics::bem::Bias;
use physics::circuit::Waveform;
use physics::lattice::Lattice;
use physics::tube::{Circuit, Emitter, OhmicContact, Tube};
use physics::zinv::{Electrode, Electrodes, Section};
use rayon::prelude::*;

mod common;

const A: f64 = 0.00125;
const MEAN_FREE_PATH: f64 = 0.5;

fn env(name: &str, default: f64) -> f64 {
    std::env::var(name).ok().and_then(|v| v.parse().ok()).unwrap_or(default)
}

/// A strip l × w: left contact grounded (ideal), right contact a floating node with a lumped
/// capacitance `c_l`, optionally kicked by charge pulses through `r_s`; a grounded frame.
fn noisy(l: f64, w: f64, c_l: f64, kick: Option<(f64, f64, f64, f64)>, carriers: usize, dt: f64, seed: u64) -> Tube {
    let contact = |x: f64| Electrode {
        section: Section::Rect {
            center: DVec3::new(x, 0.5 * w, 0.0),
            angle: 0.0,
            half_length: 0.5,
            half_thickness: 0.5 * w,
        },
        bias: Bias::Grounded,
    };
    let margin = env("MARGIN", 3.0);
    let frame = Electrode {
        section: Section::Frame {
            center: DVec3::new(0.5 * l, 0.5 * w, 0.0),
            half_width: 0.5 * l + 1.0 + margin,
            half_height: 0.5 * w + margin,
            segment: env("SEG", 0.25),
        },
        bias: Bias::Grounded,
    };
    let electrodes = Electrodes::new(vec![contact(-0.5), contact(l + 0.5), frame], &[], env("SEG", 0.25));
    let (lo, hi) = (DVec3::ZERO, DVec3::new(l, w, 0.0));
    let count = libm::round(l * w / (2.0 * A * MEAN_FREE_PATH)) as usize;
    let cores = Lattice::random_cores(lo, hi, count, A, seed);
    let lattice = Lattice::new(lo, hi, common::metal(A), cores, carriers, seed + 4);
    let weight = lattice.carrier_charge().abs();
    let emitter = |electrode: usize, toward: f64| Emitter {
        electrode,
        charge_per_mass: -1.0,
        weight,
        emit_toward: Some(DVec3::new(toward, 0.0, 0.0)),
        saturation: None,
        thermal_speed: 1.0,
        ohmic: Some(OhmicContact {
            density: 1.0,
            // RESERVOIR 1: Poisson; 2: a regular train.
            depth: match env("RESERVOIR", 0.0) as i32 {
                1 => -1.0,
                2 => -2.0,
                3 => -3.0,
                4 => -4.0,
                _ => 0.5 * MEAN_FREE_PATH,
            },
        }),
    };
    let ground = Circuit {
        electrode: 0,
        source: Waveform::Dc(0.0),
        resistance: 0.0,
        capacitance: 0.0,
        inductance: 0.0,
        to: None,
    };
    let fixed = env("FIXED", 0.0);
    let node = match kick {
        None if fixed != 0.0 => Circuit {
            electrode: 1,
            source: Waveform::Dc(fixed),
            resistance: 0.0,
            capacitance: 0.0,
            inductance: 0.0,
            to: None,
        },
        None => Circuit {
            electrode: 1,
            source: Waveform::Dc(0.0),
            resistance: f64::INFINITY,
            capacitance: c_l,
            inductance: 0.0,
            to: None,
        },
        Some((v_s, r_s, delay, period)) => Circuit {
            electrode: 1,
            source: Waveform::Pulse {
                low: 0.0,
                high: v_s,
                delay,
                rise: 0.0,
                width: 0.5 * dt,
                fall: 0.0,
                period,
            },
            resistance: r_s,
            capacitance: c_l,
            inductance: 0.0,
            to: None,
        },
    };
    Tube {
        emitters: vec![emitter(0, 1.0), emitter(1, -1.0)],
        arena: Some((DVec3::new(-2.0, -2.0, 0.0), DVec3::new(l + 2.0, w + 2.0, 0.0))),
        circuits: vec![ground, node],
        lattices: vec![lattice],
        ..Tube::new(electrodes, env("SOFT", 0.05), dt)
    }
}

fn acf(v: &[f64], max_lag: usize) -> Vec<f64> {
    let n = v.len();
    let m = v.iter().sum::<f64>() / n as f64;
    let d: Vec<f64> = v.iter().map(|x| x - m).collect();
    let c0 = d.iter().map(|x| x * x).sum::<f64>() / n as f64;
    (0..max_lag)
        .map(|k| d[..n - k].iter().zip(&d[k..]).map(|(a, b)| a * b).sum::<f64>() / ((n - k) as f64 * c0))
        .collect()
}

fn one_e(r: &[f64], dt: f64) -> f64 {
    let level = (-1.0f64).exp();
    for k in 1..r.len() {
        if r[k] < level {
            return ((k - 1) as f64 + (r[k - 1] - level) / (r[k - 1] - r[k])) * dt;
        }
    }
    f64::NAN
}

#[test]
#[ignore]
fn explore_equilibrium() {
    let (l, w) = (env("L", 4.0), env("W", 4.0));
    let c_l = env("CL", 0.0);
    let per_area = env("CPA", 31.25);
    let carriers = (per_area * l * w) as usize;
    let dt = env("DT", 0.02);
    let (t0, t1) = (env("T0", 10.0), env("T1", 210.0));
    let seeds = env("SEEDS", 4.0) as u64;
    let started = std::time::Instant::now();
    let runs: Vec<(f64, f64, f64, Vec<f64>)> = (0..seeds)
        .into_par_iter()
        .map(|k| {
            let tube = noisy(l, w, c_l, None, carriers, dt, 101 + 17 * k);
            let kt = tube.lattices[0].carrier_charge().abs();
            let c_node = tube.electrodes.capacitance[1][1] + c_l;
            let mut s = tube.state();
            let mut v = Vec::new();
            while s.t < t1 {
                tube.step(&mut s);
                if s.t >= t0 {
                    v.push(s.potentials[1]);
                }
            }
            (kt, c_node, tube.electrodes.capacitance[1][0], v)
        })
        .collect();
    let elapsed = started.elapsed().as_secs_f64();
    let max_lag = (env("LAG", 40.0) / dt) as usize;
    let mut mean_acf = vec![0.0; max_lag];
    for (kt, c_node, c10, v) in &runs {
        let n = v.len() as f64;
        let m = v.iter().sum::<f64>() / n;
        let var = v.iter().map(|x| (x - m) * (x - m)).sum::<f64>() / n;
        let r = acf(v, max_lag);
        println!(
            "kT {kt:.5} C_node {c_node:.4} (C10 {c10:.4}) mean {m:+.4} var {var:.5e} var*C/kT {:.4} tau(1/e) {:.3}",
            var * c_node / kt,
            one_e(&r, dt)
        );
        for (a, b) in mean_acf.iter_mut().zip(&r) {
            *a += b / runs.len() as f64;
        }
    }
    let fs: Vec<f64> = runs
        .iter()
        .map(|(kt, c, _, v)| {
            let n = v.len() as f64;
            let m = v.iter().sum::<f64>() / n;
            v.iter().map(|x| (x - m) * (x - m)).sum::<f64>() / n * c / kt
        })
        .collect();
    let m = fs.len() as f64;
    let fm = fs.iter().sum::<f64>() / m;
    let se = (fs.iter().map(|f| (f - fm) * (f - fm)).sum::<f64>() / (m * (m - 1.0))).sqrt();
    println!("SUMMARY var*C/kT {fm:.4} +- {se:.4} (C_node {:.4})", runs[0].1);
    let probes = [0.02, 0.1, 0.2, 0.5, 1.0, 2.0, 3.0, 5.0, 8.0, 12.0, 20.0, 30.0];
    println!(
        "mean acf: tau(1/e) {:.3}; {}",
        one_e(&mean_acf, dt),
        probes
            .iter()
            .filter(|&&p| ((p / dt) as usize) < max_lag)
            .map(|&p| format!("{p}:{:.3}", mean_acf[(p / dt) as usize]))
            .collect::<Vec<_>>()
            .join(" ")
    );
    println!("{carriers} carriers, {seeds} seeds x {} units in {elapsed:.1} s", t1 - t0 + t0);
}

/// D7 on a short strip: the charge carried's spread slope against 2kT/R.
#[test]
#[ignore]
fn explore_nyquist() {
    let (l, w) = (env("L", 2.0), env("W", 4.0));
    let v = env("FIXED", 0.2);
    let per_area = env("CPA", 31.25);
    let carriers = (per_area * l * w) as usize;
    let dt = env("DT", 0.02);
    let (t0, t1) = (env("T0", 10.0), env("T1", 410.0));
    let seeds = env("SEEDS", 16.0) as u64;
    let runs: Vec<(f64, f64)> = (0..seeds)
        .into_par_iter()
        .map(|k| {
            let tube = noisy(l, w, 0.0, None, carriers, dt, 101 + 17 * k);
            let kt = tube.lattices[0].carrier_charge().abs();
            let mut s = tube.state();
            let mut q = 0.0;
            let mut series = Vec::new();
            while s.t < t1 {
                let before: std::collections::HashMap<u32, f64> =
                    s.id.iter().zip(&s.x).map(|(id, x)| (*id, x.x)).collect();
                tube.step(&mut s);
                if s.t < t0 {
                    continue;
                }
                q += s.id.iter().zip(&s.x).zip(&s.charge)
                    .filter_map(|((id, x), c)| before.get(id).map(|x0| c * (x.x - x0)))
                    .sum::<f64>() / l;
                series.push(q);
            }
            let current = q / (t1 - t0);
            let r = v / current.abs();
            let detrended: Vec<f64> = series.iter().enumerate().map(|(i, x)| x - current * dt * (i + 1) as f64).collect();
            let spread = |lag: f64| {
                let n = (lag / dt).round() as usize;
                let m = detrended.len() - n;
                (0..m).map(|k| (detrended[k + n] - detrended[k]).powi(2)).sum::<f64>() / m as f64
            };
            let slope = (spread(4.0) - spread(1.0)) / 3.0;
            (r, slope * r / (2.0 * kt))
        })
        .collect();
    let m = runs.len() as f64;
    let mean = runs.iter().map(|r| r.1).sum::<f64>() / m;
    let se = (runs.iter().map(|r| (r.1 - mean).powi(2)).sum::<f64>() / (m * (m - 1.0))).sqrt();
    let r = runs.iter().map(|r| r.0).sum::<f64>() / m;
    println!("NYQUIST slope*R/2kT {mean:.4} +- {se:.4}, R {r:.3}");
}
