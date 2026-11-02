//! SCRATCH (not committed; excluded in .git/info/exclude): V-D10's dense rings, alone and
//! with a finer self-field cell, against the exact concentric conductors.

use num_complex::Complex64;
use physics::DVec3;
use physics::circuit::Waveform;
use physics::darwin2d::Darwin2d;
use physics::lattice::{Lattice, Material};
use physics::tube::{SolenoidWave, Tube};
use std::f64::consts::PI;

mod common;

const A: f64 = 0.00125;
const MEAN_FREE_PATH: f64 = 0.5;

fn concentric(
    rings: &[(f64, f64)],
    sigma: Complex64,
    omega: f64,
    mu0: f64,
    b0: f64,
) -> Vec<(f64, f64)> {
    let zero = Complex64::new(0.0, 0.0);
    let mut pieces = Vec::new();
    let mut reached = 0.0;
    for (a, b) in rings {
        if *a > reached {
            pieces.push((reached, *a, false));
        }
        pieces.push((*a, *b, true));
        reached = *b;
    }
    let (mut phi, mut bz) = (zero, Complex64::new(1.0, 0.0));
    let mut heats = Vec::new();
    for (a, b, conductor) in pieces {
        let s = if conductor { sigma } else { zero };
        let deriv = |r: f64, phi: Complex64, bz: Complex64| -> (Complex64, Complex64) {
            let db = if r == 0.0 {
                zero
            } else {
                Complex64::new(0.0, omega * mu0) * s * phi / (2.0 * PI * r)
            };
            (bz * (2.0 * PI * r), db)
        };
        let n = 2 * (libm::ceil(2000.0 * (b - a)) as u32).max(1);
        let h = (b - a) / f64::from(n);
        let (mut integral, mut drive) = (0.0, zero);
        for i in 0..=n {
            let r = a + h * f64::from(i);
            let weight = if i == 0 || i == n {
                1.0
            } else if i % 2 == 1 {
                4.0
            } else {
                2.0
            };
            if r > 0.0 {
                let e = phi * Complex64::new(0.0, -omega / (2.0 * PI * r));
                integral += weight * e.norm_sqr() * 2.0 * PI * r;
                drive += Complex64::new(0.0, -0.5 * omega * r) * (s * e).conj() * (weight * 2.0 * PI * r);
            }
            if i == n {
                break;
            }
            let k1 = deriv(r, phi, bz);
            let k2 = deriv(r + 0.5 * h, phi + k1.0 * (0.5 * h), bz + k1.1 * (0.5 * h));
            let k3 = deriv(r + 0.5 * h, phi + k2.0 * (0.5 * h), bz + k2.1 * (0.5 * h));
            let k4 = deriv(r + h, phi + k3.0 * h, bz + k3.1 * h);
            phi += (k1.0 + k2.0 * 2.0 + k3.0 * 2.0 + k4.0) * (h / 6.0);
            bz += (k1.1 + k2.1 * 2.0 + k3.1 * 2.0 + k4.1) * (h / 6.0);
        }
        if conductor {
            heats.push((0.5 * s.re * integral * h / 3.0, drive * (h / 3.0)));
        }
    }
    let scale = b0 * b0 / bz.norm_sqr();
    let conj = (Complex64::new(1.0, 0.0) / bz).conj();
    heats
        .into_iter()
        .map(|(p, d)| (p * scale, 0.5 * b0 * b0 * (conj * d).re))
        .collect()
}

fn annulus(inner: f64, outer: f64, material: Material, carriers: usize, seed: u64) -> Lattice {
    let half = DVec3::new(outer * 1.001, outer * 1.001, 0.0);
    let (lo, hi) = (-half, half);
    let square = (hi.x - lo.x) * (hi.y - lo.y);
    let count = libm::round(square / (2.0 * A * MEAN_FREE_PATH)) as usize;
    let cores = Lattice::random_cores(lo, hi, count, A, seed ^ 31);
    let lattice =
        Lattice::new(lo, hi, material, cores, carriers, seed).with_disc(DVec3::ZERO, outer);
    if inner > 0.0 {
        lattice.with_round_hole(DVec3::ZERO, inner)
    } else {
        lattice
    }
}

fn dense(name: &str, rings: &[(f64, f64)], cell: f64, nodes: usize, seed: u64) {
    let (c, omega, b0, dt, density, per_area) = (24.0, 1.0, 0.2, 0.0072, 60.0, 35.0);
    let mu0 = 4.0 * PI / (c * c);
    let sigma = common::mobility(1.0 / MEAN_FREE_PATH, omega) * density;
    let exact = concentric(rings, sigma, omega, mu0, b0);
    let material = Material {
        background: density,
        ..common::metal(A)
    };
    let lattices: Vec<Lattice> = rings
        .iter()
        .enumerate()
        .map(|(j, (a, b))| {
            let carriers = libm::round(per_area * PI * (b * b - a * a)) as usize;
            annulus(*a, *b, material, carriers, seed + 100 * j as u64)
        })
        .collect();
    let outer = rings[rings.len() - 1].1;
    let electrodes = common::small_floating_electrode(DVec3::new(-40.0, -40.0, 0.0));
    let half = 0.5 * cell * (nodes + 1) as f64;
    let tube = Tube {
        arena: Some((
            DVec3::new(-outer - 2.0, -outer - 2.0, 0.0),
            DVec3::new(outer + 2.0, outer + 2.0, 0.0),
        )),
        b_wave: Some(SolenoidWave {
            waveform: Waveform::Sine {
                offset: 0.0,
                amplitude: b0,
                omega,
                phase: 0.0,
            },
            axis: DVec3::ZERO,
        }),
        lattices,
        darwin: Some(Darwin2d::new(nodes, nodes, cell, (-half, -half), c)),
        ..Tube::new(electrodes, 0.05, dt)
    };
    let period = 2.0 * PI / omega;
    let (start, end) = (2.0 * period, 14.0 * period);
    let mut s = tube.state();
    let mut first = None;
    while s.t < end - 0.5 * dt {
        tube.step(&mut s);
        if first.is_none() && s.t >= start - 0.5 * dt {
            first = Some((s.t, s.induced_work.clone(), s.mutual_work.clone()));
        }
    }
    let (t0, d0, m0) = first.expect("start");
    let span = s.t - t0;
    for j in 0..rings.len() {
        let drive = (s.induced_work[j] - d0[j]) / span;
        let heat = drive + (s.mutual_work[j] - m0[j]) / span;
        println!(
            "SCRATCH {name} seed {seed} ring {:?}: heat {heat:.4} vs exact {:.4} ({:+.1} %); drive {drive:.4} vs exact {:.4} ({:+.1} %)",
            rings[j],
            exact[j].0,
            100.0 * (heat / exact[j].0 - 1.0),
            exact[j].1,
            100.0 * (drive / exact[j].1 - 1.0)
        );
    }
}

#[test]
#[ignore = "scratch"]
fn scratch_inner_alone() {
    dense("inner alone", &[(2.5, 4.5)], 0.25, 127, 47);
}

#[test]
#[ignore = "scratch"]
fn scratch_outer_alone() {
    dense("outer alone", &[(5.0, 6.0)], 0.25, 127, 47);
}

#[test]
#[ignore = "scratch"]
fn scratch_both_fine_cell() {
    dense("both, cell 0.125", &[(2.5, 4.5), (5.0, 6.0)], 0.125, 255, 47);
}

/// B(r)/B0 of the exact concentric conductors at radius r (e^{iωt}).
fn exact_b(rings: &[(f64, f64)], sigma: Complex64, omega: f64, mu0: f64, r_probe: f64) -> Complex64 {
    let zero = Complex64::new(0.0, 0.0);
    let outer = rings[rings.len() - 1].1;
    let n = 400_000;
    let h = outer / f64::from(n);
    let (mut phi, mut bz) = (zero, Complex64::new(1.0, 0.0));
    let mut at = Complex64::new(1.0, 0.0);
    let s_at = |r: f64| if rings.iter().any(|(a, b)| r >= *a && r <= *b) { sigma } else { zero };
    for i in 0..n {
        let r = h * f64::from(i);
        if r <= r_probe {
            at = bz;
        }
        let s = s_at(r + 0.5 * h);
        let f = |rr: f64, p: Complex64, b: Complex64| -> (Complex64, Complex64) {
            let db = if rr == 0.0 { zero } else { Complex64::new(0.0, omega * mu0) * s * p / (2.0 * PI * rr) };
            (b * (2.0 * PI * rr), db)
        };
        let k1 = f(r, phi, bz);
        let k2 = f(r + 0.5 * h, phi + k1.0 * (0.5 * h), bz + k1.1 * (0.5 * h));
        let k3 = f(r + 0.5 * h, phi + k2.0 * (0.5 * h), bz + k2.1 * (0.5 * h));
        let k4 = f(r + h, phi + k3.0 * h, bz + k3.1 * h);
        phi += (k1.0 + k2.0 * 2.0 + k3.0 * 2.0 + k4.0) * (h / 6.0);
        bz += (k1.1 + k2.1 * 2.0 + k3.1 * 2.0 + k4.1) * (h / 6.0);
    }
    at / bz
}

#[test]
#[ignore = "scratch"]
fn scratch_both_fields() {
    let rings = [(2.5, 4.5), (5.0, 6.0)];
    let (c, omega, b0, dt, density, per_area, seed) = (24.0, 1.0, 0.2, 0.0072, 60.0, 35.0, 47_u64);
    let mu0 = 4.0 * PI / (c * c);
    let sigma = common::mobility(1.0 / MEAN_FREE_PATH, omega) * density;
    let material = Material { background: density, ..common::metal(A) };
    let lattices: Vec<Lattice> = rings
        .iter()
        .enumerate()
        .map(|(j, (a, b))| {
            let carriers = libm::round(per_area * PI * (b * b - a * a)) as usize;
            annulus(*a, *b, material, carriers, seed + 100 * j as u64)
        })
        .collect();
    let electrodes = common::small_floating_electrode(DVec3::new(-40.0, -40.0, 0.0));
    let tube = Tube {
        arena: Some((DVec3::new(-8.0, -8.0, 0.0), DVec3::new(8.0, 8.0, 0.0))),
        b_wave: Some(SolenoidWave {
            waveform: Waveform::Sine { offset: 0.0, amplitude: b0, omega, phase: 0.0 },
            axis: DVec3::ZERO,
        }),
        lattices,
        darwin: Some(Darwin2d::new(127, 127, 0.25, (-16.0, -16.0), c)),
        ..Tube::new(electrodes, 0.05, dt)
    };
    let period = 2.0 * PI / omega;
    let (start, end) = (2.0 * period, 14.0 * period);
    let probes = [0.0, 1.5, 4.75, 7.0];
    let turns = [(1.0, 0.0), (0.0, 1.0), (-1.0, 0.0), (0.0, -1.0)];
    let mut sums = vec![Complex64::new(0.0, 0.0); probes.len()];
    let mut ext = Complex64::new(0.0, 0.0);
    let mut n = 0.0;
    let mut s = tube.state();
    while s.t < end - 0.5 * dt {
        tube.step(&mut s);
        if s.t >= start {
            let d = s.darwin.as_ref().expect("field");
            let b = d.b_z();
            let e = Complex64::new(libm::cos(omega * s.t), -libm::sin(omega * s.t));
            let applied = b0 * libm::sin(omega * s.t);
            ext += e * applied;
            for (k, r) in probes.iter().enumerate() {
                let mean: f64 = turns.iter().map(|(cx, cy)| d.b_at(&b, r * cx, r * cy)).sum::<f64>() / 4.0;
                sums[k] += e * (applied + mean);
            }
            n += 1.0;
        }
    }
    let ext = ext * (2.0 / n);
    for (k, r) in probes.iter().enumerate() {
        let measured = sums[k] * (2.0 / n) / ext;
        let exact = exact_b(&rings, sigma, omega, mu0, *r);
        println!(
            "SCRATCH fields r {r}: measured |B/B0| {:.4} arg {:.1} deg; exact {:.4} arg {:.1} deg",
            measured.norm(), measured.arg().to_degrees(), exact.norm(), exact.arg().to_degrees()
        );
    }
}

#[test]
#[ignore = "scratch"]
fn scratch_wide_gap() {
    dense("wide gap", &[(2.5, 4.0), (5.0, 6.0)], 0.25, 127, 47);
}
