//! Symmetries of the 3D flights (PHYSICS.md §9, Y1–Y6): the same scene transformed by a
//! symmetry of the equations of motion gives the transformed flight. A scene of three fixed
//! charges, a magnetic dipole and a current loop, all off-axis and tilted, a charge flown
//! through it for 12 time units at tolerance 1e-12, Newtonian (c = ∞) and relativistic
//! (c = 5). Run with `cargo test -p physics --test symmetries -- --nocapture`.

#![allow(clippy::disallowed_methods)] // references; the flights use libm (clippy.toml)
#![allow(clippy::float_cmp)] // bit-for-bit symmetries are compared exactly
#![allow(clippy::cast_precision_loss)] // small grid indices
#![allow(clippy::needless_range_loop)] // a matrix against its transpose

mod common;

use common::{Rng, cube};
use physics::DVec3;
use physics::antenna::OscillatingDipole;
use physics::beam::{BeamParticle, BeamScenario, Fates, run_beam};
use physics::bem::{BoxElectrode, Electrodes};
use physics::conductor::{Bias, Conductors, Resolution as SphereResolution, SphereConductor};
use physics::dynamics::Particle;
use physics::fdtd2d::{Fdtd, Grid};
use physics::field::{Coulomb, FixedCharge, LevelField};
use physics::geometry::{Shape, Sphere};
use physics::inductance::{around_circle, mutual_circles, mutual_polygon_circle, mutual_polygons};
use physics::magnetic::{CircularLoop, MagneticDipole, PolygonCoil};
use physics::trajectory::{RunSettings, Scenario, Trajectory, run};

const TOL: f64 = 1e-12;
const T_MAX: f64 = 12.0;

/// The scene's sources and the flown charge's start, as plain data to transform.
#[derive(Clone)]
struct Scene {
    charges: Vec<(DVec3, f64)>,
    dipoles: Vec<(DVec3, DVec3)>,
    /// Loops: centre, normal, strength κ.
    loops: Vec<(DVec3, DVec3, f64)>,
    q: f64,
    mass: f64,
    c: f64,
    x0: DVec3,
    p0: DVec3,
    t_max: f64,
}

fn scene(c: f64, with_loop: bool) -> Scene {
    Scene {
        charges: vec![
            (DVec3::new(3.0, 1.5, 0.5), 4.0),
            (DVec3::new(-2.5, 3.0, -1.0), -3.0),
            (DVec3::new(1.0, -3.5, 1.5), 2.5),
        ],
        dipoles: vec![(DVec3::new(-1.0, -1.5, 2.0), DVec3::new(0.7, -1.2, 2.5))],
        loops: if with_loop {
            vec![(
                DVec3::new(0.5, 0.75, -2.0),
                DVec3::new(0.2, 0.3, 1.0).normalize(),
                1.3,
            )]
        } else {
            Vec::new()
        },
        q: 1.0,
        mass: 1.0,
        c,
        x0: DVec3::new(-6.0, -0.5, 0.25),
        p0: DVec3::new(1.4, 0.3, -0.1),
        t_max: T_MAX,
    }
}

fn fly(s: &Scene) -> Trajectory {
    let field = LevelField {
        coulomb: Coulomb::new(
            &s.charges
                .iter()
                .map(|&(position, charge)| FixedCharge {
                    position,
                    charge,
                    radius: 0.3,
                })
                .collect::<Vec<_>>(),
        ),
        dipoles: s
            .dipoles
            .iter()
            .map(|&(position, moment)| MagneticDipole {
                position,
                moment,
                radius: 0.3,
            })
            .collect(),
        loops: s
            .loops
            .iter()
            .map(|&(center, normal, kappa)| CircularLoop {
                center,
                normal,
                radius: 2.5,
                kappa,
                wire_radius: 0.05,
                rate: 0.0,
            })
            .collect(),
        ..LevelField::default()
    };
    let scn = Scenario {
        field,
        obstacles: s
            .charges
            .iter()
            .map(|&(center, _)| {
                Shape::Sphere(Sphere {
                    center,
                    radius: 0.3,
                })
            })
            .collect(),
        particle: Particle {
            charge: s.q,
            mass: s.mass,
            radius: 0.0,
            moment: 0.0,
        },
        c: s.c,
        x0: s.x0,
        p0: s.p0,
        detector: None,
        bounds: Some(cube(60.0)),
        t_max: s.t_max,
        radiation_reaction: false,
        acceptance: None,
        gates: Vec::new(),
    };
    run(&scn, &RunSettings::with_tolerance(TOL))
}

/// A linear map of space with its determinant's sign: points and polar vectors (momenta)
/// map by it, axial vectors (magnetic moments, a loop's normal times its current) by it
/// times the sign (a mirror turns a current's sense).
#[derive(Clone, Copy)]
struct Map {
    m: [[f64; 3]; 3],
    shift: DVec3,
}

impl Map {
    fn apply(&self, v: DVec3) -> DVec3 {
        let r = |row: [f64; 3]| row[0] * v.x + row[1] * v.y + row[2] * v.z;
        DVec3::new(r(self.m[0]), r(self.m[1]), r(self.m[2]))
    }
    fn point(&self, v: DVec3) -> DVec3 {
        self.apply(v) + self.shift
    }
    fn det(&self) -> f64 {
        let m = self.m;
        m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
            - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
    }
    fn scene(&self, s: &Scene) -> Scene {
        let sign = self.det().signum();
        Scene {
            charges: s.charges.iter().map(|&(x, q)| (self.point(x), q)).collect(),
            dipoles: s
                .dipoles
                .iter()
                .map(|&(x, m)| (self.point(x), self.apply(m) * sign))
                .collect(),
            loops: s
                .loops
                .iter()
                .map(|&(x, n, k)| (self.point(x), self.apply(n), k * sign))
                .collect(),
            x0: self.point(s.x0),
            p0: self.apply(s.p0),
            ..s.clone()
        }
    }
}

/// The largest difference between a transformed flight's end and the flight's end
/// transformed, relative to the scene's size (positions) and momentum. The flight must
/// last its whole time (no collision cutting it short).
fn mismatch(map: Map, a: &Trajectory, b: &Trajectory) -> f64 {
    assert_eq!(a.outcome, b.outcome, "the outcomes differ");
    assert_eq!(
        a.outcome,
        physics::trajectory::Outcome::Timeout,
        "cut short at {}",
        a.end.t
    );
    let dx = (map.point(a.end.x) - b.end.x).length() / 10.0;
    let dp = (map.apply(a.end.p) - b.end.p).length() / a.end.p.length().max(1e-300);
    let dt = (a.end.t - b.end.t).abs() / a.end.t;
    dx.max(dp).max(dt)
}

const IDENTITY: [[f64; 3]; 3] = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];

/// Y1, translation: the scene moved by (3.37, −1.21, 0.58) flies the same path moved, to
/// 1e-8 of the scene's size (rounding of the shifted positions, amplified along the flight).
#[test]
fn y1_translation() {
    for c in [f64::INFINITY, 5.0] {
        let s = scene(c, true);
        let map = Map {
            m: IDENTITY,
            shift: DVec3::new(3.37, -1.21, 0.58),
        };
        let e = mismatch(map, &fly(&s), &fly(&map.scene(&s)));
        println!("Y1 translation, c = {c}: {e:.2e}");
        assert!(e < 1e-8, "c = {c}: {e}");
    }
}

/// Y2, rotation: a quarter turn about z (exact in floating point: a permutation with a
/// sign) gives the rotated flight bit for bit without the loop (whose field builds its own
/// basis), to 1e-8 with it; a general rotation (40° about (1, 2, 2)/3) to 1e-8.
#[test]
fn y2_rotation() {
    let quarter = Map {
        m: [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]],
        shift: DVec3::ZERO,
    };
    let (sin, cos) = (40.0_f64.to_radians().sin(), 40.0_f64.to_radians().cos());
    let n = DVec3::new(1.0, 2.0, 2.0) / 3.0;
    let k = |i: usize, j: usize| {
        let (ni, nj) = (n[i], n[j]);
        let delta = if i == j { 1.0 } else { 0.0 };
        let eps = match (i, j) {
            (0, 1) => -n.z,
            (1, 0) => n.z,
            (0, 2) => n.y,
            (2, 0) => -n.y,
            (1, 2) => -n.x,
            (2, 1) => n.x,
            _ => 0.0,
        };
        cos * delta + (1.0 - cos) * ni * nj + sin * eps
    };
    let general = Map {
        m: [
            [k(0, 0), k(0, 1), k(0, 2)],
            [k(1, 0), k(1, 1), k(1, 2)],
            [k(2, 0), k(2, 1), k(2, 2)],
        ],
        shift: DVec3::ZERO,
    };
    for c in [f64::INFINITY, 5.0] {
        let s = scene(c, false);
        let e = mismatch(quarter, &fly(&s), &fly(&quarter.scene(&s)));
        println!("Y2 quarter turn, no loop, c = {c}: {e:.2e}");
        assert!(e == 0.0, "c = {c}: {e}");
        let s = scene(c, true);
        for (name, map) in [("quarter turn", quarter), ("general rotation", general)] {
            let e = mismatch(map, &fly(&s), &fly(&map.scene(&s)));
            println!("Y2 {name}, with the loop, c = {c}: {e:.2e}");
            assert!(e < 1e-8, "{name}, c = {c}: {e}");
        }
    }
}

/// Y3, parity: the mirror x → −x, magnetic moments and currents as axial vectors (their
/// sense reversed), gives the mirrored flight: bit for bit without the loop, to 1e-8 with
/// it. A wrong treatment of B as a polar vector would turn the magnetic force the wrong way.
#[test]
fn y3_parity() {
    let mirror = Map {
        m: [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        shift: DVec3::ZERO,
    };
    for c in [f64::INFINITY, 5.0] {
        let s = scene(c, false);
        let e = mismatch(mirror, &fly(&s), &fly(&mirror.scene(&s)));
        println!("Y3 mirror, no loop, c = {c}: {e:.2e}");
        assert!(e == 0.0, "c = {c}: {e}");
        let s = scene(c, true);
        let e = mismatch(mirror, &fly(&s), &fly(&mirror.scene(&s)));
        println!("Y3 mirror, with the loop, c = {c}: {e:.2e}");
        assert!(e < 1e-8, "c = {c}: {e}");
    }
}

/// Y4, time reversal: flown for 12 time units, then from its end with the momentum and the
/// magnetic sources reversed (B is odd under T) for 12 more, the charge comes back to its
/// start with the opposite momentum, to 1e-7 of the scene's size (the adaptive integrator
/// is not reversible step by step: its local errors, ~1e-12 per step, add along the way
/// there and back).
#[test]
fn y4_time_reversal() {
    for c in [f64::INFINITY, 5.0] {
        let s = scene(c, true);
        let there = fly(&s);
        let back = Scene {
            dipoles: s.dipoles.iter().map(|&(x, m)| (x, -m)).collect(),
            loops: s.loops.iter().map(|&(x, n, k)| (x, n, -k)).collect(),
            x0: there.end.x,
            p0: -there.end.p,
            ..s.clone()
        };
        let back = fly(&back);
        let dx = (back.end.x - s.x0).length() / 10.0;
        let dp = (back.end.p + s.p0).length() / s.p0.length();
        println!("Y4 time reversal, c = {c}: position {dx:.2e}, momentum {dp:.2e}");
        assert!(dx.max(dp) < 1e-7, "c = {c}: {dx} {dp}");
    }
}

/// Y5, charge conjugation: every charge's sign and every magnetic source's sense reversed
/// (currents are moving charge) leaves the flight bit for bit the same (each force is a
/// product of two reversed signs).
#[test]
fn y5_charge_conjugation() {
    for c in [f64::INFINITY, 5.0] {
        let s = scene(c, true);
        let conj = Scene {
            charges: s.charges.iter().map(|&(x, q)| (x, -q)).collect(),
            dipoles: s.dipoles.iter().map(|&(x, m)| (x, -m)).collect(),
            loops: s.loops.iter().map(|&(x, n, k)| (x, n, -k)).collect(),
            q: -s.q,
            ..s.clone()
        };
        let (a, b) = (fly(&s), fly(&conj));
        let same = a.end.x == b.end.x && a.end.p == b.end.p && a.end.t == b.end.t;
        println!("Y5 charge conjugation, c = {c}: identical {same}");
        assert!(same, "c = {c}: {:?} against {:?}", a.end, b.end);
    }
}

/// Y6, scaling: the flown charge ×2, the fixed charges ×3, the mass ×1.5, speeds (and c)
/// ×2, the magnetic sources ×1.5 (a Coulomb force qQ/r² against m v²: 6/(1.5·4) = 1; a
/// magnetic one qvB against m v²/r: 2·2·1.5/(1.5·4) = 1) flies the same path in half the
/// time: to 1e-8.
#[test]
fn y6_scaling() {
    for c in [f64::INFINITY, 5.0] {
        let s = scene(c, true);
        let scaled = Scene {
            charges: s.charges.iter().map(|&(x, q)| (x, 3.0 * q)).collect(),
            dipoles: s.dipoles.iter().map(|&(x, m)| (x, 1.5 * m)).collect(),
            loops: s.loops.iter().map(|&(x, n, k)| (x, n, 1.5 * k)).collect(),
            q: 2.0 * s.q,
            mass: 1.5 * s.mass,
            c: 2.0 * s.c,
            p0: s.p0 * 3.0,
            t_max: 0.5 * s.t_max,
            ..s.clone()
        };
        let (a, b) = (fly(&s), fly(&scaled));
        assert_eq!(a.outcome, b.outcome);
        let dx = (a.end.x - b.end.x).length() / 10.0;
        let dp = (a.end.p * 3.0 - b.end.p).length() / (3.0 * a.end.p.length());
        println!("Y6 scaling, c = {c}: position {dx:.2e}, momentum {dp:.2e}");
        assert!(dx.max(dp) < 1e-8, "c = {c}: {dx} {dp}");
    }
}

// --- The field solver (`fdtd2d`): Y7–Y9 ------------------------------------------------------

/// Cells across the box of Y7–Y9.
const N: usize = 40;

/// A closed box N × N (h = 1, c = 1, 0.9 of the stability limit, no absorbing layer) with
/// an L of metal off its centre (`metal(i, j)` for node (i, j)), and B_z started as two
/// off-centre Gaussian bumps of opposite sign (E zero).
fn field_box(metal: &dyn Fn(usize, usize) -> bool, bump: &dyn Fn(f64, f64) -> f64) -> Fdtd {
    let grid = Grid {
        nx: N,
        ny: N,
        h: 1.0,
        origin: (0.0, 0.0),
    };
    let nodes = (0..grid.nodes())
        .map(|k| metal(k % (N + 1), k / (N + 1)))
        .collect();
    let mut f = Fdtd::new(grid, 1.0, 0.9, 0, nodes);
    for j in 0..N {
        for i in 0..N {
            #[allow(clippy::cast_precision_loss)]
            let (x, y) = (i as f64 + 0.5, j as f64 + 0.5);
            f.bz[i + N * j] = bump(x, y);
        }
    }
    f
}

/// The L of metal: nodes with 8 ≤ i ≤ 14 and 9 ≤ j ≤ 26, or 8 ≤ i ≤ 22 and 9 ≤ j ≤ 13.
fn l_metal(i: usize, j: usize) -> bool {
    ((8..=14).contains(&i) && (9..=26).contains(&j))
        || ((8..=22).contains(&i) && (9..=13).contains(&j))
}

fn bumps(x: f64, y: f64) -> f64 {
    let g = |x0: f64, y0: f64, w: f64| libm::exp(-((x - x0).powi(2) + (y - y0).powi(2)) / (w * w));
    g(28.3, 24.7, 2.5) - 0.6 * g(18.9, 31.2, 3.0)
}

/// The quarter turn (x, y) → (N − y, x) of a field: B_z unchanged in value (a proper
/// rotation), E turned: the new E_y on the y-edge (N − j, i) is the old E_x on the x-edge
/// (i, j), the new E_x on the x-edge (N − 1 − j, i) minus the old E_y on (i, j).
fn quarter(f: &Fdtd) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    let (mut ex, mut ey, mut bz) = (
        vec![0.0; f.ex.len()],
        vec![0.0; f.ey.len()],
        vec![0.0; f.bz.len()],
    );
    for j in 0..N {
        for i in 0..N {
            bz[(N - 1 - j) + N * i] = f.bz[i + N * j];
        }
    }
    for j in 0..=N {
        for i in 0..N {
            ey[(N - j) + (N + 1) * i] = f.ex[i + N * j];
        }
    }
    for j in 0..N {
        for i in 0..=N {
            ex[(N - 1 - j) + N * i] = -f.ey[i + (N + 1) * j];
        }
    }
    (ex, ey, bz)
}

/// The mirror x → N − x of a field: E_x reversed (a polar vector's x), E_y kept, B_z
/// reversed (an axial vector's z).
fn mirrored(f: &Fdtd) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    let (mut ex, mut ey, mut bz) = (
        vec![0.0; f.ex.len()],
        vec![0.0; f.ey.len()],
        vec![0.0; f.bz.len()],
    );
    for j in 0..N {
        for i in 0..N {
            bz[(N - 1 - i) + N * j] = -f.bz[i + N * j];
        }
    }
    for j in 0..=N {
        for i in 0..N {
            ex[(N - 1 - i) + N * j] = -f.ex[i + N * j];
        }
    }
    for j in 0..N {
        for i in 0..=N {
            ey[(N - i) + (N + 1) * j] = f.ey[i + (N + 1) * j];
        }
    }
    (ex, ey, bz)
}

fn max_diff(a: &(Vec<f64>, Vec<f64>, Vec<f64>), f: &Fdtd) -> f64 {
    a.0.iter()
        .zip(&f.ex)
        .chain(a.1.iter().zip(&f.ey))
        .chain(a.2.iter().zip(&f.bz))
        .map(|(u, v)| (u - v).abs())
        .fold(0.0, f64::max)
}

/// Y7, the field solver's quarter turn and mirror: a box with an L of metal and two bumps,
/// run 300 steps, against the same box turned (mirrored) and run: the turned (mirrored)
/// field bit for bit (the scheme's arithmetic is the same up to signs).
#[test]
fn y7_field_rotation_and_mirror() {
    let f0 = field_box(&l_metal, &bumps);
    let mut f = f0.clone();
    for _ in 0..300 {
        f.step();
    }
    // The turned box: its metal and its bumps turned with it.
    let turned_metal = |i: usize, j: usize| l_metal(j, N - i);
    let turned_bumps = |x: f64, y: f64| bumps(y, N as f64 - x);
    let mut g = field_box(&turned_metal, &turned_bumps);
    let start = quarter(&f0);
    let start_diff = max_diff(&start, &g);
    for _ in 0..300 {
        g.step();
    }
    let turn = max_diff(&quarter(&f), &g);
    let mirror_metal = |i: usize, j: usize| l_metal(N - i, j);
    let mirror_bumps = |x: f64, y: f64| -bumps(N as f64 - x, y);
    let mut m = field_box(&mirror_metal, &mirror_bumps);
    let mstart = max_diff(&mirrored(&f0), &m);
    for _ in 0..300 {
        m.step();
    }
    let mirror = max_diff(&mirrored(&f), &m);
    let scale =
        f.bz.iter()
            .chain(&f.ex)
            .map(|v| v.abs())
            .fold(0.0, f64::max);
    println!(
        "Y7 field: quarter turn {turn:.2e} (at the start {start_diff:.2e}), mirror {mirror:.2e} \
         (at the start {mstart:.2e}), of a largest field {scale:.3}"
    );
    assert!(
        start_diff == 0.0 && mstart == 0.0,
        "the transformed starts differ"
    );
    assert!(
        turn == 0.0 && mirror == 0.0,
        "quarter turn {turn}, mirror {mirror}"
    );
}

/// Y8, the field solver's time reversal: run 600 steps, then from E (at step 600) with B_z
/// (at step 600½) reversed, 600 more: back to the start (E zero, B_z half a step on,
/// reversed) to 1e-12 of the largest field (rounding, a few units in the last place per
/// step: the leapfrog is exactly reversible).
#[test]
fn y8_field_time_reversal() {
    let mut f = field_box(&l_metal, &bumps);
    let mut first = f.clone();
    first.step();
    let b_half = first.bz.clone();
    for _ in 0..600 {
        f.step();
    }
    let (ex, ey) = (f.ex.clone(), f.ey.clone());
    f.step();
    f.ex = ex;
    f.ey = ey;
    for b in &mut f.bz {
        *b = -*b;
    }
    for _ in 0..600 {
        f.step();
    }
    let scale = b_half.iter().map(|v| v.abs()).fold(0.0, f64::max);
    let de =
        f.ex.iter()
            .chain(&f.ey)
            .map(|v| v.abs())
            .fold(0.0, f64::max)
            / scale;
    let db =
        f.bz.iter()
            .zip(&b_half)
            .map(|(u, v)| (u + v).abs())
            .fold(0.0, f64::max)
            / scale;
    println!("Y8 field time reversal over 600 steps: E {de:.2e}, B {db:.2e} of the largest field");
    assert!(de.max(db) < 1e-12, "{de} {db}");
}

/// Y9, the field solver's energy in a closed lossless box: `energy()`, the scheme's own
/// invariant `Σ E·E (step n) + Σ B·B (steps n − ½ and n + ½)`, over 4000 steps (15 crossings
/// of the box), to 1e-12 of itself, for bumps 2.5, 5 and 10 cells wide; the plain sum
/// `Σ (E² + B²)` (E and B half a step apart) swings about it as energy passes between E and
/// B, by `(ω Δt)²`: the less resolved the field, the more (printed).
#[test]
fn y9_field_energy() {
    for width in [2.5, 5.0, 10.0] {
        let bump = |x: f64, y: f64| {
            libm::exp(-((x - 24.3).powi(2) + (y - 22.7).powi(2)) / (width * width))
        };
        let mut f = field_box(&l_metal, &bump);
        let w0 = f.energy();
        let plain = |f: &Fdtd| {
            let s: f64 = f.ex.iter().chain(&f.ey).chain(&f.bz).map(|v| v * v).sum();
            s / (8.0 * std::f64::consts::PI)
        };
        let (mut worst, mut swing) = (0.0_f64, (f64::INFINITY, f64::NEG_INFINITY));
        for _ in 0..4000 {
            f.step();
            worst = worst.max((f.energy() / w0 - 1.0).abs());
            let p = plain(&f) / w0;
            swing = (swing.0.min(p), swing.1.max(p));
        }
        println!(
            "Y9 field energy, a bump {width} cells wide, 4000 steps: energy() {worst:.2e}; \
             the plain sum within {:.4}..{:.4} of it",
            swing.0, swing.1
        );
        assert!(worst < 1e-12, "{width}: {worst}");
    }
}

// --- Electrostatics (`bem`) and inductances: Y12, Y13 ------------------------------------------

fn plate(center: DVec3, angle: f64, half: (f64, f64, f64)) -> BoxElectrode {
    BoxElectrode {
        center,
        angle,
        half_length: half.0,
        half_thickness: half.1,
        half_height: half.2,
        bias: Bias::Grounded,
    }
}

/// Y12, the electrodes' symmetries (boxes standing on the board's plane, the solver's
/// geometry: it uses their symmetry about the plane), at the verification's panels (0.35):
/// Maxwell's capacitance matrix of three unequal boxes, turned and off-centre, symmetric
/// (reciprocity, `C_ij = C_ji`) within 1e-4 of its largest entry (the collocation method is
/// not symmetric by construction; Z1 measured its reciprocation at 4e-5); the same when
/// everything moves by (2.37, −1.41) in the plane within 1e-9; two mirror images of a
/// grounded plate with a charge on the mirror plane take equal charges within 1e-4, along
/// the axes or tilted: the triangulation's diagonals all run one way, so a box's mesh is not
/// its own mirror image, and mirror images agree to the method's accuracy (Z1), not to
/// rounding.
#[test]
fn y12_electrode_symmetries() {
    let size = physics::bem::Resolution::Verify.panel_size();
    let three = |shift: DVec3| {
        vec![
            plate(DVec3::new(0.0, 0.0, 0.0) + shift, 0.3, (2.0, 0.25, 1.5)),
            plate(DVec3::new(4.0, 1.5, 0.0) + shift, 1.1, (1.5, 0.5, 1.0)),
            plate(DVec3::new(1.5, 4.5, 0.0) + shift, -0.6, (2.5, 0.25, 2.0)),
        ]
    };
    let matrix = |e: &Electrodes| e.capacitance().expect("electrodes").to_vec();
    let c = matrix(&Electrodes::with_panel_size(three(DVec3::ZERO), &[], size));
    let largest = c.iter().flatten().map(|v| v.abs()).fold(0.0, f64::max);
    let mut asym = 0.0_f64;
    for i in 0..3 {
        for j in 0..3 {
            asym = asym.max((c[i][j] - c[j][i]).abs() / largest);
        }
    }
    let moved = matrix(&Electrodes::with_panel_size(
        three(DVec3::new(2.37, -1.41, 0.0)),
        &[],
        size,
    ));
    let shift = c
        .iter()
        .flatten()
        .zip(moved.iter().flatten())
        .map(|(a, b)| (a - b).abs() / largest)
        .fold(0.0, f64::max);
    let pair = |angle: f64| {
        let e = Electrodes::with_panel_size(
            vec![
                plate(DVec3::new(-3.0, 0.0, 0.0), angle, (2.0, 0.25, 1.5)),
                plate(
                    DVec3::new(3.0, 0.0, 0.0),
                    std::f64::consts::PI - angle,
                    (2.0, 0.25, 1.5),
                ),
            ],
            &[(DVec3::new(0.0, 1.0, 0.0), 2.0)],
            size,
        );
        let (l, r) = (e.charge(0), e.charge(1));
        (l, r, (l - r).abs() / l.abs())
    };
    let (aligned, tilted) = (pair(std::f64::consts::FRAC_PI_2), pair(1.4));
    println!(
        "Y12 electrodes: C asymmetry {asym:.2e}, moved {shift:.2e}; mirror plates along y \
         {:.6} and {:.6} ({:.2e}), tilted {:.6} and {:.6} ({:.2e})",
        aligned.0, aligned.1, aligned.2, tilted.0, tilted.1, tilted.2
    );
    assert!(asym < 1e-4, "reciprocity {asym}");
    assert!(shift < 1e-9, "translation {shift}");
    assert!(aligned.2 < 1e-4, "mirror, aligned {}", aligned.2);
    assert!(tilted.2 < 1e-4, "mirror, tilted {}", tilted.2);
}

/// Y13, reciprocity of the mutual inductances (`inductance.rs`, Neumann's formula is
/// symmetric): two tilted circles, two polygons and a circle with a polygon, the flux of one
/// through the other against the reverse, within 1e-8 (the integrals converge to 1e-10).
#[test]
fn y13_inductance_reciprocity() {
    let a = CircularLoop {
        center: DVec3::new(0.0, 0.0, 0.0),
        normal: DVec3::new(0.1, 0.2, 1.0).normalize(),
        radius: 2.0,
        kappa: 1.0,
        wire_radius: 0.05,
        rate: 0.0,
    };
    let b = CircularLoop {
        center: DVec3::new(1.5, -0.5, 1.2),
        normal: DVec3::new(-0.4, 0.3, 1.0).normalize(),
        radius: 1.3,
        kappa: 1.0,
        wire_radius: 0.05,
        rate: 0.0,
    };
    let p = PolygonCoil {
        vertices: vec![
            DVec3::new(-3.0, -2.0, 0.5),
            DVec3::new(3.5, -2.5, 0.5),
            DVec3::new(2.0, 3.0, 0.8),
            DVec3::new(-2.5, 2.5, 0.2),
        ],
        kappa: 1.0,
        wire_radius: 0.05,
        rate: 0.0,
    };
    let q = PolygonCoil {
        vertices: vec![
            DVec3::new(-1.0, -1.0, 2.0),
            DVec3::new(1.5, -0.8, 2.3),
            DVec3::new(0.2, 1.7, 1.9),
        ],
        kappa: 1.0,
        wire_radius: 0.05,
        rate: 0.0,
    };
    let rel = |x: f64, y: f64| (x - y).abs() / x.abs().max(y.abs());
    let circles = rel(mutual_circles(&a, &b), mutual_circles(&b, &a));
    let polygons = rel(mutual_polygons(&p, &q), mutual_polygons(&q, &p));
    let through_circle = around_circle(
        |x| p.unit_vector_potential(x),
        a.center,
        a.normal,
        a.radius,
        1 << 14,
    );
    let mixed = rel(mutual_polygon_circle(&p, &a), through_circle);
    println!(
        "Y13 inductances: circles {circles:.2e}, polygons {polygons:.2e}, circle-polygon {mixed:.2e}"
    );
    assert!(
        circles < 1e-8 && polygons < 1e-8 && mixed < 1e-8,
        "{circles} {polygons} {mixed}"
    );
}

// --- The many-particle runner (`beam`): Y14, Y15 -----------------------------------------------

/// B2's Coulomb explosion: 12 like charges (0.3), masses 0.5–1.5, in a 4 × 4 square with small
/// random velocities, Newtonian (c = ∞), interacting, flown for 15.
fn explosion(boost: DVec3, reversed: Option<&[(DVec3, DVec3)]>) -> Vec<(DVec3, DVec3, f64)> {
    let mut rng = Rng::new(7);
    let particles: Vec<BeamParticle> = (0..12)
        .map(|k| {
            let x = DVec3::new(rng.next_f64() * 4.0 - 2.0, rng.next_f64() * 4.0 - 2.0, 0.0);
            let v = DVec3::new(rng.next_f64() - 0.5, rng.next_f64() - 0.5, 0.0) * 0.2;
            let m = 0.5 + rng.next_f64();
            let (x0, p0) = match reversed {
                Some(ends) => (ends[k].0, -ends[k].1),
                None => (x, (v + boost) * m),
            };
            BeamParticle {
                particle: Particle {
                    charge: 0.3,
                    mass: m,
                    radius: 0.0,
                    moment: 0.0,
                },
                x0,
                p0,
                detector: None,
                acceptance: None,
            }
        })
        .collect();
    let masses: Vec<f64> = particles.iter().map(|b| b.particle.mass).collect();
    let scn = BeamScenario {
        field: Coulomb::new(&[]),
        obstacles: vec![],
        particles,
        c: f64::INFINITY,
        bounds: None,
        t_max: 15.0,
        interact: true,
        gates: Vec::new(),
        radiation_reaction: false,
        retarded: false,
        fates: Fates::default(),
    };
    let r = run_beam(&scn, &RunSettings::with_tolerance(TOL));
    r.trajectories
        .iter()
        .zip(masses)
        .map(|(t, m)| (t.end.x, t.end.p, m))
        .collect()
}

/// Y14, Galilean invariance of the interacting runner (c = ∞): every particle given the same
/// extra velocity u = (0.37, −0.21, 0.05) ends moved by u t with its momentum shifted by
/// m u, within 1e-9 of the spread (the steps differ: the moving frame's positions round
/// differently).
#[test]
fn y14_galilean_boost() {
    let u = DVec3::new(0.37, -0.21, 0.05);
    let (rest, moving) = (explosion(DVec3::ZERO, None), explosion(u, None));
    let size = rest.iter().map(|e| e.0.length()).fold(0.0, f64::max);
    let mut worst = 0.0_f64;
    for ((x, p, m), (xb, pb, _)) in rest.iter().zip(&moving) {
        let dx = (*xb - (*x + u * 15.0)).length() / size;
        let dp = (*pb - (*p + u * *m)).length() / p.length();
        worst = worst.max(dx).max(dp);
    }
    println!("Y14 Galilean boost of 12 interacting charges: {worst:.2e}");
    assert!(worst < 1e-9, "{worst}");
}

/// Y15, time reversal of the interacting runner: after 15 the 12 charges' momenta reversed,
/// 15 more bring each back to its start with the opposite momentum, within 1e-7 of the
/// spread (the adaptive steps are not reversible one by one).
#[test]
fn y15_many_body_time_reversal() {
    let there = explosion(DVec3::ZERO, None);
    let ends: Vec<(DVec3, DVec3)> = there.iter().map(|e| (e.0, e.1)).collect();
    let back = explosion(DVec3::ZERO, Some(&ends));
    // The starts, drawn again as `explosion` drew them.
    let mut rng = Rng::new(7);
    let starts: Vec<(DVec3, DVec3)> = (0..12)
        .map(|_| {
            let x0 = DVec3::new(rng.next_f64() * 4.0 - 2.0, rng.next_f64() * 4.0 - 2.0, 0.0);
            let v0 = DVec3::new(rng.next_f64() - 0.5, rng.next_f64() - 0.5, 0.0) * 0.2;
            let m = 0.5 + rng.next_f64();
            (x0, v0 * m)
        })
        .collect();
    let size = there.iter().map(|e| e.0.length()).fold(0.0, f64::max);
    let p_scale = starts.iter().map(|s| s.1.length()).fold(0.0, f64::max);
    let mut worst = 0.0_f64;
    for ((x, p, _), (x0, p0)) in back.iter().zip(&starts) {
        worst = worst
            .max((*x - *x0).length() / size)
            .max((*p + *p0).length() / p_scale);
    }
    println!("Y15 time reversal of 12 interacting charges over 15 and back: {worst:.2e}");
    assert!(worst < 1e-7, "{worst}");
}

// --- Metal spheres (image series, `conductor`): Y16 ----------------------------------------------

/// The net charge the fixed sources induce on each sphere: the images inside it.
fn sphere_charges(c: &Conductors) -> Vec<f64> {
    c.spheres
        .iter()
        .map(|s| {
            c.induced
                .charges()
                .filter(|(x, _)| (*x - s.center).length() < s.radius)
                .map(|(_, q)| q)
                .sum()
        })
        .collect()
}

/// Y16, the metal spheres' symmetries (the image series, verification's resolution): two
/// equal grounded spheres mirroring each other about x = 0 with a charge on the mirror plane
/// take equal charges within 1e-10, and three unequal spheres with two charges induce the
/// same charges when everything moves by (2.37, −1.41, 0) within 1e-10 (rounding: the image
/// series' geometry is the same).
#[test]
fn y16_sphere_symmetries() {
    let sphere = |x: f64, y: f64, r: f64| SphereConductor {
        center: DVec3::new(x, y, 0.0),
        radius: r,
        bias: Bias::Grounded,
    };
    let pair = Conductors::new(
        vec![sphere(-3.0, 0.5, 1.2), sphere(3.0, 0.5, 1.2)],
        &[(DVec3::new(0.0, 1.7, 0.4), 2.0)],
        SphereResolution::Verify,
    );
    let q = sphere_charges(&pair);
    let mirror = (q[0] - q[1]).abs() / q[0].abs();
    let three = |shift: DVec3| {
        let mut spheres = vec![
            sphere(0.0, 0.0, 1.0),
            sphere(4.0, 1.5, 0.7),
            sphere(1.0, 4.0, 1.3),
        ];
        for s in &mut spheres {
            s.center += shift;
        }
        let sources = [
            (DVec3::new(2.0, 2.0, 0.3) + shift, 1.5),
            (DVec3::new(-2.0, 3.0, -0.5) + shift, -0.8),
        ];
        sphere_charges(&Conductors::new(
            spheres,
            &sources,
            SphereResolution::Verify,
        ))
    };
    let (here, there) = (three(DVec3::ZERO), three(DVec3::new(2.37, -1.41, 0.0)));
    let largest = here.iter().map(|v| v.abs()).fold(0.0, f64::max);
    let moved = here
        .iter()
        .zip(&there)
        .map(|(a, b)| (a - b).abs() / largest)
        .fold(0.0, f64::max);
    println!(
        "Y16 spheres: mirror pair {:.8} and {:.8} ({mirror:.2e}); moved {moved:.2e}",
        q[0], q[1]
    );
    assert!(mirror < 1e-10, "mirror {mirror}");
    assert!(moved < 1e-10, "translation {moved}");
}

// --- Retarded fields (`antenna`): Y17 ----------------------------------------------------------

/// Y17, the antenna's exact retarded fields (near, intermediate and radiation zones, c = 3)
/// under the mirror x → −x and a quarter turn about z: the transformed dipole (its moment a
/// polar vector) gives at the transformed point the transformed field, E as a polar vector
/// and B as an axial one (the mirror reverses B's y and z), bit for bit at 64 points and
/// times. (Time reversal is not a symmetry of these fields: they are retarded, radiation
/// leaving its source; the advanced solution would be its mirror in time.)
#[test]
fn y17_antenna_fields() {
    let dipole = OscillatingDipole {
        position: DVec3::new(0.7, -0.4, 0.3),
        amplitude: DVec3::new(0.6, 1.1, -0.4),
        omega: 1.7,
        phase: 0.4,
        c: 3.0,
        radius: 0.3,
    };
    let mirror = |v: DVec3| DVec3::new(-v.x, v.y, v.z);
    let axial_mirror = |v: DVec3| DVec3::new(v.x, -v.y, -v.z);
    let turn = |v: DVec3| DVec3::new(-v.y, v.x, v.z);
    let mirrored = OscillatingDipole {
        position: mirror(dipole.position),
        amplitude: mirror(dipole.amplitude),
        ..dipole
    };
    let turned = OscillatingDipole {
        position: turn(dipole.position),
        amplitude: turn(dipole.amplitude),
        ..dipole
    };
    let mut rng = Rng::new(17);
    let (mut worst_m, mut worst_t) = (0.0_f64, 0.0_f64);
    for _ in 0..64 {
        let x = DVec3::new(
            rng.range(-12.0, 12.0),
            rng.range(-12.0, 12.0),
            rng.range(-6.0, 6.0),
        );
        let t = rng.range(0.0, 20.0);
        let f = dipole.fields(x, t);
        let fm = mirrored.fields(mirror(x), t);
        let ft = turned.fields(turn(x), t);
        worst_m = worst_m
            .max((fm.e - mirror(f.e)).length())
            .max((fm.b - axial_mirror(f.b)).length());
        worst_t = worst_t
            .max((ft.e - turn(f.e)).length())
            .max((ft.b - turn(f.b)).length());
    }
    println!("Y17 antenna fields: mirror {worst_m:.2e}, quarter turn {worst_t:.2e}");
    assert!(
        worst_m == 0.0 && worst_t == 0.0,
        "mirror {worst_m}, quarter turn {worst_t}"
    );
}

/// A Lorentz boost by `u` (light's speed `c`; the force `q(E + v × B)`): events, a mass `m`'s
/// momenta and the fields seen from the frame moving at `u`.
struct Boost {
    u: DVec3,
    c: f64,
}

impl Boost {
    fn gamma(&self) -> f64 {
        1.0 / (1.0 - self.u.length_squared() / (self.c * self.c)).sqrt()
    }
    fn unit(&self) -> DVec3 {
        self.u.normalize_or_zero()
    }
    fn event(&self, t: f64, x: DVec3) -> (f64, DVec3) {
        let (g, n) = (self.gamma(), self.unit());
        (
            g * (t - self.u.dot(x) / (self.c * self.c)),
            x + n * ((g - 1.0) * x.dot(n)) - self.u * (g * t),
        )
    }
    fn momentum(&self, p: DVec3, m: f64) -> DVec3 {
        let (g, n, c) = (self.gamma(), self.unit(), self.c);
        let energy = (m * m * c.powi(4) + p.length_squared() * c * c).sqrt();
        p + n * ((g - 1.0) * p.dot(n)) - self.u * (g * energy / (c * c))
    }
    /// A plane wave seen from the moving frame: its phase is invariant, its wave four-vector
    /// `(ω/c, k)` transforms as an event, its amplitude and polarization as its fields at any
    /// event.
    fn wave(&self, w: &physics::external::PlaneWave) -> physics::external::PlaneWave {
        let (g, n, c) = (self.gamma(), self.unit(), self.c);
        let k = w.direction * (w.omega / c);
        let omega = g * (w.omega - self.u.dot(k));
        let k2 = k + n * ((g - 1.0) * k.dot(n)) - self.u * (g * w.omega / (c * c));
        let e = w.polarization * w.amplitude;
        let (e2, _) = self.fields(e, w.direction.cross(e) / c);
        physics::external::PlaneWave {
            amplitude: e2.length(),
            direction: k2.normalize(),
            polarization: e2.normalize(),
            omega,
            phase: w.phase,
            c,
        }
    }
    fn fields(&self, e: DVec3, b: DVec3) -> (DVec3, DVec3) {
        let (g, n, c) = (self.gamma(), self.unit(), self.c);
        (
            (e + self.u.cross(b)) * g - n * ((g - 1.0) * e.dot(n)),
            (b - self.u.cross(e) / (c * c)) * g - n * ((g - 1.0) * b.dot(n)),
        )
    }
}

/// A charge (q = 1, m = 1) flown in uniform fields from `x0`, `p0` for `t_max` (c = `c`).
fn fly_uniform(e: DVec3, b: DVec3, c: f64, x0: DVec3, p0: DVec3, t_max: f64) -> Trajectory {
    fly_external(
        vec![physics::external::External::Uniform { e, b }],
        c,
        x0,
        p0,
        t_max,
    )
}

/// A charge (q = 1, m = 1) flown in external fields from `x0`, `p0` for `t_max`.
fn fly_external(
    external: Vec<physics::external::External>,
    c: f64,
    x0: DVec3,
    p0: DVec3,
    t_max: f64,
) -> Trajectory {
    let scn = Scenario {
        field: LevelField {
            external,
            ..LevelField::default()
        },
        obstacles: Vec::new(),
        particle: Particle {
            charge: 1.0,
            mass: 1.0,
            radius: 0.0,
            moment: 0.0,
        },
        c,
        x0,
        p0,
        detector: None,
        bounds: Some(cube(1e4)),
        t_max,
        radiation_reaction: false,
        acceptance: None,
        gates: Vec::new(),
    };
    run(&scn, &RunSettings::with_tolerance(TOL))
}

/// Y20, Lorentz invariance of the relativistic engine (PHYSICS.md §9): a charge in uniform
/// crossed fields (E = 2 ŷ, B = ẑ, c = 5), flown for 12 from the origin with p₀ = (0.6, 1.1,
/// 0.4), against the same charge flown in the fields seen from a moving frame (`Boost`) from
/// the boosted start for the boosted time: the boosted end event and momentum, to 1e-8 of
/// the path's size and of the momentum (the integrator's tolerance 1e-12, amplified along the
/// flight; Y1's bound). Boosts: the E × B drift frame (u = 2 x̂, β = 0.4), where E′ vanishes
/// and B′ = B/γ, so that the boosted flight is also a gyration in closed form (checked to the
/// same bound); and an oblique one, u = (1.2, −0.8, 0.9). The same oblique boost of a plane
/// wave (ω = 1.3, amplitude 1.5, along (1, 0.3, 0)) over a uniform B = 0.8 ẑ: the wave seen
/// Doppler-shifted (`Boost::wave`), its B′ that of a plane wave (`k̂′ × E′/c`, to 1e-12), the
/// flights the same to 1e-8.
#[test]
fn y20_lorentz_boost() {
    let (e, b, c) = (DVec3::new(0.0, 2.0, 0.0), DVec3::new(0.0, 0.0, 1.0), 5.0);
    let (x0, p0, t_max) = (DVec3::ZERO, DVec3::new(0.6, 1.1, 0.4), 12.0);
    let lab = fly_uniform(e, b, c, x0, p0, t_max);
    assert_eq!(lab.outcome, physics::trajectory::Outcome::Timeout);
    for u in [DVec3::new(2.0, 0.0, 0.0), DVec3::new(1.2, -0.8, 0.9)] {
        let boost = Boost { u, c };
        let (e2, b2) = boost.fields(e, b);
        let (t0, y0) = boost.event(0.0, x0);
        let (t1, y1) = boost.event(lab.end.t, lab.end.x);
        let (q0, q1) = (boost.momentum(p0, 1.0), boost.momentum(lab.end.p, 1.0));
        let moving = fly_uniform(e2, b2, c, y0, q0, t1 - t0);
        assert_eq!(moving.outcome, physics::trajectory::Outcome::Timeout);
        let size = (y1 - y0).length().max(1.0);
        let dx = (moving.end.x - y1).length() / size;
        let dp = (moving.end.p - q1).length() / q1.length();
        println!(
            "Y20 boost u = {u:?} (γ = {:.4}): E′ = {e2:.4?}, B′ = {b2:.4?}; the end {dx:.2e} of \
             the path, the momentum {dp:.2e}",
            boost.gamma()
        );
        assert!(dx < 1e-8 && dp < 1e-8, "{u:?}: {dx} {dp}");
        if u.y == 0.0 && u.z == 0.0 {
            // The drift frame: a gyration about z at ω = B′/(γ_p m), the momentum turning.
            assert!(e2.length() < 1e-14, "E′ {e2:?}");
            let gamma_p = (1.0 + q0.length_squared() / (c * c)).sqrt();
            let omega = b2.z / gamma_p;
            let dt = t1 - t0;
            let (sn, cs) = ((-omega * dt).sin(), (-omega * dt).cos());
            let turned = DVec3::new(q0.x * cs - q0.y * sn, q0.x * sn + q0.y * cs, q0.z);
            // The position: the guiding centre plus the turned radius vector.
            let centre = y0 + DVec3::new(q0.y, -q0.x, 0.0) / b2.z;
            let at = centre - DVec3::new(turned.y, -turned.x, 0.0) / b2.z
                + DVec3::new(0.0, 0.0, q0.z / gamma_p * dt);
            let (ax, ap) = (
                (moving.end.x - at).length() / size,
                (moving.end.p - turned).length() / q1.length(),
            );
            println!(
                "Y20 the drift frame's gyration in closed form: the end {ax:.2e}, the momentum {ap:.2e}"
            );
            assert!(ax < 1e-8 && ap < 1e-8, "gyration {ax} {ap}");
        }
    }
    // A plane wave over a uniform field.
    let wave = physics::external::PlaneWave {
        amplitude: 1.5,
        direction: DVec3::new(1.0, 0.3, 0.0).normalize(),
        polarization: DVec3::new(-0.3, 1.0, 0.4).normalize(),
        omega: 1.3,
        phase: 0.4,
        c,
    };
    // The polarization made perpendicular to the direction.
    let wave = physics::external::PlaneWave {
        polarization: (wave.polarization - wave.direction * wave.polarization.dot(wave.direction))
            .normalize(),
        ..wave
    };
    let b_uniform = DVec3::new(0.0, 0.0, 0.8);
    let lab = fly_external(
        vec![
            physics::external::External::Wave(wave),
            physics::external::External::Uniform {
                e: DVec3::ZERO,
                b: b_uniform,
            },
        ],
        c,
        x0,
        p0,
        t_max,
    );
    assert_eq!(lab.outcome, physics::trajectory::Outcome::Timeout);
    let boost = Boost {
        u: DVec3::new(1.2, -0.8, 0.9),
        c,
    };
    let seen = boost.wave(&wave);
    // Its fields at an event, both ways: transformed, and the boosted wave's own.
    let (t, x) = (3.7, DVec3::new(1.0, -2.0, 0.5));
    let field = |w: &physics::external::PlaneWave, t: f64, x: DVec3| {
        let e = w.polarization * (w.amplitude * w.phase_at(x, t).cos());
        (e, w.direction.cross(e) / w.c)
    };
    let (e, b) = field(&wave, t, x);
    let (e_t, b_t) = boost.fields(e, b);
    let (t2, x2) = boost.event(t, x);
    let (e_w, b_w) = field(&seen, t2, x2);
    let consistency = (e_t - e_w).length().max((b_t - b_w).length() * c) / wave.amplitude;
    let (_, b2) = boost.fields(DVec3::ZERO, b_uniform);
    let (e2u, _) = boost.fields(DVec3::ZERO, b_uniform);
    let (t0, y0) = boost.event(0.0, x0);
    let (t1, y1) = boost.event(lab.end.t, lab.end.x);
    let (q0, q1) = (boost.momentum(p0, 1.0), boost.momentum(lab.end.p, 1.0));
    let moving = fly_external(
        vec![
            physics::external::External::Wave(seen),
            physics::external::External::Uniform { e: e2u, b: b2 },
        ],
        c,
        y0,
        q0,
        t1 - t0,
    );
    assert_eq!(moving.outcome, physics::trajectory::Outcome::Timeout);
    let size = (y1 - y0).length().max(1.0);
    let dx = (moving.end.x - y1).length() / size;
    let dp = (moving.end.p - q1).length() / q1.length();
    println!(
        "Y20 a plane wave boosted: ω {} → {:.5}, its fields at an event both ways within \
         {consistency:.2e}; the end {dx:.2e} of the path, the momentum {dp:.2e}",
        wave.omega, seen.omega
    );
    assert!(consistency < 1e-12, "the boosted wave {consistency}");
    assert!(dx < 1e-8 && dp < 1e-8, "wave: {dx} {dp}");
}
