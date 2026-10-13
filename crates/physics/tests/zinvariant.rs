//! Validation tests V1–V3 of the z-invariant world's electrostatics (`zinv.rs`, PHYSICS.md
//! §2.11; the tubes' groundwork, docs/TUBES.md). Run with
//! `cargo test --release -p physics --test zinvariant -- --nocapture --test-threads=1`.
//!
//! Only enclosed setups, whose results do not depend on the logarithm's reference.

#![allow(clippy::disallowed_methods)] // references; the engine uses libm (clippy.toml)

use physics::DVec3;
use physics::bem::Bias;
use physics::zinv::{Electrode, Electrodes, Section, segment_integrals};

/// V1: the segment integrals against Gauss–Legendre quadrature of `−ln r²` and
/// `2 (x − s)/r²` along the segment (2000 panels of 8 points, far from the segment and
/// near it, on its line beyond an end, and on it at its midpoint, where the field's
/// normal part is the average of its two sides, 0). Required: 1e-11 relative (the
/// quadrature's accuracy at these distances).
#[test]
fn v1_segment_integrals_against_quadrature() {
    let (a, b) = (DVec3::new(0.3, -0.2, 0.0), DVec3::new(1.9, 0.7, 0.0));
    let nodes = [
        (0.183_434_642_495_649_8, 0.362_683_783_378_362),
        (0.525_532_409_916_329, 0.313_706_645_877_887_3),
        (0.796_666_477_413_626_7, 0.222_381_034_453_374_5),
        (0.960_289_856_497_536_3, 0.101_228_536_290_376_3),
    ];
    let quad = |x: DVec3| {
        let panels = 2000;
        let (mut phi, mut e) = (0.0, DVec3::ZERO);
        let d = (b - a) / f64::from(panels);
        let h = d.length();
        for k in 0..panels {
            let mid = a + d * (f64::from(k) + 0.5);
            for &(t, w) in &nodes {
                for s in [mid + d * (0.5 * t), mid - d * (0.5 * t)] {
                    let r = x - s;
                    let r2 = r.x * r.x + r.y * r.y;
                    phi += -0.5 * w * h * r2.ln();
                    e += r * (0.5 * w * h * 2.0 / r2);
                }
            }
        }
        (phi, e)
    };
    let mut worst: f64 = 0.0;
    for x in [
        DVec3::new(4.0, 3.0, 0.0),
        DVec3::new(1.0, 0.6, 0.0),
        DVec3::new(-1.0, -0.4, 0.0) + (b - a) * 0.0 - (b - a) * 0.5,
    ] {
        let (p1, e1) = segment_integrals(a, b, x);
        let (p2, e2) = quad(x);
        let err = ((p1 - p2) / p2).abs().max((e1 - e2).length() / e2.length());
        println!("V1 at {x}: φ {p1:.14} (quadrature {p2:.14}), E {e1} ({e2}): {err:.1e}");
        worst = worst.max(err);
    }
    // On the segment, at its midpoint: φ = −L(ln(L²/4) − 2), E_normal = 0.
    let mid = (a + b) * 0.5;
    let l = (b - a).length();
    let (p, e) = segment_integrals(a, b, mid);
    let want = -l * ((l * l / 4.0).ln() - 2.0);
    let n = DVec3::new(-(b - a).y, (b - a).x, 0.0) / l;
    println!(
        "V1 on the midpoint: φ {p:.14} (closed form {want:.14}), E·n {:.1e}",
        e.dot(n)
    );
    worst = worst.max(((p - want) / want).abs()).max(e.dot(n).abs());
    assert!(worst < 1e-11, "V1 {worst:.3e}");
}

/// V2: coaxial cylinders (radii 1 and 3; the inner at potential 1, the outer grounded):
/// the capacitance per unit length `1/(2 ln(b/a))` (k = 1; Jackson §1.11), with polygons
/// of segments 0.2, 0.1, 0.05 long. Required: the error shrinking with every refinement
/// at an observed order of at least 1 (inscribed polygons: O(h²) expected), and the two
/// cylinders' net charge (0 in the continuum: the outer one encloses the inner, and a net
/// charge would put the grounded outer one at `−2Q ln 3`, not 0) shrinking with every
/// refinement too (the collocation enforces it only approximately; a first version of
/// this test expected it exact).
#[test]
fn v2_coaxial_capacitance() {
    let exact = 1.0 / (2.0 * 3.0_f64.ln());
    let mut errors = Vec::new();
    let mut nets: Vec<f64> = Vec::new();
    for size in [0.2, 0.1, 0.05] {
        let e = Electrodes::new(
            vec![
                Electrode {
                    section: Section::Circle {
                        center: DVec3::ZERO,
                        radius: 1.0,
                    },
                    bias: Bias::Potential(1.0),
                },
                Electrode {
                    section: Section::Circle {
                        center: DVec3::ZERO,
                        radius: 3.0,
                    },
                    bias: Bias::Grounded,
                },
            ],
            &[],
            size,
        );
        let (q_in, q_out) = (e.charge(0), e.charge(1));
        let err = q_in / exact - 1.0;
        println!(
            "V2 segments {size}: {} of them, C' = {q_in:.10} (exact {exact:.10}), {err:.2e}; outer {q_out:.10}",
            e.segments.len()
        );
        nets.push(((q_in + q_out) / q_in).abs());
        errors.push(err);
    }
    println!(
        "V2 net charge relative {:.2e} {:.2e} {:.2e}",
        nets[0], nets[1], nets[2]
    );
    assert!(nets[0] > nets[1] && nets[1] > nets[2], "V2 net {nets:?}");
    let order = (errors[1].abs() / errors[2].abs()).ln() / 2f64.ln();
    println!("V2 observed order {order:.2}");
    assert!(
        errors[0].abs() > errors[1].abs() && errors[1].abs() > errors[2].abs() && order >= 1.0,
        "V2 {errors:?}"
    );
}

/// V3: a line charge λ = 1 inside a grounded cylinder (radius 2), 0.8 off its axis: the
/// wall's induced charge is the image −λ at R²/e on the same ray, so the force per unit
/// length on the line charge is `2λ²/(R²/e − e)` towards the wall (exact). Required as V2.
#[test]
fn v3_line_charge_in_a_grounded_cylinder() {
    let (r, e) = (2.0, 0.8);
    let exact = 2.0 / (r * r / e - e);
    let x = DVec3::new(e, 0.0, 0.0);
    let mut errors = Vec::new();
    for size in [0.2, 0.1, 0.05] {
        let wall = Electrodes::new(
            vec![Electrode {
                section: Section::Circle {
                    center: DVec3::ZERO,
                    radius: r,
                },
                bias: Bias::Grounded,
            }],
            &[(x, 1.0)],
            size,
        );
        let f = wall.sample(x).1.x;
        let err = f / exact - 1.0;
        println!(
            "V3 segments {size}: force {f:.10} (image {exact:.10}), {err:.2e}; wall charge {:.10}",
            wall.charge(0)
        );
        errors.push(err);
    }
    let order = (errors[1].abs() / errors[2].abs()).ln() / 2f64.ln();
    println!("V3 observed order {order:.2}");
    assert!(
        errors[0].abs() > errors[1].abs() && errors[1].abs() > errors[2].abs() && order >= 1.0,
        "V3 {errors:?}"
    );
}

/// V4: the space-charge-limited coaxial diode (cathode radius 1, grounded; anode radius
/// 2 at potential 1; electrons, q/m = −1) against Langmuir & Blodgett: the current per
/// unit length `(2/9) √(2|q/m|) V^{3/2} / (b β²)` with β² = 0.27926716071059 at b/a = 2
/// (`scripts/wolfram/v4_langmuir_blodgett.wls`, the Langmuir–Blodgett equation solved at
/// 30 digits; its series and Langmuir & Blodgett's table agree). Measured as the anode's
/// collected charge per time over t = 4 to 8, after the flow has settled (two transit
/// times), at a base resolution and a refined one (half the macroparticle weight, the
/// step, the segments and the softening). The release half a segment out (tube.rs,
/// `EMISSION_OFFSET`) makes the error first order in the resolution: required, the error
/// shrinking at an observed order of at least 0.7 and the Richardson extrapolation
/// (2 I_fine − I_coarse) within 2 % (measured 15 %, 7.8 %, extrapolated 0.5 %). The 2 %:
/// the fine run collects ~2400 macroparticles in the window, so a Poisson bound on its
/// counting noise is 1/√2400 = 2 % (the emission is deterministic, so the real noise is
/// smaller).
#[test]
fn v4_coaxial_diode_langmuir_blodgett() {
    use physics::tube::Tube;
    let beta2 = 0.279_267_160_710_59;
    let exact = 2.0 / 9.0 * 2f64.sqrt() / (2.0 * beta2);
    let mut currents = Vec::new();
    for (weight, dt, size) in [(2e-3, 0.01, 0.1), (1e-3, 0.005, 0.05)] {
        let electrodes = Electrodes::new(
            vec![
                Electrode {
                    section: Section::Circle {
                        center: DVec3::ZERO,
                        radius: 1.0,
                    },
                    bias: Bias::Grounded,
                },
                Electrode {
                    section: Section::Circle {
                        center: DVec3::ZERO,
                        radius: 2.0,
                    },
                    bias: Bias::Potential(1.0),
                },
            ],
            &[],
            size,
        );
        let tube = Tube {
            electrodes,
            emitters: vec![physics::tube::Emitter {
                electrode: 0,
                charge_per_mass: -1.0,
                weight,
                emit_toward: None,
                saturation: None,
            }],
            projectiles: Vec::new(),
            softening: 0.5 * size,
            dt,
            arena: None,
            b_z: 0.0,
            lattices: Vec::new(),
            circuits: Vec::new(),
        };
        let mut s = tube.state();
        let start = std::time::Instant::now();
        let mut at4 = None;
        while s.t < 8.0 {
            tube.step(&mut s);
            if at4.is_none() && s.t >= 4.0 {
                at4 = Some((s.t, s.collected[1][0]));
            }
        }
        let (t4, q4) = at4.expect("reached t = 4");
        let current = (s.collected[1][0] - q4) / (s.t - t4);
        println!(
            "V4 weight {weight}, dt {dt}, segments {size}: I/L = {:.6} (Langmuir–Blodgett {exact:.6}), {:.2e}; {} particles in flight, {:.1} s",
            -current,
            -current / exact - 1.0,
            s.x.len(),
            start.elapsed().as_secs_f64()
        );
        currents.push(-current);
    }
    let (e0, e1) = (currents[0] / exact - 1.0, currents[1] / exact - 1.0);
    let order = (e0.abs() / e1.abs()).ln() / 2f64.ln();
    let extrapolated = 2.0 * currents[1] - currents[0];
    println!(
        "V4 observed order {order:.2}; extrapolated {extrapolated:.6}, {:.2e}",
        extrapolated / exact - 1.0
    );
    assert!(e0.abs() > e1.abs() && order >= 0.7, "V4 {currents:?}");
    assert!(
        (extrapolated / exact - 1.0).abs() < 0.02,
        "V4 extrapolated {extrapolated}"
    );
}

/// V5: Ramo's theorem. A line charge λ = 1 between grounded coaxial cylinders (radii 1 and
/// 3) at radius r induces `Q_a = −λ ln(b/r)/ln(b/a)` on the inner one (exact, from the
/// potential's logarithmic profile), so moving outward at speed v it drives the current
/// `dQ_a/dt = λ v/(r ln(b/a))` into it (the same for any path: only the radial speed
/// counts). Measured: the induced charge (`induced_by`, `charges_of`) at r = 1.7 and its
/// derivative by central differences (h = 1e-4: truncation ~h² ≈ 1e-8, far below the
/// mesh error), with segments 0.1 and 0.05. Required as V2: second-order convergence
/// (order ≥ 1) and the two cylinders' induced charges summing to −λ within the error.
#[test]
fn v5_ramo_current_in_a_coaxial_gap() {
    let (a, b, r): (f64, f64, f64) = (1.0, 3.0, 1.7);
    let ln_ba = (b / a).ln();
    let q_exact = -(b / r).ln() / ln_ba;
    let i_exact = 1.0 / (r * ln_ba);
    let mut errors = Vec::new();
    for size in [0.1, 0.05, 0.025] {
        let e = Electrodes::new(
            vec![
                Electrode {
                    section: Section::Circle {
                        center: DVec3::ZERO,
                        radius: a,
                    },
                    bias: Bias::Grounded,
                },
                Electrode {
                    section: Section::Circle {
                        center: DVec3::ZERO,
                        radius: b,
                    },
                    bias: Bias::Grounded,
                },
            ],
            &[],
            size,
        );
        // Off a symmetry axis of the polygons, so the mesh's discreteness shows.
        let dir = DVec3::new(0.6, 0.8, 0.0);
        let q_at = |r: f64| {
            let p = dir * r;
            e.charges_of(&e.induced_by(|y| physics::zinv::line_charge(p, 1.0, y).0))
        };
        let q = q_at(r);
        let h = 1e-4;
        let i = (q_at(r + h)[0] - q_at(r - h)[0]) / (2.0 * h);
        let (eq, ei) = (q[0] / q_exact - 1.0, i / i_exact - 1.0);
        println!(
            "V5 segments {size}: Q_a {:.10} (exact {q_exact:.10}) {eq:.2e}, I_a {i:.10} (exact {i_exact:.10}) {ei:.2e}, total {:.10}",
            q[0],
            q[0] + q[1]
        );
        assert!(
            (q[0] + q[1] + 1.0).abs() < 10.0 * eq.abs().max(ei.abs()),
            "V5 total"
        );
        errors.push(eq.abs().max(ei.abs()));
    }
    let order = (errors[1] / errors[2]).ln() / 2f64.ln();
    println!("V5 observed order {order:.2}");
    assert!(
        errors[0] > errors[1] && errors[1] > errors[2] && order >= 1.0,
        "V5 {errors:?}"
    );
}

/// V6: the space-charge-limited bipolar coaxial diode (as V4, and the outer cylinder
/// emits ions, q/m = +1) against the steady radial bipolar flow
/// (`scripts/wolfram/v6_bipolar_coax.wls`: shot from both emitters with the planar
/// start, matched in the middle, 25 digits; its unipolar case reproduces V4's
/// Langmuir–Blodgett current): I_e = 0.9647391802, I_i = 0.7785212006 (the ions lift
/// the electron current 1.7146 times). Measured over t = 6 to 12 at the coarse and fine
/// resolutions of V4. Required, as the game's verdicts assume (`level::tube`,
/// `ERROR_FACTOR`): each reference within 1.5 times the two runs' difference of the fine
/// run, and the error shrinking with the refinement.
#[test]
#[ignore = "2 minutes in release: run with --ignored (measured values in PHYSICS.md §2.11)"]
fn v6_bipolar_coaxial_diode() {
    use physics::tube::{Emitter, Tube};
    let (ie_ref, ii_ref) = (0.964_739_180_2, 0.778_521_200_6);
    let mut runs = Vec::new();
    for (weight, dt, size) in [(2e-3, 0.01, 0.1), (1e-3, 0.005, 0.05)] {
        let electrodes = Electrodes::new(
            vec![
                Electrode {
                    section: Section::Circle {
                        center: DVec3::ZERO,
                        radius: 1.0,
                    },
                    bias: Bias::Grounded,
                },
                // A tube wall (its inner face at radius 2): it emits the ions into
                // the gap. As a solid `Circle` it emitted them outward, away from the
                // gap (the first run: I_i = 0).
                Electrode {
                    section: Section::Ring {
                        center: DVec3::ZERO,
                        radius: 2.25,
                        thickness: 0.5,
                    },
                    bias: Bias::Potential(1.0),
                },
            ],
            &[],
            size,
        );
        let tube = Tube {
            electrodes,
            emitters: vec![
                Emitter {
                    electrode: 0,
                    charge_per_mass: -1.0,
                    weight,
                    emit_toward: None,
                    saturation: None,
                },
                Emitter {
                    electrode: 1,
                    charge_per_mass: 1.0,
                    weight,
                    emit_toward: None,
                    saturation: None,
                },
            ],
            projectiles: Vec::new(),
            softening: 0.5 * size,
            dt,
            arena: None,
            b_z: 0.0,
            lattices: Vec::new(),
            circuits: Vec::new(),
        };
        let mut s = tube.state();
        let start = std::time::Instant::now();
        let mut at6 = None;
        while s.t < 12.0 {
            tube.step(&mut s);
            if at6.is_none() && s.t >= 6.0 {
                at6 = Some((s.t, s.collected[1][0], s.collected[0][1]));
            }
        }
        let (t6, qe, qi) = at6.expect("reached t = 6");
        // Electrons collected by the anode (negative charge), ions by the cathode.
        let ie = -(s.collected[1][0] - qe) / (s.t - t6);
        let ii = (s.collected[0][1] - qi) / (s.t - t6);
        println!(
            "V6 weight {weight}, dt {dt}, segments {size}: I_e {ie:.5} ({:.2e}), I_i {ii:.5} ({:.2e}); {} in flight, {:.1} s",
            ie / ie_ref - 1.0,
            ii / ii_ref - 1.0,
            s.x.len(),
            start.elapsed().as_secs_f64()
        );
        runs.push((ie, ii));
    }
    let ((ec, ic), (ef, i_f)) = (runs[0], runs[1]);
    for (what, coarse, fine, exact) in [("I_e", ec, ef, ie_ref), ("I_i", ic, i_f, ii_ref)] {
        let bar = 1.5 * (fine - coarse).abs();
        println!(
            "V6 {what}: fine {fine:.5} ± {bar:.5}, reference {exact:.5}: off by {:.2} of the bar",
            (fine - exact).abs() / bar
        );
        assert!(
            (fine - exact).abs() <= bar,
            "V6 {what}: {fine} ± {bar} vs {exact}"
        );
        assert!(
            (fine - exact).abs() < (coarse - exact).abs(),
            "V6 {what} does not converge"
        );
    }
}

/// V7: a projectile (one line charge λ = −0.05, mass 1) orbiting inside a coaxial pair
/// (inner radius 1 at potential 1, outer 3 grounded), no emitters: its energy
/// `½ m v² + λ φ_static(x) + ½ λ φ_induced(x)` (the last its own induced charge's, the
/// image energy) is conserved by the Boris/leapfrog push. Required: drift below 1e-3 of
/// the energy scale over 2000 steps (far below the macroparticle method's 5 %; without
/// the image force in the dynamics the image term alone would drift by its own size).
#[test]
fn v7_projectile_conserves_energy() {
    use physics::tube::{Projectile, Tube};
    use physics::zinv::line_charge;
    let electrodes = Electrodes::new(
        vec![
            Electrode {
                section: Section::Circle {
                    center: DVec3::ZERO,
                    radius: 1.0,
                },
                bias: Bias::Potential(1.0),
            },
            Electrode {
                section: Section::Circle {
                    center: DVec3::ZERO,
                    radius: 3.0,
                },
                bias: Bias::Grounded,
            },
        ],
        &[],
        0.1,
    );
    let (lambda, mass) = (-0.05, 1.0);
    let tube = Tube {
        electrodes,
        emitters: Vec::new(),
        projectiles: vec![Projectile {
            x0: DVec3::new(2.0, 0.0, 0.0),
            // About the circular speed at r = 2 (√(λE r/m) = 0.213): with 0.25 the apogee
            // came to r = 2.8, where the image force pulled it into the wall.
            v0: DVec3::new(0.0, 0.21, 0.0),
            launch: 0.0,
            charge: lambda,
            mass,
        }],
        softening: 0.05,
        dt: 0.01,
        arena: None,
        b_z: 0.0,
        lattices: Vec::new(),
        circuits: Vec::new(),
    };
    let energy = |s: &physics::tube::TubeState| {
        let (x, v) = (s.x[0], s.v[0]);
        let phi_static = tube.electrodes.sample(x).0;
        let induced = tube.electrodes.induced_by(|y| line_charge(x, lambda, y).0);
        let phi_induced: f64 = tube
            .electrodes
            .segments
            .iter()
            .zip(&induced)
            .map(|((a, b), sg)| sg * segment_integrals(*a, *b, x).0)
            .sum();
        (
            0.5 * mass * v.length_squared(),
            lambda * phi_static,
            0.5 * lambda * phi_induced,
        )
    };
    let mut s = tube.state();
    tube.step(&mut s);
    let e0 = energy(&s);
    let total0 = e0.0 + e0.1 + e0.2;
    let (mut worst, mut r_min, mut r_max) = (0.0_f64, f64::INFINITY, 0.0_f64);
    for _ in 0..2000 {
        let before = (s.t, s.x.first().copied(), s.v.first().copied());
        tube.step(&mut s);
        assert_eq!(
            s.x.len(),
            1,
            "V7: the projectile hit an electrode after {before:?}: {:?}",
            s.fates
        );
        let e = energy(&s);
        worst = worst.max((e.0 + e.1 + e.2 - total0).abs());
        let r = s.x[0].length();
        r_min = r_min.min(r);
        r_max = r_max.max(r);
    }
    let scale = e0.0.max(e0.1.abs());
    println!(
        "V7: energy {total0:.6} (kinetic {:.4}, static {:.4}, image {:.2e}); drift {:.2e} of {scale:.3}; orbit r {r_min:.3}..{r_max:.3}",
        e0.0,
        e0.1,
        e0.2,
        worst / scale
    );
    assert!(worst / scale < 1e-3, "V7 {:.3e}", worst / scale);
}

/// V8: electrons in crossed fields (a uniform B along z in a coaxial diode: the
/// cylindrical magnetron), against canonical angular momentum. An electron leaving the
/// cathode (radius a) at rest keeps `m r² θ̇ − (eB/2)(r² − a²) = 0`, so with energy
/// conservation its largest radius solves `(B/2)² (r² − a²)²/r² = 2 |q/m| φ(r)`, with
/// `φ(r) = V ln(r/a)/ln(b/a)` (no space charge), and it reaches the anode (radius b) only
/// below Hull's cutoff `B_H = (2b/(b² − a²)) √(2V/|q/m|)` (Hull 1921). Here a = 1, b = 2,
/// V = 1, |q/m| = 1: B_H = 1.88562. Test electrons: projectiles of negligible charge
/// (1e-9) launched at rest 0.05 outside the cathode at 8 angles (the exact root uses
/// that start; at 1e-3 outside, the cathode's polygon edges, 0.1 long, spread the largest
/// radii by ±0.007 with the angle). Required: below cutoff (0.9 B_H) all reach the
/// anode; above it (1.1 B_H) none, and every largest radius within 1e-3 of the root (the
/// field of the 0.05-cell mesh is exact to ~1e-4, V2; Boris conserves the angular
/// momentum to O(dt²)).
///
/// With space charge (an emitting cathode) the trapped cloud above cutoff is unstable
/// (the diocotron instability) and carries current across the field: 14 % of the
/// field-free current at 1.05 B_H in a first version of this test, which required none.
/// That is real magnetron physics, not a check of the push; PHYSICS.md §2.11.
#[test]
fn v8_crossed_fields_hull() {
    use physics::tube::{Projectile, Tube};
    let (a, b, v) = (1.0_f64, 2.0_f64, 1.0_f64);
    let b_h = 2.0 * b / (b * b - a * a) * (2.0 * v).sqrt();
    let r0 = a + 0.05;
    let phi = |r: f64| v * (r / a).ln() / (b / a).ln();
    // Largest radius: the root of g(r) = (B/2)²(r² − r0²)²/r² − 2(φ(r) − φ(r0)) above r0.
    let r_max = |bz: f64| {
        let g = |r: f64| {
            (bz / 2.0).powi(2) * (r * r - r0 * r0).powi(2) / (r * r) - 2.0 * (phi(r) - phi(r0))
        };
        let (mut lo, mut hi) = (r0 + 1e-6, b);
        if g(hi) < 0.0 {
            return f64::INFINITY; // reaches the anode
        }
        for _ in 0..200 {
            let mid = 0.5 * (lo + hi);
            if g(mid) < 0.0 {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        0.5 * (lo + hi)
    };
    for (factor, reaches) in [(0.9, true), (1.1, false)] {
        let bz = factor * b_h;
        let electrodes = Electrodes::new(
            vec![
                Electrode {
                    section: Section::Circle {
                        center: DVec3::ZERO,
                        radius: a,
                    },
                    bias: Bias::Grounded,
                },
                Electrode {
                    section: Section::Circle {
                        center: DVec3::ZERO,
                        radius: b,
                    },
                    bias: Bias::Potential(v),
                },
            ],
            &[],
            0.05,
        );
        let projectiles: Vec<Projectile> = (0..8)
            .map(|k| {
                let ang = std::f64::consts::TAU * (f64::from(k) + 0.3) / 8.0;
                Projectile {
                    x0: DVec3::new(ang.cos(), ang.sin(), 0.0) * r0,
                    v0: DVec3::ZERO,
                    launch: 0.0,
                    charge: -1e-9,
                    mass: 1e-9,
                }
            })
            .collect();
        let tube = Tube {
            electrodes,
            emitters: Vec::new(),
            projectiles,
            softening: 0.05,
            dt: 0.002,
            arena: None,
            b_z: bz,
            lattices: Vec::new(),
            circuits: Vec::new(),
        };
        let mut s = tube.state();
        let mut largest = [0.0_f64; 8];
        while s.t < 12.0 {
            tube.step(&mut s);
            for (j, l) in largest.iter_mut().enumerate() {
                if let Some(x) = s.projectile(j) {
                    *l = l.max(x.length());
                }
            }
        }
        let anode = s
            .fates
            .iter()
            .filter(|f| matches!(f, Some(physics::tube::Fate::Absorbed { electrode: 1, .. })))
            .count();
        let exact = r_max(bz);
        println!(
            "V8 at {factor} B_H: exact largest radius {exact:.5}; measured {largest:.5?}; {anode} of 8 reach the anode"
        );
        if reaches {
            assert_eq!(anode, 8, "V8 below cutoff");
        } else {
            assert_eq!(anode, 0, "V8 above cutoff");
            for l in largest {
                assert!(
                    (l - exact).abs() < 1e-3,
                    "V8: largest radius {l} vs {exact}"
                );
            }
        }
    }
}

/// V9: temperature-limited emission. The coaxial diode of V4 with the cathode's emission
/// capped at j = 0.02 per unit length of its surface (Richardson's limit of a cooler
/// cathode): j times the perimeter (0.126) is far below the space-charge limit (0.56,
/// V4), so the field keeps pulling every electron off and the anode collects exactly
/// `j × perimeter` (the polygon's), whatever the voltage (here 1 and 4). Required: within
/// 2 % over t = 4–12 (the only error is the macroparticles' quantization at the window's
/// ends: at most one macroparticle of 2e-4 per segment, 1.3 %).
#[test]
fn v9_temperature_limited_emission() {
    use physics::tube::{Emitter, Tube};
    let j_sat = 0.02;
    for v in [1.0, 4.0] {
        let electrodes = Electrodes::new(
            vec![
                Electrode {
                    section: Section::Circle {
                        center: DVec3::ZERO,
                        radius: 1.0,
                    },
                    bias: Bias::Grounded,
                },
                Electrode {
                    section: Section::Circle {
                        center: DVec3::ZERO,
                        radius: 2.0,
                    },
                    bias: Bias::Potential(v),
                },
            ],
            &[],
            0.1,
        );
        let perimeter: f64 = electrodes
            .segments
            .iter()
            .enumerate()
            .filter(|(k, _)| electrodes.owner_of(*k) == 0)
            .map(|(_, (a, b))| (*b - *a).length())
            .sum();
        let tube = Tube {
            electrodes,
            emitters: vec![Emitter {
                electrode: 0,
                charge_per_mass: -1.0,
                weight: 2e-4,
                emit_toward: None,
                saturation: Some(j_sat),
            }],
            projectiles: Vec::new(),
            softening: 0.05,
            dt: 0.01,
            arena: None,
            b_z: 0.0,
            lattices: Vec::new(),
            circuits: Vec::new(),
        };
        let mut s = tube.state();
        let mut at4 = None;
        while s.t < 12.0 {
            tube.step(&mut s);
            if at4.is_none() && s.t >= 4.0 {
                at4 = Some((s.t, s.collected[1][0]));
            }
        }
        let (t4, q4) = at4.expect("reached t = 4");
        let current = -(s.collected[1][0] - q4) / (s.t - t4);
        let exact = j_sat * perimeter;
        println!(
            "V9 at V = {v}: current {current:.6}, j × perimeter {exact:.6} ({:.2e})",
            current / exact - 1.0
        );
        assert!(
            (current / exact - 1.0).abs() < 0.02,
            "V9 at {v}: {current} vs {exact}"
        );
    }
}

/// The coaxial pair of the circuit tests: inner radius 1 at `v_inner`, outer radius 2 at
/// `v_outer`, segments 0.1.
fn coax(v_inner: f64, v_outer: f64) -> Electrodes {
    let bias = |v: f64| {
        if v == 0.0 {
            Bias::Grounded
        } else {
            Bias::Potential(v)
        }
    };
    Electrodes::new(
        vec![
            Electrode {
                section: Section::Circle {
                    center: DVec3::ZERO,
                    radius: 1.0,
                },
                bias: bias(v_inner),
            },
            Electrode {
                section: Section::Circle {
                    center: DVec3::ZERO,
                    radius: 2.0,
                },
                bias: bias(v_outer),
            },
        ],
        &[],
        0.1,
    )
}

/// V10: a circuit node discharging. The inner cylinder of a coaxial pair (outer grounded),
/// at potential 1, discharges through R = 2 into a source at 0, with C = 0.5 to ground:
/// `τ = R (C₁₁ + C)`, `C₁₁` its capacitance per unit length (the boundary elements'). The
/// node's backward-Euler step gives exactly `V_n = (1 + Δt/τ)^(−n)`: required to 1e-12
/// (rounding); and the exponential `e^(−t/τ)` within the method's error, which is below
/// `(t/τ)(Δt/τ)` relative for these steps.
#[test]
fn v10_circuit_node_discharges() {
    use physics::circuit::Waveform;
    use physics::tube::{Circuit, Tube};
    let (r, c_ext, dt) = (2.0, 0.5, 0.01);
    let electrodes = coax(1.0, 0.0);
    let c11 = electrodes.capacitance[0][0];
    let tau = r * (c11 + c_ext);
    let tube = Tube {
        electrodes,
        emitters: Vec::new(),
        projectiles: Vec::new(),
        softening: 0.05,
        dt,
        arena: None,
        b_z: 0.0,
        lattices: Vec::new(),
        circuits: vec![Circuit {
            electrode: 0,
            source: Waveform::Dc(0.0),
            resistance: r,
            capacitance: c_ext,
            inductance: 0.0,
            to: None,
        }],
    };
    let mut s = tube.state();
    let (mut worst_discrete, mut worst_exp) = (0.0_f64, 0.0_f64);
    for n in 1..=400 {
        tube.step(&mut s);
        let v = s.potentials[0];
        let discrete = (1.0 + dt / tau).powi(n);
        worst_discrete = worst_discrete.max((v * discrete - 1.0).abs());
        let t = f64::from(n) * dt;
        let exact = (-t / tau).exp();
        worst_exp = worst_exp.max(((v - exact) / exact).abs() / ((t / tau) * (dt / tau)));
    }
    println!(
        "V10: tau = {tau:.6} (C11 = {c11:.6}); V(4) = {:.6}; discrete sequence to {worst_discrete:.1e}, exponential within {worst_exp:.3} of the bound",
        s.potentials[0]
    );
    assert!(worst_discrete < 1e-12, "V10 discrete {worst_discrete:.3e}");
    assert!(worst_exp < 1.0, "V10 exponential {worst_exp:.3}");
}

/// V11: a beam charging a floating collector. The temperature-limited cathode of V9
/// (inner, grounded; j = 0.02: a constant current I = j × perimeter) and the outer cylinder
/// a floating node at first at 4 with C = 10 to ground: the electrons it collects lower its
/// potential at `dV/dt = −I / (C₂₂ + C)` while the flow is steady (its induced charge then
/// constant), `C₂₂` its own capacitance. Measured over t = 4–12 against the collected
/// current. Required: 2 % (V9's quantization of the current; the potential follows the
/// node's charge exactly).
#[test]
fn v11_beam_charges_a_floating_collector() {
    use physics::circuit::Waveform;
    use physics::tube::{Circuit, Emitter, Tube};
    let c_ext = 10.0;
    let electrodes = coax(0.0, 4.0);
    let c22 = electrodes.capacitance[1][1];
    let tube = Tube {
        electrodes,
        emitters: vec![Emitter {
            electrode: 0,
            charge_per_mass: -1.0,
            weight: 2e-4,
            emit_toward: None,
            saturation: Some(0.02),
        }],
        projectiles: Vec::new(),
        softening: 0.05,
        dt: 0.01,
        arena: None,
        b_z: 0.0,
        lattices: Vec::new(),
        circuits: vec![Circuit {
            electrode: 1,
            source: Waveform::Dc(0.0),
            resistance: f64::INFINITY,
            capacitance: c_ext,
            inductance: 0.0,
            to: None,
        }],
    };
    let mut s = tube.state();
    let mut at4 = None;
    while s.t < 12.0 {
        tube.step(&mut s);
        if at4.is_none() && s.t >= 4.0 {
            at4 = Some((s.t, s.collected[1][0], s.potentials[1]));
        }
    }
    let (t4, q4, v4) = at4.expect("reached t = 4");
    let current = (s.collected[1][0] - q4) / (s.t - t4);
    let slope = (s.potentials[1] - v4) / (s.t - t4);
    let expected = current / (c22 + c_ext);
    println!(
        "V11: collected current {current:.6}; dV/dt {slope:.6e}, expected {expected:.6e} ({:.2e}); V from {v4:.4} to {:.4}",
        slope / expected - 1.0,
        s.potentials[1]
    );
    assert!(
        (slope / expected - 1.0).abs() < 0.02,
        "V11 {slope} vs {expected}"
    );
}

/// V12: a carrier coming back to its emitter is absorbed there (`RETURN_ABSORB`). The
/// coaxial magnetron of V8 above Hull's cutoff (B = 1.1 B_H), now emitting: a
/// temperature-limited cathode (j = 1e-4 per unit length: its electrons' field is 1e-4 of
/// the electrodes') releases electrons at rest δ = half a segment out; each flies out to
/// r_max and back, and is absorbed on crossing 1 + 1.5 δ. Little's law then fixes the
/// charge in flight: the emission rate (j × the polygon's perimeter) times the time of
/// flight, exact from energy and canonical angular momentum,
/// `v_r² = 2(φ(r) − φ(r₀)) − (B/2)²(r² − r₀²)²/r²`, out from the release radius r₀ to
/// r_max and back to 1 + 1.5 δ (substitutions remove the integrable endpoint
/// singularities; 20 000 midpoints). Before the rule the electrons skipped around the
/// cathode for ever, back at δ with zero speed every time, and the charge in flight grew
/// without bound. Required: the mean charge in flight over t = 6–12 divided by the rate
/// within 5e-3 of the exact time (the error budget: the step's quantization of release and
/// absorption, dt/T ~ 1e-3; the macroparticles' discreteness; the mesh's field, 1e-4), and
/// nothing reaching the anode. (The current back into the cathode is printed: its equal
/// segments emit together, in bursts 0.55 apart, so a 6-long window counts it to ±9 %.)
#[test]
fn v12_returning_carriers_are_absorbed() {
    use physics::tube::{EMISSION_OFFSET, Emitter, RETURN_ABSORB, Tube};
    use std::f64::consts::PI;
    let (a, b, v) = (1.0_f64, 2.0_f64, 1.0_f64);
    let bz = 1.1 * 2.0 * b / (b * b - a * a) * (2.0 * v).sqrt();
    let size = 0.1;
    let electrodes = Electrodes::new(
        vec![
            Electrode {
                section: Section::Circle {
                    center: DVec3::ZERO,
                    radius: a,
                },
                bias: Bias::Grounded,
            },
            Electrode {
                section: Section::Circle {
                    center: DVec3::ZERO,
                    radius: b,
                },
                bias: Bias::Potential(v),
            },
        ],
        &[],
        size,
    );
    let lengths: Vec<f64> = (0..electrodes.segments.len())
        .filter(|&j| electrodes.owner_of(j) == 0)
        .map(|j| {
            let (p, q) = electrodes.segments[j];
            (q - p).length()
        })
        .collect();
    let perimeter: f64 = lengths.iter().sum();
    let seg = lengths[0];
    let delta = EMISSION_OFFSET * seg;
    let j_sat = 1e-4;
    let rate = j_sat * perimeter;
    // Released along a chord h = √(a² − L²/4) from the axis, δ out along its normal: on
    // average r₀ = (h + δ) + L²/(24 (h + δ)) (uniform points along the chord).
    let h = (a * a - seg * seg / 4.0).sqrt();
    let r0 = (h + delta) + seg * seg / (24.0 * (h + delta));
    let r_abs = a + RETURN_ABSORB * delta;
    let phi = |r: f64| v * (r / a).ln() / (b / a).ln();
    let vr2 = |r: f64| {
        2.0 * (phi(r) - phi(r0)) - (bz / 2.0).powi(2) * (r * r - r0 * r0).powi(2) / (r * r)
    };
    let (mut lo, mut hi) = (r0 + 1e-9, b);
    assert!(vr2(hi) < 0.0, "V12: the orbit must turn before the anode");
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        if vr2(mid) > 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let r_max = 0.5 * (lo + hi);
    let m = 20_000;
    let mf = f64::from(m);
    // Out: r = r₀ + (r_max − r₀)(1 − cos θ)/2; the stretch r₀ → r_abs not flown back:
    // r = r₀ + (r_abs − r₀) u².
    let out: f64 = (0..m)
        .map(|i| {
            let th = PI * (f64::from(i) + 0.5) / mf;
            let r = r0 + (r_max - r0) * (1.0 - th.cos()) / 2.0;
            (r_max - r0) / 2.0 * th.sin() / vr2(r).sqrt()
        })
        .sum::<f64>()
        * PI
        / mf;
    let skipped: f64 = (0..m)
        .map(|i| {
            let u = (f64::from(i) + 0.5) / mf;
            let r = r0 + (r_abs - r0) * u * u;
            2.0 * (r_abs - r0) * u / vr2(r).sqrt()
        })
        .sum::<f64>()
        / mf;
    let exact = 2.0 * out - skipped;
    let tube = Tube {
        electrodes,
        emitters: vec![Emitter {
            electrode: 0,
            charge_per_mass: -1.0,
            weight: rate * exact / 400.0,
            emit_toward: None,
            saturation: Some(j_sat),
        }],
        projectiles: Vec::new(),
        softening: 0.5 * size,
        dt: 0.002,
        arena: None,
        b_z: bz,
        lattices: Vec::new(),
        circuits: Vec::new(),
    };
    let mut s = tube.state();
    let (mut sum, mut samples) = (0.0, 0.0);
    let mut at6 = None;
    while s.t < 12.0 {
        tube.step(&mut s);
        if s.t >= 6.0 {
            if at6.is_none() {
                at6 = Some((s.t, s.collected[0][0]));
            }
            sum += -s.charge.iter().sum::<f64>();
            samples += 1.0;
        }
    }
    let (t6, back6) = at6.expect("reached t = 6");
    let measured = sum / samples / rate;
    let back = -(s.collected[0][0] - back6) / (s.t - t6);
    println!(
        "V12: δ {delta:.5}, r₀ {r0:.5}, r_max {r_max:.5}; time of flight {measured:.6}, exact {exact:.6} ({:.2e}); \
         back to the cathode {back:.6e} vs emitted {rate:.6e} ({:.2e}); anode {:.3e}; {} in flight",
        measured / exact - 1.0,
        back / rate - 1.0,
        s.collected[1][0],
        s.x.len()
    );
    assert!(
        (measured / exact - 1.0).abs() < 5e-3,
        "V12 time of flight {measured} vs {exact}"
    );
    assert!(s.collected[1][0] == 0.0, "V12: electrons reached the anode");
}

/// A measurement, not a validation (PHYSICS.md §2.11, crossed fields): the coaxial
/// magnetron of V8 with a space-charge-limited cathode, at 1.05 B_H and without field (V4's
/// coarse resolution). Above Hull's cutoff the space-charge cloud is unstable (the
/// diocotron instability, seeded by the macroparticles) and carries some current across
/// the field; Hull's theorem covers steady axisymmetric flow only, so there is no
/// reference.
#[test]
#[ignore = "a measurement (about a minute in release): run with --ignored"]
fn magnetron_with_space_charge_measurement() {
    use physics::tube::{Emitter, Tube};
    let (a, b, v) = (1.0_f64, 2.0_f64, 1.0_f64);
    let b_h = 2.0 * b / (b * b - a * a) * (2.0 * v).sqrt();
    let mut currents = Vec::new();
    for factor in [0.0, 1.05] {
        let size = 0.1;
        let electrodes = Electrodes::new(
            vec![
                Electrode {
                    section: Section::Circle {
                        center: DVec3::ZERO,
                        radius: a,
                    },
                    bias: Bias::Grounded,
                },
                Electrode {
                    section: Section::Circle {
                        center: DVec3::ZERO,
                        radius: b,
                    },
                    bias: Bias::Potential(v),
                },
            ],
            &[],
            size,
        );
        let tube = Tube {
            electrodes,
            emitters: vec![Emitter {
                electrode: 0,
                charge_per_mass: -1.0,
                weight: 2e-3,
                emit_toward: None,
                saturation: None,
            }],
            projectiles: Vec::new(),
            softening: 0.5 * size,
            dt: 0.01,
            arena: None,
            b_z: factor * b_h,
            lattices: Vec::new(),
            circuits: Vec::new(),
        };
        let mut s = tube.state();
        let mut at6 = None;
        while s.t < 12.0 {
            tube.step(&mut s);
            if at6.is_none() && s.t >= 6.0 {
                at6 = Some((s.t, s.collected[1][0], s.collected[0][0]));
            }
        }
        let (t6, q6, back6) = at6.expect("reached t = 6");
        let current = -(s.collected[1][0] - q6) / (s.t - t6);
        let back = -(s.collected[0][0] - back6) / (s.t - t6);
        println!(
            "magnetron at {factor} B_H: anode current {current:.5}, back into the cathode {back:.5}; {} in flight",
            s.x.len()
        );
        currents.push(current);
    }
    println!(
        "magnetron: at 1.05 B_H {:.3} of the field-free current",
        currents[1] / currents[0]
    );
}

/// V13: the envelope (`Section::Frame`, meshed on its inner face only). A line charge
/// λ = 1 at (2.2, 1.7) inside a grounded rectangle [0, 6] × [0, 4], against the exact
/// Dirichlet Green's function of the rectangle (Gaussian units, ∇²φ = −4πλδ: a sine
/// series in x whose coefficients jump at the source's y),
/// `φ = Σ_n (8λ/n) sin(k_n x₀) sin(k_n x) sinh(k_n y<) sinh(k_n (b − y>)) / sinh(k_n b)`,
/// `k_n = nπ/a`, at three field points off the source's row, where it converges
/// exponentially (2000 terms, in exponentials). The BEM's potential there is the line
/// charge's plus the frame's surface charge. With segments of 0.5 (the game's envelope at
/// both resolutions, `ENVELOPE_SEGMENT`), 0.25 and 0.125. Required: second order as V2 and
/// V3 (observed order ≥ 1.5) and at the game's 0.5 within 3e-3 of the potential's scale
/// (the largest |φ| of the points): well below the macroparticle method's few per cent.
/// (Segments of 1 gave 1.2e-2: the game's envelope does not use them.)
#[test]
fn v13_line_charge_in_a_grounded_envelope() {
    use physics::zinv::line_charge;
    let (a, b) = (6.0_f64, 4.0_f64);
    let (x0, y0, lambda) = (2.2_f64, 1.7_f64, 1.0_f64);
    let points = [(4.1, 2.9), (1.0, 3.2), (3.5, 0.6)];
    let exact = |x: f64, y: f64| -> f64 {
        let (lo, hi) = if y < y0 { (y, y0) } else { (y0, y) };
        (1..=2000)
            .map(|n| {
                let k = f64::from(n) * std::f64::consts::PI / a;
                let ratio = (-k * (hi - lo)).exp()
                    * (1.0 - (-2.0 * k * lo).exp())
                    * (1.0 - (-2.0 * k * (b - hi)).exp())
                    / (2.0 * (1.0 - (-2.0 * k * b).exp()));
                8.0 * lambda / f64::from(n) * (k * x0).sin() * (k * x).sin() * ratio
            })
            .sum()
    };
    let reference: Vec<f64> = points.iter().map(|&(x, y)| exact(x, y)).collect();
    let scale = reference.iter().fold(0.0_f64, |m, v| m.max(v.abs()));
    let mut errors = Vec::new();
    for size in [0.5, 0.25, 0.125] {
        let electrodes = Electrodes::new(
            vec![Electrode {
                section: Section::Frame {
                    center: DVec3::new(0.5 * a, 0.5 * b, 0.0),
                    half_width: 0.5 * a,
                    half_height: 0.5 * b,
                    segment: size,
                },
                bias: Bias::Grounded,
            }],
            &[(DVec3::new(x0, y0, 0.0), lambda)],
            size,
        );
        let source = DVec3::new(x0, y0, 0.0);
        let worst = points
            .iter()
            .zip(&reference)
            .map(|(&(x, y), r)| {
                let p = DVec3::new(x, y, 0.0);
                let phi = electrodes.sample(p).0 + line_charge(source, lambda, p).0;
                (phi - r).abs() / scale
            })
            .fold(0.0_f64, f64::max);
        println!(
            "V13 segments {size}: {} of them, worst {worst:.2e} of the scale {scale:.4} (points {reference:.6?})",
            electrodes.segments.len()
        );
        errors.push(worst);
    }
    let order = (errors[1] / errors[2]).ln() / 2f64.ln();
    println!("V13 observed order {order:.2}");
    assert!(order >= 1.5, "V13 order {order}");
    assert!(errors[0] <= 3e-3, "V13 at 0.5: {}", errors[0]);
}

/// V14: the far field of a segment (`segment_field_fast`, the monopole and quadrupole of
/// its multipole series beyond `FAR_SEGMENTS` = 4 lengths) against the exact integrals, at
/// 20 000 points from 4 to 100 lengths in every direction about segments of every
/// orientation and of lengths 0.05 to 2. The first omitted term is `(L/r)⁴/80` of the
/// monopole: 4.9e-5 at 4 lengths. Required: below 1e-4 of the field everywhere (far below
/// the boundary elements' discretization, V2: 2e-4 at the finest mesh), and the exact
/// integrals used inside 4 lengths (bit-identical there).
#[test]
fn v14_far_field_of_a_segment() {
    use physics::zinv::{FAR_SEGMENTS, segment_field, segment_field_fast};
    let mut worst: f64 = 0.0;
    let mut worst_at = 0.0;
    let mut seed: u64 = 0x2545_f491_4f6c_dd1d;
    let mut next = move || {
        seed ^= seed << 13;
        seed ^= seed >> 7;
        seed ^= seed << 17;
        #[allow(clippy::cast_precision_loss)]
        let f = (seed >> 11) as f64 / (1u64 << 53) as f64;
        f
    };
    for _ in 0..20_000 {
        let length = 0.05 * 40.0_f64.powf(next());
        let angle = std::f64::consts::TAU * next();
        let a = DVec3::new(10.0 * next() - 5.0, 10.0 * next() - 5.0, 0.0);
        let b = a + DVec3::new(angle.cos(), angle.sin(), 0.0) * length;
        let mid = (a + b) * 0.5;
        let distance = length * FAR_SEGMENTS * 25.0_f64.powf(next());
        let dir = std::f64::consts::TAU * next();
        let x = mid + DVec3::new(dir.cos(), dir.sin(), 0.0) * distance;
        let exact = segment_field(a, b, x);
        let fast = segment_field_fast(a, b, x);
        let err = (fast - exact).length() / exact.length();
        if err > worst {
            worst = err;
            worst_at = distance / length;
        }
    }
    // Inside the switch the exact integrals are used.
    let (a, b) = (DVec3::ZERO, DVec3::new(1.0, 0.0, 0.0));
    let near = DVec3::new(0.5, 3.9, 0.0);
    assert_eq!(segment_field_fast(a, b, near), segment_field(a, b, near));
    println!("V14: worst relative error {worst:.3e} at {worst_at:.2} lengths");
    assert!(worst < 1e-4, "V14 {worst}");
}

/// V15: the emitter meshed finer than the rest (the game's `EMITTER_MESH`, a quarter). V4's
/// coaxial diode at V4's two resolutions with only the cathode's segments at a quarter: its
/// particles are released half of its own segment out, so the skipped space-charge layer,
/// V4's first-order error, shrinks with it. Measured (PHYSICS.md §2.11): the cathode at the
/// others' size, a half and a quarter gave +15.1 %, +9.2 %, +5.4 % at V4's coarse resolution
/// and +7.8 %, +4.9 %, +1.0 % at its fine one. Required: at both resolutions below V4's
/// uniform runs (the refinement helps), and the fine run within the game's verdict bar
/// (1.5 × the two runs' difference) of the exact current, the rule the verdicts rest on.
#[test]
fn v15_finer_emitter_mesh() {
    use physics::tube::{Emitter, Tube};
    let beta2 = 0.279_267_160_710_59;
    let exact = 2.0 / 9.0 * 2f64.sqrt() / (2.0 * beta2);
    let mut currents = Vec::new();
    for (weight, dt, size) in [(2e-3, 0.01, 0.1), (1e-3, 0.005, 0.05)] {
        let electrodes = Electrodes::with_sizes(
            vec![
                Electrode {
                    section: Section::Circle {
                        center: DVec3::ZERO,
                        radius: 1.0,
                    },
                    bias: Bias::Grounded,
                },
                Electrode {
                    section: Section::Circle {
                        center: DVec3::ZERO,
                        radius: 2.0,
                    },
                    bias: Bias::Potential(1.0),
                },
            ],
            &[],
            &[0.25 * size, size],
        );
        let tube = Tube {
            electrodes,
            emitters: vec![Emitter {
                electrode: 0,
                charge_per_mass: -1.0,
                weight,
                emit_toward: None,
                saturation: None,
            }],
            projectiles: Vec::new(),
            softening: 0.5 * size,
            dt,
            arena: None,
            b_z: 0.0,
            lattices: Vec::new(),
            circuits: Vec::new(),
        };
        let mut s = tube.state();
        let start = std::time::Instant::now();
        let mut at4 = None;
        while s.t < 8.0 {
            tube.step(&mut s);
            if at4.is_none() && s.t >= 4.0 {
                at4 = Some((s.t, s.collected[1][0]));
            }
        }
        let (t4, q4) = at4.expect("reached t = 4");
        let current = -(s.collected[1][0] - q4) / (s.t - t4);
        println!(
            "V15 segments {size}, cathode {:.4}: I/L = {current:.6} ({:+.2e}); {} in flight, {:.1} s",
            0.25 * size,
            current / exact - 1.0,
            s.x.len(),
            start.elapsed().as_secs_f64()
        );
        currents.push(current);
    }
    let (coarse, fine) = (currents[0], currents[1]);
    let bar = 1.5 * (fine - coarse).abs();
    println!(
        "V15: fine {fine:.6} ± {bar:.6}, exact {exact:.6}: off by {:.2} of the bar",
        (fine - exact).abs() / bar
    );
    // V4's uniform runs: +15.1 % and +7.8 %.
    assert!(coarse / exact - 1.0 < 0.151, "V15 coarse {coarse}");
    assert!(fine / exact - 1.0 < 0.078, "V15 fine {fine}");
    assert!(
        (fine - exact).abs() <= bar,
        "V15 fine {fine} ± {bar} vs {exact}"
    );
}

/// V16: the energy balance of a particle coupled to a circuit (the tubes' plan, docs/TUBES.md:
/// the sources' work = the Joule heat + the particles' energy + the field's). A projectile
/// (λ = −0.05, mass 1) orbits the inner conductor of a coaxial pair (radii 1 and 3; the inner
/// held at 1 by an ideal source, the outer a circuit node through R to ground with C = 1).
/// Its eccentric orbit changes the charge it induces on the outer conductor, which drives a
/// current through R: the orbit is braked by Joule heating. The field energy of conductors
/// and a line charge (its self-energy aside) is `W = ½ Σ Q_i V_i + ½ λ φ_c(x)`, with the
/// conductors' charges and their potential at the particle from the model's own surface
/// charge; so `E = ½ m v² + W + ½ C V² + ∫ V²/R dt − ∫ V_in dQ_in` is conserved, the last term
/// the inner source's work. The circuit is integrated by backward Euler (first order), the
/// particle by leapfrog (its velocity half a step from its position). Required: the
/// residual (the drift of E over t = 0–40) below 2 % of the Joule heat at dt = 0.002, and
/// first order: at least 1.6 times smaller than at dt = 0.004.
#[test]
fn v16_circuit_energy_balance() {
    energy_balance(0.0);
}

/// V16 with a coil (L = 5) in series with the outer conductor's resistance: the circuit is
/// a damped LC tank (its period about 20, Q ≈ 0.75), integrated by the trapezoidal rule,
/// and the balance gains the coil's energy `½ L I²`. The same requirements as V16.
#[test]
fn v16_circuit_energy_balance_with_a_coil() {
    energy_balance(5.0);
}

/// V16's balance with inductance `l` in the outer conductor's branch.
#[allow(clippy::too_many_lines)]
fn energy_balance(l: f64) {
    use physics::circuit::Waveform;
    use physics::tube::{Circuit, Projectile, Tube};
    use physics::zinv::segment_integrals;
    let (lambda, mass, r_load) = (-0.05_f64, 1.0_f64, 2.0_f64);
    let mut residuals = Vec::new();
    for dt in [0.004, 0.002] {
        let electrodes = Electrodes::new(
            vec![
                Electrode {
                    section: Section::Circle {
                        center: DVec3::ZERO,
                        radius: 1.0,
                    },
                    bias: Bias::Potential(1.0),
                },
                Electrode {
                    section: Section::Circle {
                        center: DVec3::ZERO,
                        radius: 3.0,
                    },
                    bias: Bias::Grounded,
                },
            ],
            &[],
            0.1,
        );
        let tube = Tube {
            electrodes,
            emitters: Vec::new(),
            // Slower than circular at r = 1.8 (0.213): an eccentric orbit between ~1.5 and
            // 1.8, clear of the conductors.
            projectiles: vec![Projectile {
                x0: DVec3::new(1.8, 0.0, 0.0),
                v0: DVec3::new(0.0, 0.195, 0.0),
                launch: 0.0,
                charge: lambda,
                mass,
            }],
            softening: 0.05,
            dt,
            arena: None,
            b_z: 0.0,
            lattices: Vec::new(),
            circuits: vec![Circuit {
                electrode: 1,
                source: Waveform::Dc(0.0),
                resistance: r_load,
                capacitance: 1.0,
                inductance: l,
                to: None,
            }],
        };
        let mut s = tube.state();
        // The energy now: kinetic, field, capacitor, coil.
        let energy = |s: &physics::tube::TubeState| -> (f64, f64) {
            let sigma = tube.surface_density(s);
            let q = tube.electrodes.charges_of(&sigma);
            let field: f64 = q
                .iter()
                .zip(&s.potentials)
                .map(|(q, v)| 0.5 * q * v)
                .sum::<f64>()
                + s.x
                    .iter()
                    .zip(&s.charge)
                    .map(|(x, l)| {
                        let phi: f64 = tube
                            .electrodes
                            .segments
                            .iter()
                            .zip(&sigma)
                            .map(|((a, b), d)| d * segment_integrals(*a, *b, *x).0)
                            .sum();
                        0.5 * l * phi
                    })
                    .sum::<f64>();
            let kinetic: f64 =
                s.v.iter()
                    .zip(s.charge.iter().zip(&s.qm))
                    .map(|(v, (l, qm))| 0.5 * (l / qm) * v.length_squared())
                    .sum();
            let coil = 0.5 * l * s.inductor_current[0] * s.inductor_current[0];
            (
                kinetic + field + 0.5 * 1.0 * s.potentials[1] * s.potentials[1] + coil,
                q[0],
            )
        };
        tube.step(&mut s);
        let (e0, mut q_in) = energy(&s);
        let (mut joule, mut work) = (0.0, 0.0);
        let mut r_range = (f64::INFINITY, 0.0_f64);
        while s.t < 40.0 {
            tube.step(&mut s);
            let i = s.circuit_current[0];
            joule += i * i * r_load * dt;
            let (_, q) = energy(&s);
            work += 1.0 * (q - q_in);
            q_in = q;
            let Some(x) = s.x.first() else { break };
            r_range = (r_range.0.min(x.length()), r_range.1.max(x.length()));
        }
        let (e1, _) = energy(&s);
        let residual = (e1 + joule - work) - e0;
        println!(
            "V16 L {l} dt {dt}: E {e0:.8} -> {e1:.8}, Joule {joule:.6e}, source work {work:.6e}, residual {residual:.3e} ({:.2e} of the Joule heat); r in {:.3}..{:.3}",
            residual / joule,
            r_range.0,
            r_range.1
        );
        assert!(s.fates[0].is_none(), "V16: the projectile hit an electrode");
        residuals.push((residual, joule));
    }
    let ((r0, _), (r1, j1)) = (residuals[0], residuals[1]);
    println!("V16 L {l}: residual ratio {:.2}", r0.abs() / r1.abs());
    assert!(r1.abs() < 0.02 * j1, "V16 residual {r1} vs Joule {j1}");
    assert!(
        r0.abs() >= 1.6 * r1.abs(),
        "V16 not first order: {r0} then {r1}"
    );
}

/// V17: an LC tank. The inner cylinder of the coaxial pair (outer grounded) starts at
/// potential 1 with no current, through a coil L = 1 (no resistance) to ground, with
/// C = 0.5 to ground: `C = C₁₁ + 0.5`, `ω = 1/√(LC)`. The trapezoidal rule on the coil's
/// branch with the node's charge balance is the implicit midpoint rule, which maps
/// `(√C V, √L I)` by a rotation through `θ = 2 atan(ω Δt/2)` each step: the run must give
/// `V_n = cos(nθ)` and keep `½ C V² + ½ L I²` to rounding (required 1e-11 over 4000
/// steps, 6 periods). Its frequency `θ/Δt` is then low by `(ω Δt)²/12` (printed): an LC
/// tank's only error is a small detuning, never a numerical loss (backward Euler would
/// damp it by `e^{−ω² Δt t/2}`, 15 % of the amplitude over these 6 periods).
#[test]
fn v17_lc_tank_keeps_its_energy() {
    use physics::circuit::Waveform;
    use physics::tube::{Circuit, Tube};
    let (l, c_ext, dt) = (1.0, 0.5, 0.01);
    let electrodes = coax(1.0, 0.0);
    let c = electrodes.capacitance[0][0] + c_ext;
    let omega = 1.0 / (l * c).sqrt();
    let theta = 2.0 * (0.5 * omega * dt).atan();
    let tube = Tube {
        electrodes,
        emitters: Vec::new(),
        projectiles: Vec::new(),
        softening: 0.05,
        dt,
        arena: None,
        b_z: 0.0,
        lattices: Vec::new(),
        circuits: vec![Circuit {
            electrode: 0,
            source: Waveform::Dc(0.0),
            resistance: 0.0,
            capacitance: c_ext,
            inductance: l,
            to: None,
        }],
    };
    let mut s = tube.state();
    let e0 = 0.5 * c;
    let (mut worst_v, mut worst_e) = (0.0_f64, 0.0_f64);
    for n in 1..=4000 {
        tube.step(&mut s);
        let (v, i) = (s.potentials[0], s.inductor_current[0]);
        worst_v = worst_v.max((v - (f64::from(n) * theta).cos()).abs());
        worst_e = worst_e.max(((0.5 * c * v * v + 0.5 * l * i * i) / e0 - 1.0).abs());
    }
    println!(
        "V17: C = {c:.6}, omega = {omega:.6}; V against cos(n theta) to {worst_v:.1e}, energy to {worst_e:.1e}; frequency {:+.3e} (-(omega dt)^2/12 = {:+.3e})",
        theta / dt / omega - 1.0,
        -(omega * dt).powi(2) / 12.0
    );
    assert!(worst_v < 1e-11, "V17 V {worst_v:.3e}");
    assert!(worst_e < 1e-11, "V17 energy {worst_e:.3e}");
}

/// V18: a tank between two electrodes (a branch returning to another electrode, as between
/// a magnetron's vanes). Two wires (radius 0.5, centres 2.4 apart) in a grounded cylinder
/// (radius 3): wire A starts at potential 1, wire B floats with no charge, and a coil
/// L = 2 joins them. Charge moves between them only, so `D = V_B − V_A` obeys
/// `L C_eff D'' = −D` with `1/C_eff = eᵀ C⁻¹ e`, `e = (1, −1)`, `C` the wires' 2 × 2 block of
/// the capacitance matrix (the boundary elements'): the run must give `D_n = D₀ cos(nθ)`,
/// `θ = 2 atan(ω Δt/2)`, and keep `Q_A + Q_B` (required 1e-11 over 4000 steps).
#[test]
fn v18_tank_between_two_electrodes() {
    use physics::circuit::Waveform;
    use physics::tube::{Circuit, Tube};
    let (l, dt) = (2.0, 0.01);
    let wire = |x: f64, bias: Bias| Electrode {
        section: Section::Circle {
            center: DVec3::new(x, 0.0, 0.0),
            radius: 0.5,
        },
        bias,
    };
    let electrodes = Electrodes::new(
        vec![
            wire(-1.2, Bias::Potential(1.0)),
            wire(1.2, Bias::Charge(0.0)),
            Electrode {
                section: Section::Circle {
                    center: DVec3::ZERO,
                    radius: 3.0,
                },
                bias: Bias::Grounded,
            },
        ],
        &[],
        0.1,
    );
    let m = &electrodes.capacitance;
    let det = m[0][0] * m[1][1] - m[0][1] * m[1][0];
    let inv_c_eff = (m[1][1] + m[0][1] + m[1][0] + m[0][0]) / det;
    let omega = (inv_c_eff / l).sqrt();
    let theta = 2.0 * (0.5 * omega * dt).atan();
    let tube = Tube {
        electrodes,
        emitters: Vec::new(),
        projectiles: Vec::new(),
        softening: 0.05,
        dt,
        arena: None,
        b_z: 0.0,
        lattices: Vec::new(),
        circuits: vec![Circuit {
            electrode: 0,
            source: Waveform::Dc(0.0),
            resistance: 0.0,
            capacitance: 0.0,
            inductance: l,
            to: Some(1),
        }],
    };
    let mut s = tube.state();
    let d0 = s.potentials[1] - s.potentials[0];
    let total = |s: &physics::tube::TubeState| {
        let q = tube.electrodes.charges_of(&tube.surface_density(s));
        q[0] + q[1]
    };
    let q0 = total(&s);
    let (mut worst_d, mut worst_q) = (0.0_f64, 0.0_f64);
    for n in 1..=4000 {
        tube.step(&mut s);
        let d = s.potentials[1] - s.potentials[0];
        worst_d = worst_d.max((d - d0 * (f64::from(n) * theta).cos()).abs() / d0.abs());
        worst_q = worst_q.max((total(&s) - q0).abs());
    }
    println!(
        "V18: D0 = {d0:.6}, 1/C_eff = {inv_c_eff:.6}, omega = {omega:.6}; D against D0 cos(n theta) to {worst_d:.1e}, Q_A + Q_B to {worst_q:.1e}"
    );
    assert!(worst_d < 1e-11, "V18 D {worst_d:.3e}");
    assert!(worst_q < 1e-11, "V18 charge {worst_q:.3e}");
}

/// V19: a strap (a branch with no resistance and no coil to another electrode: a wire).
/// Four wires (radius 0.5) at 1.5 from the centre of a grounded cylinder (radius 3), at
/// 0°, 90°, 180°, 270°: wires 0 and 2 strapped together (a magnetron's alternate vanes) and
/// through a coil L = 1 to ground, wires 1 and 3 grounded; 0 and 2 start at potential 1. The
/// strapped pair is one node of capacitance `C = C₀₀ + C₀₂ + C₂₀ + C₂₂`, so (as V17) the run
/// must give both wires `V_n = cos(nθ)`, `θ = 2 atan(ω Δt/2)`, `ω = 1/√(LC)`, the same to
/// the last bit, and keep `½ C V² + ½ L I²` (required 1e-11 over 4000 steps).
#[test]
fn v19_strapped_electrodes_are_one_node() {
    use physics::circuit::Waveform;
    use physics::tube::{Circuit, Tube};
    let (l, dt) = (1.0, 0.01);
    let wire = |deg: f64, bias: Bias| {
        let a = deg.to_radians();
        Electrode {
            section: Section::Circle {
                center: DVec3::new(1.5 * a.cos(), 1.5 * a.sin(), 0.0),
                radius: 0.5,
            },
            bias,
        }
    };
    let electrodes = Electrodes::new(
        vec![
            wire(0.0, Bias::Potential(1.0)),
            wire(90.0, Bias::Grounded),
            wire(180.0, Bias::Potential(1.0)),
            wire(270.0, Bias::Grounded),
            Electrode {
                section: Section::Circle {
                    center: DVec3::ZERO,
                    radius: 3.0,
                },
                bias: Bias::Grounded,
            },
        ],
        &[],
        0.1,
    );
    let m = &electrodes.capacitance;
    let (c00, c02) = (m[0][0], m[0][2]);
    let c = m[0][0] + m[0][2] + m[2][0] + m[2][2];
    let omega = 1.0 / (l * c).sqrt();
    let theta = 2.0 * (0.5 * omega * dt).atan();
    let tube = Tube {
        electrodes,
        emitters: Vec::new(),
        projectiles: Vec::new(),
        softening: 0.05,
        dt,
        arena: None,
        b_z: 0.0,
        lattices: Vec::new(),
        circuits: vec![
            Circuit {
                electrode: 0,
                source: Waveform::Dc(0.0),
                resistance: 0.0,
                capacitance: 0.0,
                inductance: l,
                to: None,
            },
            Circuit {
                electrode: 2,
                source: Waveform::Dc(0.0),
                resistance: 0.0,
                capacitance: 0.0,
                inductance: 0.0,
                to: Some(0),
            },
        ],
    };
    let mut s = tube.state();
    let e0 = 0.5 * c;
    let (mut worst_v, mut worst_e) = (0.0_f64, 0.0_f64);
    for n in 1..=4000 {
        tube.step(&mut s);
        let (v, i) = (s.potentials[0], s.inductor_current[0]);
        assert_eq!(
            v.to_bits(),
            s.potentials[2].to_bits(),
            "V19: the strapped wire"
        );
        worst_v = worst_v.max((v - (f64::from(n) * theta).cos()).abs());
        worst_e = worst_e.max(((0.5 * c * v * v + 0.5 * l * i * i) / e0 - 1.0).abs());
    }
    println!(
        "V19: C = {c:.6} (C00 {:.6}, C02 {:.6}), omega = {omega:.6}; V against cos(n theta) to {worst_v:.1e}, energy to {worst_e:.1e}",
        c00, c02
    );
    assert!(worst_v < 1e-11, "V19 V {worst_v:.3e}");
    assert!(worst_e < 1e-11, "V19 energy {worst_e:.3e}");
}
