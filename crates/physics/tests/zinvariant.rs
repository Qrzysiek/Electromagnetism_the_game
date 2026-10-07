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
            }],
            projectiles: Vec::new(),
            softening: 0.5 * size,
            dt,
            arena: None,
            b_z: 0.0,
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
#[ignore = "4.5 minutes in release: run with --ignored (measured values in PHYSICS.md §2.11)"]
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
                },
                Emitter {
                    electrode: 1,
                    charge_per_mass: 1.0,
                    weight,
                    emit_toward: None,
                },
            ],
            projectiles: Vec::new(),
            softening: 0.5 * size,
            dt,
            arena: None,
            b_z: 0.0,
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
