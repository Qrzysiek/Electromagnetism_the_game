//! Validation of the z-invariant electromagnetic field (`physics::fdtd2d`, PHYSICS.md
//! §2.11): F1 a cavity mode against the scheme's dispersion, F2 Gauss's law with moving
//! charges and a conductor, F3 a charged conductor's static field is stationary, F4 an
//! oscillating current element's radiation against the Hankel-function solution, F5 the
//! absorbing layer's reflection.

use physics::fdtd2d::{Fdtd, Grid};
use std::f64::consts::PI;

/// A line charge: position, velocity, charge.
type Particle = ((f64, f64), (f64, f64), f64);

fn empty(n: usize) -> Vec<bool> {
    vec![false; n]
}

/// F1: the TE₁₀ mode of a square perfectly conducting box (side L = 20, h = 0.5, c = 1):
/// `B_z ∝ cos(πx/L) cos ωt`. Started at its maximum (`E = 0`), it must oscillate at the
/// scheme's exact discrete frequency, `sin(ω Δt/2) = (c Δt/h) sin(k h/2)`, `k = π/L`:
/// measured from 40 periods' zero crossings, required 1e-4 (the crossings' timing). The
/// continuum's `cπ/L` differs by the scheme's dispersion (printed).
#[test]
fn f1_cavity_mode() {
    let (n, h, c) = (40, 0.5, 1.0);
    let grid = Grid {
        nx: n,
        ny: n,
        h,
        origin: (0.0, 0.0),
    };
    let mut f = Fdtd::new(grid, c, 0.9, 0, empty(grid.nodes()));
    let l = 20.0;
    let k = PI / l;
    for j in 0..n {
        for i in 0..n {
            #[allow(clippy::cast_precision_loss)]
            let x = (i as f64 + 0.5) * h;
            f.bz[i + n * j] = libm::cos(k * x);
        }
    }
    let discrete = 2.0 / f.dt * libm::asin(c * f.dt / h * libm::sin(k * h / 2.0));
    let period = 2.0 * PI / discrete;
    let mut crossings = Vec::new();
    let mut last = f.bz[0];
    while f.t < 40.5 * period {
        f.step();
        let now = f.bz[0];
        if last * now < 0.0 {
            // B_z is half a step behind E's time.
            let s = last / (last - now);
            crossings.push(f.t - 1.5 * f.dt + s * f.dt);
        }
        last = now;
    }
    #[allow(clippy::cast_precision_loss)]
    let measured =
        PI * (crossings.len() - 1) as f64 / (crossings[crossings.len() - 1] - crossings[0]);
    println!(
        "F1: omega {measured:.8} against the discrete {discrete:.8} ({:.1e}); continuum {:.8} ({:+.2e})",
        measured / discrete - 1.0,
        c * k,
        discrete / (c * k) - 1.0
    );
    assert!(
        (measured / discrete - 1.0).abs() < 1e-4,
        "F1 {measured} vs {discrete}"
    );
}

/// A metal disc (radius 5) at the centre of a 60 × 60 box (h = 1, c = 1).
fn disc_box() -> Fdtd {
    let grid = Grid {
        nx: 60,
        ny: 60,
        h: 1.0,
        origin: (0.0, 0.0),
    };
    let metal = (0..grid.nodes())
        .map(|k| {
            let (i, j) = (k % 61, k / 61);
            let (x, y) = grid.node(i, j);
            libm::hypot(x - 30.0, y - 30.0) <= 5.0
        })
        .collect();
    Fdtd::new(grid, 1.0, 0.9, 0, metal)
}

/// F2: Gauss's law. The disc at potential 1, ten line charges (±1) started in the static
/// field of the disc and of themselves, then moving on straight paths (up to 0.8 c,
/// crossing cells, reflected at the box's inner square) for 300 steps: at every node in
/// vacuum `∇·E/4π` must equal the charges' density (linear weights) to 1e-10 of its
/// largest value (rounding: the scheme conserves charge exactly).
#[test]
fn f2_gauss_law() {
    let mut f = disc_box();
    #[allow(clippy::cast_precision_loss)]
    let mut particles: Vec<Particle> = (0..10)
        .map(|k| {
            let a = k as f64 * 0.7;
            (
                (10.0 + 3.0 * k as f64, 8.0 + (k % 3) as f64 * 1.3),
                (0.8 * libm::cos(a), 0.8 * libm::sin(a)),
                if k % 2 == 0 { 1.0 } else { -1.0 },
            )
        })
        .collect();
    let density = |f: &Fdtd, ps: &[Particle]| {
        let mut rho = vec![0.0; f.grid.nodes()];
        for (x, _, q) in ps {
            f.add_charge(&mut rho, *q, x.0, x.1);
        }
        rho
    };
    let rho = density(&f, &particles);
    let phi = f.solve_static(
        |i, j| {
            if i == 0 || j == 0 || i == 60 || j == 60 {
                Some(0.0)
            } else if f.is_metal(i, j) {
                Some(1.0)
            } else {
                None
            }
        },
        &rho,
        1e-14,
    );
    f.set_static(&phi);
    let mut worst = 0.0_f64;
    for _ in 0..300 {
        for (x, v, q) in &mut particles {
            let mut next = (x.0 + v.0 * f.dt, x.1 + v.1 * f.dt);
            // Reflect inside [4, 56], outside the disc's neighbourhood.
            if !(4.0..=56.0).contains(&next.0) {
                v.0 = -v.0;
                next.0 = x.0 + v.0 * f.dt;
            }
            if !(4.0..=56.0).contains(&next.1) {
                v.1 = -v.1;
                next.1 = x.1 + v.1 * f.dt;
            }
            if libm::hypot(next.0 - 30.0, next.1 - 30.0) < 8.0 {
                *v = (-v.0, -v.1);
                next = (x.0 + v.0 * f.dt, x.1 + v.1 * f.dt);
            }
            f.deposit(*q, *x, next);
            *x = next;
        }
        f.step();
        let rho = density(&f, &particles);
        let peak = rho.iter().fold(0.0_f64, |m, r| m.max(r.abs()));
        for j in 1..60 {
            for i in 1..60 {
                let near_metal = [(0, 0), (1, 0), (0, 1), (2, 1), (1, 2), (1, 1)]
                    .iter()
                    .any(|(di, dj)| f.is_metal(i + di - 1, j + dj - 1));
                if near_metal {
                    continue;
                }
                worst = worst.max((f.gauss(i, j) - rho[i + 61 * j]).abs() / peak);
            }
        }
    }
    println!("F2: Gauss's law to {worst:.2e} of the peak density");
    assert!(worst < 1e-10, "F2 {worst:e}");
}

/// F3: the disc at potential 1 in the grounded box, no charges: the static field solved
/// on the grid is stationary under the scheme (its discrete curl vanishes): after 500
/// steps every E unchanged and B_z zero to 1e-12 of the field's largest value.
#[test]
fn f3_static_field_is_stationary() {
    let mut f = disc_box();
    let phi = f.solve_static(
        |i, j| {
            if i == 0 || j == 0 || i == 60 || j == 60 {
                Some(0.0)
            } else if f.is_metal(i, j) {
                Some(1.0)
            } else {
                None
            }
        },
        &vec![0.0; f.grid.nodes()],
        1e-15,
    );
    f.set_static(&phi);
    let (ex0, ey0) = (f.ex.clone(), f.ey.clone());
    let scale = ex0.iter().chain(&ey0).fold(0.0_f64, |m, v| m.max(v.abs()));
    for _ in 0..500 {
        f.step();
    }
    let de =
        f.ex.iter()
            .zip(&ex0)
            .chain(f.ey.iter().zip(&ey0))
            .fold(0.0_f64, |m, (a, b)| m.max((a - b).abs()));
    let b = f.bz.iter().fold(0.0_f64, |m, v| m.max(v.abs()));
    println!(
        "F3: E changed by {:.1e}, B_z {:.1e} of {scale:.3}",
        de / scale,
        b / scale
    );
    assert!(de < 1e-12 * scale && b < 1e-12 * scale, "F3 {de:e} {b:e}");
}

/// The field at a probe of a current element `K(t)` along x at the centre of a free
/// grid (`n` cells of h, a PML of `pml` cells): `B_z` there after every step.
fn current_element(
    n: usize,
    h: f64,
    pml: usize,
    k: impl Fn(f64) -> f64,
    probes: &[(usize, usize)],
    t_end: f64,
) -> (f64, Vec<Vec<(f64, f64)>>) {
    let grid = Grid {
        nx: n,
        ny: n,
        h,
        origin: (0.0, 0.0),
    };
    let mut f = Fdtd::new(grid, 1.0, 0.9, pml, empty(grid.nodes()));
    let mut out = vec![Vec::new(); probes.len()];
    let centre = n / 2;
    while f.t < t_end {
        // The element: a charge K dt / h moved by h along the central edge each step
        // carries K h... as a current density K/h² over one edge.
        let kk = k(f.t + 0.5 * f.dt);
        let (x0, y0) = grid.node(centre, centre);
        f.deposit(kk * f.dt / h, (x0, y0), (x0 + h, y0));
        f.step();
        for (o, &(i, j)) in out.iter_mut().zip(probes) {
            o.push((f.t - 0.5 * f.dt, f.bz[i + n * j]));
        }
    }
    (f.dt, out)
}

/// F4: an oscillating current element (moment `K(t) = K₀ sin ωt` per unit length, along
/// x, switched on smoothly over three periods) in free space (h = 0.25, 40 cells per
/// wavelength λ = 10, c = 1, a PML of 24 cells): after 20 periods `B_z` at distances
/// 1.5 λ and 2.5 λ along y (sin φ = 1) against the steady solution
/// `B_z = Re[(iπ K̃ k/c) H₁⁽¹⁾(kr) sin φ e^{−iωt}]`, `K̃ = i K₀`: its amplitude and phase
/// fitted over the last three periods. Required 3 % and 3° (the scheme's dispersion,
/// ~0.1 % of the phase per wavelength here, the element's one-edge size, and the 2D
/// tail of the switch-on).
#[test]
fn f4_current_element_radiates() {
    let (n, h, pml) = (320, 0.25, 24);
    let omega = 2.0 * PI / 10.0;
    let k_wave = omega;
    let ramp = 3.0 * 2.0 * PI / omega;
    let k = move |t: f64| {
        let s = if t < ramp {
            libm::sin(0.5 * PI * t / ramp).powi(2)
        } else {
            1.0
        };
        s * libm::sin(omega * t)
    };
    let centre = n / 2;
    // B_z cells straight above the element (cell centres at x = centre + 0.5).
    let rs = [15.0, 25.0];
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    let probes: Vec<(usize, usize)> = rs
        .iter()
        .map(|r| (centre, centre + (r / h) as usize))
        .collect();
    let t_end = 20.0 * 2.0 * PI / omega;
    let (_, series) = current_element(n, h, pml, k, &probes, t_end);
    for ((i, j), s) in probes.iter().zip(&series) {
        #[allow(clippy::cast_precision_loss)]
        let (dx, dy) = (
            (*i as f64 + 0.5 - centre as f64 - 0.5) * h,
            (*j as f64 + 0.5 - centre as f64) * h,
        );
        let r = libm::hypot(dx, dy);
        let sin_phi = dy / r;
        // Fit a cos ωt + b sin ωt over the last three periods.
        let window: Vec<&(f64, f64)> = s
            .iter()
            .filter(|(t, _)| *t > t_end - 3.0 * 2.0 * PI / omega)
            .collect();
        let (mut cc, mut ss, mut cs, mut yc, mut ys) = (0.0, 0.0, 0.0, 0.0, 0.0);
        for (t, y) in &window {
            let (c, sn) = (libm::cos(omega * t), libm::sin(omega * t));
            cc += c * c;
            ss += sn * sn;
            cs += c * sn;
            yc += y * c;
            ys += y * sn;
        }
        let det = cc * ss - cs * cs;
        let a = (yc * ss - ys * cs) / det;
        let b = (ys * cc - yc * cs) / det;
        // Exact: Re[A e^{−iωt}] = Re(A) cos ωt + Im(A) sin ωt, A = (iπ (iK₀) k) H₁(kr) sin φ.
        let (j1, y1) = (libm::j1(k_wave * r), libm::y1(k_wave * r));
        // iπ·i·k = −πk; times (j1 + i y1).
        let (are, aim) = (-PI * k_wave * j1 * sin_phi, -PI * k_wave * y1 * sin_phi);
        let (amp, amp_exact) = (libm::hypot(a, b), libm::hypot(are, aim));
        let phase = libm::atan2(b, a).to_degrees();
        let phase_exact = libm::atan2(aim, are).to_degrees();
        let mut dphi = phase - phase_exact;
        while dphi > 180.0 {
            dphi -= 360.0;
        }
        while dphi < -180.0 {
            dphi += 360.0;
        }
        println!(
            "F4: r = {r:.2}: amplitude {amp:.5e} against {amp_exact:.5e} ({:+.2} %), phase off by {dphi:+.2} degrees",
            100.0 * (amp / amp_exact - 1.0)
        );
        assert!(
            (amp / amp_exact - 1.0).abs() < 0.03,
            "F4 amplitude at r = {r}"
        );
        assert!(dphi.abs() < 3.0, "F4 phase at r = {r}");
    }
}

/// F5: the absorbing layer. A smooth current pulse at the centre of a 96-cell grid with a
/// 16-cell PML, against the same in a 480-cell grid (whose edges the pulse never reaches
/// in time): `B_z` at a probe 24 cells from the centre (8 from the layer) must agree to
/// 1e-3 of the pulse's peak there at all times (what returns from the layer).
#[test]
fn f5_absorbing_layer() {
    let h = 0.5;
    let width = 4.0;
    let pulse = move |t: f64| {
        let s = (t - 4.0 * width) / width;
        -s * libm::exp(-0.5 * s * s)
    };
    let t_end = 80.0;
    let small = current_element(96, h, 16, pulse, &[(48, 72)], t_end).1;
    let large = current_element(480, h, 16, pulse, &[(240, 264)], t_end).1;
    let peak = large[0].iter().fold(0.0_f64, |m, (_, v)| m.max(v.abs()));
    let worst = small[0]
        .iter()
        .zip(&large[0])
        .fold(0.0_f64, |m, ((_, a), (_, b))| m.max((a - b).abs()));
    println!("F5: reflection {:.2e} of the peak", worst / peak);
    assert!(worst < 1e-3 * peak, "F5 {worst:e} of {peak:e}");
}

/// F6: the energy the current element of F4 radiates, measured as the solver's
/// `outflow` through a square of side 30 around it (inside the PML), averaged over the
/// last three periods, against the cylindrical wave's `P̄ = π ω K₀²/(4c²)` per unit
/// length. Required 3 % (as F4's amplitude).
#[test]
fn f6_radiated_power() {
    let (n, h, pml) = (320, 0.25, 24);
    let omega = 2.0 * PI / 10.0;
    let ramp = 3.0 * 2.0 * PI / omega;
    let grid = Grid {
        nx: n,
        ny: n,
        h,
        origin: (0.0, 0.0),
    };
    let mut f = Fdtd::new(grid, 1.0, 0.9, pml, empty(grid.nodes()));
    let centre = n / 2;
    let (x0, y0) = grid.node(centre, centre);
    let t_end = 20.0 * 2.0 * PI / omega;
    let (mut energy, mut since) = (0.0, None);
    while f.t < t_end {
        let t = f.t + 0.5 * f.dt;
        let s = if t < ramp {
            libm::sin(0.5 * PI * t / ramp).powi(2)
        } else {
            1.0
        };
        let k = s * libm::sin(omega * t);
        f.deposit(k * f.dt / h, (x0, y0), (x0 + h, y0));
        f.step();
        if f.t > t_end - 3.0 * 2.0 * PI / omega {
            energy += f.outflow(x0 - 15.0, y0 - 15.0, x0 + 15.0, y0 + 15.0) * f.dt;
            since.get_or_insert(f.t);
        }
    }
    let mean = energy / (f.t - since.unwrap_or(0.0));
    let exact = PI * omega / 4.0;
    println!(
        "F6: radiated {mean:.6e} against {exact:.6e} ({:+.2} %)",
        100.0 * (mean / exact - 1.0)
    );
    assert!((mean / exact - 1.0).abs() < 0.03, "F6 {mean} vs {exact}");
}
