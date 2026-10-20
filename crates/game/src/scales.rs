//! What a level's numbers mean (the owner, 2026-10-20): for each measurable quantity the
//! real-world order of magnitude where its phenomena occur, and the dimensionless ratios that
//! decide its physics, computed for the setup shown: the particle levels (with the radiation
//! reaction's strength where it is modelled, a beam's own charge), the vacuum tubes and the
//! metals (with the real material each resembles).

use level::{Element, ElementKind, Level};

/// One line of a level's scales: the quantity, its value here, the real world's, and the
/// explanation shown on hover.
pub struct Scale {
    pub quantity: &'static str,
    pub here: String,
    pub real: &'static str,
    pub detail: &'static str,
}

/// A level's scales: the few numbers worth a glance (the section's header) and the lines.
pub struct Scales {
    pub headline: String,
    pub lines: Vec<Scale>,
}

/// The scales of the level shown with the player's `placement`; None for the families not
/// covered yet.
pub fn of(level: &Level, placement: &[Element]) -> Option<Scales> {
    match &level.tube {
        Some(spec) if !spec.metals.is_empty() => metals(level, spec, placement),
        Some(spec) => tubes(level, spec, placement),
        None => particles(level, placement),
    }
}

/// The electron's mass (kg), the vacuum's permittivity (F/m), kT at 300 K (J), ħ (J s): for the
/// materials' ratios (`resemblance`).
const ME_KG: f64 = 9.109e-31;
const EPS0: f64 = 8.854e-12;
const KT_300: f64 = 4.14e-21;
const HBAR: f64 = 1.0546e-34;

/// A real conductor at 300 K (Drude): its carriers' density (per m³), effective mass (of the
/// electron's), mobility (m²/Vs), the background's permittivity and whether its carriers are
/// degenerate (a metal's Fermi gas) or not (a doped semiconductor's).
struct Material {
    name: &'static str,
    short: &'static str,
    /// 1 everyday, 2 common, 3 specialised.
    tier: u8,
    n: f64,
    mass: f64,
    mobility: f64,
    eps: f64,
    degenerate: bool,
}

/// The conductors a metal is compared with (`resemblance`): textbook values, a few digits.
const MATERIALS: [Material; 14] = [
    Material {
        name: "copper",
        short: "copper",
        tier: 1,
        n: 8.47e28,
        mass: 1.0,
        mobility: 4.3e-3,
        eps: 1.0,
        degenerate: true,
    },
    Material {
        name: "aluminium",
        short: "aluminium",
        tier: 1,
        n: 1.81e29,
        mass: 1.0,
        mobility: 1.2e-3,
        eps: 1.0,
        degenerate: true,
    },
    Material {
        name: "iron",
        short: "iron",
        tier: 1,
        n: 1.7e29,
        mass: 1.0,
        mobility: 3.7e-4,
        eps: 1.0,
        degenerate: true,
    },
    Material {
        name: "stainless steel",
        short: "steel",
        tier: 1,
        n: 1.7e29,
        mass: 1.0,
        mobility: 5.1e-5,
        eps: 1.0,
        degenerate: true,
    },
    Material {
        name: "n-type silicon, 10¹⁶ per cm³",
        short: "silicon",
        tier: 1,
        n: 1e22,
        mass: 0.26,
        mobility: 0.135,
        eps: 11.7,
        degenerate: false,
    },
    Material {
        name: "n-type silicon, 10¹⁷ per cm³",
        short: "silicon",
        tier: 1,
        n: 1e23,
        mass: 0.26,
        mobility: 0.08,
        eps: 11.7,
        degenerate: false,
    },
    Material {
        name: "n-type silicon, 10¹⁸ per cm³",
        short: "silicon",
        tier: 1,
        n: 1e24,
        mass: 0.26,
        mobility: 0.03,
        eps: 11.7,
        degenerate: false,
    },
    Material {
        name: "silver",
        short: "silver",
        tier: 2,
        n: 5.86e28,
        mass: 1.0,
        mobility: 6.7e-3,
        eps: 1.0,
        degenerate: true,
    },
    Material {
        name: "gold",
        short: "gold",
        tier: 2,
        n: 5.9e28,
        mass: 1.0,
        mobility: 4.8e-3,
        eps: 1.0,
        degenerate: true,
    },
    Material {
        name: "nichrome",
        short: "nichrome",
        tier: 2,
        n: 1e29,
        mass: 1.0,
        mobility: 5.7e-5,
        eps: 1.0,
        degenerate: true,
    },
    Material {
        name: "n-type germanium, 10¹⁶ per cm³",
        short: "germanium",
        tier: 2,
        n: 1e22,
        mass: 0.12,
        mobility: 0.39,
        eps: 16.0,
        degenerate: false,
    },
    Material {
        name: "n-type GaAs, 10¹⁷ per cm³",
        short: "GaAs",
        tier: 2,
        n: 1e23,
        mass: 0.067,
        mobility: 0.6,
        eps: 12.9,
        degenerate: false,
    },
    Material {
        name: "n-type InAs, 10¹⁷ per cm³",
        short: "InAs",
        tier: 3,
        n: 1e23,
        mass: 0.023,
        mobility: 3.0,
        eps: 15.0,
        degenerate: false,
    },
    Material {
        name: "InSb (intrinsic)",
        short: "InSb",
        tier: 3,
        n: 2e22,
        mass: 0.014,
        mobility: 7.7,
        eps: 16.8,
        degenerate: false,
    },
];

/// A material's ω_p τ, screening length per mean free path, collision time (s) and mean free
/// path (m) (Material::ratios).
type Ratios = (f64, f64, f64, f64);

impl Material {
    /// Its plasma oscillations per collision, its screening length (Debye's, or a Fermi gas's
    /// Thomas-Fermi length) per mean free path, its collision time (s) and mean free path (m).
    fn ratios(&self) -> Ratios {
        let m = self.mass * ME_KG;
        let tau = self.mobility * m / E_SI;
        let wp = (self.n * E_SI * E_SI / (EPS0 * self.eps * m)).sqrt();
        let (screen, speed) = if self.degenerate {
            let fermi = HBAR
                * HBAR
                * libm::pow(
                    3.0 * std::f64::consts::PI * std::f64::consts::PI * self.n,
                    2.0 / 3.0,
                )
                / (2.0 * m);
            (
                (2.0 * EPS0 * self.eps * fermi / (3.0 * self.n * E_SI * E_SI)).sqrt(),
                (2.0 * fermi / m).sqrt(),
            )
        } else {
            (
                (EPS0 * self.eps * KT_300 / (self.n * E_SI * E_SI)).sqrt(),
                (KT_300 / m).sqrt() * 1.253_314_137_315_500_3,
            )
        };
        (wp * tau, screen / (speed * tau), tau, speed * tau)
    }
}

/// The real material a level's metal resembles (`resemblance`): its name, its ratios, how far
/// the metal's are from them (a factor), and its collision time (s), mean free path (m) and
/// mobility (m²/Vs).
pub struct Resemblance {
    pub name: &'static str,
    pub short: &'static str,
    pub wp_tau: f64,
    pub screen: f64,
    pub factor: f64,
    pub tau: f64,
    pub ell: f64,
    pub mobility: f64,
}

/// The real material a level's metal resembles (the owner, 2026-10-20): its own ratios, ω_p τ
/// and the screening length per mean free path (the same in any units), within ×3 (close) or
/// ×10 (roughly) of the material's (the geometric mean of the two factors), with the level's
/// strongest field reachable in the material (ω_c τ = μB at most 30 T). The most common
/// material first (`Material::tier`), close before rough, then the closest.
pub fn resemblance(level: &Level, m: &level::tube::MetalSpec) -> Option<Resemblance> {
    let qm = m.charge_per_mass.abs();
    let tau = m.mean_free_path / (m.thermal_speed * 1.253_314_137_315_500_3);
    let wp_tau = (4.0 * std::f64::consts::PI * m.density * qm).sqrt() * tau;
    let debye = (m.thermal_speed * m.thermal_speed / qm / (4.0 * std::f64::consts::PI * m.density))
        .sqrt()
        / m.mean_free_path;
    let field = level.tube.as_ref().map_or(0.0, |t| t.b_z.abs());
    let strongest = level
        .limits
        .solenoid_fields
        .iter()
        .fold(field, |a, v| a.max(v.abs()));
    let wct = qm * strongest * tau;
    let candidates: Vec<(&Material, f64, Ratios)> = MATERIALS
        .iter()
        .filter_map(|mat| {
            let r = mat.ratios();
            let factor =
                libm::exp(0.5 * (libm::log(wp_tau / r.0).abs() + libm::log(debye / r.1).abs()));
            (wct / mat.mobility <= 30.0).then_some((mat, factor, r))
        })
        .collect();
    for tier in 1..=3 {
        for limit in [3.0, 10.0] {
            if let Some((mat, factor, r)) = candidates
                .iter()
                .filter(|(mat, f, _)| mat.tier == tier && *f <= limit)
                .min_by(|a, b| a.1.total_cmp(&b.1))
            {
                return Some(Resemblance {
                    name: mat.name,
                    short: mat.short,
                    wp_tau: r.0,
                    screen: r.1,
                    factor: *factor,
                    tau: r.2,
                    ell: r.3,
                    mobility: mat.mobility,
                });
            }
        }
    }
    None
}

/// Elementary charge (C), electron rest energy (J), Coulomb constant (N m²/C²), light's speed
/// (m/s): for the real-world readings only (a few digits).
const E_SI: f64 = 1.602e-19;
const ME_C2_J: f64 = 8.187e-14;
const K_SI: f64 = 8.988e9;
const C_SI: f64 = 2.998e8;

/// The particle levels (PHYSICS.md §1: internal units, k = 1): the projectile read as an
/// electron and a cell as 1 mm. What decides a flight is dimensionless: its speed against
/// light's, the strongest charge's Coulomb energy at one cell against its kinetic energy, the
/// cell against its gyroradius in the strongest magnet's field there.
#[allow(clippy::similar_names)]
fn particles(level: &Level, placement: &[Element]) -> Option<Scales> {
    let shot = level.shots.first()?;
    let (q, mass, ke) = (
        shot.particle.charge.abs(),
        shot.particle.mass,
        shot.launch.kinetic_energy,
    );
    if !(q > 0.0 && mass > 0.0 && ke > 0.0) {
        return None;
    }
    let c = level.c();
    let gamma = if c.is_finite() {
        1.0 + ke / (mass * c * c)
    } else {
        1.0
    };
    let speed = if c.is_finite() {
        c * (1.0 - 1.0 / (gamma * gamma)).sqrt()
    } else {
        (2.0 * ke / mass).sqrt()
    };
    let momentum = gamma * mass * speed;
    let strongest = |kind: ElementKind, listed: &[f64]| {
        level
            .elements
            .iter()
            .chain(placement)
            .filter(|e| e.kind == kind)
            .map(|e| e.value.abs())
            .chain(listed.iter().map(|v| v.abs()))
            .fold(0.0, f64::max)
    };
    let charge = strongest(ElementKind::Charge, &level.limits.magnitudes);
    let magnet = strongest(ElementKind::Magnet, &level.limits.magnet_strengths);
    let coulomb = q * charge / ke;
    let bending = q * magnet / momentum;
    // The real reading: an electron at the level's v/c (Newtonian levels: at 10 keV), a cell
    // 1 mm.
    let kev = if c.is_finite() {
        (gamma - 1.0) * 511.0
    } else {
        10.0
    };
    let e_j = kev * 1.602e-16;
    let p_si = (e_j * e_j + 2.0 * e_j * ME_C2_J).sqrt() / C_SI;
    let v_si = p_si * C_SI * C_SI / (e_j + ME_C2_J);
    let num = crate::ui::fmt_num;
    let headline = format!(
        "{}{} · a cell ≈ 1 mm",
        if c.is_finite() {
            format!("v/c = {}", num(speed / c))
        } else {
            "Newtonian".into()
        },
        if charge > 0.0 {
            format!(" · Coulomb {}", num(coulomb))
        } else {
            String::new()
        }
    );
    let mut lines = vec![Scale {
        quantity: "the speed against light's, v/c",
        here: if c.is_finite() {
            format!(
                "{} (c = {} here); as an electron ≈ {} keV",
                num(speed / c),
                num(c),
                num(kev)
            )
        } else {
            "Newtonian (c = ∞): read as an electron at 10 keV".into()
        },
        real: "a TV tube's electrons 25 keV (v/c 0.3); an electron microscope's 100-300 keV \
                   (0.55-0.78); thermal electrons 0.03 eV",
        detail: "How fast the projectile goes against light. Below about 0.1 Newton's \
                     mechanics holds; near 1 its mass grows (γ), magnetic forces become as \
                     strong as electric ones and fields spread at a finite speed (retarded). \
                     The game's c is a level's choice: what matters is v/c.",
    }];
    if charge > 0.0 {
        lines.push(Scale {
            quantity: "the Coulomb push: k|qQ|/(L E_kin) at one cell",
            here: format!(
                "{} (the strongest charge, |Q| = {})",
                num(coulomb),
                num(charge)
            ),
            real: "about 1: a charge of a few nC 1 mm from a 10 keV electron (a charged bead, a \
                   deflector plate's charge)",
            detail: "The strongest charge's potential energy one cell away against the \
                     projectile's kinetic energy. Well below 1 the path barely bends; about 1 it \
                     turns hard; well above, the projectile cannot climb near a like charge (it \
                     turns back) and swings round an opposite one. Only this ratio matters: all \
                     charges and the energy scaled together give the same flight.",
        });
    }
    if magnet > 0.0 {
        let tesla = bending * p_si / (E_SI * 1e-3);
        lines.push(Scale {
            quantity: "the magnetic bending: a cell against the gyroradius",
            here: format!(
                "{} (the strongest magnet one cell away); ≈ {} T",
                num(bending),
                num(tesla)
            ),
            real: "a fridge magnet 5 mT; an MRI 1.5-3 T; a TV tube's deflection coils ~10 mT",
            detail: "How tightly the strongest magnet's field one cell away would curl the \
                     projectile: the cell over its gyroradius p/(qB). Well below 1 the path \
                     barely curves there; above 1 it circles. Read for that electron at 1 mm.",
        });
    }
    // A beam's own charge (PHYSICS.md §3.3): the energy the bunch's charge gives a particle at
    // its edge as it blows apart, k q² N/R (R the bunch's half-size), against the kinetic,
    // as the angle it turns the particle's path by (its square root).
    if let Some(b) = shot.beam.filter(|_| level.physics.beam_interaction) {
        let half = 0.5 * b.width.max(b.length).max(1e-9);
        let angle = (q * q * f64::from(b.count) / (half * ke)).sqrt();
        lines.push(Scale {
            quantity: "the beam's own charge: the angle it spreads the bunch by",
            here: format!(
                "{} rad, {}° ({} particles in a bunch about {} across)",
                num(angle),
                num(angle.to_degrees()),
                b.count,
                num(2.0 * half)
            ),
            real: "well below a milliradian in a TV tube's or a microscope's beam; tens of \
                   milliradians in a high-current gun or a dense ion bunch, which focusing \
                   must hold",
            detail: "How hard the beam's own charge pushes its particles apart: √(k q²N/(R \
                     E_kin)), the energy a particle at the bunch's edge gains as the bunch blows \
                     up against its kinetic energy, as the angle its path turns by. Over a long \
                     flight even a few degrees move the landing points by cells (measured on \
                     the beam levels: switching the interaction off moves them by 0.3-1.1 cells \
                     on average).",
        });
    }
    // The radiation reaction (PHYSICS.md §3.1): τ₀ = 2q²/3mc³ against the time to cross a
    // cell, here and for a real electron at the same speed (τ₀ = 6.26e-24 s, a cell 1 mm).
    if level.physics.radiation_reaction && c.is_finite() {
        let tau0 = 2.0 * q * q / (3.0 * mass * c * c * c);
        let ratio = tau0 * speed;
        let real = 6.26e-24 * v_si / 1e-3;
        lines.push(Scale {
            quantity: "the radiation reaction: τ₀ = 2q²/3mc³ against a cell's crossing",
            here: format!(
                "{} (a real electron at this speed: {}; the game's charges radiate {} times \
                 harder)",
                num(ratio),
                num(real),
                num(ratio / real)
            ),
            real: "an electron's τ₀ = 6.3e-24 s, the time light takes to cross two thirds of its \
                   classical radius: against a nanosecond's flight, a few parts in 10¹²",
            detail: "How strongly the projectile's own radiation brakes it: τ₀ against the time it \
                     takes to cross a cell. Real charges radiate very weakly on such scales (the \
                     reaction shows over many turns, as in a synchrotron, or over an atom's \
                     size); the game raises the charge against the mass so that it shows within \
                     one flight, and keeps the ratio small enough for the Landau–Lifshitz form \
                     of the reaction to hold (PHYSICS.md §3.1).",
        });
    }
    lines.push(Scale {
        quantity: "units: this level as a real one",
        here: format!(
            "a cell 1 mm; the electron crosses it in {} ns; a charge of 1 here ≈ {} nC",
            num(1e-3 / v_si * 1e9),
            num(1e-3 * e_j / ke / (K_SI * E_SI) * q * 1e9)
        ),
        real: "an electron beam steered by charged beads and magnets a few mm apart: a cathode-ray \
               tube, an electron lens",
        detail: "One consistent reading of the game's units: the projectile an electron at the \
                 level's v/c, a cell 1 mm. Every other unit follows (k = 1 inside the game); the \
                 flight itself depends only on the ratios above.",
    });
    Some(Scales { headline, lines })
}

/// The vacuum tubes (tube levels without metals, PHYSICS.md §2.11): electrons (q/m = −1) read
/// with a cell of 1 mm and a potential unit of 10 V, so a speed unit of `√((e/m) 10 V)` =
/// 1.33e6 m/s, a time unit of 0.75 ns and a field unit of 7.5 mT.
fn tubes(level: &Level, spec: &level::tube::TubeSpec, placement: &[Element]) -> Option<Scales> {
    const SPEED_SI: f64 = 1.326e6;
    const TIME_NS: f64 = 0.754;
    const FIELD_MT: f64 = 7.54;
    let qm = spec.charge_per_mass.abs();
    // The largest potential: the player's supplies and the level's held electrodes.
    let potential = placement
        .iter()
        .filter(|e| e.kind == ElementKind::Supply && Some(e.node) != spec.solenoid)
        .map(|e| e.value.abs())
        .chain(level.electrodes.iter().map(|e| match e.bias {
            level::ConductorBias::Potential(v) => v.abs(),
            _ => 0.0,
        }))
        .fold(0.0, f64::max);
    let b = spec
        .solenoid
        .and_then(|n| {
            placement
                .iter()
                .find(|e| e.kind == ElementKind::Supply && e.node == n)
        })
        .map_or(spec.b_z, |e| e.value)
        .abs();
    let speed = (2.0 * qm * potential).sqrt();
    let omega = spec
        .circuits
        .iter()
        .map(|c| c.source)
        .chain(spec.solenoid_drive)
        .filter_map(|src| match src {
            level::Source::Sine { omega, .. } => Some(omega),
            _ => None,
        })
        .fold(0.0, f64::max);
    let num = crate::ui::fmt_num;
    let mut headline = format!("{} V ≈ {} V real", num(potential), num(10.0 * potential));
    let mut lines = vec![Scale {
        quantity: "the electrons' speed from the largest potential",
        here: format!(
            "v = √(2|q/m|V) = {} ≈ {} m/s",
            num(speed),
            num(speed * SPEED_SI)
        ),
        real: "a radio tube's 100-300 V: 6-10e6 m/s; a TV tube's 25 kV: 9e7 m/s (0.3 c)",
        detail: "How fast an electron gets falling through the largest potential in the tube. \
                 The game's tubes are Newtonian: fine below about a tenth of light's speed \
                 (a few kV).",
    }];
    if b > 0.0 && speed > 0.0 {
        let radius = speed / (qm * b);
        headline = format!("{headline} · gyration {} cells", num(radius));
        lines.push(Scale {
            quantity: "the gyration in the field: r = v/(|q/m| B)",
            here: format!(
                "{} cells (B = {} ≈ {} mT)",
                num(radius),
                num(b),
                num(b * FIELD_MT)
            ),
            real: "a magnetron 0.1-0.3 T; a cyclotron 1-2 T; the Earth 0.05 mT",
            detail: "How tightly the field curls an electron moving at that speed. Much larger \
                     than the tube: the beam bends gently; about the gap: electrons turn back \
                     before the anode (a magnetron's cutoff); much smaller: they circle near \
                     the cathode.",
        });
    }
    if omega > 0.0 {
        let ghz = omega / (std::f64::consts::TAU * TIME_NS);
        headline = format!("{headline} · {} GHz", num(ghz));
        lines.push(Scale {
            quantity: "the radio frequency",
            here: format!("ω = {} ≈ {} GHz", num(omega), num(ghz)),
            real: "a microwave oven's magnetron 2.45 GHz; a radar klystron 1-10 GHz; FM radio \
                   0.1 GHz",
            detail: "The fastest source's frequency, read with the time unit 0.75 ns. Against \
                     the electrons' transit time across the tube it decides whether they see a \
                     steady field or a changing one (bunching, a klystron's and a magnetron's \
                     principle).",
        });
    }
    lines.push(Scale {
        quantity: "units: this tube as a real one",
        here: format!(
            "a cell 1 mm; a potential unit 10 V; time {TIME_NS} ns; B {FIELD_MT} mT; speed \
             1.33e6 m/s"
        ),
        real: "a vacuum tube or a magnetron on a millimetre scale",
        detail: "One consistent reading of the game's units: the electrons' charge-to-mass ratio \
                 as a real electron's, a cell 1 mm and a potential unit 10 V; the other units \
                 follow. Only the ratios above decide what happens.",
    });
    Some(Scales { headline, lines })
}

/// The metals (tube levels with Drude metals, PHYSICS.md §2.11).
fn metals(level: &Level, spec: &level::tube::TubeSpec, placement: &[Element]) -> Option<Scales> {
    let m = spec.metals.first()?;
    // The field: the player's solenoid, else the level's own; a swinging one's amplitude.
    let b = spec
        .solenoid
        .and_then(|n| {
            placement
                .iter()
                .find(|e| e.kind == ElementKind::Supply && e.node == n)
        })
        .map_or(spec.b_z, |e| e.value)
        .abs();
    let qm = m.charge_per_mass.abs();
    let sigma = m.thermal_speed;
    let ell = m.mean_free_path;
    // The mean speed of a two-dimensional Maxwellian, σ √(π/2), and the mean free time.
    let mean_speed = sigma * 1.253_314_137_315_500_3;
    let tau = ell / mean_speed;
    let wct = qm * b * tau;
    let kt_over_q = sigma * sigma / qm;
    let debye = (kt_over_q / (4.0 * std::f64::consts::PI * m.density)).sqrt();
    let wp_tau = (4.0 * std::f64::consts::PI * m.density * qm).sqrt() * tau;
    // The strongest driving field: the largest supply over the longest metal (a strip's
    // length between its contacts), or a swinging field's induced one at a metal's edge.
    let (lo, hi) = level.metal_rect(m);
    let length = (hi.x - lo.x).max(hi.y - lo.y);
    let width = (hi.x - lo.x).min(hi.y - lo.y);
    let supply = placement
        .iter()
        .filter(|e| e.kind == ElementKind::Supply && Some(e.node) != spec.solenoid)
        .map(|e| e.value.abs())
        .fold(0.0, f64::max);
    let field = match spec.solenoid_drive {
        Some(d) => {
            let omega = match d {
                level::Source::Sine { omega, .. } => omega,
                _ => 0.0,
            };
            0.5 * b * omega * 0.5 * length
        }
        None => supply / length,
    };
    let heating = field * ell / kt_over_q;
    let gyration = if b > 0.0 {
        sigma / (qm * b)
    } else {
        f64::INFINITY
    };
    let num = crate::ui::fmt_num;
    let like = resemblance(level, m);
    let headline = format!(
        "ω_c τ = {} · qEℓ/kT = {}{}",
        num(wct),
        num(heating),
        like.as_ref()
            .map_or(String::new(), |r| format!(" · like {}", r.short))
    );
    let mut lines = vec![
        Scale {
            quantity: "ω_c τ, the field's turn between collisions",
            here: format!(
                "{} (B = {}: circles {} across, the mean free path {})",
                num(wct),
                num(b),
                if gyration.is_finite() {
                    num(2.0 * gyration)
                } else {
                    "∞".into()
                },
                num(ell)
            ),
            real: "copper at room temperature 0.004 per tesla; GaAs 0.85 per T; InSb 7.7 per T; \
                   a clean 2D electron gas at 4 K ~1000 per T",
            detail: "How far a carrier's path turns in the magnetic field between two collisions \
                     (ω_c τ = μB, μ the mobility). Below 0.1 the paths are straight and the Hall \
                     voltage grows in proportion to B; near 1 the current turns aside (the Hall \
                     angle, 45° at 1); above a few the carriers circle on the spot between \
                     collisions. Real metals stay far below 1 even at the strongest laboratory \
                     fields; semiconductors reach it at a fraction of a tesla.",
        },
        Scale {
            quantity: "qEℓ/kT, the field's push against the thermal motion",
            here: num(heating),
            real: "copper at 1 V/m ~1e-6 (always linear); semiconductors heat up near 1 at \
                   ~1 kV/cm",
            detail: "The energy a carrier gains from the electric field over one mean free path, \
                     against its thermal energy. Well below 1, Ohm's law holds (the current in \
                     proportion to the field); near 1 the carriers warm and the conductance \
                     falls (hot carriers). The game's strips sit at 0.05-0.6.",
        },
        Scale {
            quantity: "the strip's width in mean free paths",
            here: num(width / ell),
            real: "a copper wire 1 mm thick: 25 000; a thin film or nanowire: ~1",
            detail: "How many mean free paths fit across the metal. Many: its bulk conducts as \
                     Ohm's law says; a few: its walls matter (size effects, the heat of eddy \
                     currents a few per cent below Ohm's law).",
        },
        Scale {
            quantity: "the screening (Debye) length in mean free paths",
            here: num(debye / ell),
            real: "copper 0.001 (screening within an atom); doped GaAs (1e23 per m³) ~0.1",
            detail: "How deep a charge piled on a surface reaches into the metal: the thin layer \
                     that holds a conductor's surface charge, the Hall voltage's charge, a \
                     contact's boundary layer.",
        },
        Scale {
            quantity: "ω_p τ, plasma oscillations per collision",
            here: num(wp_tau),
            real: "copper ~250; doped semiconductors ~1-10",
            detail: "How many times the carriers' charge could oscillate (the plasma frequency) \
                     before a collision damps it: a charge imbalance relaxes in a fraction of a \
                     mean free time, so the metal stays neutral inside.",
        },
    ];
    let omega = match spec.solenoid_drive {
        Some(level::Source::Sine { omega, .. }) => omega,
        _ => 0.0,
    };
    lines.push(match &like {
        Some(r) => {
            // The model's τ and ℓ read as the material's; B from ω_c τ = μB.
            let (t_unit, l_unit) = (r.tau / tau, r.ell / ell);
            Scale {
                quantity: "the real material it resembles",
                here: format!(
                    "{} ({} of its own ratios: ω_p τ {}, the screening length {} mean free paths): \
                     a cell {} µm; time {} ps; B {} T{}",
                    r.name,
                    if r.factor <= 3.0 { "within ×3" } else { "within ×10" },
                    num(r.wp_tau),
                    num(r.screen),
                    num(l_unit * 1e6),
                    num(t_unit * 1e12),
                    num(b * qm * tau / r.mobility),
                    if omega > 0.0 {
                        format!(
                            "; the swing {} GHz",
                            num(omega / t_unit / (2.0 * std::f64::consts::PI) * 1e-9)
                        )
                    } else {
                        String::new()
                    }
                ),
                real: "compared: copper, aluminium, iron, stainless steel, silicon (1e16-1e18 per \
                       cm³), silver, gold, nichrome, germanium, GaAs, InAs, InSb at 300 K. Real \
                       metals have ω_p τ 50-500, beyond the model at playable cost (about 1: its \
                       step must follow the plasma oscillation)",
                detail: "The metal's own ratios, ω_p τ (plasma oscillations per collision) and its \
                         screening length in mean free paths, against real conductors' at room \
                         temperature (Drude, from their density, effective mass and mobility); \
                         a material counts only if the level's strongest field, read through its \
                         mobility (ω_c τ = μB), is a laboratory's (at most 30 T). Named: the most \
                         common such material within ×3 of both ratios, else within ×10. The rest \
                         of the line reads the game's units as that material's: its collision \
                         time and mean free path for the model's.",
            }
        }
        None => Scale {
            quantity: "the real material it resembles",
            here: "none of the materials compared within ×10 of its ratios at a laboratory's field"
                .into(),
            real: "compared: copper, aluminium, iron, stainless steel, silicon (1e16-1e18 per cm³), \
                   silver, gold, nichrome, germanium, GaAs, InAs, InSb at 300 K",
            detail: "No real conductor among those compared has this metal's own ratios (ω_p τ and \
                     the screening length in mean free paths, within ×10) and reaches the level's \
                     strongest field at most 30 T.",
        },
    });
    Some(Scales { headline, lines })
}
