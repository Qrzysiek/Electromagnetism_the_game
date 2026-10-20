//! What a level's numbers mean (the owner, 2026-10-20): for each measurable quantity the
//! real-world order of magnitude where its phenomena occur, and the dimensionless ratios that
//! decide its physics, computed for the setup shown. One family of levels at a time: the
//! metals (tube levels with Drude metals, PHYSICS.md §2.11) first.

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
        "{} · Coulomb {} · a cell ≈ 1 mm",
        if c.is_finite() {
            format!("v/c = {}", num(speed / c))
        } else {
            "Newtonian".into()
        },
        num(coulomb)
    );
    let mut lines = vec![
        Scale {
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
        },
        Scale {
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
        },
    ];
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
    let headline = format!(
        "ω_c τ = {} · qEℓ/kT = {} · a cell ≈ 0.1 µm",
        num(wct),
        num(heating)
    );
    let lines = vec![
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
        Scale {
            quantity: "units: this model's metal as a real one",
            here: "a cell 0.1 µm; time 0.8 ps; B 0.47 T; potential 6.5 mV; speed 1.3e5 m/s".into(),
            real: "GaAs at room temperature (mobility 0.85 m²/Vs, 1e23 carriers per m³)",
            detail: "One consistent reading of the game's units, from its carriers' thermal speed, \
                     mean free path and charge-to-mass ratio set equal to GaAs's: then B = 4 is \
                     about 2 tesla, a strip 16 cells long 1.6 µm, a supply of 30 about 0.2 volt. \
                     The model is a scaled classical plasma closer to a doped semiconductor than \
                     to copper, whose numbers are in the line above each.",
        },
    ];
    Some(Scales { headline, lines })
}
