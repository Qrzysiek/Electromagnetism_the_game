//! Inspection of what lies under the pointer (the owner: every placed object, the
//! player's or the level's, shows its kind and every parameter that defines it, and a
//! time-dependent source its signal): the tooltip's content.

use level::{Coil, ConductorBias, Element, ElementKind, Level, Node, Source, tube::TubeCircuit};
use physics::DVec3;

use crate::ui::{fmt_num, fmt_si};

/// What the tooltip shows.
pub struct Inspection {
    pub title: String,
    pub lines: Vec<String>,
    /// A time-dependent source: what is plotted, and its samples `(t, value)`.
    pub chart: Option<(String, Vec<(f64, f64)>)>,
}

impl Inspection {
    fn new(title: impl Into<String>) -> Self {
        Self {
            title: title.into(),
            lines: Vec::new(),
            chart: None,
        }
    }

    fn line(mut self, s: impl Into<String>) -> Self {
        self.lines.push(s.into());
        self
    }
}

/// How close to a point element the pointer must be (cells).
const POINT_RADIUS: f64 = 0.4;
/// How close to a coil's wire (cells).
const WIRE_RADIUS: f64 = 0.3;
/// Samples of a plotted signal.
const CHART_SAMPLES: usize = 240;

/// The object under the pointer at `at` (board coordinates), the player's `placement`
/// included: point elements first, then bodies, coils, electrodes, boxes and detectors.
#[must_use]
pub fn inspect_at(level: &Level, placement: &[Element], at: DVec3) -> Option<Inspection> {
    let p = DVec3::new(at.x, at.y, 0.0);
    let pos = |n: Node| {
        let v = level.grid.position(n);
        DVec3::new(v.x, v.y, 0.0)
    };
    let near = |n: Node, r: f64| (pos(n) - p).length() <= r;
    // The player's cuts, anywhere along their slits.
    for e in placement.iter().filter(|e| e.kind == ElementKind::Cut) {
        if let Some((a, b)) = level.cut_segment(e) {
            let (d, r) = (b - a, p - DVec3::new(a.x, a.y, 0.0));
            let u = (r.dot(d) / d.length_squared()).clamp(0.0, 1.0);
            if (r - d * u).length() <= POINT_RADIUS {
                return Some(
                    Inspection::new("Cut (yours): an insulating slit")
                        .line(format!(
                            "{} long, at {}° from x",
                            fmt_num(e.value),
                            fmt_num(e.angle_deg)
                        ))
                        .line("the metal's carriers cannot cross it")
                        .line("a current loop must go round it: eddy currents weaken"),
                );
            }
        }
    }
    // Point elements: the level's, then the player's (a supply sits on its electrode or
    // its coil's terminal: described with it below).
    for (e, mine) in level
        .elements
        .iter()
        .map(|e| (e, false))
        .chain(placement.iter().map(|e| (e, true)))
    {
        if matches!(
            e.kind,
            ElementKind::Plate | ElementKind::Supply | ElementKind::Cut
        ) {
            continue;
        }
        if near(e.node, POINT_RADIUS) {
            return Some(element(level, e, mine));
        }
    }
    // Launch points of the shots, the level's free particles and the tube's projectiles.
    for (k, s) in level.shots.iter().enumerate() {
        if near(s.launch.node, POINT_RADIUS) {
            let d = s.launch.direction;
            let mut i = Inspection::new(format!("Shot {}: launch point", k + 1))
                .line(particle(&s.particle))
                .line(format!(
                    "kinetic energy {}, toward {}° (in the plane)",
                    fmt_si(s.launch.kinetic_energy),
                    fmt_num(d[1].atan2(d[0]).to_degrees())
                ));
            if s.launch.time != 0.0 {
                i = i.line(format!("launched at t = {}", fmt_num(s.launch.time)));
            }
            if let Some(b) = &s.beam {
                i = i.line(format!("a beam: {b:?}"));
            }
            return Some(i);
        }
    }
    for (k, f) in level.free_particles.iter().enumerate() {
        if near(f.node, POINT_RADIUS) {
            return Some(
                Inspection::new(format!("Free particle {} (the level's)", k + 1))
                    .line(particle(&f.particle))
                    .line(format!(
                        "launch velocity ({}, {})",
                        fmt_num(f.velocity[0]),
                        fmt_num(f.velocity[1])
                    )),
            );
        }
    }
    if let Some(t) = &level.tube {
        for (k, pr) in t.projectiles.iter().enumerate() {
            if near(pr.node, POINT_RADIUS) {
                return Some(
                    Inspection::new(format!("Projectile {} (a line charge)", k + 1))
                        .line(format!(
                            "charge {}, mass {} (per unit length)",
                            fmt_num(pr.charge),
                            fmt_num(pr.mass)
                        ))
                        .line(format!(
                            "launched at t = {} with speed {} toward {}°",
                            fmt_num(pr.launch),
                            fmt_num(pr.speed),
                            fmt_num(pr.direction_deg)
                        )),
                );
            }
        }
        if let Some(n) = t.solenoid
            && near(n, POINT_RADIUS)
        {
            let supply = placement
                .iter()
                .find(|e| e.kind == ElementKind::Supply && e.node == n);
            return Some(
                Inspection::new("Solenoid around the tube (its terminal)")
                    .line("a uniform magnetic field B along z, out of the board when positive")
                    .line(match supply {
                        Some(s) => format!("your supply: B = {}", fmt_num(s.value)),
                        None => format!("off: B = {} (the level's)", fmt_num(t.b_z)),
                    }),
            );
        }
    }
    // Spheres and clouds.
    for (k, c) in level.conductors.iter().enumerate() {
        if (pos(c.center) - p).length() <= c.radius {
            return Some(
                Inspection::new(format!("Metal sphere {}", k + 1))
                    .line(format!("radius {}", fmt_num(c.radius)))
                    .line(match c.bias {
                        ConductorBias::Grounded => "grounded".to_string(),
                        ConductorBias::Charge(q) => format!("isolated, net charge {}", fmt_si(q)),
                        ConductorBias::Potential(v) => format!("held at potential {}", fmt_si(v)),
                    }),
            );
        }
    }
    for (k, d) in level.dielectric_spheres.iter().enumerate() {
        if (pos(d.center) - p).length() <= d.radius {
            return Some(
                Inspection::new(format!("Dielectric sphere {}", k + 1))
                    .line(format!("radius {}", fmt_num(d.radius)))
                    .line(format!("relative permittivity {}", fmt_num(d.permittivity))),
            );
        }
    }
    for (k, c) in level.clouds.iter().enumerate() {
        if (pos(c.center) - p).length() <= c.radius {
            return Some(
                Inspection::new(format!("Charge cloud {}", k + 1))
                    .line(format!(
                        "radius {}, total charge {}",
                        fmt_num(c.radius),
                        fmt_si(c.charge)
                    ))
                    .line("uniform charge that particles fly through"),
            );
        }
    }
    // Coils: on the wire.
    for (k, c) in level.coils.iter().enumerate() {
        let on_wire = match c {
            Coil::Circle { center, radius, .. } => {
                ((pos(*center) - p).length() - radius).abs() <= WIRE_RADIUS
            }
            Coil::Polygon { vertices, .. } => vertices.iter().enumerate().any(|(i, v)| {
                let (a, b) = (pos(*v), pos(vertices[(i + 1) % vertices.len()]));
                segment_distance(a, b, p) <= WIRE_RADIUS
            }),
        };
        if on_wire {
            return Some(coil(level, placement, k, c));
        }
    }
    // Electrodes (the level's, then the player's plates).
    let outlines = level.electrode_outlines(placement);
    let n = level.electrodes.len();
    if let Some(i) = outlines.iter().position(|o| o.distance_to(p) == 0.0) {
        if i < n {
            return Some(electrode(level, placement, i));
        }
        let plate = placement
            .iter()
            .filter(|e| e.kind == ElementKind::Plate)
            .nth(i - n)?;
        let s = level.limits.plate;
        return Some(
            Inspection::new("Plate (yours)")
                .line(format!(
                    "{} long, {} thick, {} high, at {}°",
                    fmt_num(s.length),
                    fmt_num(s.thickness),
                    fmt_num(s.height),
                    fmt_num(plate.angle_deg)
                ))
                .line(if plate.value == 0.0 {
                    "grounded".to_string()
                } else {
                    format!("held at potential {}", fmt_si(plate.value))
                }),
        );
    }
    // Dielectric and ferrite boxes.
    for (k, d) in level.dielectrics.iter().enumerate() {
        if in_box(pos(d.center), d.length, d.thickness, d.angle_deg, p) {
            return Some(
                Inspection::new(format!("Dielectric slab {}", k + 1))
                    .line(size(d.length, d.thickness, d.height, d.angle_deg))
                    .line(format!("relative permittivity {}", fmt_num(d.permittivity))),
            );
        }
    }
    for (k, f) in level.ferrites.iter().enumerate() {
        if in_box(pos(f.center), f.length, f.thickness, f.angle_deg, p) {
            return Some(
                Inspection::new(format!("Ferrite slab {}", k + 1))
                    .line(size(f.length, f.thickness, f.height, f.angle_deg))
                    .line(format!("relative permeability {}", fmt_num(f.permeability))),
            );
        }
    }
    // A tube's metals (a contact is described with its electrode above).
    if let Some(t) = &level.tube {
        for (k, m) in t.metals.iter().enumerate() {
            let (lo, hi) = level.metal_rect(m);
            if p.x >= lo.x && p.x <= hi.x && p.y >= lo.y && p.y <= hi.y {
                return Some(metal(level, k, m, lo, hi));
            }
        }
    }
    // Detectors and gates.
    let inside = |d: &level::Detector| {
        let (a, b) = (pos(d.min), pos(d.max));
        p.x >= a.x.min(b.x) && p.x <= a.x.max(b.x) && p.y >= a.y.min(b.y) && p.y <= a.y.max(b.y)
    };
    for (k, s) in level.shots.iter().enumerate() {
        if inside(&s.detector) {
            return Some(detector(format!("Detector of shot {}", k + 1), &s.detector));
        }
    }
    for (k, g) in level.gates.iter().enumerate() {
        if inside(g) {
            return Some(detector(format!("Gate {}", k + 1), g));
        }
    }
    None
}

/// A metal (the Drude model): its rectangle, carriers, cores and contacts.
fn metal(level: &Level, k: usize, m: &level::tube::MetalSpec, lo: DVec3, hi: DVec3) -> Inspection {
    let area = (hi.x - lo.x) * (hi.y - lo.y);
    let contacts: Vec<String> = m.contacts.iter().map(|e| (e + 1).to_string()).collect();
    Inspection::new(format!("Metal {} (the Drude model)", k + 1))
        .line(format!(
            "{} × {}, from ({}, {}) to ({}, {})",
            fmt_num(hi.x - lo.x),
            fmt_num(hi.y - lo.y),
            fmt_num(lo.x),
            fmt_num(lo.y),
            fmt_num(hi.x),
            fmt_num(hi.y)
        ))
        .line(format!(
            "carriers: charge density {}, q/m = {}, thermal speed {}",
            fmt_num(m.density),
            fmt_num(m.charge_per_mass),
            fmt_num(m.thermal_speed)
        ))
        .line(format!(
            "ion cores: radius {}, {} per area ({:.0} in all): mean free path {}",
            fmt_num(m.core_radius),
            fmt_num(m.core_density()),
            area * m.core_density(),
            fmt_num(m.mean_free_path)
        ))
        .line(if contacts.is_empty() {
            "no contacts: isolated (its charge stays in it)".into()
        } else {
            format!(
                "contacts: electrode{} {} (each keeps a layer {} deep neutral)",
                if contacts.len() == 1 { "" } else { "s" },
                contacts.join(", "),
                fmt_num(level.contact_depth(m))
            )
        })
        .line(format!(
            "{} macroparticles per area in the preview; plasma frequency {}",
            fmt_num(m.carriers_per_area),
            fmt_num(m.plasma_frequency())
        ))
        .line(match crate::scales::resemblance(level, m) {
            Some(r) => format!(
                "resembles {} ({} of its ratios; the Scales section)",
                r.name,
                if r.factor <= 3.0 {
                    "within ×3"
                } else {
                    "within ×10"
                }
            ),
            None => "resembles none of the conductors compared (the Scales section)".into(),
        })
}

fn particle(s: &level::ParticleSpec) -> String {
    let mut out = format!(
        "charge {}, mass {}, radius {}",
        fmt_si(s.charge),
        fmt_si(s.mass),
        fmt_num(s.radius)
    );
    if s.moment != 0.0 {
        out.push_str(&format!(", magnetic moment {}", fmt_si(s.moment)));
    }
    out
}

fn size(length: f64, thickness: f64, height: f64, angle_deg: f64) -> String {
    format!(
        "{} long, {} thick, {} high, at {}°",
        fmt_num(length),
        fmt_num(thickness),
        fmt_num(height),
        fmt_num(angle_deg)
    )
}

fn detector(title: String, d: &level::Detector) -> Inspection {
    let mut i = Inspection::new(title).line(format!(
        "box from ({}, {}) to ({}, {})",
        d.min[0], d.min[1], d.max[0], d.max[1]
    ));
    if let Some(a) = &d.acceptance {
        if let Some([axis, half]) = a.direction {
            i = i.line(format!(
                "accepts directions within {}° of {}°",
                fmt_num(half),
                fmt_num(axis)
            ));
        }
        if let Some([lo, hi]) = a.kinetic {
            i = i.line(format!(
                "accepts kinetic energies {} to {}",
                fmt_si(lo),
                fmt_si(hi)
            ));
        }
        if a.radiation.is_some() {
            i = i.line("and a radiation goal (see the goals)");
        }
    }
    i
}

/// A point element: its kind and parameters; an antenna's dipole moment over time.
fn element(level: &Level, e: &Element, mine: bool) -> Inspection {
    let whose = if mine { "yours" } else { "the level's" };
    match e.kind {
        ElementKind::Charge => {
            Inspection::new(format!("Charge ({whose})")).line(format!("Q = {}", fmt_si(e.value)))
        }
        ElementKind::Magnet => Inspection::new(format!("Magnet ({whose})"))
            .line(format!("moment μ = {} along +z", fmt_si(e.value)))
            .line("its field in the plane: B_z = −μ/r³"),
        ElementKind::Antenna => {
            let omega = e.omega.unwrap_or(level.physics.rf_omega);
            let mut i = Inspection::new(format!("Antenna ({whose}): an oscillating dipole"))
                .line(format!(
                    "amplitude p₀ = {}, along {}°",
                    fmt_si(e.value),
                    fmt_num(e.angle_deg)
                ))
                .line(format!(
                    "ω = {}{}",
                    fmt_num(omega),
                    if e.omega.is_none() {
                        " (the level's generator)"
                    } else {
                        ""
                    }
                ));
            if omega > 0.0 {
                let end = 3.0 * std::f64::consts::TAU / omega;
                i.chart = Some((
                    "p(t) = p₀ cos ωt".into(),
                    sample(end, |t| e.value * (omega * t).cos()),
                ));
            }
            i
        }
        ElementKind::Free => Inspection::new(format!("Free charge ({whose})"))
            .line(format!(
                "charge {}, mass {}, radius {}",
                fmt_si(e.value),
                fmt_si(level.limits.free_mass),
                fmt_num(level.limits.free_radius)
            ))
            .line(format!(
                "thrown toward {}° at speed {}",
                fmt_num(e.angle_deg),
                fmt_num(e.speed.unwrap_or(0.0))
            )),
        ElementKind::Plate | ElementKind::Supply | ElementKind::Cut => Inspection::new("Element"),
    }
}

/// A coil: its shape, strength, ramp or circuit, and the player's supply on it.
fn coil(level: &Level, placement: &[Element], k: usize, c: &Coil) -> Inspection {
    let supply = c.supply_node().and_then(|n| {
        placement
            .iter()
            .find(|e| e.kind == ElementKind::Supply && e.node == n)
    });
    let t_max = level.physics.t_max;
    let (mut i, kappa, rate, drive) = match c {
        Coil::Circle {
            radius,
            kappa,
            rate,
            drive,
            tunable,
            ..
        } => {
            let mut i = Inspection::new(format!("Coil {} (a circle)", k + 1))
                .line(format!("radius {}", fmt_num(*radius)));
            if *tunable {
                i = i.line(match supply {
                    Some(s) => format!("your supply on its terminal: {}", fmt_num(s.value)),
                    None => "tunable: no supply on its terminal".to_string(),
                });
            }
            (i, *kappa, *rate, *drive)
        }
        Coil::Polygon {
            vertices,
            kappa,
            rate,
            drive,
        } => (
            Inspection::new(format!("Coil {} ({} corners)", k + 1, vertices.len())),
            *kappa,
            *rate,
            *drive,
        ),
    };
    i = i.line(format!(
        "strength κ = μ₀I/4π = {} at t = 0 (counter-clockwise when positive)",
        fmt_si(kappa)
    ));
    // A supply sets the ramp rate, or a driven coil's source level.
    let rate = match (supply, drive) {
        (Some(s), None) => s.value,
        _ => rate,
    };
    if let Some(d) = drive {
        let source = match supply {
            Some(s) => d.source.with_level(s.value),
            None => d.source,
        };
        i = i.line(format!("driven: {}", source_text(&source)));
        i = i.line(drive_text(&d, false));
        i.chart = chart(&source, t_max);
    } else if rate != 0.0 {
        i = i.line(format!("ramped: κ + {} t", fmt_si(rate)));
        i.chart = Some(("κ(t)".into(), sample(t_max, |t| kappa + rate * t)));
    }
    i
}

/// A level electrode: shape, bias, the player's supply, its circuit (a 3D level's drive
/// or a tube's circuit) and the signal of a time-dependent source.
fn electrode(level: &Level, placement: &[Element], k: usize) -> Inspection {
    let e = &level.electrodes[k];
    let shape = match e.shape {
        level::ElectrodeShape::Box => format!(
            "a box {} long, {} thick, {} high, at {}°",
            fmt_num(e.length),
            fmt_num(e.thickness),
            fmt_num(e.height),
            fmt_num(e.angle_deg)
        ),
        level::ElectrodeShape::Rod => format!("a rod of diameter {}", fmt_num(e.length)),
        level::ElectrodeShape::Ring => format!(
            "a ring of mid-diameter {}, wall {}",
            fmt_num(e.length),
            fmt_num(e.thickness)
        ),
    };
    let mut i = Inspection::new(format!("Electrode {}", k + 1)).line(shape);
    let supply = placement
        .iter()
        .find(|s| s.kind == ElementKind::Supply && e.tunable && s.node == e.center);
    let circuit = e.drive.is_some()
        || level.tube.as_ref().is_some_and(|t| {
            t.circuits
                .iter()
                .any(|c| c.electrode == k || c.linked.contains(&k))
        });
    let without = if circuit || e.tunable {
        " (without a supply or circuit)"
    } else {
        ""
    };
    i = i.line(match e.bias {
        ConductorBias::Grounded => format!("grounded{without}"),
        ConductorBias::Charge(q) => format!("isolated, net charge {}{without}", fmt_si(q)),
        ConductorBias::Potential(v) => format!("at potential {}{without}", fmt_si(v)),
    });
    if e.tunable {
        i = i.line(match supply {
            Some(s) => format!("your supply: {}", fmt_num(s.value)),
            None => "tunable: no supply placed".to_string(),
        });
    }
    let t_max = level.physics.t_max;
    if let Some(d) = e.drive {
        let source = match supply {
            Some(s) => d.source.with_level(s.value),
            None => d.source,
        };
        i = i.line(format!("driven: {}", source_text(&source)));
        i = i.line(drive_text(&d, true));
        i.chart = chart(&source, t_max);
    }
    if let Some(t) = &level.tube {
        if k == t.cathode && t.metal_of(k).is_none() {
            i = i.line(format!(
                "the cathode: emits carriers of q/m = {}{}",
                fmt_num(t.charge_per_mass),
                t.saturation.map_or(String::new(), |j| format!(
                    ", at most {} per unit length",
                    fmt_num(j)
                ))
            ));
        }
        if let Some(m) = t.metal_of(k) {
            i = i.line(format!(
                "a contact of metal {}: takes in the carriers that reach it, gives off what \
                 keeps the metal next to it neutral",
                m + 1
            ));
        }
        if let Some(ions) = t.ions.filter(|ions| ions.electrode == k) {
            i = i.line(format!(
                "emits ions of q/m = {}",
                fmt_num(ions.charge_per_mass)
            ));
        }
        for c in t
            .circuits
            .iter()
            .filter(|c| c.electrode == k || c.linked.contains(&k))
        {
            let source = tube_source(c, supply.map(|s| s.value));
            if c.electrode != k {
                i = i.line(format!("wired to electrode {}'s circuit:", c.electrode + 1));
            }
            i = i.line(format!("circuit: {}", source_text(&source)));
            let mut parts = Vec::new();
            match c.resistance {
                None => parts.push("floating".to_string()),
                Some(0.0) => parts.push("ideal source".to_string()),
                Some(r) => parts.push(format!("R = {}", fmt_num(r))),
            }
            if c.inductance > 0.0 {
                parts.push(format!("L = {}", fmt_num(c.inductance)));
            }
            if c.capacitance > 0.0 {
                parts.push(format!("C = {} to ground", fmt_num(c.capacitance)));
            }
            if let Some(r) = c.return_to {
                parts.push(format!("returns to electrode {}", r + 1));
            }
            i = i.line(parts.join(", "));
            if i.chart.is_none() {
                i.chart = chart(&source, t_max);
            }
        }
    }
    i
}

/// A tube circuit's source with the player's supply applied (as `Level::tube_circuits`).
fn tube_source(c: &TubeCircuit, supply: Option<f64>) -> Source {
    match (supply, c.source) {
        (
            Some(v),
            Source::Sine {
                amplitude,
                omega,
                phase,
                ..
            },
        ) if c.supply_offset => Source::Sine {
            offset: v,
            amplitude,
            omega,
            phase,
        },
        (Some(v), s) => s.with_level(v),
        (None, s) => s,
    }
}

fn source_text(s: &Source) -> String {
    match *s {
        Source::Dc { value } => format!("constant {}", fmt_num(value)),
        Source::Sine {
            offset,
            amplitude,
            omega,
            phase,
        } => {
            let offset = if offset == 0.0 {
                String::new()
            } else {
                format!("{} + ", fmt_num(offset))
            };
            let phase = if phase == 0.0 {
                String::new()
            } else {
                format!(" + {}", fmt_num(phase))
            };
            format!(
                "{offset}{} sin({} t{phase}), period {}",
                fmt_num(amplitude),
                fmt_num(omega),
                fmt_num(if omega > 0.0 {
                    std::f64::consts::TAU / omega
                } else {
                    f64::INFINITY
                })
            )
        }
        Source::Pulse {
            low,
            high,
            delay,
            rise,
            width,
            fall,
            period,
        } => format!(
            "pulse {} → {} after {} (rise {}, width {}, fall {}){}",
            fmt_num(low),
            fmt_num(high),
            fmt_num(delay),
            fmt_num(rise),
            fmt_num(width),
            fmt_num(fall),
            if period > 0.0 {
                format!(", every {}", fmt_num(period))
            } else {
                String::new()
            }
        ),
    }
}

fn drive_text(d: &level::Drive, electrode: bool) -> String {
    let mut parts = vec![format!("R = {}", fmt_num(d.resistance))];
    if electrode && d.inductance > 0.0 {
        parts.push(format!("L = {}", fmt_num(d.inductance)));
    }
    if !electrode && d.capacitance > 0.0 {
        parts.push(format!("C = {}", fmt_num(d.capacitance)));
    }
    if let Some(s) = d.switch {
        parts.push(format!(
            "a switch {} at first, toggling at t = {}",
            if s.closed { "closed" } else { "open" },
            fmt_num(s.at)
        ));
    }
    format!("in series: {}", parts.join(", "))
}

/// A non-constant source's signal over the run (a sine over a few periods at most).
fn chart(s: &Source, t_max: f64) -> Option<(String, Vec<(f64, f64)>)> {
    let end = match *s {
        Source::Dc { .. } => return None,
        Source::Sine { omega, .. } if omega > 0.0 => {
            (6.0 * std::f64::consts::TAU / omega).min(t_max.max(1e-9))
        }
        Source::Sine { .. } => return None,
        Source::Pulse { .. } => t_max.max(1e-9),
    };
    let w = s.waveform();
    Some((
        "the source's voltage over time".into(),
        sample(end, |t| w.value(t)),
    ))
}

#[allow(clippy::cast_precision_loss)]
fn sample(end: f64, f: impl Fn(f64) -> f64) -> Vec<(f64, f64)> {
    (0..=CHART_SAMPLES)
        .map(|k| {
            let t = end * k as f64 / CHART_SAMPLES as f64;
            (t, f(t))
        })
        .collect()
}

fn segment_distance(a: DVec3, b: DVec3, p: DVec3) -> f64 {
    let d = b - a;
    let l2 = d.length_squared();
    let s = if l2 > 0.0 {
        ((p - a).dot(d) / l2).clamp(0.0, 1.0)
    } else {
        0.0
    };
    (a + d * s - p).length()
}

/// Whether `p` lies in the box of `length` along `angle_deg` and `thickness` across it,
/// centred at `c`.
fn in_box(c: DVec3, length: f64, thickness: f64, angle_deg: f64, p: DVec3) -> bool {
    let (s, co) = angle_deg.to_radians().sin_cos();
    let r = p - c;
    let (u, v) = (r.x * co + r.y * s, -r.x * s + r.y * co);
    u.abs() <= 0.5 * length && v.abs() <= 0.5 * thickness
}
