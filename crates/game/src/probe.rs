//! The board probe (the owner, 2026-10-20): while I is held, everything known at the pointer,
//! anywhere on the board, from every source registered in `PROVIDERS`: the point and its time,
//! the fields, the particle under the pointer, a metal's carriers there, a flight passing it.
//! A new feature adds its provider to the list: a function from the game and the point to a
//! section of the tooltip, or nothing where it has nothing to say.

use level::tube::{KIND_NEGATIVE, KIND_POSITIVE, KIND_PROJECTILE};
use physics::DVec3;
use physics::field::FieldSolver;

use crate::Game;
use crate::ui::fmt_num;

/// One source's part of the probe's tooltip.
pub struct Section {
    pub title: String,
    pub lines: Vec<String>,
}

/// A source of the probe: its section at a board point.
type Provider = fn(&Game, DVec3) -> Option<Section>;

/// Every source, in the tooltip's order.
const PROVIDERS: &[Provider] = &[
    point,
    tube_fields,
    tube_particle,
    metal_here,
    slice_fields,
    flight_here,
];

/// Every section the registered sources have for the board point `at`.
#[must_use]
pub fn probe(game: &Game, at: DVec3) -> Vec<Section> {
    PROVIDERS.iter().filter_map(|p| p(game, at)).collect()
}

/// How close a particle or a flight's point must be to the pointer (cells).
const NEAR: f64 = 0.35;
/// The disc a metal's local carrier density is counted over (cells).
const DENSITY_RADIUS: f64 = 0.5;

fn distance2(a: DVec3, b: DVec3) -> f64 {
    let (dx, dy) = (a.x - b.x, a.y - b.y);
    dx * dx + dy * dy
}

/// An in-plane vector: its components, its size and its direction.
fn vector(v: DVec3) -> String {
    let size = (v.x * v.x + v.y * v.y).sqrt();
    if size == 0.0 {
        return "0".into();
    }
    format!(
        "({}, {}): {} toward {}°",
        fmt_num(v.x),
        fmt_num(v.y),
        fmt_num(size),
        fmt_num(libm::atan2(v.y, v.x).to_degrees())
    )
}

fn point(game: &Game, at: DVec3) -> Option<Section> {
    let mut lines = Vec::new();
    if game.editor.level.is_tube()
        && let Some(f) = game.tube_display_frame()
    {
        lines.push(format!("at t = {} (the frame shown)", fmt_num(f.t)));
    }
    Some(Section {
        title: format!("Point ({}, {})", fmt_num(at.x), fmt_num(at.y)),
        lines,
    })
}

/// Tube levels: the potential and the field of the frame shown (the electrodes' charge, every
/// particle's, the metals' backgrounds, a changing field's induced part) and the solenoid's B.
fn tube_fields(game: &Game, at: DVec3) -> Option<Section> {
    let view = game.tube_preview.as_ref()?;
    let frame = game.tube_display_frame()?;
    let (phi, e) = view.run.field_at(&frame, at);
    let mut lines = vec![
        format!("potential φ = {}", fmt_num(phi)),
        format!("E = {}", vector(e)),
    ];
    if let Some(w) = &view.run.b_wave {
        lines.push(format!(
            "of which induced by the changing B: {}",
            vector(w.induced(frame.t, at))
        ));
    }
    if view.run.has_field() {
        lines.push(format!(
            "B = {} along z (the solenoid's, uniform)",
            fmt_num(view.run.b_at(frame.t))
        ));
    }
    Some(Section {
        title: "Fields".into(),
        lines,
    })
}

/// Tube levels: the particle nearest the pointer in the frame shown: its kind and charge, its
/// mean velocity since the frame before, and the forces on it.
fn tube_particle(game: &Game, at: DVec3) -> Option<Section> {
    let view = game.tube_preview.as_ref()?;
    let frame = game.tube_display_frame()?;
    let position = |p: &[f32; 2]| DVec3::new(f64::from(p[0]), f64::from(p[1]), 0.0);
    let (k, d2) = frame
        .x
        .iter()
        .enumerate()
        .map(|(k, p)| (k, distance2(position(p), at)))
        .min_by(|a, b| a.1.total_cmp(&b.1))?;
    if d2 > NEAR * NEAR {
        return None;
    }
    let x = position(&frame.x[k]);
    let q = f64::from(frame.q[k]);
    let kind = match frame.kind[k] {
        KIND_NEGATIVE => "a negative carrier (electrons)",
        KIND_POSITIVE => "a positive carrier (ions, holes)",
        KIND_PROJECTILE => "a projectile",
        _ => "a particle",
    };
    let mut lines = vec![
        format!("{kind}, a line of charge {} per unit length", fmt_num(q)),
        format!("at ({}, {})", fmt_num(x.x), fmt_num(x.y)),
    ];
    // Its move since the frame before (a carrier's mean over that time: its collisions
    // within it are not seen).
    let frames = &view.run.frames;
    let before = frames.partition_point(|f| f.t < frame.t - 1e-9);
    let velocity = before
        .checked_sub(1)
        .and_then(|i| frames.get(i))
        .and_then(|f| {
            let j = f.id.iter().position(|id| *id == frame.id[k])?;
            let dt = frame.t - f.t;
            (dt > 0.0).then(|| (x - position(&f.x[j])) / dt)
        });
    let (_, e) = view.run.field_at(&frame, x);
    lines.push(format!("the field's push on it, qE = {}", vector(e * q)));
    if let Some(v) = velocity {
        lines.push(format!(
            "mean velocity since the frame before: {}",
            vector(v)
        ));
        if view.run.has_field() {
            let b = view.run.b_at(frame.t);
            lines.push(format!(
                "the magnetic push, q v × B = {}",
                vector(DVec3::new(v.y * b, -v.x * b, 0.0) * q)
            ));
        }
    }
    Some(Section {
        title: "The particle here".into(),
        lines,
    })
}

/// Tube levels: inside a metal, its carriers' charge density about the point (a disc of
/// radius 0.5, its part inside the metal) against the background's.
fn metal_here(game: &Game, at: DVec3) -> Option<Section> {
    let level = &game.editor.level;
    let spec = level.tube.as_ref()?;
    let (k, m) = spec.metals.iter().enumerate().find(|(_, m)| {
        let (lo, hi) = level.metal_rect(m);
        (lo.x..=hi.x).contains(&at.x) && (lo.y..=hi.y).contains(&at.y)
    })?;
    let frame = game.tube_display_frame()?;
    let (lo, hi) = level.metal_rect(m);
    let kind = if m.charge_per_mass < 0.0 {
        KIND_NEGATIVE
    } else {
        KIND_POSITIVE
    };
    let r2 = DENSITY_RADIUS * DENSITY_RADIUS;
    let charge: f64 = frame
        .x
        .iter()
        .zip(&frame.q)
        .zip(&frame.kind)
        .filter(|((p, _), kd)| {
            **kd == kind && distance2(DVec3::new(f64::from(p[0]), f64::from(p[1]), 0.0), at) <= r2
        })
        .map(|((_, q), _)| f64::from(*q).abs())
        .sum();
    // The disc's area inside the metal, by a grid of points over it.
    const N: i32 = 24;
    let mut inside = 0_u32;
    let mut all = 0_u32;
    for i in -N..=N {
        for j in -N..=N {
            let (u, w) = (
                f64::from(i) / f64::from(N) * DENSITY_RADIUS,
                f64::from(j) / f64::from(N) * DENSITY_RADIUS,
            );
            if u * u + w * w > r2 {
                continue;
            }
            all += 1;
            let (x, y) = (at.x + u, at.y + w);
            if (lo.x..=hi.x).contains(&x) && (lo.y..=hi.y).contains(&y) {
                inside += 1;
            }
        }
    }
    let area = std::f64::consts::PI * r2 * f64::from(inside) / f64::from(all.max(1));
    let density = charge / area.max(1e-12);
    Some(Section {
        title: format!("Metal {}", k + 1),
        lines: vec![
            format!(
                "its carriers' charge density here: {} (within {} of the point)",
                fmt_num(density),
                fmt_num(DENSITY_RADIUS)
            ),
            format!(
                "{} % of the background's {}: {}",
                fmt_num(100.0 * density / m.density),
                fmt_num(m.density),
                if density < 0.8 * m.density {
                    "depleted"
                } else if density > 1.2 * m.density {
                    "accumulated"
                } else {
                    "about neutral"
                }
            ),
        ],
    })
}

/// Particle levels: the static field at the point (at t = 0 for driven sources).
fn slice_fields(game: &Game, at: DVec3) -> Option<Section> {
    if game.editor.level.is_tube() {
        return None;
    }
    let (revision, field) = game.probe_field.as_ref()?;
    if *revision != game.sent_revision {
        return None;
    }
    let s = field.sample(at, 0.0);
    Some(Section {
        title: "Fields (at t = 0)".into(),
        lines: vec![
            format!("E = {}", vector(s.e)),
            format!("B = {} along z", fmt_num(s.b.z)),
            format!("potential φ = {} (the static sources')", fmt_num(s.phi)),
        ],
    })
}

/// Particle levels: a flight passing the point: when, how fast, its energies there.
fn flight_here(game: &Game, at: DVec3) -> Option<Section> {
    for (i, f) in game.flights.iter().enumerate() {
        let Some(preview) = &f.preview else {
            continue;
        };
        let Some(p) = preview
            .path
            .iter()
            .map(|p| (p, distance2(p.x, at)))
            .filter(|(_, d)| *d <= NEAR * NEAR)
            .min_by(|a, b| a.1.total_cmp(&b.1))
            .map(|(p, _)| p)
        else {
            continue;
        };
        let (shot, d) = game.editor.level.flight_of(i);
        let flight = if game.editor.level.disturbances.is_empty() {
            format!("Shot {}", shot + 1)
        } else {
            format!("Shot {} (disturbance {})", shot + 1, d + 1)
        };
        return Some(Section {
            title: format!("{flight} passes here"),
            lines: vec![
                format!(
                    "at t = {}, at ({}, {})",
                    fmt_num(p.t),
                    fmt_num(p.x.x),
                    fmt_num(p.x.y)
                ),
                format!("momentum {}", vector(p.p)),
                format!(
                    "kinetic energy {}, potential energy {} (from the launch point)",
                    fmt_num(p.kinetic),
                    fmt_num(p.potential)
                ),
            ],
        });
    }
    None
}
