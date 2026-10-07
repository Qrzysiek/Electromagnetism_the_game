//! Drawing of the scene with gizmos: grid, shots (launch, detector, trajectory), elements,
//! coils, regions, field lines, cursor.

use std::f32::consts::TAU;

use bevy::prelude::*;
use level::{Coil, Element, ElementKind, Grid};
use physics::DVec3;
use physics::trajectory::Outcome;
use physics::verify::Status;

use crate::{FieldLineGizmos, Game, PotentialQuad, ui};

#[allow(clippy::cast_possible_truncation)]
pub fn to_vec2(v: DVec3) -> Vec2 {
    Vec2::new(v.x as f32, v.y as f32)
}

/// Colour of a shot: the colour of its destination. Shots aimed at the same detector share
/// its hue (the palette cycles over the distinct detectors, in order of first use) and
/// differ in lightness only, so that a ray's colour tells where it must go.
pub fn shot_color(level: &level::Level, i: usize) -> Color {
    let Some(shot) = level.shots.get(i) else {
        return palette_color(i);
    };
    let base = detector_color(level, &shot.detector).to_srgba();
    let members: Vec<usize> = (0..level.shots.len())
        .filter(|&j| level.shots[j].detector == shot.detector)
        .collect();
    let n = members.len();
    if n < 2 {
        return base.into();
    }
    let rank = members.iter().position(|&j| j == i).unwrap_or(0);
    // Very slightly darker (−1) through the detector's colour (0) to slightly lighter (+1):
    // the rays of one detector must read as its colour.
    #[allow(clippy::cast_precision_loss)]
    let t = rank as f32 / (n - 1) as f32 * 2.0 - 1.0;
    let shade = |c: f32| {
        if t < 0.0 {
            c * (1.0 + 0.12 * t)
        } else {
            c + (1.0 - c) * 0.15 * t
        }
    };
    Color::srgb(shade(base.red), shade(base.green), shade(base.blue))
}

/// Colour of a detector: the palette's colour of its number among the level's distinct
/// detectors, in order of first use.
fn detector_color(level: &level::Level, detector: &level::Detector) -> Color {
    let mut targets: Vec<&level::Detector> = Vec::new();
    for s in &level.shots {
        if !targets.contains(&&s.detector) {
            targets.push(&s.detector);
        }
    }
    palette_color(targets.iter().position(|d| *d == detector).unwrap_or(0))
}

/// Colour of goal particle `k` (a level's free particle with a detector): after the
/// shots' detector colours.
pub fn free_goal_color(k: usize) -> Color {
    palette_color(4 + k)
}

/// The palette's colour number `i` (cycled). No yellow, orange, white or violet: those mark
/// unverified, lost and computing flights and gates; blue (the placement region, negative
/// charges) comes last.
fn palette_color(i: usize) -> Color {
    const PALETTE: [(f32, f32, f32); 6] = [
        (0.3, 1.0, 0.4),
        (1.0, 0.45, 0.8),
        (0.1, 0.9, 0.95),
        (0.75, 1.0, 0.2),
        (1.0, 0.35, 0.45),
        (0.45, 0.55, 1.0),
    ];
    let (r, g, b) = PALETTE[i % PALETTE.len()];
    Color::srgb(r, g, b)
}

fn draw_grid(gizmos: &mut Gizmos, grid: Grid) {
    let m = grid.max_node();
    let top_right = to_vec2(grid.position(m));
    let cell = Color::srgba(1.0, 1.0, 1.0, 0.035);
    for i in 0..=grid.nx {
        #[allow(clippy::cast_precision_loss)]
        let x = i as f32;
        gizmos.line_2d(Vec2::new(x, 0.0), Vec2::new(x, top_right.y), cell);
    }
    for j in 0..=grid.ny {
        #[allow(clippy::cast_precision_loss)]
        let y = j as f32;
        gizmos.line_2d(Vec2::new(0.0, y), Vec2::new(top_right.x, y), cell);
    }
    if grid.subdivision > 1 {
        let fine = Color::srgba(1.0, 1.0, 1.0, 0.015);
        for i in 0..=m[0] {
            let x = to_vec2(grid.position([i, 0, 0])).x;
            gizmos.line_2d(Vec2::new(x, 0.0), Vec2::new(x, top_right.y), fine);
        }
        for j in 0..=m[1] {
            let y = to_vec2(grid.position([0, j, 0])).y;
            gizmos.line_2d(Vec2::new(0.0, y), Vec2::new(top_right.x, y), fine);
        }
    }
}

/// A segment drawn as a line strip. Plain gizmo lines are rendered in a separate pass
/// before line strips (which include circles), so marks drawn with `line_2d` would end up
/// underneath the discs.
fn seg(gizmos: &mut Gizmos, a: Vec2, b: Vec2, color: Color) {
    gizmos.linestrip_2d([a, b], color);
}

/// Filled disc from concentric circles.
fn disc(gizmos: &mut Gizmos, p: Vec2, r: f32, color: Color) {
    for k in 1..=6 {
        #[allow(clippy::cast_precision_loss)]
        gizmos.circle_2d(p, r * k as f32 / 6.0, color);
    }
}

fn draw_element(
    gizmos: &mut Gizmos,
    grid: Grid,
    e: &Element,
    radius: f32,
    scale: f64,
    player: bool,
) {
    let p = to_vec2(grid.position(e.node));
    let ink = Color::srgb(0.05, 0.05, 0.05);
    let ring = if player {
        Color::srgb(1.0, 1.0, 1.0)
    } else {
        Color::srgb(0.5, 0.5, 0.5)
    };
    #[allow(clippy::cast_possible_truncation)]
    let halo = radius * (1.0 + 0.8 * (e.value.abs() / scale).min(1.0) as f32);
    match e.kind {
        ElementKind::Charge => {
            let color = if e.value > 0.0 {
                Color::srgb(1.0, 0.35, 0.3)
            } else {
                Color::srgb(0.35, 0.6, 1.0)
            };
            disc(gizmos, p, radius, color);
            gizmos.circle_2d(p, halo, color.with_alpha(0.4));
            let s = radius * 0.6;
            seg(gizmos, p - Vec2::X * s, p + Vec2::X * s, ink);
            if e.value > 0.0 {
                seg(gizmos, p - Vec2::Y * s, p + Vec2::Y * s, ink);
            }
        }
        ElementKind::Magnet => {
            // Moment out of the plane: ⊙; into the plane: ⊗.
            let color = Color::srgb(0.85, 0.55, 1.0);
            disc(gizmos, p, radius, color);
            gizmos.circle_2d(p, halo, color.with_alpha(0.4));
            if e.value > 0.0 {
                disc(gizmos, p, radius * 0.25, ink);
            } else {
                let s = radius * 0.55;
                seg(gizmos, p - Vec2::splat(s), p + Vec2::splat(s), ink);
                seg(gizmos, p + Vec2::new(-s, s), p + Vec2::new(s, -s), ink);
            }
        }
        ElementKind::Antenna => {
            // Arrow along the dipole amplitude p₀ (reversed for the opposite phase).
            let color = Color::srgb(0.45, 0.95, 0.75);
            disc(gizmos, p, radius, color);
            gizmos.circle_2d(p, halo, color.with_alpha(0.4));
            #[allow(clippy::cast_possible_truncation)]
            let dir = Vec2::from_angle(e.angle_deg.to_radians() as f32) * e.value.signum() as f32;
            let tip = p + dir * radius * 0.8;
            seg(gizmos, p - dir * radius * 0.8, tip, ink);
            let back = -dir * radius * 0.45;
            seg(gizmos, tip, tip + Vec2::from_angle(0.6).rotate(back), ink);
            seg(gizmos, tip, tip + Vec2::from_angle(-0.6).rotate(back), ink);
        }
        // Drawn with the electrodes.
        ElementKind::Plate | ElementKind::Supply => return,
        // Drawn with their velocity arrows (`draw_free`).
        ElementKind::Free => return,
    }
    gizmos.circle_2d(p, radius * 1.15, ring);
}

/// A dynamic particle: a hollow ring in its sign's colour (it moves; fixed charges are
/// filled discs), its sign, and its launch velocity as an arrow in the direction of
/// motion (`tip` at the particle: at rest). The player's have a handle behind them, joined
/// by a band: pull it back to launch the particle forward (a slingshot).
fn draw_free(gizmos: &mut Gizmos, p: Vec2, q: f64, radius: f32, tip: Vec2, handle: Option<Vec2>) {
    let color = if q > 0.0 {
        Color::srgb(1.0, 0.45, 0.4)
    } else if q < 0.0 {
        Color::srgb(0.45, 0.65, 1.0)
    } else {
        Color::srgb(0.8, 0.8, 0.8)
    };
    for k in 0..3 {
        #[allow(clippy::cast_precision_loss)]
        gizmos.circle_2d(p, radius * (1.0 - 0.08 * k as f32), color);
    }
    let s = radius * 0.5;
    if q != 0.0 {
        seg(gizmos, p - Vec2::X * s, p + Vec2::X * s, color);
    }
    if q > 0.0 {
        seg(gizmos, p - Vec2::Y * s, p + Vec2::Y * s, color);
    }
    let d = tip - p;
    let arrow = Color::srgb(1.0, 0.95, 0.6);
    if d.length() > radius {
        let start = p + d.normalize() * radius;
        gizmos.arrow_2d(start, tip, arrow);
    }
    if let Some(h) = handle {
        // The band from the particle back to the handle, and the handle.
        let band = h - p;
        if band.length() > radius {
            gizmos.line_2d(p + band.normalize() * radius, h, arrow.with_alpha(0.35));
        }
        gizmos.circle_2d(h, 0.18, arrow);
        gizmos.circle_2d(h, 0.1, arrow.with_alpha(0.6));
    }
}

/// In-plane outline of an electrode box, closed.
/// A radiation goal's receiver (far away): a band hugging the arena's boundary (`lo`–`hi`)
/// over its directions as seen from the arena's centre, closed at the ends, with a wave
/// mark on the axis.
fn draw_radiation_goal(
    gizmos: &mut Gizmos,
    lo: Vec2,
    hi: Vec2,
    goal: level::RadiationGoal,
    color: Color,
) {
    let c = (lo + hi) * 0.5;
    let half_size = (hi - lo) * 0.5;
    // Distance from the centre to the boundary along direction d.
    let edge = |d: Vec2| {
        let tx = if d.x.abs() > 1e-6 {
            half_size.x / d.x.abs()
        } else {
            f32::INFINITY
        };
        let ty = if d.y.abs() > 1e-6 {
            half_size.y / d.y.abs()
        } else {
            f32::INFINITY
        };
        tx.min(ty)
    };
    #[allow(clippy::cast_possible_truncation)]
    let (axis, half) = (
        goal.direction[0].to_radians() as f32,
        goal.direction[1].to_radians() as f32,
    );
    let n = 32;
    #[allow(clippy::cast_precision_loss)]
    let at = move |i: usize, gap: f32| {
        let d = Vec2::from_angle(axis - half + 2.0 * half * i as f32 / n as f32);
        c + d * (edge(d) + gap)
    };
    gizmos.linestrip_2d((0..=n).map(|i| at(i, 0.4)), color);
    gizmos.linestrip_2d((0..=n).map(|i| at(i, 1.2)), color);
    gizmos.line_2d(at(0, 0.4), at(0, 1.2), color);
    gizmos.line_2d(at(n, 0.4), at(n, 1.2), color);
    // A short wave travelling out along the axis, inside the band.
    let d = Vec2::from_angle(axis);
    let side = d.perp();
    let r0 = edge(d) + 0.45;
    #[allow(clippy::cast_precision_loss)]
    let wave = (0..=16).map(|i| {
        let u = i as f32 / 16.0;
        c + d * (r0 + 0.7 * u) + side * 0.2 * (u * std::f32::consts::TAU * 2.0).sin()
    });
    gizmos.linestrip_2d(wave, color);
}

fn box_outline(b: &physics::bem::BoxElectrode, grow: f32) -> [Vec2; 5] {
    let c = to_vec2(b.center);
    #[allow(clippy::cast_possible_truncation)]
    let (a, hl, ht) = (
        b.angle as f32,
        b.half_length as f32 + grow,
        b.half_thickness as f32 + grow,
    );
    let u = Vec2::from_angle(a);
    let v = u.perp();
    [
        c - u * hl - v * ht,
        c + u * hl - v * ht,
        c + u * hl + v * ht,
        c - u * hl + v * ht,
        c - u * hl - v * ht,
    ]
}

/// An electrode: filled metallic rectangle (its cross-section in the plane), rim
/// coloured by its potential (grey: grounded, red/blue: positive/negative).
fn draw_box(gizmos: &mut Gizmos, b: &physics::bem::BoxElectrode) {
    let c = to_vec2(b.center);
    #[allow(clippy::cast_possible_truncation)]
    let (a, hl, ht) = (
        b.angle as f32,
        b.half_length as f32,
        b.half_thickness as f32,
    );
    let u = Vec2::from_angle(a);
    let v = u.perp();
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    let strips = ((2.0 * ht / 0.03).ceil() as u32).max(2);
    for k in 0..=strips {
        #[allow(clippy::cast_precision_loss)]
        let t = -ht + 2.0 * ht * k as f32 / strips as f32;
        seg(
            gizmos,
            c - u * hl + v * t,
            c + u * hl + v * t,
            Color::srgb(0.42, 0.45, 0.5),
        );
    }
    let rim = match b.bias {
        physics::conductor::Bias::Charge(x) | physics::conductor::Bias::Potential(x) if x > 0.0 => {
            Color::srgb(1.0, 0.45, 0.35)
        }
        physics::conductor::Bias::Charge(x) | physics::conductor::Bias::Potential(x) if x < 0.0 => {
            Color::srgb(0.4, 0.65, 1.0)
        }
        _ => Color::srgb(0.85, 0.85, 0.85),
    };
    gizmos.linestrip_2d(box_outline(b, 0.0), rim);
}

/// A dielectric: its cross-section hatched in teal (an insulator, unlike the grey metal).
fn draw_dielectric(gizmos: &mut Gizmos, b: &physics::bem::BoxElectrode) {
    let c = to_vec2(b.center);
    #[allow(clippy::cast_possible_truncation)]
    let (a, hl, ht) = (
        b.angle as f32,
        b.half_length as f32,
        b.half_thickness as f32,
    );
    let u = Vec2::from_angle(a);
    let v = u.perp();
    let teal = Color::srgb(0.3, 0.75, 0.7);
    // Diagonal hatching every 0.3 cells across the box.
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    let lines = ((2.0 * (hl + ht) / 0.3).ceil() as u32).max(2);
    for k in 0..=lines {
        #[allow(clippy::cast_precision_loss)]
        let s = -(hl + ht) + 2.0 * (hl + ht) * k as f32 / lines as f32;
        // The line x − y = s (box coordinates) clipped to the box.
        let (x0, x1) = ((s - ht).max(-hl), (s + ht).min(hl));
        if x0 < x1 {
            seg(
                gizmos,
                c + u * x0 + v * (x0 - s),
                c + u * x1 + v * (x1 - s),
                teal.with_alpha(0.45),
            );
        }
    }
    gizmos.linestrip_2d(box_outline(b, 0.0), teal);
}

/// A ferrite: its cross-section filled dark slate with a violet-grey rim (a magnetic
/// insulator: neither metal grey nor dielectric teal).
fn draw_ferrite(gizmos: &mut Gizmos, b: &physics::bem::BoxElectrode) {
    let c = to_vec2(b.center);
    #[allow(clippy::cast_possible_truncation)]
    let (a, hl, ht) = (
        b.angle as f32,
        b.half_length as f32,
        b.half_thickness as f32,
    );
    let u = Vec2::from_angle(a);
    let v = u.perp();
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    let strips = ((2.0 * ht / 0.05).ceil() as u32).max(2);
    for k in 0..=strips {
        #[allow(clippy::cast_precision_loss)]
        let t = -ht + 2.0 * ht * k as f32 / strips as f32;
        seg(
            gizmos,
            c - u * hl + v * t,
            c + u * hl + v * t,
            Color::srgb(0.22, 0.2, 0.28),
        );
    }
    gizmos.linestrip_2d(box_outline(b, 0.0), Color::srgb(0.62, 0.55, 0.75));
}

/// A dielectric sphere: its cross-section hatched in teal, as a dielectric box.
fn draw_dielectric_sphere(gizmos: &mut Gizmos, c: Vec2, radius: f64) {
    #[allow(clippy::cast_possible_truncation)]
    let r = radius as f32;
    let teal = Color::srgb(0.3, 0.75, 0.7);
    // Diagonal hatching every 0.3 cells, clipped to the circle.
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    let lines = ((2.0 * std::f32::consts::SQRT_2 * r / 0.3).ceil() as u32).max(2);
    let (u, v) = (
        Vec2::new(1.0, 1.0).normalize(),
        Vec2::new(1.0, -1.0).normalize(),
    );
    for k in 1..lines {
        #[allow(clippy::cast_precision_loss)]
        let s = -r + 2.0 * r * k as f32 / lines as f32;
        let half = (r * r - s * s).max(0.0).sqrt();
        seg(
            gizmos,
            c + v * s - u * half,
            c + v * s + u * half,
            teal.with_alpha(0.45),
        );
    }
    gizmos.circle_2d(c, r, teal).resolution(96);
}

fn draw_coil(gizmos: &mut Gizmos, grid: Grid, coil: &Coil) {
    let copper = Color::srgb(0.95, 0.6, 0.3);
    let arrow = |gizmos: &mut Gizmos, at: Vec2, dir: Vec2| {
        let back = -dir * 0.45;
        gizmos.line_2d(at, at + Vec2::from_angle(0.5).rotate(back), copper);
        gizmos.line_2d(at, at + Vec2::from_angle(-0.5).rotate(back), copper);
    };
    match coil {
        Coil::Circle {
            center,
            radius,
            kappa,
            rate,
            drive,
            tunable,
        } => {
            let c = to_vec2(grid.position(*center));
            #[allow(clippy::cast_possible_truncation)]
            let r = *radius as f32;
            gizmos.circle_2d(c, r, copper).resolution(128);
            if *tunable {
                // Its power supply: a yellow ring inside (click the coil to switch it).
                gizmos
                    .circle_2d(c, (r - 0.3).max(0.1), Color::srgb(1.0, 0.85, 0.3))
                    .resolution(128);
            }
            // A ramped or driven coil: a second, dashed ring (its current changes).
            if *rate != 0.0 || drive.is_some() {
                for k in (0..64).step_by(2) {
                    #[allow(clippy::cast_precision_loss)]
                    let (a0, a1) = (TAU * k as f32 / 64.0, TAU * (k + 1) as f32 / 64.0);
                    gizmos.line_2d(
                        c + Vec2::from_angle(a0) * (r + 0.25),
                        c + Vec2::from_angle(a1) * (r + 0.25),
                        copper.with_alpha(0.7),
                    );
                }
            }
            // Current direction: counter-clockwise for positive κ.
            let sense = if *kappa >= 0.0 { 1.0 } else { -1.0 };
            for k in 0..8 {
                #[allow(clippy::cast_precision_loss)]
                let a = TAU * k as f32 / 8.0;
                let at = c + Vec2::from_angle(a) * r;
                arrow(gizmos, at, Vec2::from_angle(a + sense * TAU / 4.0));
            }
        }
        Coil::Polygon {
            vertices,
            kappa,
            rate,
            drive,
        } => {
            let v: Vec<Vec2> = vertices
                .iter()
                .map(|n| to_vec2(grid.position(*n)))
                .collect();
            // A ramped or driven coil: its sides doubled by a dashed line (its current
            // changes), as for a circle.
            let changing = *rate != 0.0 || drive.is_some();
            for i in 0..v.len() {
                let (a, b) = (v[i], v[(i + 1) % v.len()]);
                gizmos.line_2d(a, b, copper);
                if changing {
                    let n = (b - a).perp().normalize_or_zero() * 0.25;
                    for k in (0..16).step_by(2) {
                        #[allow(clippy::cast_precision_loss)]
                        let (s0, s1) = (k as f32 / 16.0, (k + 1) as f32 / 16.0);
                        gizmos.line_2d(
                            a + (b - a) * s0 + n,
                            a + (b - a) * s1 + n,
                            copper.with_alpha(0.7),
                        );
                    }
                }
                let dir = (b - a).normalize_or_zero() * if *kappa >= 0.0 { 1.0 } else { -1.0 };
                arrow(gizmos, (a + b) * 0.5, dir);
            }
        }
    }
}

pub fn draw(
    game: Res<Game>,
    mut gizmos: Gizmos,
    mut line_gizmos: Gizmos<FieldLineGizmos>,
    quad: Res<PotentialQuad>,
    radiation: Res<crate::radiation::RadiationView>,
    mut vis: Query<&mut Visibility>,
) {
    crate::radiation::draw_arrows(&mut gizmos, &radiation);
    crate::radiation::draw_tracers(&mut gizmos, &radiation);
    if let Ok(mut v) = vis.get_mut(quad.entity) {
        *v = if matches!(
            game.map,
            Some(crate::potential::MapMode::Potential | crate::potential::MapMode::Magnetic)
        ) {
            Visibility::Visible
        } else {
            Visibility::Hidden
        };
    }
    let level = &game.editor.level;
    let grid = level.grid;
    draw_grid(&mut gizmos, grid);

    // World bounds.
    let b = level.bounds();
    let (bmin, bmax) = (to_vec2(b.min), to_vec2(b.max));
    gizmos.rect_2d(
        (bmin + bmax) * 0.5,
        bmax - bmin,
        Color::srgba(1.0, 0.3, 0.3, 0.5),
    );

    // Field lines.
    if game.show_field_lines {
        let color = Color::srgba(0.95, 0.95, 0.75, game.field_line_opacity);
        let (head, spread) = (0.3, 0.45_f32);
        for (points, arrows) in &game.field_lines {
            line_gizmos.linestrip_2d(points.iter().copied(), color);
            // Arrowheads pointing along E.
            for &(p, d) in arrows {
                let back = -d * head;
                let tip = p + d * (0.5 * head);
                line_gizmos.line_2d(tip, tip + Vec2::from_angle(spread).rotate(back), color);
                line_gizmos.line_2d(tip, tip + Vec2::from_angle(-spread).rotate(back), color);
            }
        }
    }

    // Coils.
    for coil in &level.coils {
        draw_coil(&mut gizmos, grid, coil);
    }

    // Dielectrics: translucent boxes with a teal rim (insulators, not metal).
    for b in level.dielectric_boxes() {
        draw_dielectric(&mut gizmos, &b);
    }
    for b in level.ferrite_boxes() {
        draw_ferrite(&mut gizmos, &b);
    }
    for d in &level.dielectric_spheres {
        draw_dielectric_sphere(&mut gizmos, to_vec2(grid.position(d.center)), d.radius);
    }

    // Electrodes: the level's (at the potentials of the player's power supplies), then
    // the player's plates. Tunable electrodes get a yellow frame (click them to switch
    // their supply), player plates a white one.
    let boxes = level.all_box_electrodes(&game.editor.placement);
    for (i, b) in boxes.iter().enumerate() {
        draw_box(&mut gizmos, b);
        let frame = match level.electrodes.get(i) {
            Some(e) if e.tunable => Some(Color::srgb(1.0, 0.85, 0.3)),
            Some(_) => None,
            None => Some(Color::srgb(1.0, 1.0, 1.0)),
        };
        if let Some(color) = frame {
            gizmos.linestrip_2d(box_outline(b, 0.12), color);
        }
    }

    // Charge clouds: penetrable, so no solid fill: faint rings in the sign's colour
    // (red positive, blue negative), and a dashed rim.
    for c in &level.clouds {
        let p = to_vec2(grid.position(c.center));
        #[allow(clippy::cast_possible_truncation)]
        let r = c.radius as f32;
        let (red, green, blue) = if c.charge >= 0.0 {
            (1.0, 0.45, 0.35)
        } else {
            (0.4, 0.65, 1.0)
        };
        for k in 1..=8 {
            #[allow(clippy::cast_precision_loss)]
            let rk = r * k as f32 / 8.0;
            gizmos
                .circle_2d(p, rk, Color::srgba(red, green, blue, 0.10))
                .resolution(64);
        }
        for k in 0..48 {
            if k % 2 == 1 {
                continue;
            }
            #[allow(clippy::cast_precision_loss)]
            let (a0, a1) = (TAU * k as f32 / 48.0, TAU * (k + 1) as f32 / 48.0);
            gizmos.line_2d(
                p + Vec2::from_angle(a0) * r,
                p + Vec2::from_angle(a1) * r,
                Color::srgba(red, green, blue, 0.7),
            );
        }
    }

    // Metal spheres: a filled metallic disc, rim coloured by how it is held (grey:
    // grounded, red/blue: positive/negative charge or potential).
    for c in &level.conductors {
        let p = to_vec2(grid.position(c.center));
        #[allow(clippy::cast_possible_truncation)]
        let r = c.radius as f32;
        #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
        let rings = ((r / 0.04).ceil() as u32).max(6);
        for k in 1..=rings {
            #[allow(clippy::cast_precision_loss)]
            let rk = r * k as f32 / rings as f32;
            gizmos
                .circle_2d(p, rk, Color::srgb(0.42, 0.45, 0.5))
                .resolution(64);
        }
        let rim = match c.bias {
            level::ConductorBias::Grounded => Color::srgb(0.85, 0.85, 0.85),
            level::ConductorBias::Charge(v) | level::ConductorBias::Potential(v) if v > 0.0 => {
                Color::srgb(1.0, 0.45, 0.35)
            }
            level::ConductorBias::Charge(v) | level::ConductorBias::Potential(v) if v < 0.0 => {
                Color::srgb(0.4, 0.65, 1.0)
            }
            _ => Color::srgb(0.85, 0.85, 0.85),
        };
        gizmos.circle_2d(p, r, rim).resolution(96);
        if matches!(c.bias, level::ConductorBias::Grounded) {
            // Ground symbol.
            let g = p - Vec2::new(0.0, r * 0.35);
            for (k, w) in [0.35f32, 0.22, 0.1].iter().enumerate() {
                #[allow(clippy::cast_precision_loss)]
                let y = g.y - k as f32 * 0.12;
                seg(
                    &mut gizmos,
                    Vec2::new(g.x - w, y),
                    Vec2::new(g.x + w, y),
                    rim,
                );
            }
            seg(&mut gizmos, g, g + Vec2::new(0.0, 0.3), rim);
        }
    }

    // Gates: dashed violet boxes, numbered by dots at their top edge (1, 2, …), with the
    // accepted directions as a cone like detectors.
    let violet = Color::srgb(0.8, 0.55, 1.0);
    for (k, gate) in level.gates.iter().enumerate() {
        let (a, b) = (
            to_vec2(grid.position(gate.min)),
            to_vec2(grid.position(gate.max)),
        );
        let (lo, hi) = (a.min(b), a.max(b));
        let corners = [lo, Vec2::new(hi.x, lo.y), hi, Vec2::new(lo.x, hi.y), lo];
        for w in corners.windows(2) {
            let n = ((w[1] - w[0]).length() / 0.3).ceil().max(1.0);
            #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
            for i in (0..n as u32).step_by(2) {
                #[allow(clippy::cast_precision_loss)]
                let (s, e) = (i as f32 / n, ((i + 1) as f32 / n).min(1.0));
                seg(&mut gizmos, w[0].lerp(w[1], s), w[0].lerp(w[1], e), violet);
            }
        }
        for i in 0..=k {
            #[allow(clippy::cast_precision_loss)]
            let p = Vec2::new(lo.x + 0.25 + 0.25 * i as f32, hi.y + 0.25);
            disc(&mut gizmos, p, 0.07, violet);
        }
        if let Some([axis, half]) = gate.acceptance.and_then(|a| a.direction) {
            let c = (lo + hi) * 0.5;
            let len = (hi - lo).min_element().max(1.0) * 0.9;
            #[allow(clippy::cast_possible_truncation)]
            let (axis, half) = (axis.to_radians() as f32, half.to_radians() as f32);
            let dir = Vec2::from_angle(axis);
            let tail = c - dir * len * 0.5;
            gizmos.arrow_2d(tail, c + dir * len * 0.5, violet);
            for s in [-half, half] {
                gizmos.line_2d(tail, tail + Vec2::from_angle(axis + s) * len, violet);
            }
        }
    }

    // Shots: detectors and launch points. The active shot is drawn brighter.
    for (i, shot) in level.shots.iter().enumerate() {
        let active = i == game.active_shot;
        if !(active || game.show_all_shots) {
            continue;
        }
        let color = shot_color(level, i);
        // The detector in its own colour, the launch marker in the shot's shade of it.
        let target = detector_color(level, &shot.detector);
        let alpha = if active { 1.0 } else { 0.45 };
        let solved = matches!(game.verdict(i), Some((Status::Verified, Outcome::Arrived)));
        let d0 = to_vec2(grid.position(shot.detector.min));
        let d1 = to_vec2(grid.position(shot.detector.max));
        gizmos.rect_2d((d0 + d1) * 0.5, (d1 - d0).abs(), target.with_alpha(alpha));
        // Accepted directions: a cone (axis and half-angle) at the detector's centre.
        if let Some([axis, half]) = shot.detector.acceptance.and_then(|a| a.direction) {
            let c = (d0 + d1) * 0.5;
            let len = (d1 - d0).abs().min_element().max(1.0) * 0.9;
            #[allow(clippy::cast_possible_truncation)]
            let (axis, half) = (axis.to_radians() as f32, half.to_radians() as f32);
            let dir = Vec2::from_angle(axis);
            let tail = c - dir * len * 0.5;
            let cone = target.with_alpha(alpha * 0.8);
            gizmos.arrow_2d(tail, c + dir * len * 0.5, cone);
            for s in [-half, half] {
                gizmos.line_2d(tail, tail + Vec2::from_angle(axis + s) * len, cone);
            }
        }
        if solved || active {
            gizmos.rect_2d(
                (d0 + d1) * 0.5,
                (d1 - d0).abs() - Vec2::splat(0.15),
                target.with_alpha(alpha * 0.5),
            );
        }
        // Radiation goal: a far receiver covering an arc of directions, drawn as a band
        // just outside the arena, seen from its centre.
        if let Some(r) = shot.detector.acceptance.and_then(|a| a.radiation) {
            draw_radiation_goal(&mut gizmos, bmin, bmax, r, target.with_alpha(alpha));
        }
        let a = to_vec2(grid.position(shot.launch.node));
        let dir = shot.launch.direction;
        #[allow(clippy::cast_possible_truncation)]
        let dir = Vec2::new(dir[0] as f32, dir[1] as f32).normalize_or_zero();
        gizmos.circle_2d(a, 0.25, color.with_alpha(alpha));
        gizmos.arrow_2d(a, a + dir * 1.5, color.with_alpha(alpha));
    }

    // Elements.
    #[allow(clippy::cast_possible_truncation)]
    let charge_radius = level.physics.charge_radius as f32;
    #[allow(clippy::cast_possible_truncation)]
    let magnet_radius = level.physics.magnet_radius as f32;
    #[allow(clippy::cast_possible_truncation)]
    let antenna_radius = level.physics.antenna_radius as f32;
    let scale = |kind: ElementKind| {
        let list = crate::editor::magnitudes(level, kind);
        let from_level = level
            .elements
            .iter()
            .filter(|e| e.kind == kind)
            .map(|e| e.value.abs());
        list.iter()
            .copied()
            .chain(from_level)
            .fold(1e-300, f64::max)
    };
    let (q_scale, m_scale, a_scale) = (
        scale(ElementKind::Charge),
        scale(ElementKind::Magnet),
        scale(ElementKind::Antenna),
    );
    let draw_all = |gizmos: &mut Gizmos, list: &[Element], player: bool| {
        for e in list {
            let (r, s) = match e.kind {
                ElementKind::Charge => (charge_radius, q_scale),
                ElementKind::Magnet => (magnet_radius, m_scale),
                ElementKind::Antenna => (antenna_radius, a_scale),
                ElementKind::Plate | ElementKind::Supply => continue,
                ElementKind::Free => {
                    let t = crate::editor::arrow_tip(level, e, game.arrow_measure);
                    let h = crate::editor::arrow_handle(level, e, game.arrow_measure);
                    let p = to_vec2(grid.position(e.node));
                    #[allow(clippy::cast_possible_truncation)]
                    let (tip, handle, r) = (
                        if e.speed.unwrap_or(0.0) == 0.0 {
                            p
                        } else {
                            Vec2::new(t[0] as f32, t[1] as f32)
                        },
                        Vec2::new(h[0] as f32, h[1] as f32),
                        level.limits.free_radius.max(0.15) as f32,
                    );
                    draw_free(gizmos, p, e.value, r, tip, player.then_some(handle));
                    continue;
                }
            };
            draw_element(gizmos, grid, e, r, s, player);
        }
    };
    draw_all(&mut gizmos, &level.elements, false);
    draw_all(&mut gizmos, &game.editor.placement, true);
    // The level's dynamic particles, with their detectors (goals) if they have one.
    let scale_arrow = crate::editor::arrow_scale(level, game.arrow_measure);
    for (k, f) in level.free_particles.iter().enumerate() {
        let p = to_vec2(grid.position(f.node));
        let speed = f.velocity[0].hypot(f.velocity[1]);
        let len = crate::editor::arrow_length(level, speed, game.arrow_measure) * scale_arrow;
        #[allow(clippy::cast_possible_truncation)]
        let v = if speed > 0.0 {
            Vec2::new(f.velocity[0] as f32, f.velocity[1] as f32) * (len / speed) as f32
        } else {
            Vec2::ZERO
        };
        #[allow(clippy::cast_possible_truncation)]
        let r = f.particle.radius.max(0.15) as f32;
        draw_free(&mut gizmos, p, f.particle.charge, r, p + v, None);
        if let Some(d) = f.detector {
            let d0 = to_vec2(grid.position(d.min));
            let d1 = to_vec2(grid.position(d.max));
            let color = free_goal_color(k);
            gizmos.rect_2d((d0 + d1) * 0.5, (d1 - d0).abs(), color);
            gizmos.rect_2d(
                (d0 + d1) * 0.5,
                (d1 - d0).abs() - Vec2::splat(0.15),
                color.with_alpha(0.5),
            );
            gizmos.line_2d(p, (d0 + d1) * 0.5, color.with_alpha(0.15));
        }
    }

    // Region where player elements may be placed.
    if let Some(r) = level.limits.region {
        let p0 = to_vec2(grid.position(r.min)) - Vec2::splat(0.35);
        let p1 = to_vec2(grid.position(r.max)) + Vec2::splat(0.35);
        let blue = Color::srgba(0.45, 0.7, 1.0, 0.55);
        gizmos.rect_2d((p0 + p1) * 0.5, p1 - p0, blue);
        gizmos.rect_2d(
            (p0 + p1) * 0.5,
            p1 - p0 - Vec2::splat(0.12),
            blue.with_alpha(0.25),
        );
    }

    // Box or circle being dragged in the sandbox.
    if let Some(a) = game.sandbox.detector_drag {
        let p0 = to_vec2(grid.position(a));
        let p1 = to_vec2(grid.position(game.editor.cursor));
        let color = Color::srgba(0.3, 1.0, 0.4, 0.6);
        if game.sandbox.tool == crate::sandbox::Tool::CoilCircle {
            gizmos.circle_2d(p0, p0.distance(p1), color).resolution(96);
        } else {
            gizmos.rect_2d((p0 + p1) * 0.5, (p1 - p0).abs(), color);
        }
    }

    // Polygon coil being drawn in the sandbox.
    let pending = &game.sandbox.pending_polygon;
    if !pending.is_empty() {
        let color = Color::srgba(0.95, 0.6, 0.3, 0.8);
        let mut pts: Vec<Vec2> = pending.iter().map(|n| to_vec2(grid.position(*n))).collect();
        pts.push(to_vec2(grid.position(game.editor.cursor)));
        gizmos.linestrip_2d(pts, color);
    }

    // Cursor.
    let cur = to_vec2(grid.position(game.editor.cursor));
    let cursor_color = match (game.editor.kind, game.editor.selected_value() > 0.0) {
        (ElementKind::Magnet, _) => Color::srgb(0.9, 0.7, 1.0),
        (ElementKind::Antenna, _) => Color::srgb(0.6, 1.0, 0.85),
        (ElementKind::Charge, true) => Color::srgb(1.0, 0.6, 0.5),
        (ElementKind::Charge, false) => Color::srgb(0.6, 0.8, 1.0),
        (ElementKind::Plate | ElementKind::Supply, _) => Color::srgb(0.9, 0.9, 0.7),
        (ElementKind::Free, _) => Color::srgb(1.0, 0.95, 0.6),
    };
    gizmos.rect_2d(cur, Vec2::splat(charge_radius * 2.8), cursor_color);
    // Outline of a new plate at the cursor (not over one of the player's elements or a
    // tunable electrode, where a click acts on those).
    if game.editor.kind == ElementKind::Plate
        && game.editor.left(ElementKind::Plate) > 0
        && (!game.sandbox.active || game.sandbox.tool == crate::sandbox::Tool::PlayerElement)
        && game.editor.element_at_cursor().is_none()
        && game.editor.tunable_at(game.editor.cursor).is_none()
    {
        let ghost = level::Element::plate(
            game.editor.cursor,
            game.editor.selected_value(),
            game.editor.plate_angle_deg,
        );
        let b = level.plate_box(&ghost);
        let ok = level
            .check_placement(&[game.editor.placement.as_slice(), &[ghost]].concat())
            .is_ok();
        let color = if ok {
            Color::srgba(0.9, 0.9, 0.7, 0.7)
        } else {
            Color::srgba(1.0, 0.4, 0.3, 0.7)
        };
        gizmos.linestrip_2d(box_outline(&b, 0.0), color);
    }
    // Element being moved.
    if let Some((i, _)) = game.editor.grabbed
        && let Some(e) = game.editor.placement.get(i)
    {
        let p = to_vec2(grid.position(e.node));
        gizmos.circle_2d(p, charge_radius * 2.2, Color::srgb(1.0, 0.95, 0.3));
    }

    // Trajectories: every flight of the shown shots. Under disturbances a shot has one
    // flight per disturbance; the one selected in the panel is drawn brightest.
    let active_flight = game.active_flight();
    // The particle-field views show one particle's field: draw only that flight.
    let single = matches!(
        game.map,
        Some(crate::potential::MapMode::ParticleField | crate::potential::MapMode::Total)
    );
    for (f, view) in game.flights.iter().enumerate() {
        if single && f != active_flight {
            continue;
        }
        let (i, _) = level.flight_of(f);
        let active = i == game.active_shot;
        if !(active || game.show_all_shots) {
            continue;
        }
        let Some(p) = &view.preview else {
            continue;
        };
        let primary =
            f == game.flight_index(i, game.active_disturbance.min(level.flights_per_shot() - 1));
        let base = shot_color(level, i);
        let color = match view.verdict {
            None => Color::srgb(0.95, 0.95, 0.95),
            Some((Status::Verified, Outcome::Arrived)) => base,
            Some((Status::Verified, _)) => Color::srgb(1.0, 0.45, 0.2),
            Some(_) => Color::srgb(1.0, 0.9, 0.2),
        }
        .with_alpha(
            match (active, primary) {
                (true, true) => 1.0,
                (false, true) => 0.55,
                (true, false) => 0.4,
                (false, false) => 0.2,
            } * stale_fade(view.preview_revision, game.sent_revision),
        );
        let active = f == active_flight;
        gizmos.linestrip_2d(p.path.iter().map(|q| to_vec2(q.x)), color);
        let end = to_vec2(p.path.last().expect("path has points").x);
        if !matches!(p.outcome, Outcome::Arrived) {
            let s = 0.3;
            gizmos.line_2d(end - Vec2::splat(s), end + Vec2::splat(s), color);
            gizmos.line_2d(end + Vec2::new(-s, s), end + Vec2::new(s, -s), color);
        }
        if game.animate
            && let Some(pt) = ui::point_at(p, game.anim_time)
        {
            let x = to_vec2(pt.x);
            let dot = Color::srgb(1.0, 1.0, 0.6);
            gizmos.circle_2d(x, 0.18, dot);
            gizmos.circle_2d(x, 0.1, dot);
            if active {
                // The particle whose energy the panel shows: ringed in its shot's colour.
                gizmos.circle_2d(x, 0.32, base);
            }
            let f = to_vec2(pt.force);
            if active && f.length() > 0.0 {
                // Force arrow, length ∝ log(1 + |F|/F₀), direction exact. F₀ = the energy
                // unit per cell (T₀, or for a particle launched at rest the kinetic energy
                // it reaches): the force that changes the kinetic energy by it over a cell.
                #[allow(clippy::cast_possible_truncation)]
                let f0 = game.energy_unit(i).0 as f32;
                let len = (1.0 + f.length() / f0).ln() * 1.5;
                gizmos.arrow_2d(x, x + f.normalize() * len, Color::srgb(1.0, 0.8, 0.2));
            }
        }
    }

    draw_beams(&mut gizmos, &game);
}

/// Position along a sampled path at time `t` (linear between samples); `None` outside it.
fn path_at(path: &[(f64, physics::DVec3)], t: f64) -> Option<physics::DVec3> {
    let last = path.last()?;
    if t < path[0].0 || t > last.0 {
        return None;
    }
    // A flight that ended where it started (a particle launched touching metal) has one
    // sample: clamp(1, 0) below would panic (the owner's crash with a free charge next
    // to a grounded sphere, found by scripts/pairs.py --gui).
    if path.len() == 1 {
        return Some(last.1);
    }
    let k = path
        .partition_point(|(tk, _)| *tk <= t)
        .clamp(1, path.len() - 1);
    let ((t0, x0), (t1, x1)) = (path[k - 1], path[k]);
    let f = if t1 > t0 { (t - t0) / (t1 - t0) } else { 0.0 };
    Some(x0 + (x1 - x0) * f)
}

/// Beam levels: every particle's path, coloured by its shot like single flights (white
/// while computing, orange when lost, yellow when not verified); the selected
/// disturbance's flight is drawn brightest, lost particles fainter. The animation moves
/// a dot along every path.
fn draw_beams(gizmos: &mut Gizmos, game: &Game) {
    if !game.editor.level.has_beams() {
        return;
    }
    let d_active = game
        .active_disturbance
        .min(game.beams.len().saturating_sub(1));
    for (d, view) in game.beams.iter().enumerate() {
        // The verdict's own flight once it has arrived (exact at finite c).
        let Some(p) = view.shown() else {
            continue;
        };
        let primary = d == d_active;
        for (i, path) in p.paths.iter().enumerate() {
            let shot = p.shots[i];
            let n_shots = game.editor.level.shots.len();
            // Dynamic particles (after the shots) are always shown.
            if shot < n_shots && !(shot == game.active_shot || game.show_all_shots) {
                continue;
            }
            let verdict = view.verified.as_ref().and_then(|r| r.get(i).copied());
            let arrived = p.outcomes[i] == Outcome::Arrived;
            // The player's free charges (no goal): pale yellow; the level's free
            // particles without a detector: grey.
            let scenery = shot == usize::MAX
                || (shot >= n_shots
                    && game.editor.level.free_particles[shot - n_shots]
                        .detector
                        .is_none());
            let color = match verdict {
                _ if shot == usize::MAX => Color::srgb(1.0, 0.95, 0.6),
                _ if scenery => Color::srgb(0.8, 0.8, 0.8),
                None => Color::srgb(0.95, 0.95, 0.95),
                Some((Status::Verified, Outcome::Arrived)) if shot >= n_shots => {
                    free_goal_color(shot - n_shots)
                }
                Some((Status::Verified, Outcome::Arrived)) => shot_color(&game.editor.level, shot),
                Some((Status::Verified, _)) => Color::srgb(1.0, 0.45, 0.2),
                Some(_) => Color::srgb(1.0, 0.9, 0.2),
            }
            .with_alpha(
                match (primary, arrived) {
                    (true, true) => 0.75,
                    (true, false) => 0.4,
                    (false, _) => 0.12,
                } * stale_fade(view.preview_revision, game.sent_revision),
            );
            gizmos.linestrip_2d(path.iter().map(|(_, x)| to_vec2(*x)), color);
            if game.animate && primary {
                if let Some(x) = path_at(path, game.anim_time) {
                    gizmos.circle_2d(to_vec2(x), 0.08, Color::srgb(1.0, 1.0, 0.6));
                } else if matches!(p.outcomes[i], Outcome::Collided(_))
                    && let Some(&(t_end, x)) = path.last()
                    && game.anim_time > t_end
                {
                    // Absorbed by a body: its charge stays there at rest.
                    gizmos.circle_2d(to_vec2(x), 0.08, Color::srgb(0.75, 0.75, 0.6));
                }
            }
        }
    }
}

/// A faint spring from `anchor` to `x`: ten coils of fixed width that stretch with the
/// distance. Kept for a future level with a mechanical spring (a particle on a real
/// spring). It was drawn for charge clouds, but a cloud acts as a spring only on a bound
/// particle inside it; one with more energy than any bound state is on no spring at all
/// (the owner's review), so atoms are drawn without it.
#[allow(dead_code)]
fn draw_spring(gizmos: &mut Gizmos, anchor: Vec2, x: Vec2, color: Color) {
    let d = x - anchor;
    let len = d.length();
    if len < 0.15 {
        return;
    }
    let (u, side) = (d / len, d.perp() / len);
    let n = 40;
    #[allow(clippy::cast_precision_loss)]
    let points = (0..=n).map(|k| {
        let s = k as f32 / n as f32;
        let wiggle = match k % 4 {
            _ if k == 0 || k == n => 0.0,
            1 => 0.12,
            3 => -0.12,
            _ => 0.0,
        };
        anchor + u * (len * s) + side * wiggle
    });
    gizmos.linestrip_2d(points, color);
    gizmos.circle_2d(anchor, 0.06, color);
}

/// Paths of a previous setup (while the new one is computed) are drawn faded.
fn stale_fade(preview_revision: u64, current: u64) -> f32 {
    if preview_revision == current {
        1.0
    } else {
        0.3
    }
}
