//! The board's inspection tooltip: what lies at a point, a metal map's values.

use super::*;

/// What obstacle `i` of a level's scenario is (see `Level::field_at` for the order).
pub(super) fn obstacle_name(game: &Game, i: usize) -> &'static str {
    let all = game
        .editor
        .level
        .elements
        .iter()
        .chain(&game.editor.placement);
    let charges = all
        .clone()
        .filter(|e| e.kind == ElementKind::Charge)
        .count();
    let magnets = all
        .clone()
        .filter(|e| e.kind == ElementKind::Magnet)
        .count();
    let antennas = all.filter(|e| e.kind == ElementKind::Antenna).count();
    let wires: usize = game
        .editor
        .level
        .coils
        .iter()
        .map(|c| match c {
            level::Coil::Circle { .. } => 1,
            level::Coil::Polygon { vertices, .. } => vertices.len(),
        })
        .sum();
    let electrodes = game.editor.level.electrodes.len();
    let plates = game
        .editor
        .placement
        .iter()
        .filter(|e| e.kind == ElementKind::Plate)
        .count();
    // Obstacle order of `Level::field`: charges, magnets, antennas, coil wires, level
    // electrodes, player plates, metal spheres.
    if i < charges {
        "a charge"
    } else if i < charges + magnets {
        "a magnet"
    } else if i < charges + magnets + antennas {
        "an antenna"
    } else if i < charges + magnets + antennas + wires {
        "a coil wire"
    } else if i < charges + magnets + antennas + wires + electrodes {
        "an electrode"
    } else if i < charges + magnets + antennas + wires + electrodes + plates {
        "one of your plates"
    } else {
        "a metal sphere"
    }
}

/// The tooltip over the object under the pointer: its kind, every parameter that defines
/// it, and a time-dependent source's signal (`inspect.rs`).
pub(super) fn inspect_tooltip(
    ctx: &egui::Context,
    game: &Game,
    (camera, transform): (&Camera, &GlobalTransform),
) {
    let Some(world) = game.hover_world else {
        return;
    };
    // The pointer, or (a capture's EM_HOVER) the point's place on the screen.
    let Some(pointer) = ctx.pointer_hover_pos().or_else(|| {
        let v = camera
            .world_to_viewport(transform, world.extend(0.0))
            .ok()?;
        Some(egui::pos2(v.x, v.y))
    }) else {
        return;
    };
    let level = &game.editor.level;
    let at = physics::DVec3::new(f64::from(world.x), f64::from(world.y), 0.0);
    let info = crate::inspect::inspect_at(level, &game.editor.placement, at);
    // The probe (I held): every registered source's values at the point, anywhere.
    let probed = if game.probe {
        crate::probe::probe(game, at)
    } else {
        Vec::new()
    };
    if info.is_none() && probed.is_empty() {
        return;
    }
    let live = metal_map_lines(game, at);
    egui::Area::new(egui::Id::new("inspect_tooltip"))
        .order(egui::Order::Tooltip)
        .fixed_pos(pointer + egui::vec2(18.0, 18.0))
        .interactable(false)
        .show(ctx, |ui| {
            egui::Frame::popup(ui.style()).show(ui, |ui| {
                ui.set_max_width(300.0);
                if let Some(info) = &info {
                    ui.label(egui::RichText::new(&info.title).strong());
                    for line in &info.lines {
                        ui.label(egui::RichText::new(line).small());
                    }
                    for line in &live {
                        ui.label(egui::RichText::new(line).small().strong());
                    }
                    if let Some((what, samples)) = &info.chart {
                        ui.label(egui::RichText::new(what).small());
                        circuit_plot(ui, samples, &[], "");
                    }
                }
                for (k, section) in probed.iter().enumerate() {
                    if k > 0 || info.is_some() {
                        ui.separator();
                    }
                    ui.label(egui::RichText::new(&section.title).strong());
                    for line in &section.lines {
                        ui.label(egui::RichText::new(line).small());
                    }
                }
            });
        });
}

/// Under a metal map, its values at the point (`crate::MetalMap`): the current density and
/// the power per area of the cell there with their standard errors (the mean over the map's
/// span and the cells about it; under a swinging field the current's lock-in rms).
pub(super) fn metal_map_lines(game: &Game, at: physics::DVec3) -> Vec<String> {
    let (Some(map), Some(view)) = (&game.metal_map, &game.tube_preview) else {
        return Vec::new();
    };
    if !matches!(game.map, Some(MapMode::Current | MapMode::Heat)) || map.revision != view.revision
    {
        return Vec::new();
    }
    let maps = &map.maps;
    for ((m, r), e) in view.run.metals.iter().zip(&maps.values).zip(&maps.errors) {
        let Some(c) = m.cell_at([at.x, at.y]).filter(|c| *c < r.power.len()) else {
            continue;
        };
        let (j, ej) = (r.j[c], e.j[c]);
        let what = if maps.lock_in {
            "eddy current (rms, its way while B rises)"
        } else {
            "current density"
        };
        return vec![
            format!(
                "Here over t = {}–{}: {what} j = ({} ± {}, {} ± {})",
                fmt_num(maps.from),
                fmt_num(maps.to),
                fmt_num(j[0]),
                fmt_num(ej[0]),
                fmt_num(j[1]),
                fmt_num(ej[1]),
            ),
            format!(
                "power per area the field gives the carriers: {} ± {}",
                fmt_num(r.power[c]),
                fmt_num(e.power[c])
            ),
        ];
    }
    Vec::new()
}
