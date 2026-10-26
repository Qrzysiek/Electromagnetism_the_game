use super::*;

fn has(level: &Level, text: &str) -> bool {
    level.model_issues().iter().any(|i| i.contains(text))
}

/// A magnetic moment at finite c next to an electric source is rejected, whether the
/// moment is a shot's or a free particle's, and whether the source is a fixed charge
/// or a charge cloud (the Aharonov–Casher term is not modelled).
#[test]
fn moments_with_electric_sources_are_rejected() {
    let base = shipped("stern_gerlach");
    assert!(base.physics.c.is_some());
    assert!(base.model_issues().is_empty());
    let moment = "magnetic moments with electric fields";
    let mut cloud = base.clone();
    cloud.clouds.push(Cloud {
        center: [2, 2, 0],
        radius: 1.0,
        charge: 1.0,
    });
    assert!(has(&cloud, moment), "a cloud is an electric source");
    let mut free = base.clone();
    let spec = free.shots[0].particle;
    for s in &mut free.shots {
        s.particle.moment = 0.0;
    }
    free.elements.push(Element::charge([2, 2, 0], 1.0));
    assert!(!has(&free, moment));
    free.free_particles.push(FreeParticle {
        particle: spec,
        node: [3, 3, 0],
        velocity: [0.0, 0.0],
        detector: None,
    });
    assert!(has(&free, moment), "a free particle's moment counts");
}

/// Containment checks use each body's own radius (they once added 0.3 to every
/// element), and a launch point inside its own detector is listed.
#[test]
fn containment_uses_the_bodies_radii() {
    let mut l = shipped("high_voltage_dome");
    assert!(l.model_issues().is_empty());
    let sphere = l.conductors[0];
    let c = l.grid.position(sphere.center);
    // A magnet whose surface is 0.1 cells from the sphere's.
    let gap = 0.1;
    l.physics.magnet_radius = 0.6;
    let d = sphere.radius + l.physics.magnet_radius + gap;
    let s = f64::from(l.grid.subdivision);
    #[allow(clippy::cast_possible_truncation)]
    let node = [((c.x + d) * s).round() as i64, (c.y * s).round() as i64, 0];
    let at = l.grid.position(node);
    assert!(((at - c).length() - d).abs() < 0.5 / s, "node rounding");
    l.elements.push(Element::magnet(node, 1.0));
    // Inside with its own radius 0.6, outside with the 0.3 once assumed for all.
    assert!((at - c).length() < sphere.radius + CONTACT_DISTANCE + 0.6);
    assert!((at - c).length() > sphere.radius + CONTACT_DISTANCE + 0.3);
    assert!(has(&l, "inside or at a metal sphere"));
    l.physics.magnet_radius = 0.3;
    assert!(!has(&l, "inside or at a metal sphere"));
    // Now a launch inside its detector.
    let mut m = shipped("first_bend");
    assert!(m.model_issues().is_empty());
    m.shots[0].launch.node = m.shots[0].detector.min;
    assert!(has(&m, "inside its own detector"));
}

/// Radiation reaction near metal is rejected (the metal's response is electrostatic).
#[test]
fn radiation_reaction_near_metal_is_rejected() {
    let mut l = shipped("image_charge");
    assert!(l.model_issues().is_empty());
    l.physics.c = Some(50.0);
    l.physics.radiation_reaction = true;
    assert!(has(&l, "radiation reaction near metal"));
}
