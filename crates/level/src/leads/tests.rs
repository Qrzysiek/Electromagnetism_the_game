//! The wiring's router (DEC-26): wires along x and y with right-angle bends (but the stubs onto
//! their ends), clear of everything they do not end on, sharing no node, the same every time.

use super::*;

/// The shipped tube levels (`levels/*.json`), by file name.
fn shipped_tube_levels() -> Vec<(String, Level)> {
    let dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../levels");
    let mut out: Vec<(String, Level)> = std::fs::read_dir(&dir)
        .expect("levels directory")
        .filter_map(|e| e.ok().map(|e| e.path()))
        .filter(|p| crate::is_level_file(p))
        .filter_map(|p| {
            let level = Level::from_json(&std::fs::read_to_string(&p).ok()?).ok()?;
            level.tube.is_some().then(|| {
                let name = p.file_stem().map(|s| s.to_string_lossy().into_owned());
                (name.unwrap_or_default(), level)
            })
        })
        .collect();
    out.sort_by(|a, b| a.0.cmp(&b.0));
    out
}

fn v(p: (f64, f64)) -> DVec3 {
    DVec3::new(p.0, p.1, 0.0)
}

/// Whether segment `a`–`b` runs along x or y.
fn straight(a: (f64, f64), b: (f64, f64)) -> bool {
    (a.0 - b.0).abs() < 1e-9 || (a.1 - b.1).abs() < 1e-9
}

/// What is wrong with a wire ending on the electrodes `ends`: a segment off the axes (but a
/// short stub at either end), or a point of its run (the stubs aside) closer than the
/// clearance to an electrode it does not end on or to a metal.
fn faults(
    wire: &Wire,
    outlines: &[Outline],
    metals: &[(DVec3, DVec3)],
    ends: &[usize],
) -> Vec<String> {
    let mut out = Vec::new();
    if wire.len() < 2 {
        out.push(format!("{} points", wire.len()));
        return out;
    }
    let last = wire.len() - 2;
    for (k, s) in wire.windows(2).enumerate() {
        let stub = k == 0 || k == last;
        let len = (v(s[1]) - v(s[0])).length();
        if !straight(s[0], s[1]) && !(stub && len <= CLEARANCE + ROUTE_SPACING + 1e-9) {
            out.push(format!("segment {k} {:?}–{:?} off the axes", s[0], s[1]));
        }
        if stub && wire.len() > 2 {
            continue;
        }
        #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
        let n = (len / 0.05).ceil().max(1.0) as usize;
        for i in 0..=n {
            #[allow(clippy::cast_precision_loss)]
            let p = v(s[0]).lerp(v(s[1]), i as f64 / n as f64);
            for (e, o) in outlines.iter().enumerate() {
                if !ends.contains(&e) && o.distance_to(p) < CLEARANCE - 1e-9 {
                    out.push(format!(
                        "segment {k} within {CLEARANCE} of electrode {e} at {p}"
                    ));
                    return out;
                }
            }
            for (m, (a, b)) in metals.iter().enumerate() {
                let dx = (a.x - p.x).max(p.x - b.x).max(0.0);
                let dy = (a.y - p.y).max(p.y - b.y).max(0.0);
                if libm::hypot(dx, dy) < CLEARANCE - 1e-9 {
                    out.push(format!(
                        "segment {k} within {CLEARANCE} of metal {m} at {p}"
                    ));
                    return out;
                }
            }
        }
    }
    out
}

/// The grid nodes a wire's run passes (its stubs aside).
#[allow(clippy::cast_possible_truncation)]
fn nodes(wire: &Wire) -> Vec<(i64, i64)> {
    let key = |p: DVec3| {
        (
            libm::round(p.x / ROUTE_SPACING) as i64,
            libm::round(p.y / ROUTE_SPACING) as i64,
        )
    };
    let mut out = Vec::new();
    let (from, to) = (
        usize::from(wire.len() > 2),
        wire.len() - usize::from(wire.len() > 2),
    );
    for s in wire[from..to].windows(2) {
        let (a, b) = (v(s[0]), v(s[1]));
        #[allow(clippy::cast_sign_loss, clippy::cast_possible_truncation)]
        let n = libm::round((b - a).length() / ROUTE_SPACING) as usize;
        for i in 0..=n {
            #[allow(clippy::cast_precision_loss)]
            out.push(key(a.lerp(b, i as f64 / n.max(1) as f64)));
        }
    }
    out.sort_unstable();
    out.dedup();
    out
}

/// Every shipped tube level is wired: each electrode due a lead has one (all but an open
/// level's outermost envelope), every wire runs along x and y but its stubs, keeps clear of
/// what it does not end on, and no two wires share a node (but beside an electrode both end
/// on); and the same plan comes out twice.
#[test]
fn every_shipped_tube_level_is_wired() {
    let mut problems = Vec::new();
    for (name, level) in shipped_tube_levels() {
        let spec = level.tube.as_ref().expect("a tube level");
        let player = &level.reference_solution;
        let started = std::time::Instant::now();
        let plan = level.tube_leads(player);
        let took = started.elapsed().as_secs_f64();
        assert_eq!(plan, level.tube_leads(player), "{name}: not deterministic");
        let outlines = level.electrode_outlines(player);
        let metals: Vec<(DVec3, DVec3)> = spec.metals.iter().map(|m| level.metal_rect(m)).collect();
        let open = spec.output.is_some();
        // Each wire with the electrodes it ends on.
        let mut wires: Vec<(String, &Wire, Vec<usize>)> = Vec::new();
        for (e, lead) in plan.leads.iter().enumerate() {
            let due = !(plan.parent[e].is_none() && open && plan.envelope[e]);
            match lead {
                Some(w) => {
                    let ends: Vec<usize> = std::iter::once(e).chain(plan.parent[e]).collect();
                    wires.push((format!("lead {e}"), w, ends));
                }
                None if due => problems.push(format!("{name}: electrode {e} has no lead")),
                None => {}
            }
        }
        let straps: Vec<(usize, usize)> = spec
            .circuits
            .iter()
            .filter(|c| c.resistance.is_none())
            .flat_map(|c| c.linked.iter().map(move |&e| (c.electrode, e)))
            .collect();
        if straps.len() != plan.straps.len() {
            problems.push(format!(
                "{name}: {} straps routed of {}",
                plan.straps.len(),
                straps.len()
            ));
        }
        for (k, (w, (a, b))) in plan.straps.iter().zip(&straps).enumerate() {
            wires.push((format!("strap {k}"), w, vec![*a, *b]));
        }
        for (what, w, ends) in &wires {
            for f in faults(w, &outlines, &metals, ends) {
                problems.push(format!("{name}: {what}: {f}"));
            }
        }
        // Shared nodes, but beside an electrode both wires end on (a junction).
        let mut shared = 0;
        for (i, (_, wa, ea)) in wires.iter().enumerate() {
            for (_, wb, eb) in &wires[i + 1..] {
                let common: Vec<usize> = ea.iter().filter(|e| eb.contains(e)).copied().collect();
                let nb = nodes(wb);
                shared += nodes(wa)
                    .into_iter()
                    .filter(|n| nb.binary_search(n).is_ok())
                    .filter(|n| {
                        #[allow(clippy::cast_precision_loss)]
                        let p = DVec3::new(n.0 as f64, n.1 as f64, 0.0) * ROUTE_SPACING;
                        !common
                            .iter()
                            .any(|&e| outlines[e].distance_to(p) <= CLEARANCE + 2.0 * ROUTE_SPACING)
                    })
                    .count();
            }
        }
        if shared > 0 {
            problems.push(format!("{name}: {shared} nodes shared by two wires"));
        }
        println!(
            "WIRED {name}: {} leads, {} straps, {shared} shared nodes, {:.1} ms",
            plan.leads.iter().flatten().count(),
            plan.straps.len(),
            1e3 * took
        );
    }
    assert!(problems.is_empty(), "{}", problems.join("\n"));
}

/// A box with nothing between it and the border is wired straight to it: no bend.
#[test]
fn a_free_electrode_gets_a_straight_lead() {
    let rect = |x0: f64, y0: f64, x1: f64, y1: f64| {
        Outline::Rect([
            DVec3::new(x0, y0, 0.0),
            DVec3::new(x1, y0, 0.0),
            DVec3::new(x1, y1, 0.0),
            DVec3::new(x0, y1, 0.0),
        ])
    };
    let outlines = [rect(4.0, 1.0, 6.0, 2.0)];
    let mut router = Router::new(DVec3::ZERO, DVec3::new(10.0, 10.0, 0.0), &outlines, &[]);
    let wire = router
        .route(End::Border, End::Electrode(0))
        .expect("a route");
    assert!(
        wire.iter().all(|p| (p.0 - wire[0].0).abs() < 1e-9)
            || wire.iter().all(|p| (p.1 - wire[0].1).abs() < 1e-9),
        "{wire:?}"
    );
    assert!(faults(&wire, &outlines, &[], &[0]).is_empty(), "{wire:?}");
}

/// A plate between a box and the nearest border: the lead goes round it, keeping clear,
/// along the axes; a second lead from the other side does not run on the first's nodes.
#[test]
fn leads_go_round_what_is_in_the_way() {
    let rect = |x0: f64, y0: f64, x1: f64, y1: f64| {
        Outline::Rect([
            DVec3::new(x0, y0, 0.0),
            DVec3::new(x1, y0, 0.0),
            DVec3::new(x1, y1, 0.0),
            DVec3::new(x0, y1, 0.0),
        ])
    };
    // The box low on the board; the plate wide under it.
    let outlines = [
        rect(9.0, 3.0, 11.0, 4.0),
        rect(3.0, 1.5, 17.0, 2.0),
        rect(9.0, 6.0, 11.0, 7.0),
    ];
    let mut router = Router::new(DVec3::ZERO, DVec3::new(20.0, 20.0, 0.0), &outlines, &[]);
    let a = router
        .route(End::Border, End::Electrode(0))
        .expect("a route");
    let b = router
        .route(End::Border, End::Electrode(2))
        .expect("a route");
    for (w, e) in [(&a, 0), (&b, 2)] {
        let f = faults(w, &outlines, &[], &[e]);
        assert!(f.is_empty(), "{w:?}: {f:?}");
    }
    let (na, nb) = (nodes(&a), nodes(&b));
    assert!(
        na.iter().all(|n| nb.binary_search(n).is_err()),
        "{a:?} {b:?}"
    );
}

/// A return runs along the border the shorter way round, through the corners it passes, from
/// its start to the foot of its end (then to the end, if off the border).
#[test]
fn returns_go_round_the_border_the_short_way() {
    let (lo, hi) = (DVec3::ZERO, DVec3::new(20.0, 10.0, 0.0));
    let g = DVec3::new(10.0, 0.0, 0.0);
    // To the left side: by the lower left corner.
    assert_eq!(
        border_path(lo, hi, g, DVec3::new(0.0, 6.0, 0.0)),
        vec![g, DVec3::ZERO, DVec3::new(0.0, 6.0, 0.0)]
    );
    // To the upper side's right half: by the lower right and upper right corners.
    assert_eq!(
        border_path(lo, hi, g, DVec3::new(15.0, 10.0, 0.0)),
        vec![
            g,
            DVec3::new(20.0, 0.0, 0.0),
            DVec3::new(20.0, 10.0, 0.0),
            DVec3::new(15.0, 10.0, 0.0)
        ]
    );
    // Along the lower side, and to a point a hair inside it.
    assert_eq!(
        border_path(lo, hi, g, DVec3::new(3.0, 0.0, 0.0)),
        vec![g, DVec3::new(3.0, 0.0, 0.0)]
    );
    assert_eq!(
        border_path(lo, hi, g, DVec3::new(3.0, 0.1, 0.0)),
        vec![g, DVec3::new(3.0, 0.0, 0.0), DVec3::new(3.0, 0.1, 0.0)]
    );
}

/// With the carriers' own field the wiring closes every circuit (DEC-27): Ohm's law's two
/// contacts' leads, from the left and right sides, return along the grounded envelope to the
/// ground terminal below the middle of the lower side; each circuit's path runs from the terminal to
/// its contact (the engine's tube takes them: `tube::tests::self_field_tube_takes_the_paths`).
/// Without a speed of light, no terminal and no returns.
#[test]
fn self_field_circuits_return_to_the_ground_terminal() {
    let mut level = crate::shipped("ohms_law");
    let plain = level.tube_leads(&level.reference_solution);
    assert!(plain.ground.is_none() && plain.returns.iter().all(Option::is_none));
    let spec = level.tube.as_mut().expect("a tube level");
    spec.speed_of_light = Some(3.0);
    spec.step = 0.01;
    let plan = level.tube_leads(&level.reference_solution);
    let b = level.bounds();
    let g = plan.ground.expect("a ground terminal");
    let m = crate::tube::ENVELOPE_MARGIN;
    assert!((g.0 - 0.5 * (b.min.x + b.max.x)).abs() < 1e-12);
    assert!((g.1 - (b.min.y - m)).abs() < 1e-12);
    let paths = plan.circuit_paths();
    assert_eq!(paths.len(), 2);
    for ((back, lead), path) in plan.returns.iter().zip(&plan.leads).zip(&paths) {
        let back = back.as_ref().expect("a return");
        let lead = lead.as_ref().expect("a lead");
        assert_eq!(back.first(), Some(&g));
        assert_eq!(back.last(), lead.first());
        assert!(back.windows(2).all(|w| straight(w[0], w[1])), "{back:?}");
        let path = path.as_ref().expect("a path");
        assert_eq!(path.first().map(|p| (p.x, p.y)), Some(g));
        assert_eq!(path.last().map(|p| (p.x, p.y)), lead.last().copied());
    }
}
