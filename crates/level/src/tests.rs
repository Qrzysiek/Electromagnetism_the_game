use super::*;

fn sample_level() -> Level {
    Level {
        format_version: FORMAT_VERSION,
        engine_version: ENGINE_VERSION.to_string(),
        name: "test".into(),
        description: String::new(),
        reality: Reality::default(),
        grid: Grid {
            nx: 20,
            ny: 10,
            nz: 0,
            subdivision: 1,
        },
        physics: WorldPhysics {
            c: Some(5.0),
            charge_radius: 0.25,
            magnet_radius: 0.3,
            wire_radius: 0.1,
            antenna_radius: 0.3,
            rf_omega: 1.5,
            radiation_reaction: true,
            beam_interaction: false,
            beam_retarded: false,
            instant_drain: false,
            t_max: 100.0,
            tolerances: TolerancesSpec {
                preview: 1e-10,
                verify: 1e-12,
            },
        },
        shots: vec![Shot {
            particle: ParticleSpec {
                charge: 1.0,
                mass: 1.0,
                radius: 0.0,
                moment: 0.0,
            },
            launch: Launch {
                node: [0, 5, 0],
                direction: [1.0, 0.1, 0.0],
                kinetic_energy: 0.3,
                time: 0.25,
            },
            detector: Detector {
                min: [19, 3, 0],
                max: [20, 7, 0],
                acceptance: Some(DetectorAcceptance {
                    direction: Some([0.0, 20.0]),
                    kinetic: Some([0.1, 1.0]),
                    radiation: None,
                }),
            },
            beam: None,
        }],
        // 0.1 + 0.2 is not exactly representable in decimal.
        elements: vec![Element::charge([10, 5, 0], 0.1 + 0.2)],
        coils: vec![Coil::Circle {
            center: [10, 5, 0],
            radius: 8.0,
            kappa: 0.5,
            rate: 0.25,
            drive: None,
            tunable: false,
        }],
        limits: Limits {
            max_charges: 2,
            magnitudes: vec![1.0, 2.0],
            allow_positive: true,
            allow_negative: true,
            max_magnets: 1,
            magnet_strengths: vec![5.0],
            max_antennas: 1,
            antenna_amplitudes: vec![3.0],
            antenna_omegas: vec![0.5, 2.0],
            continuous: false,
            region: None,
            max_plates: 1,
            plate_voltages: vec![0.0, 5e4],
            plate: PlateSize {
                length: 3.0,
                thickness: 0.5,
                height: 2.0,
            },
            supply_voltages: vec![-1e4, 1e4],
            coil_rates: vec![],
            solenoid_fields: vec![],
            max_free: 2,
            max_cuts: 0,
            free_charges: vec![-1e-6, 2e-6],
            free_speeds: vec![0.0, 0.7],
            free_mass: 3.0,
            free_radius: 0.4,
        },
        reference_solution: vec![],
        disturbances: vec![Disturbance {
            name: "hum".into(),
            e: [0.01, -0.02],
            bz: 0.003,
            waves: vec![Wave {
                amplitude: 0.1,
                direction_deg: 30.0,
                omega: 2.0,
                phase_deg: 90.0,
            }],
        }],
        conductors: vec![],
        clouds: vec![Cloud {
            center: [18, 8, 0],
            radius: 1.0,
            charge: 1.5,
        }],
        free_particles: vec![FreeParticle {
            particle: ParticleSpec {
                charge: 1e-6,
                mass: 4.0,
                radius: 0.3,
                moment: 0.0,
            },
            node: [18, 1, 0],
            velocity: [0.1, -0.2],
            detector: Some(Detector {
                min: [1, 1, 0],
                max: [2, 3, 0],
                acceptance: None,
            }),
        }],
        electrodes: vec![],
        dielectrics: vec![],
        dielectric_spheres: vec![],
        ferrites: vec![],
        gates: vec![],
        tube: None,
    }
}

#[test]
fn json_round_trip_is_exact() {
    let l = sample_level();
    let back = Level::from_json(&l.to_json()).unwrap();
    assert_eq!(l, back);
    assert_eq!(
        l.elements[0].value.to_bits(),
        back.elements[0].value.to_bits()
    );
}

#[test]
fn version_1_files_are_migrated() {
    let v1 = r#"{
            "format_version": 1, "engine_version": "0.1.0", "name": "old",
            "grid": {"nx": 20, "ny": 10, "nz": 0, "subdivision": 1},
            "physics": {"c": 5.0, "charge_radius": 0.3, "t_max": 100.0,
                        "tolerances": {"preview": 1e-10, "verify": 1e-12}},
            "particle": {"charge": 1e-6, "mass": 1.0, "radius": 0.0},
            "launch": {"node": [0, 5, 0], "direction": [1.0, 0.0, 0.0], "kinetic_energy": 0.5},
            "detector": {"min": [18, 3, 0], "max": [20, 7, 0]},
            "level_charges": [{"node": [10, 5, 0], "charge": 2e6}],
            "limits": {"max_charges": 1, "magnitudes": [1e6], "allow_positive": true,
                       "allow_negative": true},
            "reference_solution": [{"node": [4, 4, 0], "charge": -1e6}]
        }"#;
    let l = Level::from_json(v1).unwrap();
    assert_eq!(l.format_version, FORMAT_VERSION);
    assert_eq!(l.shots.len(), 1);
    assert_eq!(l.shots[0].launch.node, [0, 5, 0]);
    assert_eq!(l.elements, vec![Element::charge([10, 5, 0], 2e6)]);
    assert_eq!(l.reference_solution, vec![Element::charge([4, 4, 0], -1e6)]);
    assert_eq!(l.limits.max_magnets, 0);
    assert!(l.disturbances.is_empty());
    assert_eq!(l.flight_count(), 1);
    // Saving writes version 2, which loads back identically.
    assert_eq!(Level::from_json(&l.to_json()).unwrap(), l);
}

/// A static level with one tunable electrode (at (5, 8), along x), room for two player
/// plates and no coils.
fn plate_level() -> Level {
    let mut l = sample_level();
    l.coils.clear();
    l.disturbances.clear();
    // Single flights: no dynamic particles (with them every shot flies at t = 0).
    l.free_particles.clear();
    l.limits.max_free = 0;
    l.limits.max_antennas = 0;
    // Metal does not reflect radiation: no radiation reaction next to it.
    l.physics.radiation_reaction = false;
    l.limits.max_plates = 2;
    l.limits.plate_voltages = vec![0.0, 1e4];
    l.limits.supply_voltages = vec![-2e4, 2e4];
    l.limits.plate = PlateSize::default();
    l.electrodes = vec![Electrode {
        center: [5, 8, 0],
        length: 4.0,
        thickness: 0.4,
        height: 4.0,
        angle_deg: 0.0,
        shape: ElectrodeShape::Box,
        bias: ConductorBias::Grounded,
        tunable: true,
        drive: None,
        label: String::new(),
    }];
    l
}

#[test]
fn plate_and_supply_placement() {
    let l = plate_level();
    let plate = |x, y, v, a| Element::plate([x, y, 0], v, a);
    let err = |p: &[Element]| l.check_placement(p);
    assert_eq!(err(&[plate(14, 2, 1e4, 0.0)]), Ok(()));
    assert_eq!(err(&[plate(14, 2, 0.0, 90.0)]), Ok(()), "grounded plate");
    assert!(matches!(
        err(&[plate(14, 2, 1e4, 45.0)]),
        Err(PlacementError::AngleNotAllowed(_))
    ));
    // Potentials are sliders: anything in the list's range, nothing beyond it.
    assert_eq!(err(&[plate(14, 2, 5e3, 0.0)]), Ok(()), "unlisted, in range");
    assert!(matches!(
        err(&[plate(14, 2, 2e4, 0.0)]),
        Err(PlacementError::MagnitudeNotAllowed(_))
    ));
    // Sticking out of the grid (x from −1 to 3).
    assert!(matches!(
        err(&[plate(1, 2, 1e4, 0.0)]),
        Err(PlacementError::OutsideGrid(_))
    ));
    // On the launch point, on the level charge, across the detector.
    for p in [
        plate(2, 5, 1e4, 0.0),
        plate(10, 5, 1e4, 0.0),
        plate(18, 5, 1e4, 0.0),
    ] {
        assert!(
            matches!(err(&[p]), Err(PlacementError::Occupied(_))),
            "{p:?}"
        );
    }
    // Clearance: gap 0.6 < 1 between two plates or to the level electrode; 1.6 is fine.
    assert!(matches!(
        err(&[plate(14, 2, 1e4, 0.0), plate(14, 3, 1e4, 0.0)]),
        Err(PlacementError::Occupied(_))
    ));
    assert_eq!(
        err(&[plate(14, 2, 1e4, 0.0), plate(14, 4, 1e4, 0.0)]),
        Ok(())
    );
    assert!(matches!(
        err(&[plate(5, 7, 1e4, 0.0)]),
        Err(PlacementError::Occupied(_))
    ));
    assert_eq!(err(&[plate(5, 6, 1e4, 0.0)]), Ok(()));
    // Crossing plates overlap without any corner inside the other.
    assert!(matches!(
        err(&[plate(14, 4, 1e4, 0.0), plate(14, 4, 0.0, 90.0)]),
        Err(PlacementError::Occupied(_))
    ));
    // A charge on or in a player plate.
    for c in [[14, 2, 0], [15, 2, 0]] {
        assert!(matches!(
            err(&[plate(14, 2, 1e4, 0.0), Element::charge(c, 1.0)]),
            Err(PlacementError::Occupied(_))
        ));
    }
    assert!(matches!(
        err(&[
            plate(14, 2, 1e4, 0.0),
            plate(14, 4, 1e4, 0.0),
            plate(14, 6, 1e4, 0.0)
        ]),
        Err(PlacementError::TooManyPlates)
    ));
    // Power supplies: on the tunable electrode's centre, once, at a potential on the
    // slider (the list's range).
    assert_eq!(err(&[Element::supply([5, 8, 0], 2e4)]), Ok(()));
    assert_eq!(err(&[Element::supply([5, 8, 0], 1e4)]), Ok(()), "unlisted");
    assert!(matches!(
        err(&[Element::supply([6, 8, 0], 2e4)]),
        Err(PlacementError::NoTunableElectrode(_))
    ));
    assert!(matches!(
        err(&[Element::supply([5, 8, 0], 3e4)]),
        Err(PlacementError::MagnitudeNotAllowed(_))
    ));
    assert!(matches!(
        err(&[
            Element::supply([5, 8, 0], 2e4),
            Element::supply([5, 8, 0], -2e4)
        ]),
        Err(PlacementError::Occupied(_))
    ));
    assert_eq!(l.model_issues(), Vec::<String>::new());
}

/// A player plate is the same electrode as a level electrode at the same place and
/// potential, and a power supply sets the potential of its electrode: identical
/// fields and obstacles.
#[test]
fn plates_and_supplies_are_electrodes() {
    use physics::field::FieldSolver;
    let l = plate_level();
    let player = [
        Element::plate([14, 2, 0], 1e4, 90.0),
        Element::supply([5, 8, 0], -2e4),
    ];
    let mut fixed = l.clone();
    fixed.electrodes[0].tunable = false;
    fixed.electrodes[0].bias = ConductorBias::Potential(-2e4);
    fixed.electrodes.push(Electrode {
        center: [14, 2, 0],
        length: 4.0,
        thickness: 0.4,
        height: 4.0,
        angle_deg: 90.0,
        shape: ElectrodeShape::Box,
        bias: ConductorBias::Potential(1e4),
        tunable: false,
        drive: None,
        label: String::new(),
    });
    let (a, oa) = l.field(&player);
    let (b, ob) = fixed.field(&[]);
    assert_eq!(oa, ob);
    for x in [DVec3::new(3.0, 3.0, 0.0), DVec3::new(12.0, 6.0, 0.5)] {
        let (fa, fb) = (a.sample(x, 0.0), b.sample(x, 0.0));
        assert_eq!(fa.phi.to_bits(), fb.phi.to_bits());
        assert_eq!(fa.e, fb.e);
    }
    // The supply matters: without it the electrode is grounded.
    let (c, _) = l.field(&player[..1]);
    let x = DVec3::new(5.0, 6.0, 0.0);
    assert!((c.sample(x, 0.0).phi - a.sample(x, 0.0).phi).abs() > 1e3);
    assert!(l.has_metal(&[]) && !sample_level().has_metal(&[]));
    assert!(sample_level().has_metal(&player[..1]));
}

/// A driven electrode (PHYSICS.md §2.10): the JSON round trip keeps the drive; the
/// circuit charges it from its bias towards the source's level, which a power supply
/// on a tunable electrode sets; and a drive's issues are listed.
#[test]
fn driven_electrodes() {
    let mut l = plate_level();
    l.electrodes[0].drive = Some(Drive {
        source: Source::Dc { value: 1e4 },
        resistance: 0.05,
        inductance: 0.0,
        capacitance: 0.0,
        switch: None,
    });
    assert!(l.has_drives());
    assert_eq!(Level::from_json(&l.to_json()).unwrap(), l);
    let end = l.drive_end();
    let potential = |player: &[Element]| {
        let (f, _) = l.field(player);
        let d = f.drives.clone().expect("driven");
        let n = d
            .electrodes
            .iter()
            .find(|&&(e, _)| e == 0)
            .expect("in it")
            .1;
        (d.solution.potential(n, 0.0), d.solution.potential(n, end))
    };
    // From the grounded bias to the source's level (RC ≪ t_max).
    let (v0, v1) = potential(&[]);
    assert!(v0.abs() < 1e-6 && (v1 - 1e4).abs() < 1e-2, "{v0} {v1}");
    // A supply sets the level; the bias stays the start.
    let (v0, v1) = potential(&[Element::supply([5, 8, 0], -2e4)]);
    assert!(v0.abs() < 1e-6 && (v1 + 2e4).abs() < 1e-2, "{v0} {v1}");
    let drive_issue = |l: &Level, text: &str| l.model_issues().iter().any(|i| i.contains(text));
    assert!(!drive_issue(&l, "drive") && !drive_issue(&l, "circuit"));
    let bound = l.drive_back_action_bound(&[]);
    assert!(bound > 0.0 && bound.is_finite(), "{bound}");
    let mut bad = l.clone();
    bad.electrodes[0].drive.as_mut().expect("driven").resistance = 0.0;
    assert!(drive_issue(&bad, "positive resistance"));
    let mut switched = l.clone();
    {
        let d = switched.electrodes[0].drive.as_mut().expect("driven");
        d.inductance = 1.0;
        d.switch = Some(Switch {
            closed: false,
            at: 1.0,
        });
    }
    assert!(drive_issue(&switched, "switch cannot open"));
    let mut floating = l.clone();
    floating.electrodes[0].bias = ConductorBias::Charge(0.0);
    assert!(drive_issue(&floating, "not floating"));
}

/// A driven coil (DC through R): from its strength κ at t = 0 towards the source's
/// current over c², as the RL circuit's closed form; it needs a finite c and no ramp.
#[test]
fn driven_coils() {
    let mut l = sample_level();
    l.free_particles.clear();
    l.limits.max_free = 0;
    let Coil::Circle { rate, drive, .. } = &mut l.coils[0] else {
        unreachable!("a circle")
    };
    *rate = 0.0;
    *drive = Some(Drive {
        source: Source::Dc { value: 1.0 },
        resistance: 1.0,
        inductance: 0.0,
        capacitance: 0.0,
        switch: None,
    });
    assert_eq!(Level::from_json(&l.to_json()).unwrap(), l);
    let c = l.physics.c.expect("finite c");
    let (f, _) = l.field(&[]);
    let d = f.drives.clone().expect("driven");
    let (_, k) = d.loops[0];
    let kappa = |t: f64| d.solution.current(k, t) * d.per_current;
    // κ∞ + (κ₀ − κ∞) e^{−t/τ}: κ₀ = 0.5, κ∞ = (V/R)/c², τ = L/R with L its
    // self-inductance over c².
    let tau = physics::inductance::ring_self(&f.loops[0]) / (c * c);
    let infinity = 1.0 / (c * c);
    let end = l.drive_end();
    for t in [0.0, 0.1 * end, 0.5 * end, end] {
        let want = infinity + (0.5 - infinity) * libm::exp(-t / tau);
        assert!((kappa(t) - want).abs() < 1e-9, "{t}: {} {want}", kappa(t));
    }
    let drive_issue = |l: &Level, text: &str| l.model_issues().iter().any(|i| i.contains(text));
    assert!(!drive_issue(&l, "drive") && !drive_issue(&l, "circuit"));
    let mut newtonian = l.clone();
    newtonian.physics.c = None;
    assert!(drive_issue(&newtonian, "finite speed of light"));
    let mut ramped = l.clone();
    if let Coil::Circle { rate, .. } = &mut ramped.coils[0] {
        *rate = 0.1;
    }
    assert!(drive_issue(&ramped, "no ramp"));
}

/// Coil supplies: a power supply on a tunable circular coil's terminal sets its ramp
/// rate (a value on `coil_rates`' slider), or, if a circuit drives it, its source's
/// level; not on an untunable coil, not beyond the slider.
#[test]
fn coil_supplies_set_the_ramp_or_the_source() {
    let mut l = sample_level();
    l.free_particles.clear();
    l.limits.max_free = 0;
    l.limits.coil_rates = vec![-0.5, 0.5];
    // The coil (centre (10, 5), radius 8): its terminal at (18, 5).
    let centre = [18, 5, 0];
    let supply = |v: f64| [Element::supply(centre, v)];
    // Not tunable yet: refused.
    assert!(matches!(
        l.check_placement(&supply(0.3)),
        Err(PlacementError::NoTunableElectrode(_))
    ));
    let Coil::Circle { tunable, .. } = &mut l.coils[0] else {
        unreachable!("a circle")
    };
    *tunable = true;
    assert_eq!(Level::from_json(&l.to_json()).unwrap(), l);
    assert_eq!(l.check_placement(&supply(0.3)), Ok(()));
    assert!(matches!(
        l.check_placement(&supply(0.7)),
        Err(PlacementError::MagnitudeNotAllowed(_))
    ));
    // The ramp: the supply's rate, else the coil's own (0.25).
    assert_eq!(
        l.field(&supply(0.3)).0.loops[0].rate.to_bits(),
        0.3f64.to_bits()
    );
    assert_eq!(l.field(&[]).0.loops[0].rate.to_bits(), 0.25f64.to_bits());
    // Driven: the supply sets the DC source's value, κ∞ = (V/R)/c².
    let Coil::Circle { rate, drive, .. } = &mut l.coils[0] else {
        unreachable!("a circle")
    };
    *rate = 0.0;
    *drive = Some(Drive {
        source: Source::Dc { value: 1.0 },
        resistance: 1.0,
        inductance: 0.0,
        capacitance: 0.0,
        switch: None,
    });
    // Its supply is now a voltage (`supply_voltages`, here ±1e4).
    assert_eq!(l.supply_targets()[0].0, centre);
    assert_eq!(
        l.supply_list(l.supply_targets()[0].1),
        &l.limits.supply_voltages[..]
    );
    let c = l.physics.c.expect("finite c");
    let end = l.drive_end();
    for (placement, v) in [(&supply(0.4)[..], 0.4), (&[][..], 1.0)] {
        let (f, _) = l.field(placement);
        let d = f.drives.clone().expect("driven");
        let (_, k) = d.loops[0];
        let tau = physics::inductance::ring_self(&f.loops[0]) / (c * c);
        let infinity = v / (c * c);
        let want = infinity + (0.5 - infinity) * libm::exp(-end / tau);
        let got = d.solution.current(k, end) * d.per_current;
        assert!((got - want).abs() < 1e-9, "{v}: {got} {want}");
    }
}

/// Refining keeps metal in place too (metal sphere and electrode centres are nodes),
/// and power supplies stay on their electrodes.
#[test]
fn refinement_keeps_metal_in_place() {
    let mut l = plate_level();
    l.conductors.push(Conductor {
        center: [15, 7, 0],
        radius: 1.0,
        bias: ConductorBias::Grounded,
    });
    let mut player = vec![Element::supply([5, 8, 0], 2e4)];
    let before = l.scenario(0, &player).obstacles;
    l.refine(2, &mut player);
    assert_eq!(before, l.scenario(0, &player).obstacles);
    assert_eq!(l.check_placement(&player), Ok(()));
}

#[test]
fn refinement_keeps_positions() {
    let mut l = sample_level();
    let before = l.scenario(0, &[]);
    let mut player = vec![Element::charge([3, 4, 0], 1.0)];
    let p_before = l.grid.position(player[0].node);
    l.refine(4, &mut player);
    let after = l.scenario(0, &[]);
    assert_eq!(before.x0, after.x0);
    assert_eq!(before.obstacles, after.obstacles);
    assert_eq!(p_before, l.grid.position(player[0].node));
}

#[test]
fn placement_limits() {
    let l = sample_level();
    let ok = Element::charge([5, 5, 0], -2.0);
    assert_eq!(l.check_placement(&[ok]), Ok(()));
    assert!(matches!(
        l.check_placement(&[Element::charge([5, 5, 0], 3.0)]),
        Err(PlacementError::MagnitudeNotAllowed(_))
    ));
    assert!(matches!(
        l.check_placement(&[Element::charge([10, 5, 0], 1.0)]),
        Err(PlacementError::Occupied(_))
    ));
    assert!(matches!(
        l.check_placement(&[ok, ok, ok]),
        Err(PlacementError::TooManyCharges)
    ));
    let magnet = Element::magnet([6, 6, 0], -5.0);
    assert_eq!(l.check_placement(&[ok, magnet]), Ok(()));
    assert!(matches!(
        l.check_placement(&[magnet, Element::magnet([7, 6, 0], 5.0)]),
        Err(PlacementError::TooManyMagnets)
    ));
    let antenna = Element::antenna([8, 2, 0], -3.0, 135.0).with_omega(2.0);
    assert_eq!(l.check_placement(&[ok, magnet, antenna]), Ok(()));
    assert!(matches!(
        l.check_placement(&[Element::antenna([8, 2, 0], 3.0, 30.0).with_omega(2.0)]),
        Err(PlacementError::AngleNotAllowed(_))
    ));
    // The level lists frequencies, so the generator frequency (None) or an
    // unlisted one is not allowed.
    assert!(matches!(
        l.check_placement(&[Element::antenna([8, 2, 0], 3.0, 0.0)]),
        Err(PlacementError::FrequencyNotAllowed(None))
    ));
    assert!(matches!(
        l.check_placement(&[Element::antenna([8, 2, 0], 3.0, 0.0).with_omega(1.0)]),
        Err(PlacementError::FrequencyNotAllowed(Some(_)))
    ));
    assert!(matches!(
        l.check_placement(&[
            antenna,
            Element::antenna([9, 2, 0], 3.0, 0.0).with_omega(0.5)
        ]),
        Err(PlacementError::TooManyAntennas)
    ));
    // Hardcore: anything within the ranges, any orientation.
    let mut hard = l.clone();
    hard.limits.continuous = true;
    let odd = Element::antenna([8, 2, 0], -3.0, 17.0).with_omega(1.3);
    assert!(l.check_placement(&[odd]).is_err());
    assert_eq!(
        hard.check_placement(&[odd, Element::charge([5, 5, 0], 1.7)]),
        Ok(())
    );
    assert!(matches!(
        hard.check_placement(&[Element::charge([5, 5, 0], 2.1)]),
        Err(PlacementError::MagnitudeNotAllowed(_))
    ));
    assert!(matches!(
        hard.check_placement(&[odd.with_omega(2.5)]),
        Err(PlacementError::FrequencyNotAllowed(_))
    ));
    // Each antenna oscillates at its own frequency.
    let (field, _) = l.field(&[antenna]);
    assert_eq!(field.antennas[0].omega.to_bits(), 2.0_f64.to_bits());
}
