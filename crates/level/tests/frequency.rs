//! A power supply's frequency slider (`Limits::supply_omegas`, `Element::omega` on a supply;
//! DRUDE.md Phase 2 item 2, DEC-30): the frequency a supply sets is the sine's own, bit for
//! bit, and the rules of its placement and setup. Run with
//! `cargo test --release -p level --test frequency`.

use level::tube::TubeRun;
use level::{Element, ElementKind, Level, PlacementError};

/// Two runs the same, bit for bit: every frame's time, particles, charges and potentials,
/// and every goal's collected charge.
fn assert_same_run(a: &TubeRun, b: &TubeRun, what: &str) {
    assert_eq!(a.frames.len(), b.frames.len(), "{what}: frames");
    for (f, g) in a.frames.iter().zip(&b.frames) {
        assert_eq!(f.t.to_bits(), g.t.to_bits(), "{what}: frame time");
        let bits = |v: &[f32]| v.iter().map(|x| x.to_bits()).collect::<Vec<_>>();
        let xs = |v: &[[f32; 2]]| {
            v.iter()
                .flat_map(|x| [x[0].to_bits(), x[1].to_bits()])
                .collect::<Vec<_>>()
        };
        assert_eq!(xs(&f.x), xs(&g.x), "{what}: particles at t = {}", f.t);
        assert_eq!(bits(&f.q), bits(&g.q), "{what}: charges at t = {}", f.t);
        assert_eq!(
            bits(&f.potentials),
            bits(&g.potentials),
            "{what}: potentials at t = {}",
            f.t
        );
        assert_eq!(bits(&f.sigma), bits(&g.sigma), "{what}: surface charges");
    }
    let series = |r: &TubeRun| {
        r.collected
            .iter()
            .flat_map(|s| s.iter().flat_map(|(t, q)| [t.to_bits(), q.to_bits()]))
            .collect::<Vec<_>>()
    };
    assert_eq!(series(a), series(b), "{what}: collected charge");
}

/// The level with a frequency slider over `[lo, hi]`, its reference's supplies on
/// `targets` set to `omega` (the rest of the level unchanged).
fn with_slider(level: &Level, lo: f64, hi: f64, omega: f64) -> (Level, Vec<Element>) {
    let mut l = level.clone();
    l.limits.supply_omegas = vec![lo, hi];
    let placement: Vec<Element> = l
        .reference_solution
        .iter()
        .map(|e| {
            let tuned = e.kind == ElementKind::Supply
                && l.supply_targets()
                    .iter()
                    .any(|&(c, t)| c == e.node && l.supply_frequencies(t).is_some());
            if tuned {
                Element {
                    omega: Some(omega),
                    ..*e
                }
            } else {
                *e
            }
        })
        .collect();
    (l, placement)
}

/// The level with every sine a supply could tune running at `omega` in its own file.
fn with_own(level: &Level, omega: f64) -> Level {
    let mut l = level.clone();
    let tunable: Vec<usize> = l
        .electrodes
        .iter()
        .enumerate()
        .filter(|(_, e)| e.tunable)
        .map(|(i, _)| i)
        .collect();
    let t = l.tube.as_mut().expect("a tube level");
    for c in &mut t.circuits {
        if tunable.contains(&c.electrode) {
            c.source = c.source.with_omega(omega);
        }
    }
    if t.solenoid.is_some() {
        t.solenoid_drive = t.solenoid_drive.map(|s| s.with_omega(omega));
    }
    l
}

/// A supply's frequency is its sine's own: a run with the frequency set on the slider is,
/// bit for bit, the run of the level with that frequency written into its sine. Rectifier
/// (two tube circuits' sines, their supplies setting the amplitudes), Triode (the grid's
/// sine, its supply setting the offset: `supply_offset`), Eddy currents (the solenoid's
/// swing, its supply setting B₀).
#[test]
fn a_supplys_frequency_is_its_sines_own() {
    for (slug, omega, end) in [
        ("rectifier", 0.5, 12.0),
        ("triode", 0.4, 12.0),
        ("eddy_currents", 0.2, 6.0),
    ] {
        let level = level::shipped(slug);
        let (slider, placement) = with_slider(&level, 0.05, 5.0, omega);
        assert!(
            placement.iter().any(|e| e.omega.is_some()),
            "{slug}: no supply with a frequency"
        );
        assert!(
            slider.model_issues().is_empty(),
            "{slug}: {:?}",
            slider.model_issues()
        );
        slider
            .check_placement(&placement)
            .unwrap_or_else(|e| panic!("{slug}: {e:?}"));
        let own = with_own(&level, omega);
        let a = slider.tube_run(&placement, 1, end);
        let b = own.tube_run(&own.reference_solution, 1, end);
        assert!(a.frames.len() > 5, "{slug}: {} frames", a.frames.len());
        assert_same_run(&a, &b, slug);
        // Without a frequency the supply leaves the sine's own.
        let plain = slider.tube_run(&slider.reference_solution, 1, end);
        let original = level.tube_run(&level.reference_solution, 1, end);
        assert_same_run(&plain, &original, slug);
        println!(
            "{slug}: ω = {omega} through the slider: {} frames the same",
            a.frames.len()
        );
    }
}

/// A supply that sets only a frequency (`Level::frequency_only`: a swinging solenoid whose
/// level lists no fields): switched on, its sine's own amplitude at the supply's frequency,
/// bit for bit the run of a supply that sets the same amplitude and frequency; off, no field
/// (`b_z` = 0); its value 0 only, without fields nothing to set.
#[test]
fn a_frequency_only_supply_sets_the_swings_frequency() {
    let level = level::shipped("eddy_currents");
    let b0 = level.reference_solution[0].value;
    let node = level.reference_solution[0].node;
    let mut only = level.clone();
    only.limits.supply_omegas = vec![0.05, 5.0];
    only.limits.solenoid_fields.clear();
    let t = only.tube.as_mut().expect("tube");
    t.solenoid_drive = t.solenoid_drive.map(|s| s.with_level(b0));
    assert!(only.model_issues().is_empty(), "{:?}", only.model_issues());
    assert!(only.frequency_only(level::SupplyTarget::Solenoid));
    let on = Element {
        omega: Some(0.2),
        ..Element::supply(node, 0.0)
    };
    only.check_placement(&[on]).expect("frequency-only supply");
    let (both, placement) = with_slider(&level, 0.05, 5.0, 0.2);
    let a = only.tube_run(&[on], 1, 6.0);
    let b = both.tube_run(&placement, 1, 6.0);
    assert_same_run(&a, &b, "frequency-only");
    // Off: no field.
    let off = only.tube_run(&[], 1, 1.0);
    let peak = off.b_wave.map_or(0.0, |w| w.peak());
    assert!(peak == 0.0, "off: B peaks at {peak}");
    // Its value is 0: nothing else to set.
    assert!(matches!(
        only.check_placement(&[Element { value: 0.5, ..on }]),
        Err(PlacementError::MagnitudeNotAllowed(_))
    ));
    // Without frequencies a player's solenoid needs fields.
    let mut none = only.clone();
    none.limits.supply_omegas.clear();
    assert!(
        none.model_issues()
            .iter()
            .any(|i| i.contains("needs solenoid fields")),
        "{:?}",
        none.model_issues()
    );
}

/// The frequency slider's rules: a frequency only on a supply whose target a sine drives,
/// within the slider's range, and only where the level lists frequencies (placement); the
/// list positive, set on a sine, which starts on the slider (setup issues).
#[test]
fn supply_frequencies_follow_their_rules() {
    let level = level::shipped("rectifier");
    let supply = level.reference_solution[0];
    assert_eq!(supply.kind, ElementKind::Supply);
    let tuned = |w: f64| Element {
        omega: Some(w),
        ..supply
    };
    // No frequencies listed: none allowed.
    assert!(matches!(
        level.check_placement(&[tuned(0.5)]),
        Err(PlacementError::FrequencyNotAllowed(_))
    ));
    let (slider, _) = with_slider(&level, 0.1, 2.0, 0.5);
    assert!(slider.check_placement(&[tuned(0.5)]).is_ok());
    for w in [0.05, 2.5, f64::NAN, f64::INFINITY] {
        assert!(
            matches!(
                slider.check_placement(&[tuned(w)]),
                Err(PlacementError::FrequencyNotAllowed(_))
            ),
            "ω = {w} allowed"
        );
    }
    // A supply on a target no sine drives (Ohm's law: the contact's potential).
    let mut ohm = level::shipped("ohms_law");
    ohm.limits.supply_omegas = vec![0.1, 2.0];
    let s = Element {
        omega: Some(0.5),
        ..ohm.reference_solution[0]
    };
    assert!(matches!(
        ohm.check_placement(&[s]),
        Err(PlacementError::FrequencyNotAllowed(_))
    ));
    let issue = |l: &Level, what: &str| {
        assert!(
            l.model_issues().iter().any(|i| i.contains(what)),
            "no issue '{what}': {:?}",
            l.model_issues()
        );
    };
    issue(&ohm, "driven by a sine");
    let mut negative = slider.clone();
    negative.limits.supply_omegas = vec![-1.0, 2.0];
    issue(&negative, "must be positive");
    // The sine's own frequency off the slider.
    let (off, _) = with_slider(&level, 1.0, 2.0, 1.5);
    issue(&off, "outside its slider");
}
