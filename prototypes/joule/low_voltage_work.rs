//! Probe (DEC-33's history): the field's work against I V at a low supply in Ohm's law.

#[test]
#[ignore = "a probe"]
fn ohms_law_low_voltage_work() {
    let level = level::shipped("ohms_law");
    let node = level.reference_solution[0].node;
    for v in [30.0, 4.0] {
        let placement = vec![level::Element::supply(node, v)];
        let run = level.tube_run(&placement, 1, 40.0);
        let work = |f: &level::tube::TubeFrame| -> f64 {
            f.cells
                .iter()
                .flat_map(|c| &c.work)
                .map(|w| f64::from(*w))
                .sum()
        };
        let f1 = run.frames.last().expect("frames");
        let k0 = run.frames.partition_point(|f| f.t < 10.0);
        let f0 = &run.frames[k0];
        println!(
            "Ohm's law at {v}: work rate {:.3}, goal current {:?} (I V = {:.3})",
            (work(f1) - work(f0)) / (f1.t - f0.t),
            run.currents,
            run.currents.first().copied().flatten().unwrap_or(f64::NAN) * v
        );
    }
}
