//! Probe (DEC-41, PHYSICS.md §2.11 *The contacts' layers*): the field's work at rest in a metal with contacts, against
//! the resolution's parts (which one the contacts' bias scales with).

use level::tube::TubeResolution;

fn rest_work(slug: &str, value: Option<f64>, res: TubeResolution) -> (f64, f64) {
    let level = level::shipped(slug);
    let mut player = level.reference_solution.clone();
    if let Some(v) = value {
        player[0].value = v;
    }
    let mut sim = level.tube_sim_at(&player, res);
    sim.advance_to(41.0);
    let r = sim.take_chunk();
    let at = |series: &[(f64, f64)], t: f64| {
        let k = series
            .partition_point(|p| p.0 <= t)
            .max(1)
            .min(series.len() - 1);
        let (a, b) = (series[k - 1], series[k]);
        a.1 + (b.1 - a.1) * (t - a.0) / (b.0 - a.0)
    };
    let work = (at(&r.work[0], 40.0) - at(&r.work[0], 10.0)) / 30.0;
    let current = r.currents.first().copied().flatten().unwrap_or(f64::NAN);
    (work, current)
}

#[test]
#[ignore = "a probe"]
fn contact_bias_scaling() {
    use rayon::prelude::*;
    let base = TubeResolution::refined(1);
    let cases: Vec<(&str, Option<f64>, &str, TubeResolution)> = vec![
        ("ohms_law", Some(0.0), "base", base),
        (
            "ohms_law",
            Some(0.0),
            "particles x2",
            TubeResolution {
                particles: 2.0,
                ..base
            },
        ),
        (
            "ohms_law",
            Some(0.0),
            "emitter x2",
            TubeResolution {
                emitter: 2.0,
                ..base
            },
        ),
        (
            "ohms_law",
            Some(0.0),
            "softening x2",
            TubeResolution {
                softening: 2.0,
                ..base
            },
        ),
        (
            "ohms_law",
            Some(0.0),
            "refine 2",
            TubeResolution::refined(2),
        ),
        ("ohms_law", None, "reference", base),
    ];
    let out: Vec<(f64, f64)> = cases
        .par_iter()
        .map(|(slug, v, _, res)| rest_work(slug, *v, *res))
        .collect();
    for ((slug, v, name, _), (w, i)) in cases.iter().zip(&out) {
        println!("{slug} V {v:?} {name}: field's work {w:.3} per time, current {i:.4}");
    }
}

/// The field's work and the carriers' density by column across Ohm's law's strip at rest
/// (with the account cells made 0.1 wide for the probe).
#[test]
#[ignore = "a probe"]
fn contact_bias_profile() {
    let level = level::shipped("ohms_law");
    let mut player = level.reference_solution.clone();
    player[0].value = 0.0;
    let r = level.tube_run(&player, 1, 41.0);
    let m = &r.metals[0];
    let (nx, ny) = m.cells;
    let k0 = r.frames.partition_point(|f| f.t < 10.0);
    let k1 = r
        .frames
        .partition_point(|f| f.t < 40.0)
        .min(r.frames.len() - 1);
    let (f0, f1) = (&r.frames[k0], &r.frames[k1]);
    let span = f1.t - f0.t;
    let dx = (m.hi.x - m.lo.x) / nx as f64;
    let mut density = vec![0.0; nx];
    let mut frames = 0.0;
    for f in &r.frames[k0..k1] {
        frames += 1.0;
        for (x, kind) in f.x.iter().zip(&f.kind) {
            let _ = kind;
            let u = (f64::from(x[0]) - m.lo.x) / dx;
            if u >= 0.0 && (u as usize) < nx && f64::from(x[1]) > m.lo.y && f64::from(x[1]) < m.hi.y
            {
                density[u as usize] += 1.0;
            }
        }
    }
    let mut total = 0.0;
    println!(
        "metal x {}..{}, {nx} x {ny} cells, frames {frames}",
        m.lo.x, m.hi.x
    );
    for i in 0..nx {
        let w: f64 = (0..ny)
            .map(|j| {
                f64::from(f1.cells[0].work[i + nx * j]) - f64::from(f0.cells[0].work[i + nx * j])
            })
            .sum::<f64>()
            / span;
        total += w;
        if i < 25 || i + 25 >= nx {
            println!(
                "x {:.2}: work {w:+.4} per time, carriers {:.1} per frame",
                m.lo.x + dx * (i as f64 + 0.5),
                density[i] / frames
            );
        }
    }
    println!("total {total:+.3}");
}

fn panel_heat(slug: &str, value: Option<f64>, res: TubeResolution) -> (f64, f64, f64) {
    let level = level::shipped(slug);
    let mut player = level.reference_solution.clone();
    if let Some(v) = value {
        player[0].value = v;
    }
    let mut sim = level.tube_sim_at(&player, res);
    sim.advance_to(41.0);
    let r = sim.take_chunk();
    let sum = |f: &level::tube::TubeFrame| -> f64 {
        f.cells
            .iter()
            .flat_map(|c| &c.work)
            .map(|w| f64::from(*w))
            .sum()
    };
    let k0 = r.frames.partition_point(|f| f.t < 10.0);
    let k1 = r
        .frames
        .partition_point(|f| f.t < 40.0)
        .min(r.frames.len() - 1);
    let (f0, f1) = (&r.frames[k0], &r.frames[k1]);
    let heat = (sum(f1) - sum(f0)) / (f1.t - f0.t);
    let current = r.currents.first().copied().flatten().unwrap_or(f64::NAN);
    // The supply's potential: the reference's first element unless set.
    let volts = value.unwrap_or(level.reference_solution[0].value);
    (heat, current, volts)
}

#[test]
#[ignore = "a probe"]
fn macro_heat_against_iv() {
    use rayon::prelude::*;
    let base = TubeResolution::refined(1);
    let cases: Vec<(&str, Option<f64>, &str, TubeResolution)> = vec![
        ("ohms_law", None, "reference", base),
        ("ohms_law", Some(0.0), "rest", base),
        (
            "ohms_law",
            Some(0.0),
            "rest particles x2",
            TubeResolution {
                particles: 2.0,
                ..base
            },
        ),
        (
            "ohms_law",
            Some(0.0),
            "rest refine 2",
            TubeResolution::refined(2),
        ),
        ("ohms_law", Some(4.0), "4 V", base),
        ("corbino_disc", None, "reference", base),
    ];
    let out: Vec<(f64, f64, f64)> = cases
        .par_iter()
        .map(|(slug, v, _, res)| panel_heat(slug, *v, *res))
        .collect();
    for ((slug, _, name, _), (h, i, v)) in cases.iter().zip(&out) {
        println!(
            "{slug} {name}: heat {h:.3} per time, current {i:.4}, I V {:.3}",
            i * v
        );
    }
}

#[test]
#[ignore = "a probe"]
fn macro_heat_profile() {
    let level = level::shipped("ohms_law");
    for v in [0.0, 30.0] {
        let mut player = level.reference_solution.clone();
        player[0].value = v;
        let r = level.tube_run(&player, 1, 41.0);
        let m = &r.metals[0];
        let (nx, ny) = m.cells;
        let k0 = r.frames.partition_point(|f| f.t < 10.0);
        let k1 = r
            .frames
            .partition_point(|f| f.t < 40.0)
            .min(r.frames.len() - 1);
        let (f0, f1) = (&r.frames[k0], &r.frames[k1]);
        let span = f1.t - f0.t;
        let cols: Vec<String> = (0..nx)
            .map(|i| {
                let w: f64 = (0..ny)
                    .map(|j| {
                        f64::from(f1.cells[0].work[i + nx * j])
                            - f64::from(f0.cells[0].work[i + nx * j])
                    })
                    .sum::<f64>()
                    / span;
                format!("{w:+.2}")
            })
            .collect();
        let rows: Vec<String> = (0..ny)
            .map(|j| {
                let w: f64 = (0..nx)
                    .map(|i| {
                        f64::from(f1.cells[0].work[i + nx * j])
                            - f64::from(f0.cells[0].work[i + nx * j])
                    })
                    .sum::<f64>()
                    / span;
                format!("{w:+.2}")
            })
            .collect();
        println!("V {v}: columns {}", cols.join(" "));
        println!("V {v}: rows {}", rows.join(" "));
    }
}
