"""Adds a temporary test listing the tube levels' leads that cross a metal."""
p = "crates/level/tests/tubes.rs"
s = open(p, encoding="utf-8").read()
anchor = "/// The metals' chevrons move with the carriers' drift"
new = '''/// TEMP exploration: leads that cross a metal.
#[test]
#[ignore = "exploration"]
fn lead_metal_explore() {
    for slug in ["ohms_law", "long_and_thin", "potentiometer", "hall_probe", "positive_carriers"] {
        let level = level::shipped(slug);
        let run = level.tube_run(&level.reference_solution, 1, 0.2);
        let em = level::tube_em::TubeEm::new(&level, &run);
        for (e, lead) in em.leads.iter().enumerate() {
            let Some(((ax, ay), (bx, by))) = lead else {
                continue;
            };
            let crosses: Vec<usize> = run
                .metals
                .iter()
                .enumerate()
                .filter(|(_, m)| {
                    (1..200).any(|k| {
                        let u = f64::from(k) / 200.0;
                        let (x, y) = (ax + (bx - ax) * u, ay + (by - ay) * u);
                        x > m.lo.x && x < m.hi.x && y > m.lo.y && y < m.hi.y
                    })
                })
                .map(|(i, _)| i)
                .collect();
            println!(
                "{slug} electrode {}: ({ax:.2}, {ay:.2}) -> ({bx:.2}, {by:.2}) crosses metals {crosses:?}",
                e + 1
            );
        }
    }
}

'''
assert s.count(anchor) == 1
s = s.replace(anchor, new + anchor)
open(p, "w", encoding="utf-8", newline="\n").write(s)
print("ok")
