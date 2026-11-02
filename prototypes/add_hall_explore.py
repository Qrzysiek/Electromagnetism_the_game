"""Adds a temporary test printing the field and carriers across each strip of a Hall level."""
p = "crates/level/tests/tubes.rs"
s = open(p, encoding="utf-8").read()
anchor = "/// The metals' chevrons move with the carriers' drift"
new = '''/// TEMP exploration: the field across the strips of the Hall levels.
#[test]
#[ignore = "exploration"]
fn hall_profile_explore() {
    use physics::DVec3;
    for (slug, metals) in [("hall_probe", 1), ("positive_carriers", 2)] {
        let level = level::shipped(slug);
        let run = level.tube_run(&level.reference_solution, 1, 41.0);
        println!("{slug}: b_z {}", run.b_z);
        for metal in 0..metals {
            let m = &run.metals[metal];
            let (lo, hi) = (m.lo, m.hi);
            let cx = 0.5 * (lo.x + hi.x);
            let len = hi.x - lo.x;
            let depths: Vec<f64> = (0..9).map(|i| lo.y + 0.5 * f64::from(i)).collect();
            let mut ey = vec![0.0; depths.len()];
            let mut phi_side = [0.0; 2];
            let mut density = vec![0.0; 8];
            let (mut ratio, mut electric, mut magnetic, mut drift) = (0.0, 0.0, 0.0, 0.0);
            let mut n = 0.0;
            for k in 0..=10 {
                let t = 20.0 + 2.0 * f64::from(k);
                let now = run.frame_between(t).expect("a frame");
                let before = run.frame_between(t - 1.0).expect("a frame");
                let h = run.hall_forces(&now, &before, metal).expect("hall");
                ratio += h.electric.length() / h.magnetic.length();
                electric += h.electric.length();
                magnetic += h.magnetic.length();
                drift += h.drift.x;
                n += 1.0;
                for (i, y) in depths.iter().enumerate() {
                    for j in 0..7 {
                        let x = cx + 0.3 * len * (f64::from(j) / 3.0 - 1.0);
                        ey[i] += run.e_at(&now, DVec3::new(x, *y, 0.0)).y / 7.0;
                    }
                }
                for (side, y) in [lo.y, hi.y].iter().enumerate() {
                    for j in 0..7 {
                        let x = cx + 0.3 * len * (f64::from(j) / 3.0 - 1.0);
                        phi_side[side] += run.field_at(&now, DVec3::new(x, *y, 0.0)).0 / 7.0;
                    }
                }
                for (x, kind) in now.x.iter().zip(&now.kind) {
                    if *kind > 1 {
                        continue;
                    }
                    let (px, py) = (f64::from(x[0]), f64::from(x[1]));
                    if (px - cx).abs() < 0.3 * len && py >= lo.y && py < hi.y {
                        let b = (((py - lo.y) / 0.5) as usize).min(7);
                        density[b] += 1.0;
                    }
                }
            }
            println!(
                "  metal {metal}: ratio {:.3}, |qE| {:.4}, |qvB| {:.4}, drift x {:.4}",
                ratio / n,
                electric / n,
                magnetic / n,
                drift / n
            );
            let row: Vec<String> = ey.iter().map(|e| format!("{:.3}", e / n)).collect();
            println!("  E_y at y = lo + 0, 0.5, .. 4: {}", row.join(" "));
            println!(
                "  potential at the sides (mid 3/5): lo {:.3}, hi {:.3}, hi - lo {:.3}",
                phi_side[0] / n,
                phi_side[1] / n,
                (phi_side[1] - phi_side[0]) / n
            );
            let total: f64 = density.iter().sum();
            let row: Vec<String> = density.iter().map(|d| format!("{:.3}", 8.0 * d / total)).collect();
            println!("  carriers per half-unit band (relative): {}", row.join(" "));
        }
    }
}

'''
assert s.count(anchor) == 1
s = s.replace(anchor, new + anchor)
open(p, "w", encoding="utf-8", newline="\n").write(s)
print("ok")
