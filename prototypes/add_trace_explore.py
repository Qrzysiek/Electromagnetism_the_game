"""Adds a temporary test replaying the traced-carrier replacement rules."""
p = "crates/level/tests/tubes.rs"
s = open(p, encoding="utf-8").read()
anchor = "/// The magnetron's antenna radiates its load"
new = '''/// TEMP exploration: traced carriers near the contacts under two replacement rules.
#[test]
#[ignore = "exploration"]
fn traced_spread_explore() {
    let level = level::shipped("ohms_law");
    let run = level.tube_run(&level.reference_solution, 1, 81.0);
    let m = &run.metals[0];
    let scatter = |id: u32| u64::from(id).wrapping_mul(0x9e37_79b9_7f4a_7c15).rotate_left(32);
    for rule in ["newest", "hashed"] {
        let present = |f: &level::tube::TubeFrame| -> Vec<(u32, f32)> {
            f.id.iter()
                .zip(&f.x)
                .filter(|(_, x)| {
                    f64::from(x[0]) >= m.lo.x
                        && f64::from(x[0]) <= m.hi.x
                        && f64::from(x[1]) >= m.lo.y
                        && f64::from(x[1]) <= m.hi.y
                })
                .map(|(id, x)| (*id, x[0]))
                .collect()
        };
        let mut p0 = present(&run.frames[0]);
        p0.sort_unstable_by_key(|(id, _)| *id);
        let every = (p0.len() / 40).max(1);
        let mut ids: Vec<u32> = p0.iter().step_by(every).map(|(id, _)| *id).collect();
        let (mut near, mut total) = (0usize, 0usize);
        for f in &run.frames[1..] {
            let p = present(f);
            let alive: std::collections::HashMap<u32, f32> = p.iter().copied().collect();
            let kept: std::collections::HashSet<u32> =
                ids.iter().copied().filter(|id| alive.contains_key(id)).collect();
            let mut candidates: Vec<u32> =
                p.iter().map(|(id, _)| *id).filter(|id| !kept.contains(id)).collect();
            if rule == "newest" {
                candidates.sort_unstable_by(|a, b| b.cmp(a));
            } else {
                candidates.sort_unstable_by_key(|id| scatter(*id));
            }
            let mut c = candidates.into_iter();
            for id in &mut ids {
                if !alive.contains_key(id)
                    && let Some(n) = c.next()
                {
                    *id = n;
                }
            }
            if f.t >= 20.0 {
                for id in &ids {
                    let x = f64::from(alive[id]);
                    total += 1;
                    if x - m.lo.x < 0.5 || m.hi.x - x < 0.5 {
                        near += 1;
                    }
                }
            }
        }
        println!("{rule}: {:.3} of the traced within 0.5 of the contacts (t = 20-80)", near as f64 / total as f64);
    }
}

'''
assert s.count(anchor) == 1
s = s.replace(anchor, new + anchor)
open(p, "w", encoding="utf-8", newline="\n").write(s)
print("ok")
