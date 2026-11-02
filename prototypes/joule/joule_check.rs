//! Scratch (not committed): a tube level's reference run, the lattice's heat (`joule`) and
//! the heat goals' drive work, to compare the panel's "Joule heat" with the goal's measure.
//! Usage: cargo run --release -p level --example joule_check -- LEVEL.json [END]

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let text = std::fs::read_to_string(&args[1]).expect("level");
    let level = level::Level::from_json(&text).expect("level json");
    let end: f64 = args.get(2).map_or(24.0, |s| s.parse().expect("end"));
    let mut sim = level.tube_sim_at(
        &level.reference_solution,
        level::tube::TubeResolution::refined(1),
    );
    let mut rows = Vec::new();
    let mut t = 0.0;
    while t < end - 1e-9 {
        t += 1.0;
        sim.advance_to(t);
        let chunk = sim.take_chunk();
        if let Some(f) = chunk.frames.last() {
            rows.push((f.t, f.heat));
            println!("t {:.3} joule {:.4}", f.t, f.heat);
        }
    }
    let rate = |a: f64, b: f64| {
        let pa = rows.iter().min_by(|x, y| (x.0 - a).abs().total_cmp(&(y.0 - a).abs()));
        let pb = rows.iter().min_by(|x, y| (x.0 - b).abs().total_cmp(&(y.0 - b).abs()));
        match (pa, pb) {
            (Some(pa), Some(pb)) if pb.0 > pa.0 => (pb.1 - pa.1) / (pb.0 - pa.0),
            _ => f64::NAN,
        }
    };
    println!("lattice heat per time 4..{end}: {:.3}", rate(4.0, end));
    println!("lattice heat per time 0..4: {:.3}", rate(0.0, 4.0));
    println!("heat goals (drive work): {:?}", sim.results().heat);
}
