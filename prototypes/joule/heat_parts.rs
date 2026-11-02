//! Scratch (not committed): a tube level's reference run at a given resolution, its metal's
//! energy accounts as rates over [t0, end]: the lattice's heat, the field's work, the drive's,
//! the induced fields' (with the self-field's part), the self-field part by carrier kind.
//! Usage: heat_parts LEVEL.json END REFINE PARTICLES SOFTENING

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let text = std::fs::read_to_string(&args[1]).expect("level");
    let level = level::Level::from_json(&text).expect("level json");
    let end: f64 = args[2].parse().expect("end");
    let refine: u32 = args[3].parse().expect("refine");
    let particles: f64 = args[4].parse().expect("particles");
    let softening: f64 = args[5].parse().expect("softening");
    let res = level::tube::TubeResolution {
        refine,
        particles,
        softening,
        emitter: 1.0,
    };
    let mut sim = level.tube_sim_at(&level.reference_solution, res);
    let t0 = 4.0;
    sim.advance_to(t0);
    let a = sim.debug_state();
    sim.advance_to(end);
    let b = sim.debug_state();
    let span = b.0 - a.0;
    let names = [
        "lattice heat",
        "field work",
        "drive work",
        "inductive",
        "self free",
        "self jumped",
        "self born",
        "jumped count",
    ];
    print!("{} r{refine} p{particles} s{softening}:", args[1]);
    for (k, n) in names.iter().enumerate() {
        print!(" {n} {:.2};", (b.1[k] - a.1[k]) / span);
    }
    println!(" carriers {}", sim.in_flight());
}
