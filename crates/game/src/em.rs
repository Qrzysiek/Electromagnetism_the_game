//! Tube levels: the electromagnetic field of the preview (`level::tube_em`, PHYSICS.md
//! §2.11), computed on a thread of its own as the animation plays, for the magnetic and
//! energy-flow maps. The thread keeps the field at the latest time asked for (rebuilding
//! it from the start when asked for an earlier one) and sends snapshots on a display grid.

use std::sync::mpsc::{Receiver, Sender, TryRecvError, channel};
use std::sync::{Arc, Mutex};

use level::Level;
use level::tube::TubeRun;
use level::tube_em::TubeEm;

/// Display samples per board cell.
const SAMPLES_PER_CELL: f64 = 3.0;
/// The thread advances in slices of this much time, looking for new requests between and
/// showing its progress. (Slices of 0.5 made a field slower than the playback jump half a
/// time unit at a time, the playback held back with it: Klystron's half a period.)
const SLICE: f64 = 0.02;

/// What the game asks of the field's thread.
pub enum EmRequest {
    /// A new run (its frames so far).
    Reset {
        revision: u64,
        level: Arc<Level>,
        run: Box<TubeRun>,
    },
    /// More of the run's frames.
    Frames { revision: u64, chunk: Arc<TubeRun> },
    /// The field at time `t`.
    Want { revision: u64, t: f64 },
}

/// The field at one time on a display grid.
#[derive(Clone, Debug)]
pub struct EmSnapshot {
    pub revision: u64,
    pub t: f64,
    pub origin: (f64, f64),
    pub spacing: f64,
    pub nx: usize,
    pub ny: usize,
    /// Per sample (row by row): `B_z` and the Poynting vector
    /// `S = (c/4π) E × B_z ẑ = (c/4π)(E_y B_z, −E_x B_z)`.
    pub samples: Vec<[f32; 3]>,
    /// The display's speed of light.
    pub c: f64,
    /// The power per unit length radiated out of the board (the energy flow through a
    /// contour inside its edge), averaged over the last two time units.
    pub radiated: f64,
    /// The field has reached the end of the run's frames so far (it can go no further).
    pub at_end: bool,
    /// A closed box's smoothing's mean delay (`TubeEm::delay`; 0: the field as it is).
    pub delay: f64,
}

/// The field's thread.
pub struct EmWorker {
    tx: Sender<EmRequest>,
    /// (A receiver is not `Sync`; the game's state must be.)
    rx: Mutex<Receiver<EmSnapshot>>,
}

impl EmWorker {
    pub fn start() -> Self {
        let (tx, requests) = channel::<EmRequest>();
        let (out, rx) = channel::<EmSnapshot>();
        // A thread pool of its own, half the cores: in the shared one each of the field's
        // small parallel steps waited behind the preview's and verification's long tasks,
        // and the field (and the playback that waits for it) crawled while they ran.
        let threads = std::thread::available_parallelism().map_or(4, |n| (n.get() / 2).clamp(2, 8));
        std::thread::Builder::new()
            .name("tube field".into())
            .spawn(move || {
                match rayon::ThreadPoolBuilder::new()
                    .num_threads(threads)
                    .thread_name(|i| format!("tube field {i}"))
                    .build()
                {
                    Ok(pool) => pool.install(move || run(&requests, &out)),
                    Err(_) => run(&requests, &out),
                }
            })
            .expect("the field's thread");
        Self {
            tx,
            rx: Mutex::new(rx),
        }
    }

    pub fn send(&self, r: EmRequest) {
        let _ = self.tx.send(r);
    }

    /// The latest snapshot sent since the last call, if any.
    pub fn latest(&self) -> Option<EmSnapshot> {
        let mut last = None;
        if let Ok(rx) = self.rx.lock() {
            while let Ok(s) = rx.try_recv() {
                last = Some(s);
            }
        }
        last
    }
}

struct State {
    revision: u64,
    level: Arc<Level>,
    run: TubeRun,
    em: Option<TubeEm>,
}

#[allow(clippy::needless_pass_by_value)]
fn apply(state: &mut Option<State>, want: &mut Option<(u64, f64)>, r: EmRequest) {
    match r {
        EmRequest::Reset {
            revision,
            level,
            run,
        } => {
            *state = Some(State {
                revision,
                level,
                run: *run,
                em: None,
            });
        }
        EmRequest::Frames { revision, chunk } => {
            if let Some(s) = state.as_mut().filter(|s| s.revision == revision) {
                s.run.extend((*chunk).clone());
            }
        }
        EmRequest::Want { revision, t } => *want = Some((revision, t)),
    }
}

/// How often a catch-up shows its progress (a display frame's time).
const PROGRESS: std::time::Duration = std::time::Duration::from_millis(16);

fn run(requests: &Receiver<EmRequest>, out: &Sender<EmSnapshot>) {
    let mut state: Option<State> = None;
    let mut want: Option<(u64, f64)> = None;
    // Work left from an interrupted request (no need to wait for another).
    let mut dirty = false;
    // When the last snapshot was sent.
    let mut shown = std::time::Instant::now();
    loop {
        if !dirty {
            match requests.recv() {
                Ok(r) => apply(&mut state, &mut want, r),
                Err(_) => return,
            }
        }
        dirty = false;
        loop {
            match requests.try_recv() {
                Ok(r) => apply(&mut state, &mut want, r),
                Err(TryRecvError::Empty) => break,
                Err(TryRecvError::Disconnected) => return,
            }
        }
        let Some(s) = state.as_mut() else { continue };
        let Some((revision, t)) = want.filter(|(r, _)| *r == s.revision) else {
            continue;
        };
        if s.run.frames.len() < 2 {
            continue;
        }
        if s.em.as_ref().is_none_or(|em| em.shown().t > t + 1e-9) {
            s.em = Some(TubeEm::new(&s.level, &s.run));
        }
        let em = s.em.as_mut().expect("the field");
        // The field shown at t is the field's own at t + its smoothing's delay.
        let target = t + em.delay();
        let last = s.run.frames.last().map_or(0.0, |f| f.t);
        // Advance in slices, leaving the work to the next round if a request arrives (one
        // comes with every displayed frame while the animation plays); while catching up,
        // the field shows its progress every `PROGRESS` of wall time, and at once when it
        // reaches the time asked for.
        let mut next = None;
        while em.fdtd.t + em.fdtd.dt <= target {
            let before = em.fdtd.t;
            em.advance_to(&s.run, (em.fdtd.t + SLICE).min(target));
            if em.fdtd.t <= before {
                break; // the frames end here
            }
            if let Ok(r) = requests.try_recv() {
                next = Some(r);
                break;
            }
            if shown.elapsed() > PROGRESS {
                let _ = out.send(snapshot(revision, Some(em), last));
                shown = std::time::Instant::now();
            }
        }
        if next.is_none() || shown.elapsed() > PROGRESS {
            let _ = out.send(snapshot(revision, Some(em), last));
            shown = std::time::Instant::now();
        }
        if let Some(r) = next {
            apply(&mut state, &mut want, r);
            dirty = true;
        }
    }
}

#[allow(
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::cast_precision_loss
)]
fn snapshot(revision: u64, em: Option<&TubeEm>, last: f64) -> EmSnapshot {
    let Some(em) = em else {
        return EmSnapshot {
            revision,
            t: 0.0,
            origin: (0.0, 0.0),
            spacing: 1.0,
            nx: 0,
            ny: 0,
            samples: Vec::new(),
            c: 1.0,
            radiated: 0.0,
            at_end: false,
            delay: 0.0,
        };
    };
    // The field the maps show (smoothed in time).
    let field = em.shown();
    let g = em.fdtd.grid;
    let spacing = 1.0 / SAMPLES_PER_CELL;
    let (w, hgt) = (g.nx as f64 * g.h, g.ny as f64 * g.h);
    let (nx, ny) = ((w / spacing) as usize, (hgt / spacing) as usize);
    let k = em.fdtd.c / (4.0 * std::f64::consts::PI);
    // At the cell centres: B_z and S, then smoothed over a board cell (the macroparticles'
    // graininess: each one's current speckles the field within a cell or two).
    let (gx, gy) = (g.nx, g.ny);
    let mut b = vec![0.0; gx * gy];
    let mut sx = vec![0.0; gx * gy];
    let mut sy = vec![0.0; gx * gy];
    for j in 0..gy {
        for i in 0..gx {
            let c = i + gx * j;
            let bz = field.bz[c];
            let ex = 0.5 * (field.ex[i + gx * j] + field.ex[i + gx * (j + 1)]);
            let ey = 0.5 * (field.ey[i + (gx + 1) * j] + field.ey[i + 1 + (gx + 1) * j]);
            b[c] = bz;
            sx[c] = k * ey * bz;
            sy[c] = -k * ex * bz;
        }
    }
    let half = (0.5 / g.h).round() as usize;
    for f in [&mut b, &mut sx, &mut sy] {
        box_filter(f, gx, gy, half);
    }
    let at = |f: &[f64], x: f64, y: f64| {
        let u = ((x - g.origin.0) / g.h - 0.5).clamp(0.0, (gx - 1) as f64 - 1e-9);
        let v = ((y - g.origin.1) / g.h - 0.5).clamp(0.0, (gy - 1) as f64 - 1e-9);
        let (i, j) = (u as usize, v as usize);
        let (fu, fv) = (u - i as f64, v - j as f64);
        let c = |i: usize, j: usize| f[i.min(gx - 1) + gx * j.min(gy - 1)];
        (1.0 - fu) * (1.0 - fv) * c(i, j)
            + fu * (1.0 - fv) * c(i + 1, j)
            + (1.0 - fu) * fv * c(i, j + 1)
            + fu * fv * c(i + 1, j + 1)
    };
    let samples = (0..nx * ny)
        .map(|n| {
            let (x, y) = (
                g.origin.0 + ((n % nx) as f64 + 0.5) * spacing,
                g.origin.1 + ((n / nx) as f64 + 0.5) * spacing,
            );
            [
                at(&b, x, y) as f32,
                at(&sx, x, y) as f32,
                at(&sy, x, y) as f32,
            ]
        })
        .collect();
    EmSnapshot {
        revision,
        t: field.t,
        origin: g.origin,
        spacing,
        nx,
        ny,
        samples,
        c: em.fdtd.c,
        radiated: em.mean_outflow(2.0),
        at_end: em.fdtd.t + em.fdtd.dt > last,
        delay: em.delay(),
    }
}

/// A box filter of half-width `half` cells along each axis of a `gx × gy` field (running
/// sums; the edges average what lies inside).
#[allow(clippy::cast_precision_loss)]
fn box_filter(f: &mut [f64], gx: usize, gy: usize, half: usize) {
    if half == 0 {
        return;
    }
    let mut line = Vec::new();
    for j in 0..gy {
        line.clear();
        line.extend_from_slice(&f[gx * j..gx * (j + 1)]);
        for i in 0..gx {
            let (a, b) = (i.saturating_sub(half), (i + half).min(gx - 1));
            f[i + gx * j] = line[a..=b].iter().sum::<f64>() / (b - a + 1) as f64;
        }
    }
    for i in 0..gx {
        line.clear();
        line.extend((0..gy).map(|j| f[i + gx * j]));
        for j in 0..gy {
            let (a, b) = (j.saturating_sub(half), (j + half).min(gy - 1));
            f[i + gx * j] = line[a..=b].iter().sum::<f64>() / (b - a + 1) as f64;
        }
    }
}
