//! The physics worker kept in step with the placement: what it is sent, what
//! comes back.

use super::*;

/// Current setup revision (changes with every edit and level switch).
pub(super) fn revision(game: &Game) -> u64 {
    game.editor.revision + ((game.level_index as u64) << 48)
}

/// Sends changed setups to the physics thread.
pub(super) fn sync_physics(mut game: ResMut<Game>, worker: Res<PhysicsWorker>) {
    let revision = revision(&game);
    if revision == game.sent_revision {
        return;
    }
    game.sent_revision = revision;
    let level = &game.editor.level;
    let n = level.shots.len();
    let flights = level.flight_count();
    let per_shot = level.flights_per_shot();
    let request_level = level.clone();
    // Keep old previews on screen until new ones arrive; verdicts are recomputed.
    game.flights.resize(flights, FlightView::default());
    for s in &mut game.flights {
        s.verdict = None;
    }
    if request_level.has_beams() {
        game.beams.resize(per_shot, BeamView::default());
        for b in &mut game.beams {
            b.verified = None;
            // The old setup's best flight stays on screen (faded) until the new preview
            // arrives.
            if let Some(f) = b.verdict_flight.take() {
                b.preview = Some(f);
            }
        }
    } else {
        game.beams.clear();
    }
    if game.active_shot >= n {
        game.active_shot = 0;
    }
    if game.active_disturbance >= per_shot {
        game.active_disturbance = 0;
    }
    let placement = game.editor.placement.clone();
    worker.0.submit(Request {
        revision,
        level: request_level,
        placement,
    });
}

pub(super) fn poll_physics(mut game: ResMut<Game>, worker: Res<PhysicsWorker>) {
    sandbox::poll(&mut game);
    let current = game.sent_revision;
    for r in worker.0.poll() {
        match r {
            Response::Preview {
                revision,
                flight: shot,
                preview,
            } if revision == current => {
                if let Some(s) = game.flights.get_mut(shot) {
                    s.preview = Some(preview);
                    s.preview_revision = revision;
                }
            }
            Response::Verified {
                revision,
                flight: shot,
                status,
                outcome,
            } if revision == current => {
                if let Some(s) = game.flights.get_mut(shot) {
                    s.verdict = Some((status, outcome));
                }
            }
            Response::BeamPreview {
                revision,
                flight,
                preview,
            } if revision == current => {
                if let Some(b) = game.beams.get_mut(flight) {
                    b.preview = Some(preview);
                    b.preview_revision = revision;
                    b.verdict_flight = None;
                }
            }
            Response::BeamVerified {
                revision,
                flight,
                results,
                run,
            } if revision == current => {
                if let Some(b) = game.beams.get_mut(flight) {
                    b.verified = Some(results);
                    b.verdict_flight = run;
                }
            }
            Response::Cost { revision, cost } if revision == current => {
                game.cost = Some((revision, cost));
            }
            Response::Circuit { revision, view } if revision == current => {
                game.circuit = Some((revision, view));
            }
            Response::Invalid { revision, issues } if revision == current => {
                game.setup_issues = Some((revision, issues));
            }
            Response::TubeChunk {
                revision,
                chunk,
                done,
            } if revision == current => {
                // The field's thread follows the same run.
                if game.em_run == Some(revision)
                    && let Some(w) = &game.em
                {
                    w.send(em::EmRequest::Frames {
                        revision,
                        chunk: chunk.clone(),
                    });
                }
                let chunk = std::sync::Arc::try_unwrap(chunk).unwrap_or_else(|a| (*a).clone());
                match &mut game.tube_preview {
                    Some(v) if v.revision == revision => {
                        v.run.extend(chunk);
                        v.done = done;
                    }
                    _ => {
                        game.tube_preview = Some(TubeView {
                            revision,
                            run: chunk,
                            done,
                        });
                    }
                }
            }
            Response::TubeVerified { revision, verdict } if revision == current => {
                game.tube_verdict = Some((revision, verdict));
            }
            _ => {}
        }
    }
}
