//! A tube level run and judged: the runs' ends, the runs, their results, records
//! and verdicts, and the objective the generator searches.

use super::*;

impl Level {
    /// The end of a tube level's verification runs: the current goal's window end, and
    /// with projectile detectors the level's `t_max` (they must arrive by then).
    pub fn tube_verify_end(&self) -> f64 {
        self.tube.as_ref().map_or(0.0, |t| {
            let goal = t
                .goals
                .iter()
                .map(|g| g.end)
                .chain(t.potential_goals.iter().map(|g| g.end))
                .chain(t.heat_goals.iter().map(|g| g.end))
                .fold(0.0, f64::max);
            if t.projectiles.iter().any(|p| p.detector.is_some()) {
                goal.max(self.physics.t_max)
            } else {
                goal
            }
        })
    }

    /// The end of a tube level's display run: its `t_max`, at least the verification's.
    pub fn tube_display_end(&self) -> f64 {
        self.tube_verify_end().max(self.physics.t_max)
    }

    /// Runs the tube at `refine` to `t_end` (whole), as one run.
    ///
    /// # Panics
    /// If the level is not a tube level.
    pub fn tube_run(&self, player: &[Element], refine: u32, t_end: f64) -> TubeRun {
        let mut sim = self.tube_sim(player, refine);
        sim.advance_to(t_end);
        sim.take_chunk()
    }

    /// The results a verdict needs from a run at `refine`: the goal's current and the
    /// projectiles' arrivals, run to `tube_verify_end`.
    ///
    /// # Panics
    /// If the level is not a tube level.
    pub fn tube_results(&self, player: &[Element], refine: u32) -> TubeResults {
        self.tube_record(player, TubeResolution::refined(refine))
            .results
    }

    /// One run at `res` to `tube_verify_end`: its results, its projectiles' closest
    /// approaches, its size and its time; from the offline tools' run cache when one is on
    /// (`run_cache`: hash-locked to the level, the placement, the resolution and the engine).
    ///
    /// # Panics
    /// If the level is not a tube level.
    pub fn tube_record(
        &self,
        player: &[Element],
        res: TubeResolution,
    ) -> crate::run_cache::TubeRecord {
        let end = self.tube_verify_end();
        let compute = || {
            let started = std::time::Instant::now();
            let mut sim = self.tube_sim_at(player, res);
            sim.advance_to(end);
            crate::run_cache::TubeRecord {
                results: sim.results(),
                closest: sim.closest().to_vec(),
                segments: sim.segment_count(),
                steps: sim.steps(),
                seconds: started.elapsed().as_secs_f64(),
            }
        };
        match crate::run_cache::get() {
            Some(cache) => {
                let key = cache.tube_key("tube-record", self, player, res, end, &[]);
                cache.tube_record(&key, compute)
            }
            None => compute(),
        }
    }

    /// The verdict of a placement: the preview and the verification run.
    ///
    /// # Panics
    /// If the level is not a tube level.
    pub fn tube_verdict(&self, player: &[Element]) -> TubeVerdict {
        let spec = self.tube.as_ref().expect("a tube level");
        TubeVerdict::new(
            spec,
            &self.tube_results(player, 1),
            &self.tube_results(player, 2),
        )
    }

    /// The search objective (`solve::objective`) of a tube level, from the preview: the
    /// current's and each potential goal's distance outside its range in units of its
    /// width (0 inside), plus each projectile's closest approach to its detector (cells),
    /// and `Arrived` when all are 0, else `Rejected` (as a particle in its detector outside
    /// its acceptance).
    ///
    /// # Panics
    /// If the level is not a tube level.
    pub fn tube_objective(&self, player: &[Element]) -> (f64, physics::trajectory::Outcome) {
        use physics::trajectory::Outcome;
        let spec = self.tube.as_ref().expect("a tube level");
        let run = self.tube_record(player, TubeResolution::refined(1));
        let mut score = 0.0;
        for (g, i) in spec.goals.iter().zip(&run.results.currents) {
            let i = i.unwrap_or(0.0);
            score += (g.min - i).max(i - g.max).max(0.0) / (g.max - g.min).max(1e-300);
        }
        let measured = &run.results.potential;
        for (k, g) in spec.potential_goals.iter().enumerate() {
            // Not measured (the run ended early): one width, finite for the annealing.
            score += measured.get(k).copied().flatten().map_or(1.0, |v| {
                let v = g.compared(v);
                (g.min - v).max(v - g.max).max(0.0) / (g.max - g.min).max(1e-300)
            });
        }
        let heat = &run.results.heat;
        for (k, g) in spec.heat_goals.iter().enumerate() {
            score += heat.get(k).copied().flatten().map_or(1.0, |v| {
                (g.min - v).max(v - g.max).max(0.0) / (g.max - g.min).max(1e-300)
            });
        }
        score += run.closest.iter().filter(|d| d.is_finite()).sum::<f64>();
        if score == 0.0 {
            (0.0, Outcome::Arrived)
        } else {
            (score, Outcome::Rejected)
        }
    }
}
