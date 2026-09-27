//! Beams in levels (SPEC §3, Beams; physics in PHYSICS.md §3.3): a shot may be fired as a
//! beam of many particles with spreads in position, direction and energy. All beam
//! shots of a level fly together as one system (one flight per disturbance), so that
//! the particles of different shots (for example two spin states) interact and can be
//! sent to different detectors.
//!
//! The particles are a fixed, deterministic sample of the beam's distribution (a
//! scrambled Halton sequence), so a level's beam is always the same and verifiable. The
//! goal of a beam shot is a verified transmission: at least the required fraction of its
//! particles must arrive, counting only arrivals that are verified.

use physics::DVec3;
use physics::beam::{BeamParticle, BeamRun, BeamScenario, run_beam};
use physics::conductor::Resolution;
use physics::dynamics::Kinematics;
use physics::field::LevelField;
use physics::trajectory::{Outcome, RunSettings};
use physics::verify::{Status, classify};
use serde::{Deserialize, Serialize};

use crate::{Element, Level};

/// Shape of the beam's spreads.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum Distribution {
    /// Gaussian with the given standard deviations, truncated at ±3σ (so that no stray
    /// particle is certain to be lost).
    #[default]
    Gaussian,
    /// Uniform within ± the given spreads.
    Uniform,
}

/// A shot fired as a beam.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct BeamSpec {
    /// Number of particles.
    pub count: u32,
    /// Relative energy spread, of `T₀`.
    pub energy_spread: f64,
    /// Spread of the launch direction, degrees.
    pub angle_spread_deg: f64,
    /// Spread of the launch position across the direction, cells.
    pub width: f64,
    /// Spread along the direction (the bunch length), cells.
    pub length: f64,
    #[serde(default)]
    pub distribution: Distribution,
    /// Fraction of the particles that must arrive, verified.
    pub transmission: f64,
    /// Scrambling of the sample (a different, equally valid sample of the same beam).
    #[serde(default)]
    pub seed: u64,
}

impl Default for BeamSpec {
    fn default() -> Self {
        Self {
            count: 16,
            energy_spread: 0.02,
            angle_spread_deg: 1.0,
            width: 0.2,
            length: 0.2,
            distribution: Distribution::Gaussian,
            transmission: 0.9,
            seed: 0,
        }
    }
}

/// Inverse of the standard normal CDF (P. J. Acklam's rational approximation, relative
/// error below 1.2e-9: the sample defines the beam, so this only has to be deterministic
/// and close to Gaussian).
fn inverse_normal(p: f64) -> f64 {
    const A: [f64; 6] = [
        -3.969_683_028_665_376e1,
        2.209_460_984_245_205e2,
        -2.759_285_104_469_687e2,
        1.383_577_518_672_69e2,
        -3.066_479_806_614_716e1,
        2.506_628_277_459_239,
    ];
    const B: [f64; 5] = [
        -5.447_609_879_822_406e1,
        1.615_858_368_580_409e2,
        -1.556_989_798_598_866e2,
        6.680_131_188_771_972e1,
        -1.328_068_155_288_572e1,
    ];
    const C: [f64; 6] = [
        -7.784_894_002_430_293e-3,
        -3.223_964_580_411_365e-1,
        -2.400_758_277_161_838,
        -2.549_732_539_343_734,
        4.374_664_141_464_968,
        2.938_163_982_698_783,
    ];
    const D: [f64; 4] = [
        7.784_695_709_041_462e-3,
        3.224_671_290_700_398e-1,
        2.445_134_137_142_996,
        3.754_408_661_907_416,
    ];
    let tail = |q: f64| {
        (((((C[0] * q + C[1]) * q + C[2]) * q + C[3]) * q + C[4]) * q + C[5])
            / ((((D[0] * q + D[1]) * q + D[2]) * q + D[3]) * q + 1.0)
    };
    const P_LOW: f64 = 0.024_25;
    if p < P_LOW {
        tail((-2.0 * p.ln()).sqrt())
    } else if p > 1.0 - P_LOW {
        -tail((-2.0 * (1.0 - p).ln()).sqrt())
    } else {
        let q = p - 0.5;
        let r = q * q;
        (((((A[0] * r + A[1]) * r + A[2]) * r + A[3]) * r + A[4]) * r + A[5]) * q
            / (((((B[0] * r + B[1]) * r + B[2]) * r + B[3]) * r + B[4]) * r + 1.0)
    }
}

/// Φ(−3): the lower cut of the truncated Gaussian.
const PHI_MINUS_3: f64 = 0.001_349_898_031_630_094_6;

/// Radical inverse of `i` in base `b` (van der Corput).
fn radical_inverse(mut i: u64, b: u64) -> f64 {
    let (mut f, mut r) = (1.0, 0.0);
    #[allow(clippy::cast_precision_loss)]
    let inv = 1.0 / b as f64;
    while i > 0 {
        f *= inv;
        #[allow(clippy::cast_precision_loss)]
        {
            r += f * (i % b) as f64;
        }
        i /= b;
    }
    r
}

fn splitmix(mut z: u64) -> u64 {
    z = z.wrapping_add(0x9E37_79B9_7F4A_7C15);
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

impl BeamSpec {
    /// The standardized offsets `(across, along, angle, energy)` of particle `i`: points
    /// of the Halton sequence (bases 2, 3, 5, 7), shifted by a seed-dependent rotation,
    /// mapped to the distribution (Gaussian: in σ, within ±3; uniform: within ±1).
    pub fn offsets(&self, i: u32) -> [f64; 4] {
        let mut out = [0.0; 4];
        for (d, base) in [2u64, 3, 5, 7].into_iter().enumerate() {
            #[allow(clippy::cast_precision_loss)]
            let shift = (splitmix(self.seed.wrapping_mul(8).wrapping_add(d as u64)) >> 11) as f64
                / (1u64 << 53) as f64;
            let u = (radical_inverse(u64::from(i) + 1, base) + shift).fract();
            out[d] = match self.distribution {
                Distribution::Gaussian => {
                    inverse_normal(PHI_MINUS_3 + u * (1.0 - 2.0 * PHI_MINUS_3))
                }
                Distribution::Uniform => 2.0 * u - 1.0,
            };
        }
        out
    }

    /// Largest standardized offset (3 for the truncated Gaussian, 1 for uniform).
    pub fn reach(&self) -> f64 {
        match self.distribution {
            Distribution::Gaussian => 3.0,
            Distribution::Uniform => 1.0,
        }
    }
}

/// Verification of a level's beam in one flight (one disturbance).
#[derive(Clone, Debug)]
pub struct BeamVerification {
    /// Shot of each particle.
    pub shot: Vec<usize>,
    pub status: Vec<Status>,
    pub preview: BeamRun,
    pub verified: BeamRun,
}

impl BeamVerification {
    /// Particles of `shot` that arrive, verified, and its particle count.
    pub fn transmitted(&self, shot: usize) -> (usize, usize) {
        let mine = || (0..self.shot.len()).filter(move |&i| self.shot[i] == shot);
        let ok = mine()
            .filter(|&i| {
                self.status[i].is_verified()
                    && self.verified.trajectories[i].outcome == Outcome::Arrived
            })
            .count();
        (ok, mine().count())
    }
}

impl Level {
    /// Whether the level's shots are beams (all of them: mixing beam and single shots
    /// is not supported, `model_issues`).
    pub fn has_beams(&self) -> bool {
        self.shots.iter().any(|s| s.beam.is_some())
    }

    /// The beam flights (one per disturbance) with the field at `resolution`.
    pub fn beam_scenarios(
        &self,
        player: &[Element],
        resolution: Resolution,
    ) -> Vec<BeamScenario<LevelField>> {
        (0..self.flights_per_shot())
            .map(|d| self.beam_scenario(player, resolution, d))
            .collect()
    }

    fn beam_scenario(
        &self,
        player: &[Element],
        resolution: Resolution,
        disturbance: usize,
    ) -> BeamScenario<LevelField> {
        let (field, obstacles) = self.field_at(player, resolution);
        let base = self.scenario_with(0, disturbance, field, obstacles);
        let mut particles = Vec::new();
        // All particles of the level are one sample sequence: shot `s` continues where
        // the previous shots stopped, so particles of different shots never coincide.
        let mut index = 0;
        for (s, shot) in self.shots.iter().enumerate() {
            let one = self.scenario_with(s, disturbance, LevelField::default(), Vec::new());
            let spec = shot.beam.unwrap_or(BeamSpec {
                count: 1,
                ..BeamSpec::default()
            });
            let d = one.p0.normalize();
            let across = DVec3::new(-d.y, d.x, 0.0);
            let kin = Kinematics::new(one.particle.mass, one.c);
            let t0 = shot.launch.kinetic_energy;
            for i in 0..spec.count {
                let [a, l, ang, e] = if shot.beam.is_some() {
                    spec.offsets(index + i)
                } else {
                    [0.0; 4]
                };
                let (sin, cos) = (ang * spec.angle_spread_deg).to_radians().sin_cos();
                let dir = d * cos + across * sin;
                particles.push(BeamParticle {
                    particle: one.particle,
                    x0: one.x0 + across * (a * spec.width) + d * (l * spec.length),
                    p0: kin.momentum_from_kinetic_energy(t0 * (1.0 + e * spec.energy_spread), dir),
                    detector: one.detector,
                    acceptance: one.acceptance,
                });
            }
            index += spec.count;
        }
        BeamScenario {
            field: base.field,
            obstacles: base.obstacles,
            particles,
            c: base.c,
            bounds: base.bounds,
            t_max: base.t_max,
            interact: self.physics.beam_interaction,
        }
    }

    /// Shot of each particle of the beam flights, in order.
    pub fn beam_shots(&self) -> Vec<usize> {
        self.shots
            .iter()
            .enumerate()
            .flat_map(|(s, shot)| std::iter::repeat_n(s, shot.beam.map_or(1, |b| b.count) as usize))
            .collect()
    }

    /// Verifies the beam in every flight: preview at preview resolution and tolerance,
    /// the tighter run with the field at verification resolution; each particle is
    /// classified as a single flight is (PHYSICS.md §7).
    pub fn verify_beams(&self, player: &[Element]) -> Vec<BeamVerification> {
        let tol = self.tolerances();
        let preview = self.beam_scenarios(player, Resolution::Preview);
        let fine = if self.has_metal(player) {
            self.beam_scenarios(player, Resolution::Verify)
        } else {
            preview.clone()
        };
        let shots = self.beam_shots();
        preview
            .iter()
            .zip(&fine)
            .map(|(a, b)| {
                let pa = run_beam(a, &RunSettings::with_tolerance(tol.preview));
                let pb = run_beam(b, &RunSettings::with_tolerance(tol.verify));
                let status = pa
                    .trajectories
                    .iter()
                    .zip(&pb.trajectories)
                    .map(|(x, y)| classify(x, y, a.t_max))
                    .collect();
                BeamVerification {
                    shot: shots.clone(),
                    status,
                    preview: pa,
                    verified: pb,
                }
            })
            .collect()
    }

    /// Whether every beam shot reaches its required verified transmission in every
    /// flight.
    pub fn beams_solved(&self, verifications: &[BeamVerification]) -> bool {
        verifications.iter().all(|v| {
            self.shots.iter().enumerate().all(|(s, shot)| {
                let need = shot.beam.map_or(1.0, |b| b.transmission);
                let (ok, n) = v.transmitted(s);
                #[allow(clippy::cast_precision_loss)]
                let enough = ok as f64 >= need * n as f64 - 1e-9;
                enough
            })
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn inverse_normal_matches_reference_values() {
        // Φ⁻¹(0.975) = 1.959963984540054, Φ⁻¹(Φ(−3)) = −3, symmetric.
        assert!((inverse_normal(0.975) - 1.959_963_984_540_054).abs() < 3e-9);
        assert!((inverse_normal(PHI_MINUS_3) + 3.0).abs() < 1e-8);
        assert!((inverse_normal(0.3) + inverse_normal(0.7)).abs() < 1e-12);
        assert!(inverse_normal(0.5).abs() < 1e-15);
    }

    fn first_bend_beam(interact: bool) -> Level {
        let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../../levels/01_first_bend.json");
        let mut l = Level::from_json(&std::fs::read_to_string(path).unwrap()).unwrap();
        l.physics.c = None;
        l.physics.beam_interaction = interact;
        // Strongly charged particles (allowed for c = ∞: no radiation), so that the
        // interaction matters: charge, mass and energy scaled together keep each
        // particle's path (same q/m and T/q), while the interaction grows with q/Q.
        let k = 5e7;
        l.shots[0].particle.charge *= k;
        l.shots[0].particle.mass *= k;
        l.shots[0].launch.kinetic_energy *= k;
        l.shots[0].beam = Some(BeamSpec {
            count: 12,
            transmission: 0.5,
            ..BeamSpec::default()
        });
        l
    }

    /// A shipped level fired as a beam: the particles are sampled around the launch,
    /// the level is consistent, the interaction changes the flights, and the verified
    /// transmission is counted per shot.
    #[test]
    fn a_level_as_an_interacting_beam() {
        let l = first_bend_beam(true);
        assert_eq!(l.model_issues(), Vec::<String>::new());
        let scn = &l.beam_scenarios(&l.reference_solution, Resolution::Preview)[0];
        assert_eq!(scn.particles.len(), 12);
        let a = l.grid.position(l.shots[0].launch.node);
        let mean: DVec3 = scn.particles.iter().map(|p| p.x0).sum::<DVec3>() / 12.0;
        assert!((mean - a).length() < 0.2);
        let v = l.verify_beams(&l.reference_solution);
        let (ok, n) = v[0].transmitted(0);
        let free = first_bend_beam(false).verify_beams(&l.reference_solution);
        let moved = v[0].verified.trajectories[0].end.x - free[0].verified.trajectories[0].end.x;
        println!(
            "first bend as a beam: {ok}/{n} arrive verified; energy drift {:.1e}; the interaction moves particle 0 by {:.3} cells",
            v[0].verified.energy_max_rel_error,
            moved.length()
        );
        assert!(moved.length() > 1e-3);
        assert!(v[0].verified.energy_max_rel_error < 1e-9);
        assert_eq!(l.beams_solved(&v), ok * 2 >= n);
        // Interacting beams at finite c are rejected.
        let mut c = l.clone();
        c.physics.c = Some(5.0);
        assert!(!c.model_issues().is_empty());
    }

    /// The sample is deterministic, within the truncation, and its moments are close to
    /// those of the distribution.
    #[test]
    fn beam_samples() {
        for distribution in [Distribution::Gaussian, Distribution::Uniform] {
            let spec = BeamSpec {
                count: 512,
                distribution,
                ..BeamSpec::default()
            };
            let pts: Vec<[f64; 4]> = (0..spec.count).map(|i| spec.offsets(i)).collect();
            assert_eq!(
                pts,
                (0..spec.count).map(|i| spec.offsets(i)).collect::<Vec<_>>()
            );
            for d in 0..4 {
                let n = f64::from(spec.count);
                let mean = pts.iter().map(|p| p[d]).sum::<f64>() / n;
                let var = pts.iter().map(|p| (p[d] - mean).powi(2)).sum::<f64>() / n;
                let max = pts.iter().map(|p| p[d].abs()).fold(0.0, f64::max);
                // Truncated at ±3σ: variance 0.9733; uniform on ±1: 1/3.
                let expected = match distribution {
                    Distribution::Gaussian => 0.973_4,
                    Distribution::Uniform => 1.0 / 3.0,
                };
                assert!(mean.abs() < 0.02, "{distribution:?} dim {d}: mean {mean}");
                assert!(
                    (var / expected - 1.0).abs() < 0.03,
                    "{distribution:?} dim {d}: {var}"
                );
                assert!(max <= spec.reach() + 1e-9);
            }
        }
        // Another seed gives another sample of the same beam.
        let a = BeamSpec::default();
        let b = BeamSpec { seed: 1, ..a };
        assert_ne!(
            a.offsets(0).map(f64::to_bits),
            b.offsets(0).map(f64::to_bits)
        );
    }
}
