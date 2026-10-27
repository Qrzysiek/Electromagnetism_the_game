//! A tube level's simulation built from the level and a placement: the engine's
//! electrodes, emitters, metals, circuits and drives.

use super::*;

impl Level {
    /// A metal's lattice (`lo`..`hi`) with `many` times its preview's carriers (`piece` 0:
    /// the metal's own seeds).
    pub(super) fn metal_lattice(
        &self,
        m: &MetalSpec,
        lo: DVec3,
        hi: DVec3,
        piece: u64,
        many: f64,
    ) -> Lattice {
        let seed = m.seed.wrapping_add(piece.wrapping_mul(0x9e37_79b9));
        let area = (hi.x - lo.x) * (hi.y - lo.y);
        #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
        let cores = libm::round(area * m.core_density()) as usize;
        #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
        let carriers = libm::round(area * m.carriers_per_area * many).max(1.0) as usize;
        Lattice::new(
            lo,
            hi,
            physics::lattice::Material {
                background: -m.charge_per_mass.signum() * m.density,
                core_radius: m.core_radius,
                charge_per_mass: m.charge_per_mass,
                thermal_speed: m.thermal_speed,
            },
            Lattice::random_cores(lo, hi, cores, m.core_radius, seed),
            carriers,
            seed.wrapping_add(METAL_SEED),
        )
    }

    /// A tube run at `refine` = 1 (preview) or 2 (verification), at t = 0.
    ///
    /// # Panics
    /// If the level is not a tube level.
    pub fn tube_sim(&self, player: &[Element], refine: u32) -> TubeSim {
        self.tube_sim_at(player, TubeResolution::refined(refine))
    }

    /// A tube run at a resolution whose parts may be varied separately (convergence
    /// studies: `generator tube-scan --particles --softening`), at t = 0.
    ///
    /// # Panics
    /// If the level is not a tube level.
    pub fn tube_sim_at(&self, player: &[Element], res: TubeResolution) -> TubeSim {
        let spec = self.tube.clone().expect("a tube level");
        let refine = res.refine;
        let r = f64::from(refine);
        let size = PREVIEW_SEGMENT / r;
        // The emitters meshed finer (EMITTER_MESH).
        let shapes = self.tube_prisms(player);
        let sizes: Vec<f64> = (0..shapes.len())
            .map(|e| {
                let emits = e == spec.cathode
                    || spec.ions.is_some_and(|i| i.electrode == e)
                    || spec.metal_of(e).is_some();
                if emits {
                    size * EMITTER_MESH / (res.emitter * spec.emitter_mesh)
                } else {
                    size
                }
            })
            .collect();
        let prisms = Prisms::with_sizes(shapes, &[], &sizes);
        let toward = |a: Option<f64>| {
            a.map(|a| {
                let (s, c) = libm::sincos(a.to_radians());
                DVec3::new(c, s, 0.0)
            })
        };
        // The macroparticles' weights are set once the tube is built (below); their floor.
        let v_ref = self.tube_reference_voltage();
        let floor = |e: usize| {
            let others: f64 = prisms.capacitance[e]
                .iter()
                .enumerate()
                .filter(|(j, _)| *j != e)
                .map(|(_, c)| c.abs())
                .sum();
            WEIGHT_FLOOR * v_ref * others
        };
        // A temperature-limited cathode emits at most its saturation over its surface:
        // the weight also gives `TL_RATE` macroparticles per unit time at that rate (the
        // whole perimeter: an upper bound of the emitting faces, so at least that many).
        let perimeter: f64 = prisms
            .segments
            .iter()
            .enumerate()
            .filter(|(k, _)| prisms.owner_of(*k) == spec.cathode)
            .map(|(_, (a, b))| (*b - *a).length())
            .sum();
        // The metals, and their contacts as emitters of their carriers (at the carriers'
        // weight, toward the metal).
        let many = r * res.particles * spec.particles;
        // A lattice per metal, the player's cuts in it as walls (one that ends on the metal's
        // edge reaching a tenth of a unit past it, so that no carrier slips round its end).
        let mut lattices: Vec<Lattice> = spec
            .metals
            .iter()
            .map(|m| {
                let (lo, hi) = self.metal_rect(m);
                self.metal_lattice(m, lo, hi, 0, many)
            })
            .collect();
        for e in player.iter().filter(|e| e.kind == ElementKind::Cut) {
            let (Some(k), Some((a, b))) = (self.cut_metal(e), self.cut_segment(e)) else {
                continue;
            };
            let (lo, hi) = (lattices[k].lo, lattices[k].hi);
            let d = (b - a) / (b - a).length();
            let tol = 1e-9 * (1.0 + hi.x.abs() + hi.y.abs());
            let edge = |p: DVec3| {
                (p.x - lo.x).abs() < tol
                    || (p.x - hi.x).abs() < tol
                    || (p.y - lo.y).abs() < tol
                    || (p.y - hi.y).abs() < tol
            };
            let a = if edge(a) { a - d * 0.1 } else { a };
            let b = if edge(b) { b + d * 0.1 } else { b };
            lattices[k].walls.push((a, b));
        }
        let lattice_metal: Vec<usize> = (0..lattices.len()).collect();
        let contact = |e: usize| {
            let k = spec.metal_of(e)?;
            let m = &spec.metals[k];
            let first = lattice_metal.iter().position(|&j| j == k)?;
            Some(Emitter {
                electrode: e,
                charge_per_mass: m.charge_per_mass,
                weight: lattices[first].carrier_charge().abs(),
                emit_toward: self.contact_side(m, e),
                saturation: None,
                thermal_speed: m.thermal_speed,
                ohmic: Some(physics::tube::OhmicContact {
                    density: m.density,
                    depth: self.contact_depth(m),
                }),
            })
        };
        let mut floors = vec![floor(spec.cathode)];
        let mut emitters = vec![contact(spec.cathode).unwrap_or(Emitter {
            electrode: spec.cathode,
            charge_per_mass: spec.charge_per_mass,
            weight: 1.0,
            emit_toward: toward(spec.emit_toward_deg),
            saturation: spec.saturation,
            thermal_speed: 0.0,
            ohmic: None,
        })];
        if let Some(ions) = spec.ions {
            floors.push(floor(ions.electrode));
            emitters.push(Emitter {
                electrode: ions.electrode,
                charge_per_mass: ions.charge_per_mass,
                weight: 1.0,
                emit_toward: toward(ions.emit_toward_deg),
                saturation: None,
                thermal_speed: 0.0,
                ohmic: None,
            });
        }
        for m in &spec.metals {
            for &e in &m.contacts {
                if e != spec.cathode && spec.ions.is_none_or(|i| i.electrode != e) {
                    emitters.extend(contact(e));
                }
            }
        }
        // The level's projectiles, then the player's free charges (thrown at t = 0).
        let mut projectiles: Vec<Projectile> = spec
            .projectiles
            .iter()
            .map(|p| {
                let (s, c) = libm::sincos(p.direction_deg.to_radians());
                Projectile {
                    x0: self.grid.position(p.node),
                    v0: DVec3::new(c, s, 0.0) * p.speed,
                    launch: p.launch,
                    charge: p.charge,
                    mass: p.mass,
                }
            })
            .collect();
        let mut detectors: Vec<_> = spec
            .projectiles
            .iter()
            .map(|p| {
                p.detector
                    .map(|[a, b]| (self.grid.position(a), self.grid.position(b)))
            })
            .collect();
        for e in player.iter().filter(|e| e.kind == ElementKind::Free) {
            let (s, c) = libm::sincos(e.angle_deg.to_radians());
            projectiles.push(Projectile {
                x0: self.grid.position(e.node),
                v0: DVec3::new(c, s, 0.0) * e.speed.unwrap_or(0.0),
                launch: 0.0,
                charge: e.value,
                mass: self.limits.free_mass,
            });
            detectors.push(None);
        }
        let n_proj = detectors.len();
        let max = self.grid.position(self.grid.max_node());
        // The solenoid's field: the player's supply on its terminal, else the level's.
        let b_z = spec
            .solenoid
            .and_then(|n| {
                player
                    .iter()
                    .find(|e| e.kind == ElementKind::Supply && e.node == n)
            })
            .map_or(spec.b_z, |e| e.value);
        let mut tube = Tube {
            emitters,
            projectiles,
            arena: Some((DVec3::ZERO, DVec3::new(max.x, max.y, 0.0))),
            b_z,
            b_wave: spec.solenoid_drive.map(|s| physics::tube::SolenoidWave {
                waveform: s.with_level(b_z).waveform(),
                // The solenoid around the board: its axis at the board's centre.
                axis: DVec3::new(0.5 * max.x, 0.5 * max.y, 0.0),
            }),
            circuits: self.tube_circuits(player),
            lattices,
            darwin: spec
                .speed_of_light
                .map(|c| super::self_field_grid(max.x, max.y, c)),
            ..Tube::new(prisms, SOFTENING * size * res.softening, spec.step / r)
        };
        // Each emitter's weight: 1/PREVIEW_PARTICLES of the largest charge its emitting face
        // holds (no particles) over 16 phases of the sources' longest period (t = 0 alone
        // without periodic sources), at least the floor; a temperature-limited cathode's at
        // most what gives TL_RATE macroparticles per unit time.
        let period = tube
            .circuits
            .iter()
            .filter_map(|c| match c.source {
                physics::circuit::Waveform::Sine { omega, .. } if omega > 0.0 => {
                    Some(std::f64::consts::TAU / omega)
                }
                physics::circuit::Waveform::Pulse { period, .. } if period > 0.0 => Some(period),
                _ => None,
            })
            .fold(0.0, f64::max);
        let phases: Vec<f64> = if period > 0.0 {
            (0..16).map(|k| period * f64::from(k) / 16.0).collect()
        } else {
            vec![0.0]
        };
        let samples: Vec<Vec<f64>> = phases
            .iter()
            .map(|&t| {
                let mut v = tube.electrodes.potentials.clone();
                for c in &tube.circuits {
                    if c.resistance.is_finite() {
                        v[c.electrode] = c.source.value(t);
                    }
                }
                v
            })
            .collect();
        for (i, fl) in floors.iter().enumerate() {
            if spec.metal_of(tube.emitters[i].electrode).is_some() {
                // A contact: its metal's carriers' weight.
                continue;
            }
            let q = samples
                .iter()
                .map(|v| tube.emitting_charge(i, v))
                .fold(0.0, f64::max);
            let mut w = q.max(*fl).max(1e-300) / (PREVIEW_PARTICLES * many);
            // (A cold cathode, saturation 0, emits nothing: its weight stays.)
            if i == 0
                && let Some(j) = spec.saturation.filter(|j| *j > 0.0)
            {
                w = w.min(j * perimeter / (spec.particle_rate * many));
            }
            tube.emitters[i].weight = w;
        }
        let state = tube.state();
        let potential_range = samples
            .iter()
            .flatten()
            .fold((0.0_f64, 0.0_f64), |(lo, hi), v| (lo.min(*v), hi.max(*v)));
        let self_b_grid = tube
            .darwin
            .as_ref()
            .map(|_| super::DisplayGrid::board(max.x, max.y, super::SELF_B_SPACING));
        let pending = TubeRun {
            self_b_grid,
            segments: tube.electrodes.segments.clone(),
            owners: (0..tube.electrodes.segments.len())
                .map(|j| tube.electrodes.owner_of(j))
                .collect(),
            potentials: tube.electrodes.potentials.clone(),
            potential_range,
            #[allow(clippy::cast_possible_truncation)]
            metals: tube
                .lattices
                .iter()
                .enumerate()
                .map(|(k, l)| MetalView {
                    lo: l.lo,
                    hi: l.hi,
                    background: l.background,
                    cores: l.cores.iter().map(|c| [c.x as f32, c.y as f32]).collect(),
                    core_radius: l.core_radius,
                    thermal_speed: l.thermal_speed,
                    cells: l.account_size(),
                    walls: l
                        .walls
                        .iter()
                        .map(|(a, b)| [[a.x, a.y], [b.x, b.y]])
                        .collect(),
                    // A lattice per metal (`lattice_metal`).
                    isolated: spec.metals.get(k).is_some_and(|m| m.contacts.is_empty()),
                    contacts: spec
                        .metals
                        .get(k)
                        .map(|m| self.contact_faces(m))
                        .unwrap_or_default(),
                })
                .collect(),
            b_z: tube.b_z,
            b_wave: tube.b_wave,
            leads: self.tube_leads(player),
            ..TubeRun::default()
        };
        let n_goals = spec.potential_goals.len();
        let n_heat = spec.heat_goals.len();
        let n_currents = spec.goals.len();
        // A metal's run keeps a frame every few steps: its particles over the whole run
        // within `FRAME_BUDGET` (other tubes: every step).
        #[allow(
            clippy::cast_possible_truncation,
            clippy::cast_sign_loss,
            clippy::cast_precision_loss
        )]
        let frame_stride = if spec.metals.is_empty() {
            1
        } else {
            let steps = self.physics.t_max / tube.dt;
            (steps * state.x.len() as f64 / FRAME_BUDGET)
                .ceil()
                .max(1.0) as u32
        };
        // Between frames, a sample of the carriers at every step (identities are dealt out
        // in order: every so many of them).
        let mut pending = pending;
        pending.dt = tube.dt;
        pending.fine_every = if frame_stride > 1 {
            u32::try_from((state.x.len() / FINE_SAMPLE).max(1)).unwrap_or(u32::MAX)
        } else {
            0
        };
        TubeSim {
            self_b_grid,
            tube,
            state,
            spec,
            detectors,
            at_start: vec![None; n_currents],
            layout: pending.clone(),
            lattice_metal,
            pending,
            steps: 0,
            frame_stride,
            fine: std::collections::BTreeMap::new(),
            currents: vec![None; n_currents],
            windows: vec![(0.0, 0.0, f64::INFINITY, f64::NEG_INFINITY); n_goals],
            series: vec![Vec::new(); n_goals],
            potential: vec![None; n_goals],
            heat_start: vec![None; n_heat],
            heat: vec![None; n_heat],
            arrivals: vec![None; n_proj],
            closest: vec![f64::INFINITY; n_proj],
        }
    }

    /// The physics circuits of a placement: each tube circuit with its source set by the
    /// player's supply on its electrode, if tunable and placed.
    pub fn tube_circuits(&self, player: &[Element]) -> Vec<physics::tube::Circuit> {
        let Some(spec) = &self.tube else {
            return Vec::new();
        };
        spec.circuits
            .iter()
            .map(|c| {
                let electrode = &self.electrodes[c.electrode];
                let supply = player.iter().find(|e| {
                    e.kind == ElementKind::Supply && electrode.tunable && e.node == electrode.center
                });
                let source = match (supply, c.source) {
                    (
                        Some(s),
                        crate::Source::Sine {
                            amplitude,
                            omega,
                            phase,
                            ..
                        },
                    ) if c.supply_offset => crate::Source::Sine {
                        offset: s.value,
                        amplitude,
                        omega,
                        phase,
                    },
                    (Some(s), source) => source.with_level(s.value),
                    (None, source) => source,
                };
                let primary = physics::tube::Circuit {
                    electrode: c.electrode,
                    source: source.waveform(),
                    resistance: c.resistance.unwrap_or(f64::INFINITY),
                    capacitance: c.capacitance,
                    inductance: c.inductance,
                    to: c.return_to,
                };
                // Linked electrodes follow the same ideal source, or are strapped to the
                // circuit's electrode.
                let ideal =
                    primary.resistance == 0.0 && primary.inductance == 0.0 && primary.to.is_none();
                std::iter::once(primary).chain(c.linked.iter().map(move |&e| {
                    if ideal {
                        physics::tube::Circuit {
                            electrode: e,
                            ..primary
                        }
                    } else {
                        physics::tube::Circuit {
                            electrode: e,
                            source: physics::circuit::Waveform::Dc(0.0),
                            resistance: 0.0,
                            capacitance: 0.0,
                            inductance: 0.0,
                            to: Some(primary.electrode),
                        }
                    }
                }))
            })
            .flat_map(|it| it.collect::<Vec<_>>())
            .collect()
    }

    /// The extremes of every circuit source of a placement (over all time).
    pub(super) fn tube_source_extremes(&self, player: &[Element]) -> (f64, f64) {
        let mut lo: f64 = 0.0;
        let mut hi: f64 = 0.0;
        for c in self.tube_circuits(player) {
            let (a, b) = match c.source {
                physics::circuit::Waveform::Dc(v) => (v, v),
                physics::circuit::Waveform::Sine {
                    offset, amplitude, ..
                } => (offset - amplitude.abs(), offset + amplitude.abs()),
                physics::circuit::Waveform::Pulse { low, high, .. } => {
                    (low.min(high), low.max(high))
                }
            };
            lo = lo.min(a);
            hi = hi.max(b);
        }
        (lo, hi)
    }

    /// The level's largest voltage: of its electrodes' biases and of the sliders' ranges
    /// (supplies, plates); at least 1. Sets the macroparticles' weight floor.
    pub fn tube_reference_voltage(&self) -> f64 {
        let biases = self.electrodes.iter().map(|e| match e.bias {
            crate::ConductorBias::Potential(v) => v.abs(),
            crate::ConductorBias::Grounded | crate::ConductorBias::Charge(_) => 0.0,
        });
        let lists = self
            .limits
            .supply_voltages
            .iter()
            .chain(&self.limits.plate_voltages)
            .map(|v| v.abs());
        let (lo, hi) = self.tube_source_extremes(&[]);
        biases
            .chain(lists)
            .chain([lo.abs(), hi.abs()])
            .fold(1.0, f64::max)
    }
}
