//! Whether a level stays within the editor's ranges (`check_editable`): every
//! shipped level must be reproducible in the sandbox.

use super::*;

/// Checks that every value of `level` can be produced with the sandbox editors (the same
/// ranges and constraints as the controls above, and the same exhaustive destructuring).
/// Returns the first problem found.
pub fn check_editable(level: &Level) -> Result<(), String> {
    let Level {
        format_version: _,
        engine_version: _,
        name: _,
        description: _,
        // Free text (the remarks' kinds: a combo box of every kind).
        reality: _,
        grid,
        physics,
        shots,
        elements,
        coils,
        limits,
        reference_solution,
        disturbances,
        conductors,
        clouds,
        free_particles,
        electrodes,
        dielectrics,
        dielectric_spheres,
        ferrites,
        gates,
        tube,
    } = level;
    if let Some(t) = tube {
        let TubeSpec {
            cathode,
            charge_per_mass,
            emit_toward_deg,
            saturation,
            ions,
            projectiles,
            b_z,
            solenoid,
            solenoid_drive,
            particle_rate,
            particles,
            emitter_mesh,
            step,
            goals,
            circuits,
            potential_goals,
            heat_goals,
            output,
            metals,
        } = t;
        let n = electrodes.len();
        let any = |v: f64| (-MAX_POSITIVE..=MAX_POSITIVE).contains(&v);
        let nonneg = |v: f64| (0.0..=MAX_POSITIVE).contains(&v);
        let angle = |a: &Option<f64>| a.is_none_or(|a| (-360.0..=360.0).contains(&a));
        let drive_ok = solenoid_drive.is_none_or(|s| s.is_valid());
        let ratio_ok = (MIN_POSITIVE..=MAX_POSITIVE).contains(&charge_per_mass.abs())
            && (1.0..=MAX_POSITIVE).contains(particle_rate)
            && (0.1..=100.0).contains(particles)
            && (1.0..=16.0).contains(emitter_mesh)
            && output.is_none_or(|o| {
                o.circuit < circuits.len() && (MIN_POSITIVE..=MAX_POSITIVE).contains(&o.c)
            });
        let goal_ok = goals.len() <= MAX_COUNT as usize
            && goals.iter().all(|g| {
                let CurrentGoal {
                    electrode,
                    min,
                    max,
                    start,
                    end,
                    // Either value is editable (a combo box); `setup_issues` checks
                    // that an ion goal has an ion source.
                    carriers: _,
                } = *g;
                electrode < n
                    && any(min)
                    && any(max)
                    && nonneg(start)
                    && (MIN_POSITIVE..=MAX_POSITIVE).contains(&end)
            });
        let ions_ok = ions.is_none_or(|i| {
            let IonSource {
                electrode,
                charge_per_mass,
                emit_toward_deg,
            } = i;
            electrode < n
                && (MIN_POSITIVE..=MAX_POSITIVE).contains(&charge_per_mass)
                && angle(&emit_toward_deg)
        });
        let projectiles_ok = projectiles.len() <= MAX_COUNT as usize
            && projectiles.iter().all(|p| {
                let ProjectileSpec {
                    node,
                    direction_deg,
                    speed,
                    launch,
                    charge,
                    mass,
                    detector,
                } = p;
                grid.contains(*node)
                    && any(*direction_deg)
                    && nonneg(*speed)
                    && nonneg(*launch)
                    && any(*charge)
                    && (MIN_POSITIVE..=MAX_POSITIVE).contains(mass)
                    && detector.is_none_or(|[a, b]| grid.contains(a) && grid.contains(b))
            });
        let circuits_ok = circuits.len() <= MAX_COUNT as usize
            && circuits.iter().all(|c| {
                let TubeCircuit {
                    electrode,
                    source,
                    resistance,
                    capacitance,
                    supply_offset: _,
                    linked,
                    inductance,
                    return_to,
                } = c;
                *electrode < n
                    && linked.iter().all(|e| *e < n)
                    && source_editable(source)
                    && resistance.is_none_or(&nonneg)
                    && nonneg(*capacitance)
                    && nonneg(*inductance)
                    && return_to.is_none_or(|e| e < n)
            });
        let potential_ok = potential_goals.len() <= MAX_COUNT as usize
            && potential_goals.iter().all(|g| {
                let PotentialGoal {
                    electrode,
                    reference,
                    measure: _,
                    omega,
                    min,
                    max,
                    start,
                    end,
                } = g;
                *electrode < n
                    && reference.is_none_or(|r| r < n)
                    && nonneg(*omega)
                    && any(*min)
                    && any(*max)
                    && nonneg(*start)
                    && (MIN_POSITIVE..=MAX_POSITIVE).contains(end)
            });
        let heat_ok = heat_goals.len() <= MAX_COUNT as usize
            && heat_goals.iter().all(|g| {
                let HeatGoal {
                    metal,
                    min,
                    max,
                    start,
                    end,
                } = g;
                *metal < metals.len()
                    && any(*min)
                    && any(*max)
                    && nonneg(*start)
                    && (MIN_POSITIVE..=MAX_POSITIVE).contains(end)
            });
        let metals_ok = metals.len() <= MAX_COUNT as usize
            && metals.iter().all(|m| {
                let MetalSpec {
                    corners,
                    contacts,
                    density,
                    charge_per_mass,
                    thermal_speed,
                    mean_free_path,
                    core_radius,
                    carriers_per_area,
                    seed: _,
                    // Free text.
                    label: _,
                } = m;
                let positive = |v: f64| (MIN_POSITIVE..=MAX_POSITIVE).contains(&v);
                corners.iter().all(|c| grid.contains(*c))
                    && contacts.len() <= MAX_COUNT as usize
                    && contacts.iter().all(|e| *e < n)
                    && positive(*density)
                    && positive(charge_per_mass.abs())
                    && positive(*thermal_speed)
                    && positive(*mean_free_path)
                    && positive(*core_radius)
                    && positive(*carriers_per_area)
            });
        let field_ok = circuits_ok
            && metals_ok
            && potential_ok
            && heat_ok
            && solenoid.is_none_or(|n| grid.contains(n))
            && drive_ok
            && any(*b_z)
            && (MIN_POSITIVE..=MAX_POSITIVE).contains(step)
            && saturation.is_none_or(|j| (0.0..=MAX_POSITIVE).contains(&j));
        if *cathode >= n
            || !ratio_ok
            || !goal_ok
            || !angle(emit_toward_deg)
            || !ions_ok
            || !projectiles_ok
            || !field_ok
        {
            return Err("tube value outside the editor's range".into());
        }
    }
    if gates.len() > MAX_COUNT as usize {
        return Err("too many gates".into());
    }
    for g in gates {
        let Detector {
            min,
            max,
            acceptance,
        } = g;
        let acc_ok = acceptance.is_none_or(|a| acceptance_in_range(&a, false));
        if !grid.contains(*min) || !grid.contains(*max) || !acc_ok {
            return Err("gate outside the editor's range".into());
        }
    }
    let Grid {
        nx,
        ny,
        nz,
        subdivision,
    } = grid;
    let fail = |what: &str| Err(what.to_string());
    if !GRID_CELLS.contains(nx) || !GRID_CELLS.contains(ny) {
        return fail("grid size outside the editor's range");
    }
    if *nz != 0 {
        return fail("3D grids cannot be edited yet");
    }
    if !SUBDIVISION.contains(subdivision) {
        return fail("subdivision outside the editor's range");
    }
    let positive = |v: f64| (MIN_POSITIVE..=MAX_POSITIVE).contains(&v);
    let radius = |v: f64| (MIN_POSITIVE..=MAX_RADIUS).contains(&v);
    let WorldPhysics {
        c,
        charge_radius,
        magnet_radius,
        wire_radius,
        antenna_radius,
        rf_omega,
        radiation_reaction: _,
        // Checkboxes.
        beam_interaction: _,
        beam_retarded: _,
        instant_drain: _,
        t_max,
        tolerances: TolerancesSpec { preview, verify },
    } = physics;
    if c.is_some_and(|c| !positive(c))
        || !radius(*charge_radius)
        || !radius(*magnet_radius)
        || !radius(*wire_radius)
        || !radius(*antenna_radius)
        || !(0.0..=MAX_POSITIVE).contains(rf_omega)
        || !positive(*t_max)
        || !TOLERANCE.contains(preview)
        || !TOLERANCE.contains(verify)
    {
        return fail("physics value outside the editor's range");
    }
    let on_grid = |n: &Node| grid.contains(*n);
    for (i, s) in shots.iter().enumerate() {
        let Shot {
            particle:
                ParticleSpec {
                    charge,
                    mass,
                    radius: r,
                    moment,
                },
            launch:
                Launch {
                    node,
                    direction,
                    kinetic_energy,
                    time,
                },
            detector:
                Detector {
                    min,
                    max,
                    acceptance,
                },
            beam,
        } = s;
        let beam_ok = beam.is_none_or(|b| {
            let BeamSpec {
                count,
                energy_spread,
                angle_spread_deg,
                width,
                length,
                distribution: _,
                transmission,
                seed: _,
            } = b;
            (1..=MAX_BEAM).contains(&count)
                && (0.0..=0.3).contains(&energy_spread)
                && (0.0..=45.0).contains(&angle_spread_deg)
                && (0.0..=5.0).contains(&width)
                && (0.0..=5.0).contains(&length)
                && (0.01..=1.0).contains(&transmission)
        });
        if !beam_ok {
            return Err(format!("shot {}: beam outside the editor's range", i + 1));
        }
        if !charge.is_finite()
            || !moment.is_finite()
            || !positive(*mass)
            || !(0.0..=MAX_RADIUS).contains(r)
            || !on_grid(node)
            || !direction.iter().all(|d| d.is_finite())
            || direction[2] != 0.0
            || !positive(*kinetic_energy)
            || !time.is_finite()
            || !on_grid(min)
            || !on_grid(max)
            || !acceptance.is_none_or(|a| acceptance_in_range(&a, true))
        {
            return Err(format!(
                "shot {} has a value outside the editor's range",
                i + 1
            ));
        }
    }
    for (e, reference) in elements
        .iter()
        .map(|e| (e, false))
        .chain(reference_solution.iter().map(|e| (e, true)))
    {
        let Element {
            node,
            kind,
            value,
            angle_deg,
            omega,
            speed,
        } = e;
        match kind {
            ElementKind::Charge | ElementKind::Magnet | ElementKind::Antenna => {}
            // Player-only kinds (the level's own electrodes are `electrodes`, its dynamic
            // particles `free_particles`).
            ElementKind::Plate | ElementKind::Supply | ElementKind::Free | ElementKind::Cut
                if reference => {}
            ElementKind::Plate | ElementKind::Supply | ElementKind::Free | ElementKind::Cut => {
                return fail("plates, power supplies, free charges and cuts are player elements");
            }
        }
        if !on_grid(node)
            || !value.is_finite()
            || !angle_deg.is_finite()
            || omega.is_some_and(|w| !(0.0..=MAX_POSITIVE).contains(&w))
            || speed.is_some_and(|v| !v.is_finite() || v < 0.0)
        {
            return fail("element outside the grid or with an invalid value");
        }
    }
    for c in coils {
        match c {
            Coil::Circle {
                center,
                radius: r,
                kappa,
                rate,
                drive,
                // A checkbox.
                tunable: _,
            } => {
                if !on_grid(center)
                    || !positive(*r)
                    || !kappa.is_finite()
                    || !rate.is_finite()
                    || !drive.as_ref().is_none_or(drive_editable)
                {
                    return fail("circular coil outside the editor's range");
                }
            }
            Coil::Polygon {
                vertices,
                kappa,
                rate,
                drive,
            } => {
                if vertices.len() < 3
                    || !vertices.iter().all(on_grid)
                    || !kappa.is_finite()
                    || !rate.is_finite()
                    || !drive.as_ref().is_none_or(drive_editable)
                {
                    return fail(
                        "polygon coil not reproducible (needs 3 or more vertices on the grid)",
                    );
                }
            }
        }
    }
    if electrodes.len() > MAX_COUNT as usize {
        return fail("too many electrodes");
    }
    for e in electrodes {
        let Electrode {
            center,
            length,
            thickness,
            height,
            angle_deg,
            // A combo box of every shape (`setup_issues` keeps rods and rings to tubes).
            shape: _,
            bias,
            // A checkbox.
            tunable: _,
            drive,
            // Free text.
            label: _,
        } = e;
        let value_ok = match bias {
            ConductorBias::Grounded => true,
            ConductorBias::Charge(v) | ConductorBias::Potential(v) => v.is_finite(),
        };
        let size_ok = [length, thickness, height]
            .iter()
            .all(|v| (MIN_POSITIVE..=MAX_POSITIVE).contains(*v));
        if !on_grid(center)
            || !size_ok
            || !angle_deg.is_finite()
            || !value_ok
            || !drive.as_ref().is_none_or(drive_editable)
        {
            return fail("electrode outside the editor's range");
        }
    }
    if dielectrics.len() > MAX_COUNT as usize {
        return fail("too many dielectrics");
    }
    for d in dielectrics {
        let Dielectric {
            center,
            length,
            thickness,
            height,
            angle_deg,
            permittivity,
        } = d;
        let size_ok = [length, thickness, height]
            .iter()
            .all(|v| (MIN_POSITIVE..=MAX_POSITIVE).contains(*v));
        if !on_grid(center)
            || !size_ok
            || !angle_deg.is_finite()
            || !(1.0..=MAX_PERMITTIVITY).contains(permittivity)
        {
            return fail("dielectric outside the editor's range");
        }
    }
    if ferrites.len() > MAX_COUNT as usize {
        return fail("too many ferrites");
    }
    for f in ferrites {
        let FerritePlate {
            center,
            length,
            thickness,
            height,
            angle_deg,
            permeability,
        } = f;
        let size_ok = [length, thickness, height]
            .iter()
            .all(|v| (MIN_POSITIVE..=MAX_POSITIVE).contains(*v));
        if !on_grid(center)
            || !size_ok
            || !angle_deg.is_finite()
            || !(1.0..=MAX_PERMITTIVITY).contains(permeability)
        {
            return fail("ferrite outside the editor's range");
        }
    }
    if dielectric_spheres.len() > MAX_COUNT as usize {
        return fail("too many dielectric spheres");
    }
    for d in dielectric_spheres {
        let DielectricSphere {
            center,
            radius,
            permittivity,
        } = d;
        if !on_grid(center)
            || !(MIN_POSITIVE..=MAX_POSITIVE).contains(radius)
            || !(1.0..=MAX_PERMITTIVITY).contains(permittivity)
        {
            return fail("dielectric sphere outside the editor's range");
        }
    }
    if conductors.len() > MAX_COUNT as usize {
        return fail("too many metal spheres");
    }
    for c in conductors {
        let Conductor {
            center,
            radius,
            bias,
        } = c;
        let value_ok = match bias {
            ConductorBias::Grounded => true,
            ConductorBias::Charge(v) | ConductorBias::Potential(v) => v.is_finite(),
        };
        if !on_grid(center) || !(MIN_POSITIVE..=MAX_POSITIVE).contains(radius) || !value_ok {
            return fail("metal sphere outside the editor's range");
        }
    }
    if free_particles.len() > MAX_COUNT as usize {
        return fail("too many free particles");
    }
    for f in free_particles {
        let FreeParticle {
            particle:
                ParticleSpec {
                    charge,
                    mass,
                    radius: r,
                    moment,
                },
            node,
            velocity,
            detector,
        } = f;
        let detector_ok = detector.is_none_or(|d| {
            let Detector {
                min,
                max,
                acceptance,
            } = d;
            on_grid(&min)
                && on_grid(&max)
                && acceptance.is_none_or(|a| acceptance_in_range(&a, false))
        });
        if !on_grid(node)
            || !charge.is_finite()
            || !(MIN_POSITIVE..=MAX_POSITIVE).contains(mass)
            || !(0.0..=MAX_RADIUS).contains(r)
            || !moment.is_finite()
            || !velocity.iter().all(|v| v.is_finite())
            || !detector_ok
        {
            return fail("free particle outside the editor's range");
        }
    }
    if clouds.len() > MAX_COUNT as usize {
        return fail("too many charge clouds");
    }
    for c in clouds {
        let Cloud {
            center,
            radius,
            charge,
        } = c;
        if !on_grid(center)
            || !(MIN_POSITIVE..=MAX_POSITIVE).contains(radius)
            || !charge.is_finite()
        {
            return fail("charge cloud outside the editor's range");
        }
    }
    if disturbances.len() > MAX_COUNT as usize {
        return fail("too many disturbances");
    }
    for d in disturbances {
        let Disturbance {
            name: _,
            e,
            bz,
            waves,
        } = d;
        let waves_ok = waves.iter().all(|w| {
            let Wave {
                amplitude,
                direction_deg,
                omega,
                phase_deg,
            } = w;
            amplitude.is_finite()
                && direction_deg.is_finite()
                && (0.0..=MAX_POSITIVE).contains(omega)
                && phase_deg.is_finite()
        });
        if !e.iter().all(|v| v.is_finite()) || !bz.is_finite() || !waves_ok {
            return fail("disturbance outside the editor's range");
        }
    }
    let Limits {
        max_charges,
        magnitudes,
        allow_positive: _,
        allow_negative: _,
        max_magnets,
        magnet_strengths,
        max_antennas,
        antenna_amplitudes,
        antenna_omegas,
        // Any value is editable (a checkbox).
        continuous: _,
        region,
        max_plates,
        plate_voltages,
        plate,
        supply_voltages,
        coil_rates,
        solenoid_fields,
        max_free,
        max_cuts,
        free_charges,
        free_speeds,
        free_mass,
        free_radius,
    } = limits;
    if *max_cuts > MAX_COUNT {
        return fail("cut limit outside the editor's range");
    }
    if *max_free > MAX_COUNT
        || !free_charges.iter().all(|v| v.is_finite())
        || !free_speeds.iter().all(|v| v.is_finite() && *v >= 0.0)
        || !(MIN_POSITIVE..=MAX_POSITIVE).contains(free_mass)
        || !(0.0..=MAX_RADIUS).contains(free_radius)
    {
        return fail("free-charge limits outside the editor's range");
    }
    let plate_ok = [plate.length, plate.thickness, plate.height]
        .iter()
        .all(|v| (MIN_POSITIVE..=MAX_POSITIVE).contains(v));
    if !plate_ok
        || !plate_voltages
            .iter()
            .chain(supply_voltages)
            .chain(coil_rates)
            .chain(solenoid_fields)
            .all(|v| v.is_finite())
    {
        return fail("plate limits outside the editor's range");
    }
    if *max_charges > MAX_COUNT
        || *max_magnets > MAX_COUNT
        || *max_antennas > MAX_COUNT
        || *max_plates > MAX_COUNT
        || !magnitudes
            .iter()
            .chain(magnet_strengths)
            .chain(antenna_amplitudes)
            .chain(antenna_omegas)
            .all(|v| *v > 0.0 && v.is_finite())
    {
        return fail("limits outside the editor's range");
    }
    if let Some(Region2 { min, max }) = region
        && !(on_grid(min) && on_grid(max))
    {
        return fail("player region outside the grid");
    }
    Ok(())
}
