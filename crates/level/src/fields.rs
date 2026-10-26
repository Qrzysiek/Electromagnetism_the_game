//! The level's static field: its electrodes, plates, dielectrics and ferrites as
//! the solvers' bodies, and a placement's field and obstacles.

use super::*;

/// Box electrodes with the surface charge the fixed charges induce on them.
pub(super) fn electrodes_for(
    boxes: Vec<physics::bem::BoxElectrode>,
    dielectrics: Vec<physics::bem::Dielectric>,
    charges: &[FixedCharge],
    resolution: Resolution,
) -> physics::bem::Electrodes {
    if boxes.is_empty() && dielectrics.is_empty() {
        return physics::bem::Electrodes::default();
    }
    let sources: Vec<(DVec3, f64)> = charges.iter().map(|c| (c.position, c.charge)).collect();
    let res = match resolution {
        Resolution::Preview => physics::bem::Resolution::Preview,
        Resolution::Verify => physics::bem::Resolution::Verify,
        Resolution::Display => physics::bem::Resolution::Display,
    };
    physics::bem::Electrodes::new_with_dielectrics(boxes, dielectrics, &sources, res)
}

impl Level {
    /// The box of a player plate.
    pub fn plate_box(&self, e: &Element) -> physics::bem::BoxElectrode {
        let s = self.limits.plate;
        physics::bem::BoxElectrode {
            center: self.grid.position(e.node),
            angle: e.angle_deg.to_radians(),
            half_length: s.length / 2.0,
            half_thickness: s.thickness / 2.0,
            half_height: s.height / 2.0,
            bias: if e.value == 0.0 {
                Bias::Grounded
            } else {
                Bias::Potential(e.value)
            },
        }
    }

    /// Every electrode of a placement: the level's (at the potentials of the player's
    /// power supplies where set), then the player's plates.
    pub fn all_box_electrodes(&self, player: &[Element]) -> Vec<physics::bem::BoxElectrode> {
        let mut boxes = self.box_electrodes();
        for (b, e) in boxes.iter_mut().zip(&self.electrodes) {
            if let Some(s) = self.supply_of(e, player)
                && e.drive.is_none()
            {
                b.bias = Bias::Potential(s.value);
            }
        }
        boxes.extend(
            player
                .iter()
                .filter(|e| e.kind == ElementKind::Plate)
                .map(|e| self.plate_box(e)),
        );
        boxes
    }

    /// Whether a placement has metal (spheres, electrodes or player plates), whose model
    /// differs between preview and verification resolution.
    pub fn has_metal(&self, player: &[Element]) -> bool {
        !self.conductors.is_empty()
            || !self.electrodes.is_empty()
            || !self.dielectrics.is_empty()
            || !self.dielectric_spheres.is_empty()
            || !self.ferrites.is_empty()
            || player.iter().any(|e| e.kind == ElementKind::Plate)
    }

    /// The static field and the obstacles for a placement at preview resolution.
    pub fn field(&self, player: &[Element]) -> (LevelField, Vec<Shape>) {
        self.field_at(player, Resolution::Preview)
    }

    /// The field and the obstacles for a placement: level elements first, then player
    /// elements, in the given order (the summation order of the field). Metal spheres
    /// are computed at the given resolution (PHYSICS.md §2.6).
    pub fn field_at(&self, player: &[Element], resolution: Resolution) -> (LevelField, Vec<Shape>) {
        let all: Vec<&Element> = self.elements.iter().chain(player).collect();
        let charges: Vec<FixedCharge> = all
            .iter()
            .filter(|e| e.kind == ElementKind::Charge)
            .map(|e| FixedCharge {
                position: self.grid.position(e.node),
                charge: e.value,
                radius: self.physics.charge_radius,
            })
            .collect();
        let dipoles: Vec<MagneticDipole> = all
            .iter()
            .filter(|e| e.kind == ElementKind::Magnet)
            .map(|e| MagneticDipole {
                position: self.grid.position(e.node),
                moment: DVec3::new(0.0, 0.0, e.value),
                radius: self.physics.magnet_radius,
            })
            .collect();
        let antennas: Vec<OscillatingDipole> = all
            .iter()
            .filter(|e| e.kind == ElementKind::Antenna)
            .map(|e| {
                let a = e.angle_deg.to_radians();
                OscillatingDipole {
                    position: self.grid.position(e.node),
                    amplitude: DVec3::new(libm::cos(a), libm::sin(a), 0.0) * e.value,
                    omega: e.omega.unwrap_or(self.physics.rf_omega),
                    phase: 0.0,
                    c: self.c(),
                    radius: self.physics.antenna_radius,
                }
            })
            .collect();
        let wire = self.physics.wire_radius;
        let mut loops = Vec::new();
        let mut polygons = Vec::new();
        let mut obstacles: Vec<Shape> = charges
            .iter()
            .map(|c| {
                Shape::Sphere(Sphere {
                    center: c.position,
                    radius: c.radius,
                })
            })
            .collect();
        obstacles.extend(dipoles.iter().map(|d| {
            Shape::Sphere(Sphere {
                center: d.position,
                radius: d.radius,
            })
        }));
        obstacles.extend(antennas.iter().map(|a| {
            Shape::Sphere(Sphere {
                center: a.position,
                radius: a.radius,
            })
        }));
        for coil in &self.coils {
            match coil {
                Coil::Circle {
                    center,
                    radius,
                    kappa,
                    rate,
                    drive,
                    tunable: _,
                } => {
                    // A tunable coil's supply sets its ramp rate (a driven one's: its
                    // source's level, `drives_for`).
                    let supplied = coil.supply_node().and_then(|n| self.supply_at(n, player));
                    let rate = match supplied {
                        Some(s) if drive.is_none() => s.value,
                        _ => *rate,
                    };
                    let l = CircularLoop {
                        center: self.grid.position(*center),
                        normal: DVec3::Z,
                        radius: *radius,
                        kappa: *kappa,
                        wire_radius: wire,
                        rate,
                    };
                    obstacles.push(Shape::Torus(Torus {
                        center: l.center,
                        normal: l.normal,
                        major: l.radius,
                        minor: wire,
                    }));
                    loops.push(l);
                }
                Coil::Polygon {
                    vertices,
                    kappa,
                    rate,
                    drive: _,
                } => {
                    let v: Vec<DVec3> = vertices.iter().map(|n| self.grid.position(*n)).collect();
                    for i in 0..v.len() {
                        obstacles.push(Shape::Capsule(Capsule {
                            a: v[i],
                            b: v[(i + 1) % v.len()],
                            radius: wire,
                        }));
                    }
                    polygons.push(PolygonCoil {
                        vertices: v,
                        kappa: *kappa,
                        wire_radius: wire,
                        rate: *rate,
                    });
                }
            }
        }
        // Electrode boxes (the level's, then the player's plates) with the contact shell.
        let boxes = self.all_box_electrodes(player);
        obstacles.extend(
            physics::bem::Electrodes::shapes_only(boxes.clone()).obstacles(CONTACT_DISTANCE),
        );
        obstacles.extend(
            physics::bem::Electrodes::shapes_only(self.body_boxes()).obstacles(CONTACT_DISTANCE),
        );
        for d in &self.dielectric_spheres {
            obstacles.push(Shape::Sphere(Sphere {
                center: self.grid.position(d.center),
                radius: d.radius + CONTACT_DISTANCE,
            }));
        }
        let dielectrics = self.physics_dielectrics();
        for c in &self.conductors {
            obstacles.push(Shape::Sphere(Sphere {
                center: self.grid.position(c.center),
                radius: c.radius + CONTACT_DISTANCE,
            }));
        }
        let clouds: Vec<ChargeCloud> = self
            .clouds
            .iter()
            .map(|c| ChargeCloud {
                position: self.grid.position(c.center),
                charge: c.charge,
                radius: c.radius,
            })
            .collect();
        // Metal and electrodes see the clouds as point charges: the validation keeps
        // every cloud clear of them, where its field is a point charge's.
        let sources: Vec<FixedCharge> = charges
            .iter()
            .copied()
            .chain(clouds.iter().map(|c| FixedCharge {
                position: c.position,
                charge: c.charge,
                radius: 0.0,
            }))
            .collect();
        let mut field = LevelField {
            coulomb: Coulomb::with_clouds(&charges, &clouds),
            dipoles,
            loops,
            polygons,
            antennas,
            external: Vec::new(),
            conductors: self.conductors_for(&sources, resolution),
            electrodes: electrodes_for(boxes, dielectrics, &sources, resolution),
            time_offset: 0.0,
            drives: None,
            ferrites: physics::bem::Ferrites::default(),
        };
        // Ferrites, magnetized by the static magnets and coils (time-dependent coils,
        // outside fields and moments are refused with them: `model_issues`).
        if !self.ferrites.is_empty() {
            let size = match resolution {
                Resolution::Preview => physics::bem::Resolution::Preview,
                Resolution::Verify => physics::bem::Resolution::Verify,
                Resolution::Display => physics::bem::Resolution::Display,
            }
            .panel_size();
            let ferrites =
                physics::bem::Ferrites::new(self.physics_ferrites(), |x| field.magnetic(x), size);
            field.ferrites = ferrites;
        }
        // A circuit that cannot be solved is a model issue (`drive_issues`): the field is
        // then without it.
        field.drives = self
            .drives_for(player, &field)
            .ok()
            .flatten()
            .map(std::sync::Arc::new);
        (field, obstacles)
    }

    /// The level's dielectric boxes as boxes (for geometry: obstacles, clearances).
    pub fn dielectric_boxes(&self) -> Vec<physics::bem::BoxElectrode> {
        self.dielectrics
            .iter()
            .map(|d| physics::bem::BoxElectrode {
                center: self.grid.position(d.center),
                angle: d.angle_deg.to_radians(),
                half_length: d.length / 2.0,
                half_thickness: d.thickness / 2.0,
                half_height: d.height / 2.0,
                bias: Bias::Charge(0.0),
            })
            .collect()
    }

    /// The level's ferrite plates as boxes.
    pub fn ferrite_boxes(&self) -> Vec<physics::bem::BoxElectrode> {
        self.ferrites
            .iter()
            .map(|f| physics::bem::BoxElectrode {
                center: self.grid.position(f.center),
                angle: f.angle_deg.to_radians(),
                half_length: f.length / 2.0,
                half_thickness: f.thickness / 2.0,
                half_height: f.height / 2.0,
                bias: Bias::Charge(0.0),
            })
            .collect()
    }

    /// The boxes of the insulating bodies, dielectric and ferrite (for geometry:
    /// obstacles, clearances, contacts).
    pub fn body_boxes(&self) -> Vec<physics::bem::BoxElectrode> {
        let mut b = self.dielectric_boxes();
        b.extend(self.ferrite_boxes());
        b
    }

    /// The level's ferrites for the physics.
    pub fn physics_ferrites(&self) -> Vec<physics::bem::Ferrite> {
        self.ferrite_boxes()
            .into_iter()
            .zip(&self.ferrites)
            .map(|(b, f)| physics::bem::Ferrite {
                shape: physics::bem::BodyShape::Box(b),
                permeability: f.permeability,
            })
            .collect()
    }

    /// The level's dielectrics for the physics: the boxes, then the spheres.
    pub fn physics_dielectrics(&self) -> Vec<physics::bem::Dielectric> {
        let boxes = self
            .dielectric_boxes()
            .into_iter()
            .zip(&self.dielectrics)
            .map(|(b, d)| physics::bem::Dielectric {
                shape: physics::bem::BodyShape::Box(b),
                permittivity: d.permittivity,
            });
        let spheres = self
            .dielectric_spheres
            .iter()
            .map(|d| physics::bem::Dielectric {
                shape: physics::bem::BodyShape::Sphere {
                    center: self.grid.position(d.center),
                    radius: d.radius,
                },
                permittivity: d.permittivity,
            });
        boxes.chain(spheres).collect()
    }

    /// The level's electrodes as physics boxes, at their own bias (a rod or a ring: its
    /// bounding square; their own outlines are `electrode_outlines`).
    pub fn box_electrodes(&self) -> Vec<physics::bem::BoxElectrode> {
        self.electrodes
            .iter()
            .map(|e| physics::bem::BoxElectrode {
                center: self.grid.position(e.center),
                angle: e.angle_deg.to_radians(),
                half_length: e.outer_radius().unwrap_or(e.length / 2.0),
                half_thickness: e.outer_radius().unwrap_or(e.thickness / 2.0),
                half_height: e.height / 2.0,
                bias: match e.bias {
                    ConductorBias::Grounded => Bias::Grounded,
                    ConductorBias::Charge(q) => Bias::Charge(q),
                    ConductorBias::Potential(v) => Bias::Potential(v),
                },
            })
            .collect()
    }
}
