//! The leads under the carriers' own magnetic field (docs/DARWIN.md §2.3, DEC-27). An
//! electrode with a lead (`Tube::leads`, from the board's border) is a circuit node whose
//! branch runs from ground through its supply (its circuit's source, resistance and
//! inductance, or its fixed potential) and along its lead. The lead's current is a source of
//! the field like any current (`darwin_step` deposits it), and the field's EMF along the lead,
//! `−dΛ/dt` with `Λ = ∫ A·dl`, is in the branch's equation: its flux is
//! `Φᵢ = Lᵢ Iᵢ + Σⱼ Mᵢⱼ Iⱼ + Λᵢ′`, `Lᵢ` its own lumped inductance, `Mᵢⱼ` the leads' inductances
//! (the flux along lead i per unit current in lead j, screened by the metals' carriers as the
//! solve screens: one solve per lead at the start), taken implicitly with the circuit, and
//! `Λᵢ′` the rest of the flux along it (the carriers' currents') from the last solve, a step
//! behind (stable: the carriers' coupling to a lead is below one, `M²/(L_lead L_carriers) < 1`).

use super::*;
use crate::darwin2d::{Darwin2d, Edges};

/// A lead's branch (`Tube::lead_branches`).
#[derive(Clone, Copy, Debug)]
pub(super) struct LeadBranch {
    pub electrode: usize,
    /// Its supply: its circuit's source, resistance, inductance and capacitance, or its fixed
    /// potential (an ideal source: no resistance, inductance or capacitance).
    pub source: Waveform,
    pub resistance: f64,
    pub inductance: f64,
    pub capacitance: f64,
    /// The circuit it stands for (None: a fixed potential).
    pub circuit: Option<usize>,
}

/// The leads' currents and fluxes under the self-field (`TubeState::leads`), per branch in
/// `Tube::lead_branches`'s order.
#[derive(Clone, Debug, Default)]
pub struct LeadField {
    /// The electrode of each branch.
    pub electrodes: Vec<usize>,
    /// Each branch's current now, from the border along its lead into its electrode.
    pub current: Vec<f64>,
    /// The leads' inductances `Mᵢⱼ` (m × m, row by row; symmetric).
    pub inductance: Vec<f64>,
    /// Per branch: its lead's flux but the leads' own (`Λᵢ − Σⱼ Mᵢⱼ Iⱼ`) after the last solve,
    /// and the one its last equation took.
    pub rest: Vec<f64>,
    pub rest_used: Vec<f64>,
    /// Per branch: its lead's quadrature (`lead_points`).
    points: Vec<Vec<(DVec3, DVec3)>>,
}

/// A lead's quadrature: the midpoints of equal pieces of each segment, at most half the
/// grid's spacing `h` long, with each piece's vector `dl`. The deposit (`I dl` at each point)
/// and the flux (`Σ A·dl`) use the same points and the solver's same weights, so that the
/// leads' inductances are symmetric and `½ Σ Iᵢ Mᵢⱼ Iⱼ` is the field's energy.
#[allow(
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::cast_precision_loss
)]
pub(super) fn lead_points(path: &[DVec3], h: f64) -> Vec<(DVec3, DVec3)> {
    let mut out = Vec::new();
    for w in path.windows(2) {
        let d = w[1] - w[0];
        let len = libm::sqrt(d.x * d.x + d.y * d.y);
        if len == 0.0 {
            continue;
        }
        let n = libm::ceil(len / (0.5 * h)).max(1.0) as usize;
        let dl = d / n as f64;
        for k in 0..n {
            out.push((w[0] + dl * (k as f64 + 0.5), dl));
        }
    }
    out
}

/// Adds a current `current` along a lead's quadrature to `j` (per area).
fn deposit_lead(d: &Darwin2d, j: &mut Edges, points: &[(DVec3, DVec3)], current: f64) {
    for (p, dl) in points {
        d.deposit_vector(j, p.x, p.y, current * dl.x, current * dl.y);
    }
}

/// The flux of `a` along a lead's quadrature, `∫ A·dl`.
fn flux(d: &Darwin2d, a: &Edges, points: &[(DVec3, DVec3)]) -> f64 {
    points
        .iter()
        .map(|(p, dl)| {
            let (ax, ay) = d.interpolate(a, p.x, p.y);
            ax * dl.x + ay * dl.y
        })
        .sum()
}

impl LeadField {
    /// Adds the branches' currents along their leads to `j`.
    pub(super) fn deposit(&self, d: &Darwin2d, j: &mut Edges) {
        for (points, i) in self.points.iter().zip(&self.current) {
            deposit_lead(d, j, points, *i);
        }
    }

    /// After a solve with these currents (`a` its `A_T`): each lead's flux but the leads' own.
    pub(super) fn update_rest(&mut self, d: &Darwin2d, a: &Edges) {
        let m = self.current.len();
        for i in 0..m {
            let own: f64 = (0..m)
                .map(|j| self.inductance[i * m + j] * self.current[j])
                .sum();
            self.rest[i] = flux(d, a, &self.points[i]) - own;
        }
    }

    /// The field's energy in the leads' own currents, `½ Σ Iᵢ Mᵢⱼ Iⱼ`.
    #[must_use]
    pub fn energy(&self) -> f64 {
        let m = self.current.len();
        let mut w = 0.0;
        for i in 0..m {
            for j in 0..m {
                w += self.current[i] * self.inductance[i * m + j] * self.current[j];
            }
        }
        0.5 * w
    }
}

impl Tube {
    /// The leads' branches: with the self-field, every electrode with a lead (of two corners
    /// at least) and a supply to ground: a circuit (not returning to another electrode, not a
    /// strap, not floating), else its fixed potential (grounded or held). A floating electrode
    /// (a fixed charge) has none; neither does a circuit returning to another electrode
    /// (the level refuses those with the self-field).
    pub(super) fn lead_branches(&self) -> Vec<LeadBranch> {
        if self.darwin.is_none() {
            return Vec::new();
        }
        let mut out = Vec::new();
        for (e, lead) in self.leads.iter().enumerate() {
            if lead.as_ref().is_none_or(|p| p.len() < 2) || e >= self.electrodes.electrodes.len() {
                continue;
            }
            match self.circuits.iter().position(|c| c.electrode == e) {
                Some(ci) => {
                    let c = &self.circuits[ci];
                    if c.to.is_some() || !c.resistance.is_finite() {
                        continue;
                    }
                    out.push(LeadBranch {
                        electrode: e,
                        source: c.source,
                        resistance: c.resistance,
                        inductance: c.inductance,
                        capacitance: c.capacitance,
                        circuit: Some(ci),
                    });
                }
                None => {
                    let v = match self.electrodes.electrodes[e].bias {
                        Bias::Grounded => 0.0,
                        Bias::Potential(v) => v,
                        Bias::Charge(_) => continue,
                    };
                    out.push(LeadBranch {
                        electrode: e,
                        source: Waveform::Dc(v),
                        resistance: 0.0,
                        inductance: 0.0,
                        capacitance: 0.0,
                        circuit: None,
                    });
                }
            }
        }
        out
    }

    /// Whether circuit `ci` is one a lead's branch stands for.
    pub(super) fn has_lead(&self, branches: &[LeadBranch], ci: usize) -> bool {
        branches.iter().any(|b| b.circuit == Some(ci))
    }

    /// The leads' field at the start (None without branches): their quadratures and
    /// inductances, each from a solve with a unit current along one lead and every particle
    /// of the state responding (the metals' carriers screening as in the steps' solves, their
    /// response `β = Σ q²/m` as it is at the start; it changes only with their density's
    /// fluctuations, and the difference falls into the lagged rest), from χ = 0. Symmetrised:
    /// the solves agree with their transposes to their tolerance.
    pub(super) fn lead_field(&self, s: &TubeState) -> Option<LeadField> {
        let d = s.darwin.as_ref()?;
        let branches = self.lead_branches();
        if branches.is_empty() {
            return None;
        }
        let points: Vec<Vec<(DVec3, DVec3)>> = branches
            .iter()
            .map(|b| lead_points(self.leads[b.electrode].as_deref().unwrap_or(&[]), d.h))
            .collect();
        let mut beta = Edges::zeros(d);
        for k in 0..s.x.len() {
            d.deposit_scalar(&mut beta, s.x[k].x, s.x[k].y, s.charge[k] * s.qm[k]);
        }
        let m = branches.len();
        let mut l = vec![0.0; m * m];
        let mut solver = d.clone();
        for k in 0..m {
            let mut j = Edges::zeros(d);
            deposit_lead(d, &mut j, &points[k], 1.0);
            solver.chi.fill(0.0);
            solver.solve(&j, &beta, DARWIN_TOLERANCE);
            let a = solver.a();
            for i in 0..m {
                l[i * m + k] = flux(d, &a, &points[i]);
            }
        }
        for i in 0..m {
            for k in 0..i {
                let mean = 0.5 * (l[i * m + k] + l[k * m + i]);
                l[i * m + k] = mean;
                l[k * m + i] = mean;
            }
        }
        Some(LeadField {
            electrodes: branches.iter().map(|b| b.electrode).collect(),
            current: vec![0.0; m],
            inductance: l,
            rest: vec![0.0; m],
            rest_used: vec![0.0; m],
            points,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::zinv::{Electrode, Electrodes, Section};

    /// A box 2 × 1 at `(x0, 10)` driven through a lead straight up from `(x0, 0)` to its
    /// lower face by a source `s0` from t = 0, through `resistance` with `capacitance` to
    /// ground, a grounded box 4 × 1 at `(x0, 13)` above it; no particles; the self-field at
    /// c = 3 on nodes 0.5 apart over −6…26 (symmetric about x = 10).
    fn lead_tube(x0: f64, s0: f64, resistance: f64, capacitance: f64) -> Tube {
        let rect = |cx: f64, cy: f64, half_length: f64| Section::Rect {
            center: DVec3::new(cx, cy, 0.0),
            angle: 0.0,
            half_length,
            half_thickness: 0.5,
        };
        let electrodes = Electrodes::new(
            vec![
                Electrode {
                    section: rect(x0, 10.0, 1.0),
                    bias: Bias::Grounded,
                },
                Electrode {
                    section: rect(x0, 13.0, 2.0),
                    bias: Bias::Grounded,
                },
            ],
            &[],
            0.25,
        );
        Tube {
            circuits: vec![Circuit {
                electrode: 0,
                source: Waveform::Dc(s0),
                resistance,
                capacitance,
                inductance: 0.0,
                to: None,
            }],
            darwin: Some(Darwin2d::new(63, 63, 0.5, (-6.0, -6.0), 3.0)),
            leads: vec![
                Some(vec![DVec3::new(x0, 0.0, 0.0), DVec3::new(x0, 9.5, 0.0)]),
                None,
            ],
            ..Tube::new(electrodes, 0.05, 0.02)
        }
    }

    /// The leads' field adds up: a lead's inductance is positive and two leads' matrix
    /// positive definite, and the field of their currents holds `½ Σ Iᵢ Mᵢⱼ Iⱼ`, the solver's
    /// own magnetic energy, to its tolerance (the deposit and the flux pair exactly). Two
    /// boxes side by side, each with its lead from below, one driven through a resistance,
    /// the other held at −0.5.
    #[test]
    fn leads_inductance_is_the_field_energy() {
        let mut tube = lead_tube(8.0, 1.0, 1.0, 0.0);
        let rect = |cx: f64| Section::Rect {
            center: DVec3::new(cx, 10.0, 0.0),
            angle: 0.0,
            half_length: 1.0,
            half_thickness: 0.5,
        };
        tube.electrodes = Electrodes::new(
            vec![
                Electrode {
                    section: rect(8.0),
                    bias: Bias::Grounded,
                },
                Electrode {
                    section: rect(12.0),
                    bias: Bias::Potential(-0.5),
                },
            ],
            &[],
            0.25,
        );
        tube.leads[1] = Some(vec![DVec3::new(12.0, 0.0, 0.0), DVec3::new(12.0, 9.5, 0.0)]);
        let mut s = tube.state();
        let lf = s.leads.as_ref().expect("two leads");
        assert_eq!(lf.electrodes, vec![0, 1]);
        let l = lf.inductance.clone();
        println!("the leads' inductances: {l:?}");
        assert!(
            l[0] > 0.0 && l[3] > 0.0 && l[0] * l[3] > l[1] * l[2],
            "{l:?}"
        );
        assert_eq!(l[1].to_bits(), l[2].to_bits());
        for _ in 0..20 {
            tube.step(&mut s);
        }
        let lf = s.leads.as_ref().expect("two leads");
        let w = s.darwin.as_ref().expect("the field").magnetic_energy();
        println!(
            "currents {:?}: field energy {w:.6e}, ½ IMI {:.6e}",
            lf.current,
            lf.energy()
        );
        assert!(lf.current[0] != 0.0 && lf.current[1] != 0.0);
        assert!((w - lf.energy()).abs() < 1e-8 * w, "{w} vs {}", lf.energy());
    }

    /// Two straight antiparallel leads, 4 apart and 24 long, have the loop inductance
    /// `L₁₁ + L₂₂ − 2L₁₂` of the infinite plane's Darwin field in closed form,
    /// `μ₀ [8dl atan(l/d) + l² ln(1 + d²/l²) − 3d² ln(1 + l²/d²)]/(4π)` (102.44 at c = 3; long
    /// leads `μ₀dl`: `scripts/wolfram/d8_two_leads.wls`, checked by summing ψ² on a fine
    /// grid), within 3 % (set before measuring: the leads on the grid's edge lines, sharp; the
    /// frame ~50 away, its share of the pair's quadrupole field below 0.5 %; the ends' grid
    /// error ~h/l = 2 %). An open lead's own inductance has no such limit (its field falls as
    /// 1/r, its energy grows as the log of the frame's distance): only a pair whose currents
    /// close does.
    #[test]
    fn two_leads_make_the_loop_inductance() {
        let (d, l, c) = (4.0, 24.0, 3.0);
        let rect = |cx: f64| Section::Rect {
            center: DVec3::new(cx, 22.5, 0.0),
            angle: 0.0,
            half_length: 1.0,
            half_thickness: 0.5,
        };
        // Nodes 0.5 apart over −54…74: the y-edges at x = −53.75 + 0.5 i.
        let (x1, x2) = (7.75, 7.75 + d);
        let electrodes = Electrodes::new(
            vec![
                Electrode {
                    section: rect(x1),
                    bias: Bias::Grounded,
                },
                Electrode {
                    section: rect(x2),
                    bias: Bias::Grounded,
                },
            ],
            &[],
            0.25,
        );
        let lead = |x: f64| Some(vec![DVec3::new(x, -2.0, 0.0), DVec3::new(x, -2.0 + l, 0.0)]);
        let tube = Tube {
            darwin: Some(Darwin2d::new(255, 255, 0.5, (-54.0, -54.0), c)),
            leads: vec![lead(x1), lead(x2)],
            ..Tube::new(electrodes, 0.05, 0.02)
        };
        let s = tube.state();
        let m = &s.leads.as_ref().expect("two leads").inductance;
        let loop_l = m[0] + m[3] - 2.0 * m[1];
        let mu0 = 4.0 * std::f64::consts::PI / (c * c);
        let exact = mu0
            * (8.0 * d * l * libm::atan(l / d) + l * l * libm::log(1.0 + d * d / (l * l))
                - 3.0 * d * d * libm::log(1.0 + l * l / (d * d)))
            / (4.0 * std::f64::consts::PI);
        println!(
            "two leads: L₁₁ {:.4}, L₁₂ {:.4}: loop {loop_l:.4} against {exact:.4} ({:+.2} %; μ₀dl {:.4})",
            m[0],
            m[1],
            100.0 * (loop_l / exact - 1.0),
            mu0 * d * l
        );
        assert!((loop_l / exact - 1.0).abs() < 0.03, "{loop_l} vs {exact}");
    }

    /// A source switched on behind a lead charges its electrode as a series RLC circuit,
    /// `L dI/dt + R I + Q/C = S`, with L the lead's inductance from the field and C the
    /// electrode's (to the grounded box) and its capacitor's: the current against the
    /// underdamped solution `I = S/(L ω_d) e^{−αt} sin ω_d t` over one and a half periods,
    /// within 1e-4 of its peak (the trapezoidal rule's phase error, `(ω dt)²/12` per radian,
    /// ~1e-6 here; set before measuring). The circuit and the field act on each other through
    /// the lead alone: no particles, so nothing is lagged.
    #[test]
    fn a_lead_makes_an_rlc_circuit() {
        let (s0, r, c_ext) = (1.0, 0.5, 2.0);
        let tube = lead_tube(10.0, s0, r, c_ext);
        let mut s = tube.state();
        let l = s.leads.as_ref().expect("a lead").inductance[0];
        let c = tube.electrodes.capacitance[0][0] + c_ext;
        let alpha = r / (2.0 * l);
        let w0 = 1.0 / libm::sqrt(l * c);
        assert!(alpha < w0, "overdamped: α {alpha}, ω₀ {w0}");
        let wd = libm::sqrt(w0 * w0 - alpha * alpha);
        let end = 1.5 * std::f64::consts::TAU / wd;
        let (mut worst, mut peak) = (0.0_f64, 0.0_f64);
        while s.t < end {
            tube.step(&mut s);
            let i = s.leads.as_ref().expect("a lead").current[0];
            let exact = s0 / (l * wd) * libm::exp(-alpha * s.t) * libm::sin(wd * s.t);
            worst = worst.max((i - exact).abs());
            peak = peak.max(exact.abs());
        }
        println!(
            "RLC: L {l:.4}, C {c:.4}, R {r}: ω_d {wd:.5}, α {alpha:.5}; worst {:.2e} of the peak {peak:.4}",
            worst / peak
        );
        assert!(worst < 1e-4 * peak, "{worst} of {peak}");
    }

    /// Without resistance the lead, a small inductor of its circuit's (1, against the lead's
    /// ~16; with neither the supply would be an ideal source, its electrode charged from the
    /// start) and the capacitances make an LC circuit whose energy is the source's work:
    /// `½ (L_c + M) I² + ½ Q²/C = S ΔQ` at every step, to rounding (1e-9 of the work; the
    /// implicit midpoint rule keeps the quadratic energy exactly), `½ M I²` the field's own
    /// energy. And run back: the current reversed, as many steps return the charge and the
    /// current to their start (1e-9 of their swing): the LC circuit with its field is
    /// symmetric under time reversal (I and B odd).
    #[test]
    fn a_lead_lc_keeps_its_energy_and_runs_back() {
        let (s0, l_c) = (1.0, 1.0);
        let mut tube = lead_tube(10.0, s0, 0.0, 2.0);
        tube.circuits[0].inductance = l_c;
        let mut s = tube.state();
        let l = s.leads.as_ref().expect("a lead").inductance[0];
        let c = tube.electrodes.capacitance[0][0] + 2.0;
        let q0 = s.node_charge[0];
        let mut worst = 0.0_f64;
        for _ in 0..400 {
            tube.step(&mut s);
            let i = s.leads.as_ref().expect("a lead").current[0];
            let q = s.node_charge[0] - q0;
            let w = s.darwin.as_ref().expect("the field").magnetic_energy();
            let work = s0 * q;
            let energy = w + 0.5 * l_c * i * i + 0.5 * q * q / c;
            worst = worst.max((energy - work).abs() / work.abs().max(1e-12));
            assert!((w - 0.5 * l * i * i).abs() < 1e-8 * w.max(1e-30));
        }
        println!("LC: energy against the source's work, worst {worst:.2e}");
        assert!(worst < 1e-9, "{worst}");
        let (q1, i1) = (
            s.node_charge[0],
            s.leads.as_ref().expect("a lead").current[0],
        );
        s.leads.as_mut().expect("a lead").current[0] = -i1;
        for _ in 0..400 {
            tube.step(&mut s);
        }
        let (q2, i2) = (
            s.node_charge[0],
            s.leads.as_ref().expect("a lead").current[0],
        );
        println!(
            "run back: charge {q2:.3e} from {q0:.3e} (at the turn {q1:.4}), current {i2:.3e} (at the turn {i1:.4})"
        );
        assert!((q2 - q0).abs() < 1e-9 * (q1 - q0).abs(), "{q2} vs {q0}");
        assert!(i2.abs() < 1e-9 * i1.abs(), "{i2}");
    }

    /// A setup symmetric under the mirror x → 20 − x (the lead on x = 10, the grid's middle)
    /// has a field odd under it (B_z a pseudoscalar in the plane), to the solver's tolerance.
    #[test]
    fn a_lead_field_is_odd_under_a_mirror() {
        let tube = lead_tube(10.0, 1.0, 0.5, 2.0);
        let mut s = tube.state();
        for _ in 0..100 {
            tube.step(&mut s);
        }
        let d = s.darwin.as_ref().expect("the field");
        let b = d.b_z();
        let top = b.iter().fold(0.0_f64, |m, v| m.max(v.abs()));
        let mut worst = 0.0_f64;
        for j in 0..d.ny {
            for i in 0..d.nx {
                let m = d.nx - 1 - i;
                worst = worst.max((b[i + d.nx * j] + b[m + d.nx * j]).abs());
            }
        }
        println!(
            "mirror: the field's even part {:.2e} of its largest {top:.3e}",
            worst / top
        );
        assert!(worst < 1e-8 * top, "{worst} of {top}");
    }
}
