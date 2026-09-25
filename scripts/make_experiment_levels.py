"""Levels modelled on real experiments and instruments (point charges only).

Writes levels/1x_*.json. Reference solutions are then found and stored with
    cargo run --release -p generator -- solve levels/<file> --write

Units are the game's dimensionless units (k = 1, cell = 1). Particles are weakly
charged and electrodes strongly charged so that neglected radiation stays below
1e-10 of the launch energy (PHYSICS.md section 8).
"""

import json
import math


def level(name, desc, particle, launch, direction_deg, ke, detector, charges, limits,
          region=None, c=5.0, grid=(30, 20), t_max=300.0):
    a = math.radians(direction_deg)
    return {
        "format_version": 1,
        "engine_version": "0.1.0",
        "name": name,
        "description": desc,
        "grid": {"nx": grid[0], "ny": grid[1], "nz": 0, "subdivision": 1},
        "physics": {"c": c, "charge_radius": 0.3, "t_max": t_max,
                    "tolerances": {"preview": 1e-10, "verify": 1e-12}},
        "particle": {"charge": particle[0], "mass": particle[1], "radius": 0.0},
        "launch": {"node": [launch[0], launch[1], 0],
                   "direction": [math.cos(a), math.sin(a), 0.0], "kinetic_energy": ke},
        "detector": {"min": [detector[0], detector[1], 0], "max": [detector[2], detector[3], 0]},
        "level_charges": [{"node": [x, y, 0], "charge": q} for (x, y, q) in charges],
        "limits": {"max_charges": limits[0], "magnitudes": limits[1],
                   "allow_positive": limits[2], "allow_negative": limits[3],
                   **({"region": {"min": [region[0], region[1], 0], "max": [region[2], region[3], 0]}}
                      if region else {})},
        "reference_solution": [],
    }


M = 1e6
levels = {
    "11_geiger_marsden": level(
        "Geiger–Marsden (1909)",
        "Alpha particles fired at gold foil sometimes bounced straight back: the atom has "
        "a tiny, heavy, positive nucleus. Nudge the alpha so it backscatters into the "
        "detector. The nucleus field is exact Coulomb; the deflection follows "
        "tan(θ/2) = qQ / (2 T₀ b).",
        particle=(2e-6, 4.0), launch=(0, 11), direction_deg=0.0, ke=4.0,
        detector=(0, 16, 3, 20), charges=[(18, 10, 8 * M)],
        limits=(1, [0.25 * M, 0.5 * M, 1 * M], True, True), region=(3, 6, 7, 14), c=20.0),
    "12_thomson_crt": level(
        "Thomson's cathode-ray tube (1897)",
        "J. J. Thomson deflected cathode rays with charged plates and showed they are "
        "negative particles: electrons. Here an electron at 0.2 c must hit a one-cell "
        "spot on the screen. Build deflection plates from point charges.",
        particle=(-1e-6, 1.0), launch=(0, 10), direction_deg=0.0, ke=0.5,
        detector=(29, 15, 30, 16), charges=[],
        limits=(2, [0.5 * M, 1 * M], True, True), region=(8, 5, 14, 15)),
    "13_einzel_lens": level(
        "Einzel lens",
        "An Einzel lens focuses charged-particle beams in electron microscopes and ion "
        "sources without changing their energy. The beam leaves the source diverging at "
        "20°; bring it back onto the axis at the small detector. Only charges of the "
        "particle's sign are available (a decelerating lens).",
        particle=(1e-6, 1.0), launch=(0, 10), direction_deg=20.0, ke=0.5,
        detector=(28, 9, 30, 11), charges=[],
        limits=(2, [0.5 * M, 1 * M, 2 * M], True, False), region=(12, 3, 16, 17)),
    "14_hemispherical_analyzer": level(
        "Hemispherical analyzer (XPS)",
        "Photoelectron spectrometers sort electrons by energy between two concentric "
        "hemispheres. Between them the field is exactly that of a point charge at the "
        "centre. Bend the particle through 180° into the detector below the source: a "
        "circular orbit of radius R needs |qQ| = p·v·R.",
        particle=(1e-6, 1.0), launch=(4, 15), direction_deg=0.0, ke=0.5,
        detector=(3, 4, 5, 6), charges=[],
        limits=(1, [2 * M, 3 * M, 4 * M, 5 * M, 6 * M], False, True), region=(4, 3, 4, 14)),
    "15_reflectron": level(
        "Reflectron",
        "Time-of-flight mass spectrometers use an electrostatic mirror, the reflectron, "
        "to turn ions around and sharpen their arrival times. Build a mirror that sends "
        "the ion back to the detector beside the source.",
        particle=(1e-6, 1.0), launch=(0, 13), direction_deg=0.0, ke=0.5,
        detector=(0, 4, 2, 7), charges=[],
        limits=(3, [1 * M, 2 * M, 4 * M], True, False), region=(20, 0, 29, 20)),
    "16_beta_spectrometer": level(
        "Beta-ray spectrometer",
        "Beta decay electrons are fast. Here v = 0.91 c (γ = 2.41). As in the "
        "hemispherical analyzer, bend the electron through 180° into the detector, but "
        "now the circular-orbit condition |qQ| = γ m v² R differs from the Newtonian "
        "2 T₀ R. Trust relativity, not intuition.",
        particle=(-1e-6, 1.0), launch=(4, 17), direction_deg=0.0, ke=math.sqrt(2.0),
        detector=(3, 2, 5, 4), charges=[],
        limits=(1, [10 * M, 12 * M, 14 * M, 17 * M, 20 * M], True, False), region=(4, 3, 4, 16),
        c=1.0),
}

for key, value in levels.items():
    with open(f"levels/{key}.json", "w", newline="\n", encoding="utf-8") as f:
        f.write(json.dumps(value, indent=2, ensure_ascii=False) + "\n")
print("wrote", len(levels), "levels")
