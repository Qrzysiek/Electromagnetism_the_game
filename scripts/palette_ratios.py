"""Read-only (`python scripts/palette_ratios.py`; docs/LEVELS.md, Showcase review): each particle
level's charge palette as `k|qQ|/E_kin`, the head-on closest approach in cells. A charge at
distance d bends the path by about `k|qQ|/(d E_kin)`: the player's choice of place (one to tens
of cells) takes every palette from gentle bends (far away) to turning the projectile back
(close up)."""
import glob
import json
import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
rows = []
for f in sorted(glob.glob(os.path.join(ROOT, "levels", "[0-9]*.json")),
                key=lambda p: int(os.path.basename(p).split("_")[0])):
    level = json.load(open(f, encoding="utf-8"))
    if level.get("tube") or not level.get("shots"):
        continue
    mags = level.get("limits", {}).get("magnitudes") or []
    if not mags:
        continue
    shot = level["shots"][0]
    q = abs(shot["particle"]["charge"])
    ke = shot["launch"]["kinetic_energy"]
    lo, hi = min(abs(m) for m in mags), max(abs(m) for m in mags)
    rows.append((os.path.basename(f).split("_")[0], level["name"], q * lo / ke, q * hi / ke))
for n, name, a, b in rows:
    print(f"{n:>3} {name[:34]:34} closest approach {a:8.3g} .. {b:8.3g} cells (x{b / a:4.1f})")
print(len(rows), "levels with charge palettes")
