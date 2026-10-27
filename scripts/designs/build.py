"""The machinery: probing and placing automatic detectors, solving and checking references, writing the level files and their texts (`--text`), rendering the formulas."""

import json
import math
import os
import re
import subprocess
import sys
import tempfile
from reality import REALITY
from designs.common import Auto, ROOT, box, run_generator, slider_problems
from designs.jackson import FINALES, MUST_FAIL_ALONE, NEEDS_BODY, check_indirect, needs_elements
from designs.curriculum import LEVELS, write_curriculum

def landing_points(lvl):
    """Where the reference flights end, from `generator check` on a probe level."""
    with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False, encoding="utf-8") as f:
        json.dump(lvl, f)
        path = f.name
    try:
        text = run_generator("check", path)
    finally:
        os.unlink(path)
    # One list of (outcome, x, y) per shot: one entry per disturbance.
    ends = {}
    for line in text.splitlines():
        m = re.search(r"shot (\d+)(?: disturbance \d+)?: outcome (\w+).*"
                      r"ends at \(([-\d.]+), ([-\d.]+)\)", line)
        if m:
            ends.setdefault(int(m.group(1)), []).append(
                (m.group(2), float(m.group(3)), float(m.group(4))))
    return [ends[k] for k in sorted(ends)]



def resolve_auto(lvl):
    autos = [s["detector"] for s in lvl["shots"]]
    if not any(isinstance(d, Auto) for d in autos):
        return lvl
    probe = json.loads(json.dumps(lvl, default=lambda d: box(d.strip)))
    if not probe["reference_solution"]:
        # No designed reference: let the solver find a setup that brings every shot into
        # its probe strip, and use it as the reference.
        with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False,
                                         encoding="utf-8") as f:
            json.dump(probe, f)
            path = f.name
        try:
            run_generator("solve", path, "--write", "--restarts", "64", "--iterations", "600")
            probe = json.load(open(path, encoding="utf-8"))
        finally:
            os.unlink(path)
        if not probe["reference_solution"]:
            sys.exit(f"{lvl['name']}: no reference setup found for the probe strips")
        lvl["reference_solution"] = probe["reference_solution"]
    ends = landing_points(probe)
    placed = []
    for s, d, flights in zip(lvl["shots"], autos, ends):
        if not isinstance(d, Auto):
            continue
        for outcome, _, _ in flights:
            if outcome != "Arrived":
                sys.exit(f"{lvl['name']}: a reference flight misses its probe strip ({outcome})")
        x0, y0, x1, y1 = d.strip
        acrosses = [y if d.axis == "x" else x for _, x, y in flights]
        a_lo, a_hi = min(acrosses), max(acrosses)
        size = d.size
        if len(acrosses) == 1:
            across = a_lo
            lo = math.floor(across - d.size / 2 + 0.5)
            # Keep the landing point at least 0.2 cells inside the detector.
            while across - lo < 0.2:
                lo -= 1
            while lo + d.size - across < 0.2:
                lo += 1
        else:
            # Under disturbances the landing points spread: the detector covers all of
            # them with 0.2 cells to spare, growing beyond `size` if it must.
            while True:
                lo = math.floor(a_lo - 0.2)
                if lo + size - a_hi >= 0.2:
                    break
                size += 1
            # Centre the spread in the detector.
            lo += max(0, math.floor(((lo + size - a_hi) - (a_lo - lo)) / 2))
        # Stay inside the arena.
        n = lvl["grid"]["ny" if d.axis == "x" else "nx"] * lvl["grid"]["subdivision"]
        lo = max(min(lo, n - size), 0)
        across = (a_lo, a_hi)
        # Shots of the same species (particle and energy, differing only in direction)
        # share one detector; different species get detectors of their own.
        screen = (x0, x1) if d.axis == "x" else (y0, y1)
        species = (d.axis, screen, json.dumps(s["particle"]), s["launch"]["kinetic_energy"],
                   id(d) if d.own else None)
        placed.append([species, [s], across[0], across[1], lo, lo + size])
    groups = {}
    for p in placed:
        g = groups.setdefault(p[0], p)
        if g is not p:
            g[1] += p[1]
            g[2], g[3] = min(g[2], p[2]), max(g[3], p[3])
            g[4], g[5] = min(g[4], p[4]), max(g[5], p[5])
    # Split overlapping detectors on the same screen at the node halfway between the
    # landing points.
    placed = sorted(groups.values(), key=lambda p: (p[0][:2], p[2]))
    for a, b in zip(placed, placed[1:]):
        if a[0][:2] == b[0][:2] and a[5] > b[4]:
            m = round((a[3] + b[2]) / 2)
            if m - a[3] < 0.2 or b[2] - m < 0.2 or m - a[4] < 1 or b[5] - m < 1:
                sys.exit(f"{lvl['name']}: two landing points are too close to separate")
            a[5], b[4] = m, m
    for (axis, (c0, c1), _, _, _), shots, _, _, lo, hi in placed:
        for s in shots:
            s["detector"] = box((c0, lo, c1, hi) if axis == "x" else (lo, c0, hi, c1))
    return lvl



def previous_build(slug):
    """The level file of `slug` from the previous build (whatever its number), or None."""
    for name in os.listdir(os.path.join(ROOT, "levels")):
        if re.fullmatch(rf"\d+_{re.escape(slug)}\.json", name):
            return os.path.join(ROOT, "levels", name)
    return None



def with_reality(slug, lvl):
    """The level with what it represents (`scripts/reality.py`): its headline and remarks
    after the description, and each electrode's and metal's function."""
    r = REALITY.get(slug)
    if r is None:
        sys.exit(f"{slug}: no entry in scripts/reality.py (what the level represents)")
    for key, items in (("electrodes", lvl.get("electrodes", [])),
                       ("metals", (lvl.get("tube") or {}).get("metals", []))):
        if len(r[key]) != len(items):
            sys.exit(f"{slug}: scripts/reality.py names {len(r[key])} {key} of {len(items)}")
        for item, label in zip(items, r[key]):
            item["label"] = label
    return with_reality_dict(lvl, {"represents": r["represents"], "remarks": [
        {"kind": kind, "line": line, "detail": detail} for kind, line, detail in r["remarks"]]})



def text_of(lvl):
    """A level's texts: what `--text` rewrites (they change no run)."""
    return {"name": lvl["name"], "description": lvl["description"],
            "reality": lvl.get("reality"),
            "electrodes": [e.get("label") for e in lvl.get("electrodes", [])],
            "metals": [m.get("label") for m in (lvl.get("tube") or {}).get("metals", [])]}



def with_text(lvl, text):
    """`lvl` with the texts `text` (`text_of`) and nothing else changed."""
    out = with_reality_dict(lvl, text["reality"])
    out["name"], out["description"] = text["name"], text["description"]
    for e, label in zip(out.get("electrodes", []), text["electrodes"]):
        e["label"] = label
    for m, label in zip((out.get("tube") or {}).get("metals", []), text["metals"]):
        m["label"] = label
    return out



def with_reality_dict(lvl, reality):
    """`lvl` with `reality` right after its description."""
    out = {}
    for k, v in lvl.items():
        if k != "reality":
            out[k] = v
        if k == "description":
            out["reality"] = reality
    return out



def same_physics(a, b):
    """Whether two levels differ at most in their texts (`text_of`), limits and reference:
    then a reference of one is a solution of the other (if within limits)."""
    automatic = any(isinstance(sh["detector"], Auto) for sh in a["shots"] + b["shots"])

    def strip(l):
        l = {k: v for k, v in l.items()
             if k not in ("name", "description", "reality", "limits", "reference_solution")}
        for key in ("electrodes",):
            if key in l:
                l[key] = [{k: v for k, v in e.items() if k != "label"} for e in l[key]]
        if l.get("tube") and l["tube"].get("metals"):
            l["tube"] = dict(l["tube"], metals=[{k: v for k, v in m.items() if k != "label"}
                                                 for m in l["tube"]["metals"]])
        if automatic:
            # Automatic detectors are placed from the reference: compare without them.
            l["shots"] = [{k: v for k, v in sh.items() if k != "detector"} for sh in l["shots"]]
        return json.dumps(numbers(l), sort_keys=True)
    return strip(a) == strip(b)



def numbers(v):
    """`v` with every number a float: a file the generator wrote holds 8 as 8.0."""
    if isinstance(v, dict):
        return {k: numbers(x) for k, x in v.items()}
    if isinstance(v, list):
        return [numbers(x) for x in v]
    if isinstance(v, int) and not isinstance(v, bool):
        return float(v)
    return v



def text_only(only):
    """`--text`: rewrite the built levels' texts from this script, refusing any level whose
    physics, limits or designed reference differ from its built file (build those in full)."""
    for i, (slug, make) in enumerate(LEVELS, start=1):
        key = f"{i:02d}_{slug}"
        if only and not (key.startswith(only) or slug in only.split(",")):
            continue
        path = os.path.join(ROOT, "levels", f"{key}.json")
        if not os.path.exists(path):
            sys.exit(f"{key}: not built yet: build it in full")
        previous = json.load(open(path, encoding="utf-8"))
        lvl = with_reality(slug, make())
        # Empty lists are absent to the game (serde's default).
        limits = [{k: v for k, v in l["limits"].items() if v != []} for l in (lvl, previous)]
        if (not same_physics(lvl, previous) or limits[0] != limits[1]
                or (lvl["reference_solution"]
                    and lvl["reference_solution"] != previous["reference_solution"])):
            sys.exit(f"{key}: more than its text changed: build it in full (--only {slug})")
        new = with_text(previous, text_of(lvl))
        if new != previous:
            with open(path, "w", newline="\n", encoding="utf-8") as f:
                f.write(json.dumps(new, indent=2, ensure_ascii=False) + "\n")
            print(f"{key}: text rewritten")
    render_math()



def main():
    only = sys.argv[sys.argv.index("--only") + 1] if "--only" in sys.argv else None
    if "--text" in sys.argv:
        return text_only(only)
    # Renumber: move each slug's previous file to its place in the curriculum first.
    wanted = {slug: f"{i:02d}_{slug}.json" for i, (slug, _) in enumerate(LEVELS, start=1)}
    for slug, name in wanted.items():
        old = previous_build(slug)
        new = os.path.join(ROOT, "levels", name)
        if old and os.path.normcase(old) != os.path.normcase(new):
            os.replace(old, new)
    for name in os.listdir(os.path.join(ROOT, "levels")):
        m = re.fullmatch(r"\d+_(\w+)\.json", name)
        if m and m.group(1) not in wanted:
            print(f"removing {name} (no longer in the curriculum)")
            os.remove(os.path.join(ROOT, "levels", name))
    write_curriculum()
    for i, (slug, make) in enumerate(LEVELS, start=1):
        key = f"{i:02d}_{slug}"
        if only and not (key.startswith(only) or slug in only.split(",")):
            continue
        path = os.path.join(ROOT, "levels", f"{key}.json")
        previous = json.load(open(path, encoding="utf-8")) if os.path.exists(path) else None
        kept = False
        lvl = with_reality(slug, make())
        if (previous and previous["reference_solution"] and not lvl["reference_solution"]
                and same_physics(lvl, previous)):
            # Only the text or the limits changed: keep the verified reference (checked
            # again below) instead of searching anew.
            lvl["reference_solution"] = previous["reference_solution"]
            print(f"{key}: reference kept from the previous build")
            kept = True
        lvl = resolve_auto(lvl)
        if slug in MUST_FAIL_ALONE:
            alone = dict(lvl, reference_solution=[])
            with open(path, "w", newline="\n", encoding="utf-8") as f:
                f.write(json.dumps(alone, indent=2, ensure_ascii=False) + "\n")
            report = run_generator("check", path)
            lines = [l for l in report.splitlines() if "shot " in l]
            works = all("Arrived" in l and "Verified" in l for l in lines)
            print(f"{key}: idealised design alone: {'WORKS (not a real iteration)' if works else 'fails, as it should'}")
            for l in lines:
                print("   ", l.strip())
            if works:
                sys.exit(f"{key}: the idealised design still works on its own")
        with open(path, "w", newline="\n", encoding="utf-8") as f:
            f.write(json.dumps(lvl, indent=2, ensure_ascii=False) + "\n")
        indirect = check_indirect(lvl)
        if indirect:
            sys.exit(f"{key}: goal particles are reached directly: " + "; ".join(indirect))
        if slug in FINALES:
            fewer = FINALES[slug]
            ok_fewer = needs_elements(lvl, fewer)
            print(f"{key}: not solved with {fewer} charges: {ok_fewer}", flush=True)
            if not ok_fewer:
                sys.exit(f"{key}: solvable with only {fewer} charges")
        if not lvl["reference_solution"]:
            run_generator("solve", path, "--write", "--restarts", "64", "--iterations", "400")
        report = run_generator("check", path)
        if "placement Ok" not in report and kept:
            # A kept reference outside the new limits (e.g. a smaller region): solve anew.
            print(f"{key}: kept reference no longer allowed; solving anew", flush=True)
            lvl["reference_solution"] = []
            with open(path, "w", newline="\n", encoding="utf-8") as f:
                f.write(json.dumps(lvl, indent=2, ensure_ascii=False) + "\n")
            run_generator("solve", path, "--write", "--restarts", "64", "--iterations", "400")
            report = run_generator("check", path)
        ok = "placement Ok" in report and all(
            "Arrived" in l and "Verified" in l for l in report.splitlines() if "shot " in l)
        print(f"{key}: reference {'verified' if ok else 'NOT VERIFIED'}", flush=True)
        if slug in NEEDS_BODY:
            built = json.load(open(path, encoding="utf-8"))
            vacuum = dict(built, dielectrics=[], dielectric_spheres=[], ferrites=[])
            with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False,
                                             encoding="utf-8") as f:
                json.dump(vacuum, f)
                probe_path = f.name
            try:
                alone = run_generator("check", probe_path)
            finally:
                os.unlink(probe_path)
            works = all("Arrived" in l and "Verified" in l
                        for l in alone.splitlines() if "shot " in l)
            print(f"{key}: reference without the dielectrics and ferrites: "
                  f"{'WORKS (the body does not matter)' if works else 'fails, as it should'}")
            if works:
                sys.exit(f"{key}: the reference works without its dielectrics and ferrites")
        if not ok:
            print(report)
        for problem in slider_problems(json.load(open(path, encoding="utf-8"))):
            sys.exit(f"{key}: {problem}")

    render_math()



def render_math():
    """Renders the LaTeX formulas ($...$) of every level's description (shipped and
    custom) to SVG with MathJax (scripts/math/render.mjs; `npm install` in scripts/math
    once), into levels/math/ with a manifest the game reads (crates/game/src/math.rs)."""
    formulas = []
    for folder in (os.path.join(ROOT, "levels"), os.path.join(ROOT, "levels", "custom")):
        if not os.path.isdir(folder):
            continue
        for name in sorted(os.listdir(folder)):
            if not name.endswith(".json") or name in ("curriculum.json", "golden_hashes.json"):
                continue
            text = json.load(open(os.path.join(folder, name), encoding="utf-8")).get(
                "description", "")
            if text.count("$") % 2:
                sys.exit(f"{name}: unmatched $ in the description")
            formulas += re.findall(r"\$([^$]+)\$", text)
    with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False,
                                     encoding="utf-8") as f:
        json.dump(sorted(set(formulas)), f)
        path = f.name
    try:
        out = subprocess.run(["node", os.path.join(ROOT, "scripts", "math", "render.mjs"), path,
                              os.path.join(ROOT, "levels", "math")],
                             capture_output=True, text=True, encoding="utf-8")
    finally:
        os.unlink(path)
    if out.returncode != 0:
        sys.exit("rendering the formulas failed: " + out.stdout + out.stderr)
    print(out.stdout.strip())
    # Formulas are not broken across lines (MathJax 3's SVG output has no line breaking);
    # the game shrinks one wider than the panel to fit. Warn well before that: 30 ex is
    # about 60 % of the side panel (340 points, 1 ex = 0.52 x the 12.5-point body font).
    manifest = json.load(open(os.path.join(ROOT, "levels", "math", "manifest.json"),
                              encoding="utf-8"))
    for tex, m in sorted(manifest.items()):
        if m["width_ex"] > MAX_FORMULA_EX:
            print(f"warning: formula {m['width_ex']:.1f} ex wide (over {MAX_FORMULA_EX}): "
                  f"${tex}$ -- split it into several $...$ at its relation signs, so that it "
                  "can wrap between them")



MAX_FORMULA_EX = 30.0



if __name__ == "__main__":
    main()
