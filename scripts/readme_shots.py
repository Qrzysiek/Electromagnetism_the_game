"""The README's screenshots, taken by the game's own capture mode (`EM_CAPTURE`).

    python scripts/readme_shots.py              # every shot
    python scripts/readme_shots.py magnetron    # these shots only

Each shot opens a level with its reference placement, waits for the verdicts (`EM_WAIT=1`),
parks the editor's cursor in the board's corner (`EM_CURSOR`), saves the window and exits;
ImageMagick (`magick`) then writes it to `docs/screenshots/<name>.webp` (quality 90, a
quarter of the PNG's size). Levels are named by slug, so the curriculum's order can change.
Needs the game built: `cargo build --release -p game` (or `CARGO_TARGET_DIR`'s). A capture
opens a window on the desktop for a few seconds to minutes (the tube levels' verdicts).
"""

import os
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TARGET = os.environ.get("CARGO_TARGET_DIR") or os.path.join(ROOT, "target")
GAME = os.path.join(TARGET, "release", "game.exe" if os.name == "nt" else "game")
PNG = os.path.join(ROOT, "target", "screenshots")
OUT = os.path.join(ROOT, "docs", "screenshots")

# (file name, level slug, capture settings)
SHOTS: list[tuple[str, str, dict[str, str]]] = [
    ("analyzer", "hemispherical_analyzer", {"EM_MAP": "potential", "EM_LINES": "1"}),
    ("calutron", "calutron_space_charge", {}),
    ("synchrotron", "synchrotron_light", {"EM_MAP": "particle"}),
    ("klystron", "klystron", {"EM_MAP": "potential"}),
    ("magnetron", "magnetron", {"EM_MAP": "magnetic"}),
    ("magnetoresistance", "magnetoresistance", {}),
]


def level_number(slug: str) -> str:
    """The 1-based number of the level `levels/<n>_<slug>.json`."""
    for f in os.listdir(os.path.join(ROOT, "levels")):
        n, _, rest = f.partition("_")
        if n.isdigit() and rest == slug + ".json":
            return str(int(n))
    raise SystemExit(f"no level {slug} in levels/")


def shoot(name: str, slug: str, env: dict[str, str]) -> bool:
    png = os.path.join(PNG, name + ".png")
    if os.path.exists(png):
        os.remove(png)
    e = dict(os.environ, EM_CAPTURE=png, EM_LEVEL=level_number(slug), EM_PLACE="reference",
             EM_WAIT="1", EM_CURSOR="0,0", **env)
    try:
        r = subprocess.run([GAME], cwd=ROOT, env=e, capture_output=True, text=True,
                           timeout=900, errors="replace")
    except subprocess.TimeoutExpired:
        print(f"{name}: timed out", flush=True)
        return False
    if not os.path.exists(png):
        print(f"{name}: exit {r.returncode}, no image\n{r.stderr[-1500:]}", flush=True)
        return False
    webp = os.path.join(OUT, name + ".webp")
    subprocess.run(["magick", png, "-quality", "90", "-define", "webp:method=6", webp],
                   check=True)
    print(f"{name}: {os.path.relpath(webp, ROOT)}", flush=True)
    return True


def main() -> None:
    if not os.path.exists(GAME):
        raise SystemExit(f"{GAME} not found: cargo build --release -p game")
    os.makedirs(PNG, exist_ok=True)
    os.makedirs(OUT, exist_ok=True)
    only = set(sys.argv[1:])
    unknown = only - {s[0] for s in SHOTS}
    if unknown:
        raise SystemExit(f"unknown shots: {', '.join(sorted(unknown))}")
    failed = [name for name, slug, env in SHOTS
              if (not only or name in only) and not shoot(name, slug, env)]
    if failed:
        raise SystemExit(f"failed: {', '.join(failed)}")


if __name__ == "__main__":
    main()
