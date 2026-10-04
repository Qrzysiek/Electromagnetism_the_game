"""Show the GPU upload cap with two captures. Does not write a level file.

Both scenes contain 1025 charges. In one the only nonzero charge is last, so a
prefix of 1024 drops it. In the other it is first, so the upload keeps it.
"""

from __future__ import annotations

import json
import os
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
GAME = ROOT / "target" / "release" / "game.exe"
LOG = ROOT / "audit" / "results" / "cap_launch.txt"


def scene(strong_first: bool, png: Path) -> str:
    zeros = [{"node": [0, 0, 0], "value": 0.0} for _ in range(1024)]
    strong = [{"node": [15, 10, 0], "value": 1.0e6}]
    items = strong + zeros if strong_first else zeros + strong
    payload = json.dumps(items, separators=(",", ":"))
    env = os.environ.copy()
    env["EM_PLACE"] = payload
    env["EM_CAPTURE"] = str(png)
    env["EM_LEVEL"] = "1"
    env["EM_MAP"] = "potential"
    env["EM_WAIT"] = "1"
    env["EM_FRAME"] = "30"
    proc = subprocess.run(
        [str(GAME)],
        cwd=ROOT,
        env=env,
        capture_output=True,
        text=True,
        timeout=90,
    )
    return (
        f"file {png.name} strong_first {strong_first} json_chars {len(payload)} "
        f"exit {proc.returncode}\n{proc.stderr}\n"
    )


def main() -> None:
    shots = ROOT / "audit" / "screenshots"
    parts = [
        scene(False, shots / "cap_strong_last.png"),
        scene(True, shots / "cap_strong_first.png"),
    ]
    LOG.write_text("".join(parts), encoding="utf-8")
    print("".join(parts))


if __name__ == "__main__":
    main()
