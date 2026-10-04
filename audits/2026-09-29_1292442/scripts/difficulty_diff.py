"""Compare generator analyze rows with the difficulty table in docs/LEVELS.md.

Reads audit/results/analyze_curriculum.txt. A row still missing from that log
is counted as not yet measured, not as a mismatch.
"""

from __future__ import annotations

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def table_rows(text: str) -> dict[str, str]:
    rows = {}
    start = text.find("Measured 2026-09-28")
    body = text[start:]
    for line in body.splitlines():
        if not line.startswith("| ") or line.startswith("| level") or line.startswith("|---"):
            continue
        cells = [c.strip() for c in line.strip("|").split("|")]
        if len(cells) < 7 or not re.match(r"\d", cells[0]):
            continue
        rows[cells[0]] = " | ".join(cells[1:])
    return rows


def measured_rows(text: str) -> dict[str, str]:
    rows = {}
    for line in text.splitlines():
        if not line.startswith("| ") or line.startswith("| level") or line.startswith("|---"):
            continue
        cells = [c.strip() for c in line.strip("|").split("|")]
        if len(cells) < 7:
            continue
        rows[cells[0]] = " | ".join(cells[1:])
    return rows


def main() -> None:
    doc = table_rows((ROOT / "docs" / "LEVELS.md").read_text(encoding="utf-8"))
    got = measured_rows((ROOT / "audit" / "results" / "analyze_curriculum.txt").read_text(encoding="utf-8"))
    mismatch = 0
    for name, published in doc.items():
        if name not in got:
            continue
        if got[name] != published:
            mismatch += 1
            print(f"DIFF {name}")
            print(f"  doc  {published}")
            print(f"  new  {got[name]}")
    present = sum(1 for name in doc if name in got)
    print(f"published rows: {len(doc)}")
    print(f"remeasured rows: {len(got)}")
    print(f"published rows present: {present}")
    print(f"mismatches among remeasured: {mismatch}")


if __name__ == "__main__":
    main()
