"""Count a few patterns in the Rust sources. Read-only.

    py -3 audit/scripts/static_scan.py

Writes ``audit/results/static_scan.txt``.
"""

from __future__ import annotations

from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "audit" / "results" / "static_scan.txt"

PATTERNS = (
    "unwrap()",
    "expect(",
    "todo!",
    "unimplemented!",
    "mul_add",
    "unsafe ",
    "f32",
)

ROOTS = (
    ROOT / "crates" / "physics" / "src",
    ROOT / "crates" / "level" / "src",
    ROOT / "crates" / "game" / "src",
    ROOT / "crates" / "generator" / "src",
)


def main() -> None:
    lines = ["Pattern counts in *.rs (comments included; this is a scan, not a verdict).", ""]
    for root in ROOTS:
        files = sorted(root.rglob("*.rs"))
        lines.append(f"## {root.relative_to(ROOT)}  ({len(files)} files)")
        totals = {p: 0 for p in PATTERNS}
        per_file: list[str] = []
        for path in files:
            text = path.read_text(encoding="utf-8", errors="replace")
            counts = {p: text.count(p) for p in PATTERNS}
            if any(counts.values()):
                shown = ", ".join(f"{p}={n}" for p, n in counts.items() if n)
                per_file.append(f"  {path.relative_to(ROOT)}: {shown}")
            for p, n in counts.items():
                totals[p] += n
        lines.append("  totals: " + ", ".join(f"{p}={n}" for p, n in totals.items()))
        lines.extend(per_file)
        lines.append("")
    text = "\n".join(lines)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(text, encoding="utf-8")
    print(text)
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
