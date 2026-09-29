#!/usr/bin/env python3
"""Typesets PHYSICS.md into PHYSICS.pdf.

PHYSICS.md stays the source. Pandoc converts it (with docs/physics/filter.lua) into
docs/physics/PHYSICS.tex, the body that docs/physics/main.tex sets with LuaLaTeX. The
compilation runs twice (for the contents) in a temporary directory, so only the PDF is
kept; missing glyphs are reported. Both outputs are generated (not committed).

Usage: python scripts/physics_pdf.py   (pandoc and LuaLaTeX from PATH, or the user-scope
installations of Pandoc and MiKTeX on Windows; override with $PANDOC and $LUALATEX).
"""

import os
import re
import shutil
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DOCS = os.path.join(ROOT, "docs", "physics")
LOCAL = os.environ.get("LOCALAPPDATA", "")


def tool(name, env, candidates):
    for path in [os.environ.get(env), shutil.which(name), *candidates]:
        if path and os.path.exists(path):
            return path
    sys.exit(f"{name} not found (set ${env})")


def main():
    pandoc = tool("pandoc", "PANDOC", [os.path.join(LOCAL, "Pandoc", "pandoc.exe")])
    lualatex = tool("lualatex", "LUALATEX", [
        os.path.join(LOCAL, "Programs", "MiKTeX", "miktex", "bin", "x64", "lualatex.exe")])
    subprocess.run([pandoc, os.path.join(ROOT, "PHYSICS.md"),
                    "-f", "gfm+tex_math_dollars+smart", "-t", "latex",
                    "--shift-heading-level-by=-1",
                    "--lua-filter", os.path.join(DOCS, "filter.lua"),
                    "-o", os.path.join(DOCS, "PHYSICS.tex")], check=True)
    with tempfile.TemporaryDirectory() as build:
        for _ in range(2):
            r = subprocess.run([lualatex, "-interaction=nonstopmode", "-halt-on-error",
                                f"-output-directory={build}", "main.tex"],
                               cwd=DOCS, capture_output=True, text=True, errors="replace")
            if r.returncode != 0:
                print(r.stdout[-3000:])
                sys.exit("LuaLaTeX failed")
        log = open(os.path.join(build, "main.log"), encoding="utf-8", errors="replace").read()
        missing = sorted(set(re.findall(r"Missing character: There is no (\S+)", log)))
        if missing:
            print("missing glyphs:", " ".join(missing))
        shutil.copy(os.path.join(build, "main.pdf"), os.path.join(ROOT, "PHYSICS.pdf"))
    pages = re.search(r"Output written on .*?\((\d+) pages", log)
    print(f"PHYSICS.pdf: {pages.group(1) if pages else '?'} pages")


if __name__ == "__main__":
    main()
