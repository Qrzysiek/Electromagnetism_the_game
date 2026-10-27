"""Measured values in the documents (DEC-22): a test records each number a document quotes
(`measured(id, value, digits)` in the physics tests, `tests/common`), to
`target/measured/<id>`, and the documents mark where they quote it:
`<!--m:<id>-->0.027<!--/m-->` (GitHub and the PDF show the number alone).

    python scripts/measured.py check     # quotes that differ from their record (exit 1)
    python scripts/measured.py update    # writes the recorded values into the documents

Run the tests first (those that measure are often `--ignored`: run them by name). A quote
without a record is listed, not an error: its test has not run in this checkout.
"""
import glob
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RECORDS = os.path.join(ROOT, "target", "measured")
MARK = re.compile(r"<!--m:([\w.-]+)-->(.*?)<!--/m-->")


def documents():
    return [os.path.join(ROOT, "PHYSICS.md")] + sorted(glob.glob(os.path.join(ROOT, "docs", "*.md")))


def recorded(key):
    """The recorded value, as the documents write it (a leading minus typeset as −)."""
    try:
        with open(os.path.join(RECORDS, key), encoding="utf-8") as f:
            text = f.read().strip()
    except FileNotFoundError:
        return None
    return "−" + text[1:] if text.startswith("-") else text


def quotes(text):
    """(id, quoted value) of every marker in a document."""
    return [(m.group(1), m.group(2)) for m in MARK.finditer(text)]


def check():
    differ, missing = [], []
    for path in documents():
        text = open(path, encoding="utf-8").read()
        for key, value in quotes(text):
            rec = recorded(key)
            if rec is None:
                missing.append(key)
            elif rec != value:
                differ.append(f"{os.path.relpath(path, ROOT)}: {key}: quotes {value}, measured {rec}")
    for line in differ:
        print(line)
    if missing:
        print(f"no record (their tests have not run here): {', '.join(sorted(set(missing)))}")
    print(f"{len(differ)} quotes differ from their record")
    return 1 if differ else 0


def update():
    changed = 0
    for path in documents():
        text = open(path, encoding="utf-8").read()

        def put(m):
            nonlocal changed
            rec = recorded(m.group(1))
            if rec is None or rec == m.group(2):
                return m.group(0)
            changed += 1
            print(f"{os.path.relpath(path, ROOT)}: {m.group(1)}: {m.group(2)} -> {rec}")
            return f"<!--m:{m.group(1)}-->{rec}<!--/m-->"
        new = MARK.sub(put, text)
        if new != text:
            with open(path, "w", encoding="utf-8", newline="\n") as f:
                f.write(new)
    print(f"{changed} quotes updated")
    return 0


if __name__ == "__main__":
    if len(sys.argv) != 2 or sys.argv[1] not in ("check", "update"):
        sys.exit(__doc__)
    sys.exit(check() if sys.argv[1] == "check" else update())
