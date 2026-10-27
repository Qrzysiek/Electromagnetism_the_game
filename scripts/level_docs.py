#!/usr/bin/env python3
"""The level documents: the physics of every level, derived from first principles.

docs/levels/README.md describes the system. A document lives in docs/levels/<slug>/
(<slug>: the level's name in levels/curriculum.json, stable when levels are renumbered):
  doc.tex          the text (\\documentclass{leveldoc}, docs/levels/leveldoc.cls)
  calc.py          its companion: computes every number the text prints from the level file
                   and compares them with the game's (`compute(doc)`, the Doc class below)
  game.json        the game's numbers, written by `game` (never edited by hand)
  placements.json  optional: further placements for `game` (the reference is always first):
                   a list of placements, or {"detailed": [...], "family": {"kind": "charge",
                   "values": [...]}}: a family is one element of each value on every node
                   (node-major: i, then j, then k, then the values), of which game.json keeps
                   one verdict letter each (`verdict_code`)
  NOTES.md         optional: questions and problems for the owner
  fig_*.pdf|png|jpg  optional figures (TikZ belongs in doc.tex)
The Foundations (docs/levels/foundations/) hold what every level builds on.

Usage:
  python scripts/level_docs.py status              every level: document, game data, review
  python scripts/level_docs.py check [slug ...]    fast checks, no TeX (all documents by default)
  python scripts/level_docs.py calc slug           run a companion: its values and comparisons
  python scripts/level_docs.py build [slug ...]    check, then typeset into docs/levels/pdf/
  python scripts/level_docs.py game slug [...]     record the game's numbers (generator scan)
  python scripts/level_docs.py new slug            start a document from docs/levels/_template

Typesetting needs LuaLaTeX (from PATH, $LUALATEX, or the user-scope MiKTeX on Windows); `game`
needs the generator (cargo build --release -p generator; $EM_GENERATOR overrides its path).
The PDFs (docs/levels/pdf/) and the build directories (target/level-docs/) are generated.
"""

import ast
import copy
import datetime
import hashlib
import importlib.util
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
import traceback

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DOCS = os.path.join(ROOT, "docs", "levels")
PDF_DIR = os.path.join(DOCS, "pdf")
BUILD = os.path.join(ROOT, "target", "level-docs")
LEVELS = os.path.join(ROOT, "levels")
TEMPLATE = os.path.join(DOCS, "_template")
FOUNDATIONS = "foundations"
GENERATOR = os.environ.get("EM_GENERATOR") or os.path.join(ROOT, "target", "release", "generator")
LOCAL = os.environ.get("LOCALAPPDATA", "")

# The sections of a level document, in this order (the Foundations are free-form).
SECTIONS = [
    "What the level asks",
    "Setup",
    "Principles and model",
    "Derivation",
    "The reference solution by hand",
    "Comparison with the game",
    "Other solutions",
    "What the game computes",
    "Further reading",
]
# Generated material each level document must show (macro: the section it belongs in).
REQUIRED_MACROS = {
    r"\leveldescription": "What the level asks",
    r"\setuptable": "Setup",
    r"\comparisontable": "Comparison with the game",
    r"\levelfileappendix": None,
}
# Not allowed in doc.tex: packages, file and Lua access, catcodes, redefinitions, and the
# class's internals (scripts/level_docs.py generates those).
FORBIDDEN = [
    r"\usepackage", r"\RequirePackage", r"\input", r"\include", r"\write", r"\immediate",
    r"\openout", r"\openin", r"\directlua", r"\luaexec", r"\luadirect", r"\catcode",
    r"\makeatletter", r"\makeatother", r"\renewcommand", r"\providecommand", r"\def",
    r"\gdef", r"\edef", r"\xdef", r"\let", r"\hypersetup", r"\externaldocument",
    r"\DocEntry", r"\DocValue", r"\ExplSyntaxOn", r"\csname", r"\expandafter",
    r"\special", r"\pdfobj", r"\shellescape",
]
LABEL = re.compile(r"^(eq|sec|fig|tab|as|st):[a-z0-9][a-z0-9-]*$")
# Modules a companion may import (deterministic numerics; no files, processes or network).
ALLOWED_IMPORTS = {"math", "cmath", "fractions", "decimal", "itertools", "functools",
                   "statistics", "numpy", "scipy", "mpmath", "sympy"}
FORBIDDEN_CALLS = {"open", "exec", "eval", "compile", "__import__", "input", "breakpoint",
                   "globals", "locals", "vars", "setattr", "delattr"}
FIGURE = re.compile(r"^fig_[a-z0-9_-]+\.(pdf|png|jpg)$")
DOC_FILES = {"doc.tex", "calc.py", "game.json", "placements.json", "NOTES.md"}

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
# Importing the companions must not leave __pycache__ folders in docs/levels/.
sys.dont_write_bytecode = True


class DocError(Exception):
    pass


# ---- The curriculum and the level files ---------------------------------------------------

def load_json(path):
    with open(path, encoding="utf-8") as f:
        return json.load(f)


def curriculum():
    """Every level in curriculum order: dicts with slug, number, title, arc, tier, path."""
    files = {}
    for name in os.listdir(LEVELS):
        m = re.match(r"^(\d+)_(.+)\.json$", name)
        if m and name != "golden_hashes.json":
            files[m.group(2)] = (int(m.group(1)), os.path.join(LEVELS, name))
    out = []
    for arc in load_json(os.path.join(LEVELS, "curriculum.json"))["arcs"]:
        for tier in arc["tiers"]:
            for slug in tier["levels"]:
                if slug not in files:
                    raise DocError(f"curriculum lists {slug}, but levels/ has no file for it")
                number, path = files[slug]
                out.append({"slug": slug, "number": number, "arc": arc["name"],
                            "tier": tier["name"], "path": path,
                            "title": load_json(path)["name"]})
    out.sort(key=lambda e: e["number"])
    return out


def level_sha(path):
    """The level file's fingerprint, 16 digits of SHA-256: everything a document uses, so
    without what the level represents (`reality`) and its electrodes' and metals' labels,
    which no document shows (rewording them leaves the documents current)."""
    level = load_json(path)
    level.pop("reality", None)
    for e in level.get("electrodes", []):
        e.pop("label", None)
    for m in (level.get("tube") or {}).get("metals", []):
        m.pop("label", None)
    data = json.dumps(level, sort_keys=True, ensure_ascii=False).encode("utf-8")
    return hashlib.sha256(data).hexdigest()[:16]


def golden_hash(path):
    stem = os.path.splitext(os.path.basename(path))[0]
    return load_json(os.path.join(LEVELS, "golden_hashes.json")).get(stem)


def documents():
    """Slugs with a document (a folder with doc.tex), Foundations first."""
    out = []
    for name in sorted(os.listdir(DOCS)):
        if name.startswith("_") or name == "pdf":
            continue
        if os.path.isfile(os.path.join(DOCS, name, "doc.tex")):
            out.append(name)
    if FOUNDATIONS in out:
        out.remove(FOUNDATIONS)
        out.insert(0, FOUNDATIONS)
    return out


# ---- TeX helpers ----------------------------------------------------------------------------

def tex_escape(s):
    rep = {"\\": r"\textbackslash{}", "{": r"\{", "}": r"\}", "$": r"\$", "&": r"\&",
           "#": r"\#", "_": r"\_", "%": r"\%", "~": r"\textasciitilde{}",
           "^": r"\textasciicircum{}"}
    return "".join(rep.get(ch, ch) for ch in s)


def tex_text_with_math(s):
    """A level's text: math between $...$ kept, the rest escaped."""
    parts = re.split(r"(\$[^$]*\$)", s)
    return "".join(p if p.startswith("$") and p.endswith("$") and len(p) > 1 else tex_escape(p)
                   for p in parts)


def strip_comments(tex):
    """The source without TeX comments (an unescaped % to the end of the line)."""
    return "\n".join(re.sub(r"(?<!\\)%.*$", "", line) for line in tex.split("\n"))


# ---- Numbers ----------------------------------------------------------------------------------

def fmt_number(x, digits=10, fixed=None, keep_zeros=False, sci=False):
    """A number as siunitx \\num{...}: `digits` significant digits (trailing zeros dropped
    unless `keep_zeros`), `fixed` decimals, or with `sci` always as a power of ten."""
    if isinstance(x, bool):
        raise DocError("a value must be a number, not a bool")
    if isinstance(x, int) and not sci:
        return rf"\num{{{x}}}"
    try:
        import mpmath
        is_mp = isinstance(x, mpmath.mpf)
    except ImportError:
        is_mp = False
    if is_mp:
        if not mpmath.isfinite(x):
            raise DocError(f"a value is not finite: {x}")
        if fixed is not None:
            x = float(x)
        elif sci:
            s = mpmath.nstr(x, digits, min_fixed=0, max_fixed=0, strip_zeros=not keep_zeros)
        else:
            s = mpmath.nstr(x, digits, min_fixed=-5, max_fixed=9, strip_zeros=not keep_zeros)
    if not is_mp or fixed is not None:
        x = float(x)
        if x != x or x in (float("inf"), float("-inf")):
            raise DocError(f"a value is not finite: {x}")
        if fixed is not None:
            s = f"{x:.{fixed}f}"
        elif sci:
            s = f"{x:.{digits - 1}e}"
        else:
            s = f"{x:#.{digits}g}" if keep_zeros else f"{x:.{digits}g}"
    mant, _, exp = s.partition("e")
    if not keep_zeros and fixed is None and "." in mant:
        mant = mant.rstrip("0").rstrip(".")
    s = mant + (f"e{int(exp)}" if exp else "")
    return rf"\num{{{s}}}"


def as_float(x):
    try:
        return float(x)
    except Exception as e:
        raise DocError(f"not a number: {x!r}") from e


class Doc:
    """What a companion (calc.py) sees: `compute(doc)` reads the level and the game's numbers
    and defines every number of the text.

    doc.level             the level file (dict)
    doc.value(name, x, digits=10, fixed=None, keep_zeros=False, sci=False)   \\val{name} in
                           the text: significant digits, fixed decimals, or a power of ten
    doc.text(name, latex)  a value that is LaTeX text (a word, a formula)
    doc.table(name, header, rows, align=None, digits=10, sci=False)  \\val{name}: a booktabs
                           tabular (digits, sci: one value or a list per column)
    doc.data(name, rows)   writes name.dat (whitespace-separated columns) for pgfplots
    doc.check(ok, message) an internal consistency check of the derivation (must hold)
    doc.result(i=0)        the game's verdict for placement i (0: the reference), game.json
    doc.n_results()        how many placements game.json holds in detail
    doc.flight(i=0, flight=0, run="verified")   one flight of it (run: verified or preview)
    doc.compare(name, what, ours, game, tol, rel=False, digits=13)   a row of the comparison
                           table; strings compare for equality. The game's numbers may enter
                           only here (checked: the companion is run again with them changed)
    doc.family()           the placements of placements.json's family, in order
    doc.compare_verdicts(name, what, ours)   ours: one letter per family placement (A, L, C,
                           T, R, S, X as in verdict_code); every placement the game verified
                           must agree
    doc.verdicts()         the game's verdict letters for the family, for a map the text
                           presents as the game's (with doc.game_value for its counts)
    doc.game_value(name, x, ...)  shows a game number in the text (not one of the document's)
    """

    def __init__(self, level, game):
        self.level = level
        self._game = game
        self.values = {}
        self.raw = {}
        self.game_values = set()
        self.rows = []
        self.datas = {}
        self.checks = []

    def _define(self, name, tex):
        if not re.match(r"^[A-Za-z][A-Za-z0-9]*$", name):
            raise DocError(f"value name {name!r}: letters and digits only (LaTeX \\val{{...}})")
        if name in self.values:
            raise DocError(f"value {name} defined twice")
        self.values[name] = tex

    def value(self, name, x, digits=10, fixed=None, keep_zeros=False, sci=False):
        self._define(name, fmt_number(x, digits, fixed, keep_zeros, sci))
        self.raw[name] = repr(float(x)) if not isinstance(x, int) else str(x)
        return x

    def text(self, name, latex):
        self._define(name, str(latex))
        return latex

    def table(self, name, header, rows, align=None, digits=10, sci=False):
        """digits and sci: one for every column, or a list per column."""
        align = align or ("l" + "r" * (len(header) - 1))
        per = lambda opt, j: opt[j] if isinstance(opt, (list, tuple)) else opt
        cell = lambda v, j: v if isinstance(v, str) else fmt_number(v, per(digits, j), sci=per(sci, j))
        body = [r"\begin{tabular}{@{}" + align + r"@{}}", r"\toprule",
                " & ".join(header) + r" \\", r"\midrule"]
        body += [" & ".join(cell(v, j) for j, v in enumerate(r)) + r" \\" for r in rows]
        body += [r"\bottomrule", r"\end{tabular}"]
        self._define(name, "\n".join(body))

    def data(self, name, rows):
        if not re.match(r"^[a-z][a-z0-9_]*$", name):
            raise DocError(f"data name {name!r}: lower-case letters, digits, underscores")
        self.datas[name] = "\n".join(" ".join(repr(float(v)) for v in r) for r in rows) + "\n"

    def check(self, ok, message):
        self.checks.append((bool(ok), message))
        if not ok:
            raise DocError(f"check failed: {message}")

    def result(self, i=0):
        res = self._game.get("results") if self._game else None
        if not res:
            raise DocError("no game numbers: run python scripts/level_docs.py game <slug>")
        if i >= len(res):
            raise DocError(f"game.json has {len(res)} placements, not {i + 1}")
        return res[i]

    def n_results(self):
        """How many placements game.json holds in detail (the reference, then placements.json's)."""
        return len((self._game or {}).get("results", []))

    def flight(self, i=0, flight=0, run="verified"):
        r = self.result(i)
        if "refused" in r:
            raise DocError(f"placement {i} was refused: {r['refused']}")
        f = r["flights"][flight]
        return f[run] | {"outcome": f["outcome"], "status": f["status"]} if run in f else f

    def family(self):
        spec = (self._game or {}).get("family")
        if not spec:
            raise DocError("game.json has no family: list one in placements.json, then run game")
        return family_placements(self.level, spec)

    def verdicts(self):
        """The game's verdict letters for the family (verdict_code): for a map that the text
        presents as the game's, when the document does not derive it."""
        v = (self._game or {}).get("verdicts")
        if v is None:
            raise DocError("game.json has no verdicts: list a family in placements.json")
        return v

    def compare_verdicts(self, name, what, ours):
        game = (self._game or {}).get("verdicts")
        if not game:
            raise DocError("game.json has no verdicts: list a family in placements.json")
        if len(ours) != len(game):
            raise DocError(f"{len(ours)} verdicts for a family of {len(game)} placements")
        verified = [k for k, g in enumerate(game) if g.isupper() or g == "X"]
        wrong = [k for k in verified if ours[k] != game[k]]
        ok = not wrong
        self.rows.append({"name": name, "what": what,
                          "ours": rf"\num{{{len(verified) - len(wrong)}}} agree",
                          "game": rf"\num{{{len(verified)}}} verified",
                          "diff": rf"\num{{{len(wrong)}}}", "tol": r"\num{0}", "ok": ok,
                          "ours_raw": ours})
        self._define(f"cmp{name}", rf"\num{{{len(verified) - len(wrong)}}}")
        self._define(f"cmp{name}Game", rf"\num{{{len(verified)}}}")
        self._define(f"cmp{name}Total", rf"\num{{{len(game)}}}")
        self._define(f"cmp{name}Diff", rf"\num{{{len(wrong)}}}")
        if wrong:
            fam = self.family()
            self.checks.append((False, "verdicts differ at " + ", ".join(
                f"{fam[k][0]['node']} {fam[k][0]['value']}: ours {ours[k]}, game {game[k]}"
                for k in wrong[:10])))
        return ok

    def game_value(self, name, x, digits=10, fixed=None, keep_zeros=False, sci=False):
        self.game_values.add(name)
        return self.value(name, x, digits, fixed, keep_zeros, sci)

    def compare(self, name, what, ours, game, tol=0.0, rel=False, digits=13):
        if isinstance(ours, str) or isinstance(game, str):
            ok = str(ours) == str(game)
            row = {"name": name, "what": what, "ours": tex_escape(str(ours)),
                   "game": tex_escape(str(game)), "diff": "", "tol": "equal", "ok": ok,
                   "ours_raw": str(ours)}
        else:
            o, g = as_float(ours), as_float(game)
            diff = abs(o - g)
            if rel:
                diff = diff / abs(g) if g != 0 else float("inf")
            ok = diff <= tol
            row = {"name": name, "what": what, "ours": fmt_number(ours, digits),
                   "game": fmt_number(g, digits),
                   "diff": fmt_number(diff, 2, sci=diff != 0) if diff == diff else "nan",
                   "tol": fmt_number(tol, 2, sci=tol != 0) + (" (rel.)" if rel else ""), "ok": ok,
                   "ours_raw": fmt_number(ours, digits)}
            self._define(f"cmp{name}", row["ours"])
            self._define(f"cmp{name}Game", row["game"])
            self._define(f"cmp{name}Diff", row["diff"])
        self.rows.append(row)
        return ok


def family_placements(level, spec):
    """The placements of a family: one element of each value on every node of the grid
    (node-major: i, then j, then k, then the values)."""
    g = level["grid"]
    s = g.get("subdivision", 1)
    out = []
    for i in range(g["nx"] * s + 1):
        for j in range(g["ny"] * s + 1):
            for k in range(g.get("nz", 0) * s + 1):
                for v in spec["values"]:
                    out.append([{"node": [i, j, k], "kind": spec["kind"], "value": v}])
    return out


# One letter per verdict of a family: upper case when the game verified it.
CODES = {"Arrived": "A", "Rejected": "R", "SkippedGate": "S", "Collided": "C",
         "LeftBounds": "L", "Timeout": "T", "Failed": "F"}


def verdict_code(result):
    """A placement's verdict as one letter (its first flight): A arrived, L left the world,
    C collided, T timed out, R rejected by the detector's conditions, S skipped a gate,
    F failed; lower case if not verified; X refused (not a legal placement)."""
    if "refused" in result:
        return "X"
    f = result["flights"][0]
    code = CODES[f["outcome"].split("(")[0]]
    return code if f["status"] == "Verified" else code.lower()


def read_placements(slug):
    """placements.json as (detailed placements, family spec or None)."""
    path = os.path.join(DOCS, slug, "placements.json")
    if not os.path.isfile(path):
        return [], None
    p = load_json(path)
    if isinstance(p, list):
        return p, None
    return p.get("detailed", []), p.get("family")


def perturbed(game):
    """game.json with every number the game computed changed (the independence check); the
    placements asked about stay as they are."""
    def walk(v):
        if isinstance(v, bool) or isinstance(v, int):
            return v
        if isinstance(v, float):
            return v * 1.001 + 1e-3
        if isinstance(v, list):
            return [walk(x) for x in v]
        if isinstance(v, dict):
            return {k: (x if k in ("placement", "family") else walk(x)) for k, x in v.items()}
        return v
    return walk(copy.deepcopy(game))


def check_calc_source(path):
    """Imports and calls a companion may not use."""
    problems = []
    tree = ast.parse(open(path, encoding="utf-8").read(), path)
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for a in node.names:
                if a.name.split(".")[0] not in ALLOWED_IMPORTS:
                    problems.append(f"calc.py imports {a.name} (allowed: {sorted(ALLOWED_IMPORTS)})")
        elif isinstance(node, ast.ImportFrom):
            if (node.module or "").split(".")[0] not in ALLOWED_IMPORTS or node.level:
                problems.append(f"calc.py imports from {node.module} (allowed: {sorted(ALLOWED_IMPORTS)})")
        elif isinstance(node, ast.Call) and isinstance(node.func, ast.Name):
            if node.func.id in FORBIDDEN_CALLS:
                problems.append(f"calc.py calls {node.func.id}()")
        elif isinstance(node, ast.Attribute) and node.attr.startswith("__") and node.attr != "__init__":
            problems.append(f"calc.py uses {node.attr}")
    if not any(isinstance(n, ast.FunctionDef) and n.name == "compute" for n in tree.body):
        problems.append("calc.py must define compute(doc)")
    return problems


def run_calc(slug, level, game):
    """Runs docs/levels/<slug>/calc.py; returns the Doc. Raises DocError on failure."""
    path = os.path.join(DOCS, slug, "calc.py")
    spec = importlib.util.spec_from_file_location(f"calc_{slug}", path)
    module = importlib.util.module_from_spec(spec)
    try:
        spec.loader.exec_module(module)
        doc = Doc(level, game)
        t0 = time.time()
        module.compute(doc)
        if time.time() - t0 > 120:
            raise DocError(f"calc.py took {time.time() - t0:.0f} s (limit 120 s)")
        return doc
    except DocError:
        raise
    except Exception:
        raise DocError("calc.py failed:\n" + traceback.format_exc(limit=6))


def calc_with_independence(slug, level, game):
    """The companion's Doc, after checking that its own numbers ignore the game's."""
    doc = run_calc(slug, level, game)
    if game:
        try:
            twin = run_calc(slug, level, perturbed(game))
        except DocError as e:
            raise DocError("calc.py fails when the game's numbers change: it must use them only "
                           "in doc.compare (or doc.game_value):\n" + str(e))
        changed = [n for n, v in doc.values.items()
                   if n not in doc.game_values and not n.startswith("cmp")
                   and twin.values.get(n) != v]
        changed += [r["name"] for r, s in zip(doc.rows, twin.rows) if r["ours_raw"] != s["ours_raw"]]
        if changed:
            raise DocError("these numbers change with the game's: compute them from the level "
                           "file, not from game.json: " + ", ".join(changed))
    return doc


# ---- Checks -----------------------------------------------------------------------------------

def parse_uses(tex):
    m = re.findall(r"\\DocUses\{([^}]*)\}", tex)
    if len(m) > 1:
        raise DocError("\\DocUses appears more than once")
    return [u.strip() for u in m[0].split(",") if u.strip()] if m else []


def check_doc(slug, levels, verbose=True):
    """All checks without TeX. Returns (errors, warnings, doc or None)."""
    errors, warnings = [], []
    by_slug = {e["slug"]: e for e in levels}
    order = {e["slug"]: i for i, e in enumerate(levels)}
    folder = os.path.join(DOCS, slug)
    is_level = slug != FOUNDATIONS
    if is_level and slug not in by_slug:
        return [f"{slug}: not a level of levels/curriculum.json"], [], None
    for name in os.listdir(folder):
        if name not in DOC_FILES and not FIGURE.match(name) and name != "__pycache__":
            errors.append(f"unexpected file docs/levels/{slug}/{name} (allowed: {sorted(DOC_FILES)}, fig_*.pdf|png|jpg)")
    raw = open(os.path.join(folder, "doc.tex"), encoding="utf-8").read()
    tex = strip_comments(raw)
    for name in ("doc.tex", "calc.py"):
        path = os.path.join(folder, name)
        if os.path.isfile(path) and "TEMPLATE" in open(path, encoding="utf-8").read():
            errors.append(f"{name} still has the template's TEMPLATE placeholders")
    first = next((l.strip() for l in tex.split("\n") if l.strip()), "")
    if first != r"\documentclass{leveldoc}":
        errors.append(r"doc.tex must begin with \documentclass{leveldoc}")
    for f in FORBIDDEN:
        if re.search(re.escape(f) + r"(?![A-Za-z@])", tex):
            errors.append(f"doc.tex uses {f} (not allowed: ask for what you need in the level's NOTES.md)")
    try:
        uses = parse_uses(tex)
    except DocError as e:
        errors.append(str(e))
        uses = []
    if is_level:
        if FOUNDATIONS not in uses:
            errors.append(r"\DocUses must list foundations")
        for u in uses:
            if u == FOUNDATIONS:
                continue
            if u not in by_slug:
                errors.append(f"\\DocUses lists {u}: not a level")
            elif order[u] >= order[slug]:
                errors.append(f"\\DocUses lists {u}: only earlier levels may be cited")
            elif not os.path.isfile(os.path.join(DOCS, u, "doc.tex")):
                errors.append(f"\\DocUses lists {u}, which has no document yet")
    elif uses:
        errors.append("the Foundations cite no other document")
    # References.
    for kind, s in re.findall(r"\\(levelref|levelnum)\{([^}]*)\}", tex):
        if s != FOUNDATIONS and s not in by_slug:
            errors.append(f"\\{kind}{{{s}}}: unknown level")
    for kind, s, _ in re.findall(r"\\(xeqref|xref|xsecref)\{([^}]*)\}\{([^}]*)\}", tex):
        if s not in uses:
            errors.append(f"\\{kind}{{{s}}}{{...}}: list {s} in \\DocUses")
    labels = re.findall(r"\\label\{([^}]*)\}", tex)
    for l in labels:
        if not LABEL.match(l):
            errors.append(f"label {l!r}: use eq:, sec:, fig:, tab:, as: or st: and lower-case words")
    dup = sorted({l for l in labels if labels.count(l) > 1})
    if dup:
        errors.append("labels defined twice: " + ", ".join(dup))
    # Structure.
    if is_level:
        found = re.findall(r"\\section\*?\{([^}]*)\}", tex)
        if found != SECTIONS:
            errors.append("sections must be exactly, in order: " + " | ".join(SECTIONS)
                          + "\n      found: " + " | ".join(found))
        for macro, section in REQUIRED_MACROS.items():
            if not re.search(re.escape(macro) + r"(?![A-Za-z])", tex):
                errors.append(f"doc.tex must use {macro}" + (f" (in '{section}')" if section else ""))
        fp = re.findall(r"\\DocFingerprint\{([^}]*)\}", tex)
        current = level_sha(by_slug[slug]["path"])
        if len(fp) != 1:
            errors.append(r"\DocFingerprint{...} must appear once (the level file's fingerprint)")
        elif fp[0] != current:
            warnings.append(f"the level changed since the text was reviewed: review it, then set "
                            f"\\DocFingerprint{{{current}}}")
    # The companion.
    doc = None
    calc = os.path.join(folder, "calc.py")
    game = None
    if is_level:
        gpath = os.path.join(folder, "game.json")
        if not os.path.isfile(gpath):
            errors.append(f"no game.json: run python scripts/level_docs.py game {slug}")
        else:
            game = load_json(gpath)
            entry = by_slug[slug]
            if game.get("level_sha256") != level_sha(entry["path"]) or game.get("golden_hash") != golden_hash(entry["path"]):
                errors.append(f"game.json is stale (the level or its physics changed): run "
                              f"python scripts/level_docs.py game {slug}")
    if os.path.isfile(calc):
        problems = check_calc_source(calc)
        errors += problems
        if not problems:
            level = load_json(by_slug[slug]["path"]) if is_level else {}
            try:
                doc = calc_with_independence(slug, level, game)
                for r in doc.rows:
                    if not r["ok"]:
                        errors.append(f"comparison {r['name']} fails: ours {r['ours']}, game {r['game']}, "
                                      f"difference {r['diff']} > tolerance {r['tol']}")
                if is_level and not doc.rows:
                    errors.append("calc.py makes no comparison with the game (doc.compare)")
            except DocError as e:
                errors.append(str(e))
    elif is_level:
        errors.append("no calc.py")
    if doc is not None:
        used = set(re.findall(r"\\val\{([^}]*)\}", tex))
        missing = sorted(used - set(doc.values))
        if missing:
            errors.append("\\val of values calc.py does not define: " + ", ".join(missing))
        for name in re.findall(r"table\s*(?:\[[^\]]*\])?\s*\{([a-z0-9_]+)\.dat\}", tex):
            if name not in doc.datas:
                errors.append(f"{name}.dat is used but calc.py writes no doc.data('{name}', ...)")
    numbers = re.findall(r"(?<![\\A-Za-z_{])\d+\.\d{4,}", re.sub(r"\\(label|ref|eqref|xeqref|xref|cite)\{[^}]*\}", "", tex))
    if numbers:
        warnings.append("numbers with 4+ decimals typed in the text (use \\val from calc.py): "
                        + ", ".join(sorted(set(numbers))[:8]))
    return errors, warnings, doc


# ---- Generated inputs ----------------------------------------------------------------------------

def levelinfo(slug, levels):
    by_slug = {e["slug"]: e for e in levels}
    have = set(documents())
    lines = []
    if slug == FOUNDATIONS:
        lines += [r"\def\DocNumber{F}", r"\def\DocTitle{Foundations}",
                  r"\def\DocArc{}", r"\def\DocTier{}", r"\def\DocLevelFile{}"]
    else:
        e = by_slug[slug]
        rel = os.path.relpath(e["path"], ROOT).replace("\\", "/")
        lines += [rf"\def\DocNumber{{{e['number']}}}",
                  rf"\def\DocTitle{{{tex_escape(e['title'])}}}",
                  rf"\def\DocArc{{{tex_escape(e['arc'])}}}",
                  rf"\def\DocTier{{{tex_escape(e['tier'])}}}",
                  rf"\def\DocLevelFile{{{tex_escape(rel)}}}"]
    lines.append(rf"\def\DocSlugText{{{tex_escape(slug)}}}")
    lines.append(rf"\DocEntry{{{FOUNDATIONS}}}{{F}}{{Foundations}}{{{int(FOUNDATIONS in have)}}}")
    for e in levels:
        lines.append(rf"\DocEntry{{{e['slug']}}}{{{e['number']}}}{{{tex_escape(e['title'])}}}{{{int(e['slug'] in have)}}}")
    return "\n".join(lines) + "\n"


def setup_tex(level):
    """The level file as a table: every field, grouped by its top-level key."""
    skip = {"name", "description", "format_version", "engine_version"}
    rows = []

    def short(v):
        return json.dumps(v, ensure_ascii=False)

    def walk(path, v):
        if isinstance(v, dict):
            if not v:
                rows.append((path, "{}"))
            for k, x in v.items():
                walk(f"{path}.{k}", x)
        elif isinstance(v, list) and any(isinstance(x, (dict, list)) for x in v) and not all(
                isinstance(x, list) and all(not isinstance(y, (dict, list)) for y in x) for x in v):
            if not v:
                rows.append((path, "[]"))
            for i, x in enumerate(v):
                walk(f"{path}[{i}]", x)
        else:
            rows.append((path, short(v)))

    groups = []
    for key, v in level.items():
        if key in skip:
            continue
        rows = []
        walk(key, v)
        groups.append(rows)
    out = [r"\begingroup\small",
           r"\begin{longtable}{@{}>{\raggedright\arraybackslash}p{0.42\linewidth}>{\raggedright\arraybackslash}p{0.54\linewidth}@{}}",
           r"\toprule Field & Value \\ \midrule \endhead",
           r"\bottomrule \endfoot"]
    for gi, rows in enumerate(groups):
        if gi:
            out.append(r"\midrule")
        for p, v in rows:
            out.append(rf"\texttt{{{tex_escape(p)}}} & \texttt{{{tex_escape(v)}}} \\")
    out += [r"\end{longtable}", r"\endgroup"]
    return "\n".join(out) + "\n"


def comparison_tex(doc, game):
    when = (game or {}).get("recorded", "?")
    out = [r"\begin{table}[htbp]", r"\centering\small",
           r"\begin{tabular}{@{}>{\raggedright\arraybackslash}p{0.27\linewidth}rrrrc@{}}",
           r"\toprule", r"Quantity & This document & The game & Difference & Tolerance & \\",
           r"\midrule"]
    for r in doc.rows:
        mark = r"{\color{green!50!black}\checkmark}" if r["ok"] else r"{\color{red!70!black}$\times$}"
        out.append(f"{r['what']} & {r['ours']} & {r['game']} & {r['diff']} & {r['tol']} & {mark} \\\\")
    out += [r"\bottomrule", r"\end{tabular}",
            rf"\caption{{This document's numbers against the game's (\code{{game.json}}, recorded "
            rf"{tex_escape(when)} with \code{{generator scan}}: the verified run). A row passes when the "
            r"difference is within its tolerance, set in \code{calc.py} from the error estimates "
            r"of this document and of the game.}", r"\label{tab:comparison}", r"\end{table}"]
    return "\n".join(out) + "\n"


def values_tex(doc):
    return ("".join(rf"\DocValue{{{n}}}{{{v}}}" + "\n" for n, v in doc.values.items())
            + "".join(rf"\DocRawValue{{{n}}}{{{v}}}" + "\n" for n, v in doc.raw.items()))


# ---- Typesetting -----------------------------------------------------------------------------------

def tool_lualatex():
    for p in [os.environ.get("LUALATEX"), shutil.which("lualatex"),
              os.path.join(LOCAL, "Programs", "MiKTeX", "miktex", "bin", "x64", "lualatex.exe")]:
        if p and os.path.exists(p):
            return p
    raise DocError("lualatex not found (install MiKTeX or TeX Live, or set $LUALATEX)")


def typeset(slug, levels, doc, force=False):
    """Compiles one document into docs/levels/pdf/<slug>.pdf. Returns a list of problems."""
    by_slug = {e["slug"]: e for e in levels}
    folder = os.path.join(DOCS, slug)
    build = os.path.join(BUILD, slug)
    os.makedirs(build, exist_ok=True)
    files = {"leveldoc.cls": open(os.path.join(DOCS, "leveldoc.cls"), encoding="utf-8").read(),
             "doc.tex": open(os.path.join(folder, "doc.tex"), encoding="utf-8").read(),
             "levelinfo.tex": levelinfo(slug, levels)}
    game = None
    if slug != FOUNDATIONS:
        e = by_slug[slug]
        level = load_json(e["path"])
        gpath = os.path.join(folder, "game.json")
        game = load_json(gpath) if os.path.isfile(gpath) else None
        files["setup.tex"] = setup_tex(level)
        files["description.tex"] = (r"\begin{quote}\itshape " + tex_text_with_math(level.get("description", ""))
                                    + r"\end{quote}" + "\n")
        files["level.json"] = open(e["path"], encoding="utf-8").read()
    if doc is not None:
        files["values.tex"] = values_tex(doc)
        if doc.rows:
            files["comparison.tex"] = comparison_tex(doc, game)
        for name, text in doc.datas.items():
            files[f"{name}.dat"] = text
    for name in os.listdir(folder):
        if FIGURE.match(name):
            shutil.copy(os.path.join(folder, name), os.path.join(build, name))
    for name, text in files.items():
        with open(os.path.join(build, name), "w", encoding="utf-8", newline="\n") as f:
            f.write(text)
    # Up to date? The inputs, the figures and the cited documents' labels.
    h = hashlib.sha256()
    for name in sorted(files):
        h.update(name.encode() + b"\0" + files[name].encode("utf-8") + b"\0")
    for name in sorted(os.listdir(folder)):
        if FIGURE.match(name):
            h.update(open(os.path.join(folder, name), "rb").read())
    for u in parse_uses(strip_comments(files["doc.tex"])):
        aux = os.path.join(BUILD, u, "doc.aux")
        if os.path.isfile(aux):
            h.update(open(aux, "rb").read())
    stamp = os.path.join(build, "stamp")
    pdf = os.path.join(PDF_DIR, f"{slug}.pdf")
    if not force and os.path.isfile(pdf) and os.path.isfile(stamp) and open(stamp).read() == h.hexdigest():
        return []
    lualatex = tool_lualatex()
    log = ""
    for _ in range(4):
        r = subprocess.run([lualatex, "-interaction=nonstopmode", "-halt-on-error", "-file-line-error",
                            "doc.tex"], cwd=build, capture_output=True, text=True, errors="replace",
                           timeout=600)
        log = open(os.path.join(build, "doc.log"), encoding="utf-8", errors="replace").read() \
            if os.path.isfile(os.path.join(build, "doc.log")) else r.stdout
        if r.returncode != 0:
            errs = [l for l in log.split("\n") if l.startswith("!") or ":error:" in l or re.match(r"^\S+\.tex:\d+:", l)]
            return [f"LuaLaTeX failed ({build}\\doc.log):\n      " + "\n      ".join(errs[:12] or log.split("\n")[-30:])]
        if not re.search(r"Rerun to get|Label\(s\) may have changed|rerunfilecheck", log):
            break
    problems = []
    for m in re.findall(r"LaTeX Warning: (Reference|Hyper reference) `([^']+)' on page \d+ undefined", log):
        problems.append(f"undefined reference {m[1]}")
    if "There were undefined references" in log and not problems:
        problems.append("undefined references (see doc.log)")
    for m in re.findall(r"Package xr(?:-hyper)? Warning: No file (\S+)", log):
        problems.append(f"a cited document is not built: {m}")
    missing = sorted(set(re.findall(r"Missing character: There is no (\S+)", log)))
    if missing:
        problems.append("missing glyphs: " + " ".join(missing))
    overfull = [float(x) for x in re.findall(r"Overfull \\hbox \(([\d.]+)pt too wide\)", log)]
    if any(x > 20 for x in overfull):
        problems.append(f"{sum(x > 20 for x in overfull)} lines stick out of the margin by more than 20 pt")
    if problems:
        return problems
    os.makedirs(PDF_DIR, exist_ok=True)
    shutil.copy(os.path.join(build, "doc.pdf"), pdf)
    with open(stamp, "w") as f:
        f.write(h.hexdigest())
    return []


# ---- Commands ----------------------------------------------------------------------------------------

def select(slugs, levels):
    have = documents()
    if not slugs:
        return have
    for s in slugs:
        if s not in have:
            raise DocError(f"{s}: no document (python scripts/level_docs.py new {s})")
    return slugs


def cmd_check(slugs, quiet=False):
    levels = curriculum()
    failed = 0
    for slug in select(slugs, levels):
        errors, warnings, _ = check_doc(slug, levels)
        status = "ok" if not errors else "FAILED"
        if errors:
            failed += 1
        if errors or warnings or not quiet:
            print(f"{slug}: {status}")
        for e in errors:
            print(f"  error: {e}")
        for w in warnings:
            print(f"  warning: {w}")
    return failed


def cmd_calc(slug):
    levels = curriculum()
    by_slug = {e["slug"]: e for e in levels}
    game = None
    gpath = os.path.join(DOCS, slug, "game.json")
    if os.path.isfile(gpath):
        game = load_json(gpath)
    level = load_json(by_slug[slug]["path"]) if slug != FOUNDATIONS else {}
    doc = run_calc(slug, level, game)
    for n, v in doc.values.items():
        print(f"\\val{{{n}}} = {v[:200]}")
    for r in doc.rows:
        print(f"compare {r['name']}: ours {r['ours']}, game {r['game']}, diff {r['diff']}, "
              f"tol {r['tol']}: {'ok' if r['ok'] else 'FAILS'}")
    for ok, m in doc.checks:
        print(f"check {'ok' if ok else 'FAILS'}: {m}")
    return 0


def dependencies_first(slugs, levels):
    """The documents to build, each after the ones it cites."""
    order = []

    def visit(s, stack=()):
        if s in order:
            return
        if s in stack:
            raise DocError(f"documents cite each other in a cycle: {' -> '.join(stack + (s,))}")
        tex = strip_comments(open(os.path.join(DOCS, s, "doc.tex"), encoding="utf-8").read())
        for u in parse_uses(tex):
            if os.path.isfile(os.path.join(DOCS, u, "doc.tex")):
                visit(u, stack + (s,))
        order.append(s)

    for s in slugs:
        visit(s)
    return order


def cmd_build(slugs, force=False):
    levels = curriculum()
    todo = dependencies_first(select(slugs, levels), levels)
    failed = 0
    for slug in todo:
        errors, warnings, doc = check_doc(slug, levels)
        for w in warnings:
            print(f"{slug}: warning: {w}")
        if errors:
            failed += 1
            print(f"{slug}: not typeset, the checks fail:")
            for e in errors:
                print(f"  error: {e}")
            continue
        t0 = time.time()
        problems = typeset(slug, levels, doc, force)
        if problems:
            failed += 1
            print(f"{slug}: FAILED")
            for p in problems:
                print(f"  {p}")
        else:
            pdf = os.path.relpath(os.path.join(PDF_DIR, f"{slug}.pdf"), ROOT)
            print(f"{slug}: ok ({pdf}, {time.time() - t0:.1f} s)")
    return failed


def cmd_game(slugs):
    levels = curriculum()
    by_slug = {e["slug"]: e for e in levels}
    exe = GENERATOR
    if not os.path.exists(exe) and os.path.exists(exe + ".exe"):
        exe += ".exe"
    if not os.path.exists(exe):
        raise DocError(f"the generator is not built ({exe}): cargo build --release -p generator")
    for slug in slugs:
        if slug not in by_slug:
            raise DocError(f"{slug}: not a level")
        e = by_slug[slug]
        level = load_json(e["path"])
        detailed, family = read_placements(slug)
        placements = [level["reference_solution"]] + detailed
        members = family_placements(level, family) if family else []
        placements += members
        with tempfile.TemporaryDirectory() as tmp:
            pfile = os.path.join(tmp, "placements.json")
            with open(pfile, "w", encoding="utf-8") as f:
                json.dump(placements, f)
            t0 = time.time()
            r = subprocess.run([exe, "scan", e["path"], pfile], capture_output=True, text=True,
                               errors="replace")
        if r.returncode != 0:
            raise DocError(f"generator scan failed for {slug}:\n{r.stderr[-2000:]}")
        results = json.loads(r.stdout)
        n = 1 + len(detailed)
        out = {
            "level_file": os.path.relpath(e["path"], ROOT).replace("\\", "/"),
            "level_sha256": level_sha(e["path"]),
            "golden_hash": golden_hash(e["path"]),
            "recorded": datetime.date.today().isoformat(),
            "command": "generator scan <level file> <the reference, then placements.json>",
            "results": results[:n],
        }
        if family:
            out["family"] = family
            out["verdicts"] = "".join(verdict_code(x) for x in results[n:])
        os.makedirs(os.path.join(DOCS, slug), exist_ok=True)
        with open(os.path.join(DOCS, slug, "game.json"), "w", encoding="utf-8", newline="\n") as f:
            json.dump(out, f, indent=1)
            f.write("\n")
        print(f"{slug}: game.json written ({n} placements in detail"
              f"{f', {len(members)} verdicts' if family else ''}, {time.time() - t0:.1f} s)")
    return 0


def cmd_new(slug):
    levels = curriculum()
    by_slug = {e["slug"]: e for e in levels}
    if slug not in by_slug:
        raise DocError(f"{slug}: not a level of levels/curriculum.json")
    folder = os.path.join(DOCS, slug)
    if os.path.exists(os.path.join(folder, "doc.tex")):
        raise DocError(f"docs/levels/{slug}/doc.tex exists already")
    os.makedirs(folder, exist_ok=True)
    e = by_slug[slug]
    for name in ("doc.tex", "calc.py"):
        text = open(os.path.join(TEMPLATE, name), encoding="utf-8").read()
        text = (text.replace("@FINGERPRINT@", level_sha(e["path"]))
                .replace("@TITLE@", e["title"]).replace("@NUMBER@", str(e["number"]))
                .replace("@SLUG@", slug))
        with open(os.path.join(folder, name), "w", encoding="utf-8", newline="\n") as f:
            f.write(text)
    print(f"docs/levels/{slug}/: doc.tex and calc.py from the template. Next:\n"
          f"  python scripts/level_docs.py game {slug}\n"
          f"  (write calc.py, then doc.tex)\n"
          f"  python scripts/level_docs.py build {slug}")
    return 0


def cmd_status():
    levels = curriculum()
    have = set(documents())
    counts = {"documents": 0, "ok": 0}
    print(f"{'#':>4}  {'level':28} document")
    if FOUNDATIONS in have:
        errors, warnings, _ = check_doc(FOUNDATIONS, levels)
        print(f"{'F':>4}  {FOUNDATIONS:28} {'ok' if not errors else 'FAILS'}"
              f"{'' if os.path.isfile(os.path.join(PDF_DIR, FOUNDATIONS + '.pdf')) else ', no PDF'}")
    for e in levels:
        s = e["slug"]
        if s not in have:
            blocked = os.path.isfile(os.path.join(DOCS, s, "NOTES.md"))
            print(f"{e['number']:>4}  {s:28} {'blocked (NOTES.md)' if blocked else '-'}")
            continue
        counts["documents"] += 1
        errors, warnings, _ = check_doc(s, levels)
        if not errors:
            counts["ok"] += 1
        state = "ok" if not errors else f"FAILS ({len(errors)})"
        if any("reviewed" in w for w in warnings):
            state += ", review"
        if not os.path.isfile(os.path.join(PDF_DIR, s + ".pdf")):
            state += ", no PDF"
        print(f"{e['number']:>4}  {s:28} {state}")
    print(f"{counts['documents']} of {len(levels)} levels have a document; {counts['ok']} pass the checks")
    return 0


def main(argv):
    if len(argv) < 2 or argv[1] in ("-h", "--help", "help"):
        print(__doc__)
        return 0
    cmd, args = argv[1], argv[2:]
    force = "--force" in args
    args = [a for a in args if a != "--force"]
    try:
        if cmd == "check":
            return 1 if cmd_check(args) else 0
        if cmd == "calc" and len(args) == 1:
            return cmd_calc(args[0])
        if cmd == "build":
            return 1 if cmd_build(args, force) else 0
        if cmd == "game" and args:
            return cmd_game(args)
        if cmd == "new" and len(args) == 1:
            return cmd_new(args[0])
        if cmd == "status":
            return cmd_status()
    except DocError as e:
        print(f"error: {e}")
        return 1
    print(__doc__)
    return 2


if __name__ == "__main__":
    sys.exit(main(sys.argv))
