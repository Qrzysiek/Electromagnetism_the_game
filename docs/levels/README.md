# The level documents

Every level gets a document that derives what happens in it from first principles (Maxwell's
equations, the Lorentz force, special relativity, and the numerical methods any textbook
gives) down to the numbers the game computes: where the particle arrives, when, how fast. An
undergraduate with unlimited time should be able to reproduce every number by hand from the
level's document, the documents it cites, and the Foundations. The game opens a level's
document from its panel ("The physics of this level", "Foundations").

Status: the Foundations and two worked examples, [`first_bend`](first_bend/doc.tex) (Level 1:
a closed form, the exact relativistic Coulomb orbit, and an exact map of every single-charge
solution) and [`slingshot`](slingshot/doc.tex) (Level 2: no closed form, a hand method, RK4
with an error estimate). The rest are to be written to [`GROK.md`](GROK.md);
`python scripts/level_docs.py status` shows where things stand.

## Commands

```
python scripts/level_docs.py status              every level: document, checks, PDF
python scripts/level_docs.py new <slug>          start a document from _template/
python scripts/level_docs.py game <slug>         record the game's numbers (generator scan)
python scripts/level_docs.py calc <slug>         run a document's companion: values, comparisons
python scripts/level_docs.py check [<slug> ...]  every check except typesetting (no TeX needed)
python scripts/level_docs.py build [<slug> ...]  check, then typeset into docs/levels/pdf/
```

`build` without arguments builds every document (each after the ones it cites; unchanged ones
are skipped). Typesetting needs LuaLaTeX with the packages of `leveldoc.cls` (MiKTeX installs
missing ones on first use); `game` needs the generator (`cargo build --release -p generator`;
`$EM_GENERATOR` points elsewhere); the companions need Python 3 with numpy, scipy, mpmath and
sympy. The PDFs (`docs/levels/pdf/`) and the build directories (`target/level-docs/`) are
generated and not committed, like `PHYSICS.pdf`.

## Layout

```
docs/levels/
  README.md            this file: the system
  GROK.md              the brief for writing the documents (scope, workflow, rules, checklist)
  leveldoc.cls         the LaTeX class: layout, title block, notation, environments, references
  _template/           doc.tex and calc.py that `new` copies
  foundations/         the Foundations: units, Maxwell's equations, mechanics, the board's
                       rules, verification, radiation, the by-hand numerical toolkit
  <slug>/              one folder per level, named by its slug in levels/curriculum.json:
    doc.tex            the text
    calc.py            its companion: every number of the text, and the comparison with the game
    game.json          the game's numbers (written by `game`, never by hand)
    placements.json    optional: further placements to compare (in detail, or a whole family)
    NOTES.md           optional: questions and problems for the owner
    fig_*.pdf|png|jpg  optional figures (TikZ is preferred, inside doc.tex)
  pdf/                 generated: <slug>.pdf, foundations.pdf
```

A document is filed under the level's **slug** (`first_bend`), not its number: when levels are
inserted and renumbered, the documents stay where they are and their numbers follow
`levels/curriculum.json` at the next build.

## How a document is made

```
levels/<NN>_<slug>.json ──► calc.py (independent of the game) ──► values.tex, comparison.tex, *.dat
        │                        ▲                                         │
        │                        │ only in doc.compare                     ▼
        └──► generator scan ──► game.json                 doc.tex + leveldoc.cls ──► LuaLaTeX ──► pdf/<slug>.pdf
```

1. **The text** (`doc.tex`) has nine sections in a fixed order: *What the level asks*, *Setup*,
   *Principles and model*, *Derivation*, *The reference solution by hand*, *Comparison with the
   game*, *Other solutions*, *What the game computes*, *Further reading*; and the level file
   verbatim as an appendix. The class typesets the title block, the level's own text from the
   game, a table of every field of the level file, the comparison table and the appendix:
   authors write only physics.
2. **Every number** in the text comes from the companion `calc.py` (`\val{name}`), which reads
   the level file and implements the document's formulas or method. A number typed into the
   text could drift from the level; a computed one cannot.
3. **The game's numbers** come from `generator scan` (a full-precision JSON of each flight's
   outcome, verification status, end state, radiated energy; for tube and beam levels their
   goals' values and error bars), recorded in `game.json` by `level_docs.py game`, for the
   reference solution and the placements in `placements.json`.
4. **The comparison**: `calc.py` compares its numbers with the game's (`doc.compare`, each with a
   tolerance set from the error estimates, with the reason in a comment); the table in the
   PDF is its output.

## What the checks guarantee

`check` (and `build`, before it typesets) fails a document when:

- its structure differs from the nine sections, or it lacks the generated material (the
  level's text, the setup table, the comparison table, the appendix);
- it loads packages, reads or writes files, runs Lua, redefines macros or touches the class's
  internals (the list is `FORBIDDEN` in `scripts/level_docs.py`); labels must look like `eq:...`,
  `sec:...`, `fig:...`, `tab:...`, `as:...`, `st:...`;
- it cites a level that is not earlier in the curriculum, or one without a document;
- a `\val{...}` has no value from `calc.py`, or a plot's data file is not written by it;
- `calc.py` imports anything but the standard library's math modules, numpy, scipy, mpmath and
  sympy, opens files or runs processes, fails, or takes over two minutes;
- **any of `calc.py`'s own numbers depends on the game's**: the companion is run a second time
  with every number of `game.json` changed, and none of its values may move (the game's
  numbers may enter only the comparison, or values the text presents as the game's,
  `doc.game_value`);
- a comparison fails its tolerance;
- `game.json` is stale: the level file or the golden hash of its reference flight
  (`levels/golden_hashes.json`, which changes with the physics) differs from the recorded one;
- the template's `TEMPLATE` placeholders remain.

It warns when the level file has changed since the text was last reviewed
(`\DocFingerprint`, the first 16 hex digits of the level file's SHA-256): the numbers update
themselves, the prose may not. Typesetting fails on any LaTeX error, undefined reference,
missing glyph, or a line sticking out of the margin by more than 20 pt.

## References between documents

Equations, figures and tables carry their document's number: (1.7) is Level 1's seventh
equation in every document, (F.3) the Foundations' third. A document lists the documents it
cites in `\DocUses{foundations,first_bend}` (the Foundations and earlier levels only; the build
compiles them first, and their labels resolve through LaTeX's `xr-hyper`), then writes
`\xeqref{first_bend}{eq:orbit}` → "(1.6)", `\xsecref{foundations}{sec:board}` → "§F.5",
`\xsecref{first_bend}{sec:time}` → "Level 1, §4.5", `\levelref{slingshot}` → "Level 2
(*Slingshot*)", `\physref{7}` → "PHYSICS.md §7". The links work between the PDFs in
`docs/levels/pdf/` in a PDF viewer.

A later level cites an earlier one instead of repeating its derivation; the Foundations hold what
every level uses. Each arc's first level (the curriculum introduces every element in an
easy level first) is the natural home of its element's general derivation: magnets and coils,
conductors, waves, circuits, beams, the tubes' z-invariant world, metals.

## Keeping the documents true

- **A level changes** (its file): `check` reports `game.json` stale; run `game <slug>`, then
  `check` and `build`. The numbers follow the level automatically; where the comparisons still
  pass, review the prose and set `\DocFingerprint` to the value `check` prints.
- **The physics changes** (the golden hashes change): the same, for every documented level whose
  hash changed. A comparison that now fails means either the document's model no longer
  describes the game (update the document) or the change broke something (fix the engine).
  Tolerances are never widened to make a row pass.
- **Levels are renumbered**: nothing to do; numbers follow the curriculum at the next build.
- **The class or the build changes**: the owner's (or Claude's) job, never the authors'; a
  missing macro or environment is requested in the level's `NOTES.md`.

## Relation to the engine

`PHYSICS.md` documents the engine: its models, methods and validation tests. The level documents
explain the levels: each derives its level's physics independently and shows that the engine
agrees. They cite `PHYSICS.md` for the engine's methods (*What the game computes*) and never
replace it. When the engine is one day separated from the game, `PHYSICS.md` and the
Foundations' physics would go with the engine, and the level documents stay with the game.
