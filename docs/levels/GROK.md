# Brief: writing the level documents

You are writing the physics documents of the levels of *Electromagnetism – the game*, a puzzle
game on accurate classical electrodynamics. Read this whole file before you start, then
[`README.md`](README.md) (how the system works), the Foundations
([`foundations/doc.tex`](foundations/doc.tex)) and the two examples:
[`first_bend`](first_bend/) (a closed form) and [`slingshot`](slingshot/) (a numerical method).
Build them once (`python scripts/level_docs.py build`) and read the PDFs in `docs/levels/pdf/`:
they are the standard your documents must meet.

## 1. Your task, and nothing else

For each level, in curriculum order, write `docs/levels/<slug>/doc.tex` and
`docs/levels/<slug>/calc.py` (plus `placements.json` and figures when useful), record the
game's numbers with the script, make every check pass, and commit. That is the whole task.

**You may create or change only these files, and only in the folders of the levels you are
documenting:**

- `docs/levels/<slug>/doc.tex`
- `docs/levels/<slug>/calc.py`
- `docs/levels/<slug>/placements.json`
- `docs/levels/<slug>/game.json`, **only** by running `python scripts/level_docs.py game <slug>`
  (never by hand)
- `docs/levels/<slug>/NOTES.md`
- `docs/levels/<slug>/fig_<name>.pdf|png|jpg` (figures you made; prefer TikZ inside doc.tex)

**You must not change anything else**, in particular: anything under `crates/`, `levels/`
(including `curriculum.json`, `golden_hashes.json`, `math/`), `scripts/`, `tools/`, `audits/`,
`Cargo.toml`, `Cargo.lock`, `PHYSICS.md`, `SPEC.md`, `README.md`, `CLAUDE.md`, the other
`docs/*.md`, `.gitignore`, and inside `docs/levels/`: `README.md`, this file, `leveldoc.cls`,
`_template/`, `foundations/`, `first_bend/`, `slingshot/`, and other levels' folders. Do not
reformat, "fix" or "improve" any of them, even if they look wrong to you. Do not delete files.
Do not add files elsewhere (no helper scripts in the repository; keep scratch files outside it
or in `target/`, which is ignored).

If something outside your scope seems wrong or missing (a macro you need, a bug in the game, a
level whose physics disagrees with your derivation, a check that seems mistaken), **do not fix
it**: describe it in `docs/levels/<slug>/NOTES.md` (what, where, the numbers, how to
reproduce) and carry on with what you can do.

## 2. Git

- Work on a branch of your own, e.g. `git switch -c level-docs/<batch name>` from the commit you
  were given. Never commit to `main`, never push to `main`, never force-push, never rewrite
  history (no `rebase -i`, `commit --amend` of pushed commits, `reset --hard`), never run
  `git stash`, `git clean` or `git checkout -- <path>` on files you did not create.
- One level per commit. Stage only that level's folder: `git add docs/levels/<slug>`, then check
  `git status` shows nothing else staged. Never commit `docs/levels/pdf/` or `target/` (they are
  generated and ignored).
- Commit message: `Level document: <number> <Title>` with a one-line summary of the method
  (closed form, numerical, statistical) and of the comparison (e.g. "12 comparisons, largest
  difference 3e-11"). If the repository's history uses a co-author line, end with your own.
- Commit only documents that pass `python scripts/level_docs.py build <slug>`. If a level is
  blocked (§8), commit its `NOTES.md` alone (and your work in progress renamed to
  `draft.tex` / `draft_calc.py`, which the checks ignore), never a failing `doc.tex`.

## 3. Tools

- Python 3 with numpy, scipy, mpmath, sympy (`python -c "import numpy, scipy, mpmath, sympy"`).
- LuaLaTeX (MiKTeX or TeX Live) with the packages of `leveldoc.cls`.
- The game's generator, built once and after every pull: `cargo build --release -p generator`
  (`$EM_GENERATOR` may point to another build). Building is allowed; editing Rust code is not.
- Useful, read-only: `cargo run --release -p generator -- check levels/<file>.json` (the verdict
  of the reference solution), `... -- trace levels/<file>.json --every 1` (the reference flight's
  steps), `... -- window levels/<file>.json` (for supplies and plates: the solving intervals).

## 4. The workflow for one level

1. **Read** the level: `levels/<NN>_<slug>.json`; its entry in `scripts/levels.py` (the design
   comments explain the intent and often the physics); its row in `docs/LEVELS.md`; the
   sections of `PHYSICS.md` for every element it contains (that is how the engine computes
   them); the Foundations; the documents of earlier levels in its arc.
2. **Start**: `python scripts/level_docs.py new <slug>` (copies the template, fills in the
   level's fingerprint).
3. **Record the game's numbers**: write `placements.json` if you want more placements than the
   reference compared (see §6), then `python scripts/level_docs.py game <slug>`.
4. **Write `calc.py` first**: the physics, computed independently of the game, from the level
   file (§6). Run `python scripts/level_docs.py calc <slug>` until its comparisons pass for the
   right reasons.
5. **Write `doc.tex`** (§5, §7), every number through `\val{...}`.
6. **Check and build**: `python scripts/level_docs.py build <slug>`. Read the PDF in full:
   layout, figures, every number in context, every reference.
7. **Review** with the checklist (§9), then commit (§2).

Work in curriculum order (`python scripts/level_docs.py status` lists the levels in order): a
document may cite only earlier levels, and the first level of each arc should carry the general
derivations its arc needs.

## 5. What a document must contain

The section titles are fixed and checked. Model every section on the two examples.

1. **What the level asks**: `\leveldescription` (the game's text, inserted for you), then the
   puzzle in physical terms; the goal **exactly** (boxes, acceptance windows, goal values and
   their time windows, tolerances); what the player controls (`limits`); what ends a flight
   (`\xsecref{foundations}{sec:board}` and anything this level adds); the reference solution and,
   in a sentence, why it works; what the document does (closed form, numerical, statistical).
2. **Setup**: `\setuptable` (inserted), then every symbol the derivation uses with its value
   (`\val`), and the geometry. A figure of the board (pgfplots/TikZ, coordinates in cells; data
   from `calc.py` via `doc.data` and `\addplot table {name.dat};`; `\rawval{...}` inside
   coordinates).
3. **Principles and model**: a `principles` list, each item with its source (an equation of the
   Foundations or an earlier level, or a textbook section you have checked); the model as
   equations specialised to this level, with initial and boundary conditions; every assumption
   or approximation in an `assumption` box: what it neglects, why, and how large that is, as a
   number from `calc.py`. Say which approximations the game makes too.
4. **Derivation**: from the principles to formulas that can be evaluated, every step shown or
   cited. No "it can be shown", no "it is well known that" without a citation, no skipped
   algebra that a student could not reconstruct. Where an earlier document derived something,
   cite it (`\xeqref`, `\xsecref`) instead of repeating it. Where no closed form exists, state the
   exact mathematical problem and a numerical method a person could carry out (the Foundations'
   toolkit: RK4, Simpson, Newton, bisection, Hermite interpolation), with its error estimate.
   Derive the independent checks the physics offers (conserved quantities, limits, symmetries).
5. **The reference solution by hand**: a `procedure` of numbered steps a student can follow with
   a calculator, every intermediate value shown (`\val`), tables where an iteration or an
   integration is long (`doc.table`), ending in a `result` box.
6. **Comparison with the game**: `\comparisontable` (generated), then a discussion: are the
   differences the size the error estimates predict? Name the largest and explain it.
7. **Other solutions**: why the reference works and what else does: alternatives checked in
   detail against the game (`placements.json` "detailed"), a map of solutions when it is cheap
   (a "family", see first_bend), the edges of a slider's solving window (`generator window`).
   A map taken from the game's verdicts rather than derived must say so (see slingshot).
8. **What the game computes**: the engine's method for this level with `\physref{...}` to
   `PHYSICS.md`, what in it is exact and what approximate, and why its numbers and the
   document's agree or how far they may differ.
9. **Further reading**: textbook sections (checked) and the earlier levels that lead here.

Then `\levelfileappendix`. Length: as long as the physics needs; the examples are 7–10 pages.

**The standard "by hand, in principle".** Given the level file, this document, the documents it
cites and the Foundations, a reader who knows undergraduate physics and mathematics and has
unlimited time must be able to compute every number in the comparison table. Nothing may rest
on "the computer says so". A numerical method counts if it is fully specified (method, step,
stopping rule, event location) and its error is estimated. For a stochastic model (the tubes'
macroparticles, the metals' Drude carriers) the hand computation is the model's exact
statement plus its analytic consequences (kinetic theory, the Child–Langmuir law, Ohm's law,
the Hall effect, ...) with their accuracy; the comparison is then statistical (§6).

## 6. `calc.py`

- Defines `compute(doc)`. Imports only `math`, `cmath`, `fractions`, `decimal`, `itertools`,
  `functools`, `statistics`, `numpy`, `scipy`, `mpmath`, `sympy`. No files, no processes, no
  network, no randomness without a fixed seed. Under two minutes (the check runs it twice).
- Reads the level from `doc.level`; never hard-code a number that the level file contains.
  Use the level file's numbers as the game does (binary doubles: `mpmath.mpf(1e-06)` is exact
  for the double, not for 10⁻⁶).
- Implements **your** formulas and methods, with the equation labels of `doc.tex` in comments.
  Never import or copy the engine's code; the point is an independent computation.
- `doc.value(name, x, digits=10, sci=False, fixed=None)` for every number the text prints
  (names: letters and digits only); `doc.text` for words; `doc.table` for tables; `doc.data`
  for plot data; `doc.check(condition, message)` for the derivation's internal checks
  (conserved quantities, two routes to the same number).
- The game's numbers: `doc.flight(i, flight=0)` (flight levels: `["end"]["t"]`, `["end"]["x"]`,
  `["end"]["p"]`, `["outcome"]`, `["status"]`, `["radiated_energy"]`, `["steps"]`), or
  `doc.result(i)` for everything `generator scan` reports (tube levels: `["tube"]["currents"]`,
  `["potentials"]` with `value`, `error`, `preview`; beam levels: `["beams"]`). Placement 0 is
  the reference; 1, 2, ... are `placements.json`'s "detailed" list. They may be used **only**
  in `doc.compare` (or shown as the game's with `doc.game_value`); the check reruns `calc.py`
  with all game numbers changed and fails if any of your values moves.
- `doc.compare(name, what, ours, game, tol, rel=False)`: one row of the comparison table. Every
  document compares at least the outcome and the quantities the level's goal depends on.
  **Tolerances are set from error estimates before you look at the differences**, with the
  reasoning in a comment: an exact result against a verified flight gets 1e-9 (ten times the
  game's verified accuracy, Foundations §6); a numerical method gets ten times its own error
  estimate or 1e-9, whichever is larger; an approximation gets its estimated size; a statistical
  game value gets its error bar (`error`) plus your model's own uncertainty. **Never widen a
  tolerance to make a comparison pass.** A failing comparison means your physics, your
  arithmetic, or the game is wrong: find out which (§8).
- `placements.json`: `{"detailed": [placement, ...], "family": {"kind": "charge", "values": [...]}}`
  (a placement is a list of elements as in the level file's `reference_solution`). A family is
  every node of the grid with each value: use one only for fast flight levels, and either
  derive its verdicts (`doc.compare_verdicts`, see first_bend) or present them as the game's
  (`doc.verdicts`, `doc.game_value`, see slingshot). Never for tube or beam levels (too slow).

## 7. LaTeX rules

- Start with `\documentclass{leveldoc}`, then `\DocUses{foundations,...}` (the documents you cite
  by equation or section: the Foundations and earlier levels) and `\DocFingerprint{...}`
  (filled in by `new`; update it only after reviewing the text against a changed level).
- Not allowed (checked): `\usepackage`, `\input`, `\include`, file or Lua access, `\def`,
  `\renewcommand`, `\let`, `\makeatletter`, `\csname`, `\hypersetup` and the class's internals.
  `\newcommand` for local notation is allowed. Everything else you need is in the class: if
  something is missing, ask in `NOTES.md`.
- Use the class's notation: `\vect{p}`, `\uvec{r}`, `\dd`, `\deriv{x}{t}`, `\pderiv{}{}`,
  `\abs{}`, `\norm{}`, `\order{}`; environments `principles`, `assumption` (label `as:...` right
  after `\begin{assumption}{Title}`), `result`, `procedure` (items labelled `st:...`), `remark`;
  `align`, `equation`, `tikzpicture`, `axis`, `groupplot`, `figure`, `table`.
- Labels: `eq:`, `sec:`, `fig:`, `tab:`, `as:`, `st:` followed by lower-case words and digits.
- References: within the document `\eqref`, `\ref`, `§\ref`; to other documents `\xeqref{slug}{eq:...}`,
  `\xsecref{slug}{sec:...}`, `\xref{slug}{...}`, `\levelref{slug}`; to the engine `\physref{5}`.
- Books: `\Jackson{§14.2}`, `\LandauII{§39}`, `\Goldstein{§3.10}`, `\Griffiths{Ch.~7}`,
  `\Brau{...}`, `\HairerI{§II.4}`. Cite a section or equation number **only if you have checked
  it** (Jackson and Brau are in `references/` on the owner's machine); otherwise cite the
  chapter or the topic. Never invent a reference.
- Numbers: only through `\val{...}` (or `\rawval` in coordinates). Small exact integers of the
  setup (a node, a face) may be typed; anything computed may not. The check warns about typed
  numbers with four or more decimals.
- Writing: English, present tense, concise and precise; define every symbol at its first use;
  units are the game's (§2 of the Foundations), no SI unless converting. No filler, no
  marketing, no emojis, no rhetorical questions.

## 8. When something does not agree

1. Recheck your own work first: the derivation, the arithmetic, the level file's numbers (the
   game reads binary doubles), the events (which boundary first? the world is the grid plus one
   cell on each side), relativity (`c` finite?), the units (`μ₀/4π = 1/c²`).
2. Read the relevant `PHYSICS.md` section: the engine may model something differently than you
   assumed (that is a difference of model, which the document must describe, not hide).
3. If the difference remains and you believe the game is wrong: do not change tolerances, do not
   change the game. Write `NOTES.md` with the numbers, your derivation, and how to reproduce
   (`python scripts/level_docs.py calc <slug>`), rename your files to `draft.tex` /
   `draft_calc.py`, commit that, and go on to the next level.
4. If a level is beyond a by-hand treatment in some part (say, a field computed by a boundary
   element method with thousands of panels), state the exact problem the engine solves, its
   method, and the analytic approximations you can derive with their accuracy; compare within
   that accuracy. Say clearly what is exact, what is approximate, and how large the difference
   is. Never present an approximation as exact.

## 9. Done: the checklist for each level

- [ ] `python scripts/level_docs.py build <slug>` passes, with no warning about typed numbers
      or the fingerprint.
- [ ] The PDF reads well from start to end; figures are right; no `??`.
- [ ] Every assumption has the size of what it neglects, as a number.
- [ ] Every step of the derivation is shown or cited; every citation checked.
- [ ] The by-hand procedure reaches every number of the comparison table.
- [ ] The comparison covers the outcome and every quantity the goal depends on; every tolerance
      has a comment with its reason, set before comparing.
- [ ] The largest difference is explained in the comparison section.
- [ ] Only `docs/levels/<slug>/` is staged; the commit message follows §2.

## 10. Notes on the arcs

- **Charges** (levels 1–14): Coulomb fields, relativistic dynamics; single centres have closed
  forms (Level 1), several need numerics (Level 2); detectors with acceptance windows and gates
  (`PHYSICS.md` §6.1–6.2); Thomson's atom: charge clouds (§2.1).
- **Metal and electrodes**: metal spheres by images (§2.6: exact series), electrodes by the
  boundary element method (§2.7: a linear system; by hand: the method and the analytic limits,
  e.g. a parallel-plate field, with the method's measured accuracy), dielectrics, power supplies.
- **Relativity and magnetism**: magnets and coils (§2.2: Biot–Savart, `μ₀/4π = 1/c²`), ferrites,
  magnetic moments (§3.2: Stern–Gerlach), crossed fields.
- **Time**: stray fields, plane waves (§2.3), ramped coils and induction, antennas (§2.4),
  radiation reaction (§3.1, the Landau–Lifshitz force) and radiation goals (§3.4).
- **Circuits** (§2.8–2.10): modified nodal analysis, the stiff integrator (§5.4).
- **Beams** (§3.3): many particles, their mutual fields; `doc.result(i)["beams"]`.
- **Jackson levels** (`docs/JACKSON.md`): each cites its section of Jackson; check it there.
- **Tubes** (levels 100–108, §2.11): the z-invariant world (logarithmic potentials, line
  charges), macroparticles, space charge; the game's numbers are statistical
  (`["tube"]["currents"][k]` with `value` and `error`); derive the analytic law the level
  shows (Child–Langmuir, the triode's amplification, ...) and compare within its accuracy.
- **Metals** (levels 109–112, §2.11 and `docs/DRUDE.md`): the Drude model of carriers among ion
  cores; Ohm's law, resistances, the Hall effect; compare the game's currents and potentials
  within their error bars and the kinetic theory's corrections (the carriers' heating, the
  contacts, the walls' surface charge are all in `PHYSICS.md`).
