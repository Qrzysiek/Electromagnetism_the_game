# Glass screen

## Resolved 2026-11-15: the disagreement was the draft's

The draft's flight (2026-10-16) left the board where the game's reference arrives. The
cause was a bug in `draft_calc.py`, not in the engine: `Plates._solve` projected each
source panel's field on the *source* panel's normal (`self._en(..., pack)`, the source
pack's kind and outward sign) where the dielectric condition needs the *receiving* panel's.
Every interaction between faces of different orientation, and between opposite faces (whose
outward signs differ), was wrong. The source term and the self term already used the
receiving panel's normal, which is why the solve looked well behaved (residual 5.5e-11) and
converged in the panel size: it converged to the solution of a wrong equation.

Fixed (the receiving panel's `kind` and `sign`). The draft's field then matches the game's
(`target/scratch` probes, the game's `LevelField` at its verification mesh against the
draft's rectangles at 0.25):

| point | game Ey/Ey_bare | draft, fixed | draft, before |
|---|---|---|---|
| (15, 10) launch line under the charge | 1.0793 | 1.0789 | 1.358 |
| (15, 12) centre of the glass | 0.392 | 0.387 | 0.233 |
| (10, 10) under the glass's left end | 1.4317 | 1.4298 | 1.118 |
| (22, 10) | 1.1600 | 1.1571 | 0.967 |
| (25, 12) right of the glass | 0.8485 | 0.8494 | 0.788 |
| (15, 13.5) between the charge and the glass | 1.0258 | 1.0261 | 1.026 |

And the flights agree (the draft's step 0.01, panels 0.5; the game's verified run): the
reference arrives at y = 8.5375 against 8.4815, the empty placement leaves the bottom at
x = 22.4779 against 22.4760.

## What it changed in the level: its text

The level's text said the glass weakened the field beyond it "by roughly the permittivity".
It does not: inside the glass the charge's field falls to 0.39, but on the beam's line it
grows (the table: +8 % under the charge, +43 % under the glass's ends), because the far
face's positive bound charge is nearer the beam than the near face's negative one, and a
slab 2 high in z does not hold the field lines in. With no charges placed the beam leaves
the bottom at x = 22.48 with the glass and at 23.73 without it: the glass throws it off
harder. Even an infinite slab would leave 0.73 of the field on the launch line (the image
series `q(1 − β²) Σ β^{2n} / (a + 2nd + d + h)²`, β = (ε − 1)/(ε + 1)). The design's
comment assumed the screening and never measured it. The description, the reality entry
(`scripts/reality.py`), docs/LEVELS.md and docs/CURRICULUM.md were corrected on 2026-11-15;
the level's physics, reference and golden hash are unchanged.

## Left for the document

- `calc.py` from `draft_calc.py`: under the 120 s limit (it takes 136 s: vectorize the
  rectangles over all panels and corners at once), and an error estimate that means
  something: doubling the panel edge to 1.0 moves the reference's end by 0.16, so the
  tolerance (ten times that) is 1.6, which passes anything. Panels 0.5 against 0.25 instead.
- `doc.tex` from `draft.tex` (GROK.md §5), with the field table above as the physics the
  level shows.
