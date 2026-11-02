# Self-inductance introduction: design notes (2026-11-02)

Ring 20 x 20 at (2,2)-(22,22), hole 16 x 16 (walls 2), density 2, thermal speed 2, mfp 0.5,
cores 0.005, 50 per area (7200 carriers), step 0.012, B0 = 0.5, heat over t = 10-40.
Estimates: sigma = rho mu0 = 2 x (3/8)sqrt(pi/2)/(kappa s = 4) = 0.235; R = P/(sigma w) =
72/(0.235 x 2) = 153; L = 4 pi A / c^2 = 4 pi 324 / 4 = 1018; L/R = 6.6, bend at 0.15.

| omega | c = 2 | c = 100 (no self-field) | ratio | 1/(1+(omega L/R)^2) |
|---|---|---|---|---|
| 0.05 | 0.106 | | | |
| 0.1 | 0.643 | 0.920 | 0.70 | 0.70 |
| 0.2 | 1.50 | | | |
| 0.4 | 2.42 | 13.6 | 0.18 | 0.125 |
| 0.8 | 6.50 | | | |
| 1.6 | 22.5 | 218.6 | 0.10 | 0.009 |
| 3.2 | 73.2 | | | |

Unscreened heat goes as omega^2 (0.92 -> 13.6 -> 219). The ring's loop screening matches
L/R = 6.55 at omega = 0.1. At high omega ~10 % of the unscreened heat remains: the wall's
own eddy currents (skin depth c/sqrt(2 pi sigma omega) = 1.3 at 1.6, below the wall 2).
Run time: ~7 min per 40 time units with the verdicts running, ~4 min alone (c = 2);
c = 100: 70 s (coarser self-field grid).

Next: wall 1, mfp 0.25 (sigma halves, R x4): the walls' own share ~4x smaller; c ~1.2 to keep
the bend near 0.2 (L x2.8). Then the frequency-only supply (B0 fixed), the slider range
(log, e.g. 0.03-1), the heat goal at the plateau's start.

## Walls 1, mfp 0.25, c = 1.2 (sigma 0.1175, R 647, L 3150, L/R 4.9), B0 0.5

| omega | c = 1.2 | c = 100 | ratio | loop theory |
|---|---|---|---|---|
| 0.05 | 0.0127 | | | |
| 0.1 | 0.225 | 0.331 | 0.68 | 0.81 |
| 0.2 | 0.671 | | | |
| 0.4 | 1.058 | 3.89 | 0.27 | 0.21 |
| 0.8 | 3.08 | | | |
| 1.6 | 10.6 | 64.8 | 0.16 | 0.016 |

The residual at high omega grew (16 % against 10 %): not the wall's own eddy currents
(those fall as w^2). What grew is omega/(4 pi sigma) (1.5 here, 2.95 before): the square
ring's corners need surface charges to steer the current, which the metal cannot move
faster than 4 pi sigma. Cure: a better conductor. Density 20 (omega_p dt = 0.19): sigma
2.35, 4 pi sigma 30; walls 2, mfp 0.5, c = 7.3 keeps L/R = 5 (bend at 0.2), and c/omega =
36 > the board at the bend (the Darwin limit's condition holds there).

## Density 20 (sigma 2.35, 4 pi sigma 30), walls 2, mfp 0.5, c = 7.3 (L/R 5), B0 0.15

| omega | 0.05 | 0.1 | 0.2 | 0.4 | 0.8 | 1.6 | 3.2 |
|---|---|---|---|---|---|---|---|
| c = 7.3 | 0.099 | 0.761 | 0.982 | 2.64 | 7.29 | 17.4 | 63.7 |
| c = 100 | | 1.054 | | 12.5 | | 184 | |

Ratio 0.72, 0.21, 0.095 at 0.1, 0.4, 1.6; the loop theory 0.80, 0.20, 0.015. Above the loop's
prediction a constant ~8-9.5 % of the unscreened heat remains at 0.8, 1.6, 3.2 (the c = 100
heat itself falls short of omega^2 by 26-32 % at 0.4-1.6).

## Density 20, walls 1, mfp 0.25, c = 3.65 (L/R 5.3 estimated), B0 0.15

0.4: 1.056 against 2.658 at c = 100 (0.40; theory 0.18); 1.6: 9.03 against 59.4 (0.152;
theory 0.015). Thinner walls screen less, not more: the wall-flux hypothesis is refuted.

## Summary of the residual

Beyond the loop's screening ~9-10 % of the unscreened heat remains with walls 2, ~14-16 %
with walls 1, whatever the density (2 or 20) or c (1.2-7.3). Not the wall's own eddy
currents (they fall as w^2), not the corners' steering charges (4 pi sigma 30 changes
nothing), not the wall's flux. Untested: the walls' effective resistance above P/(sigma w)
(the thin ring's 0.40 at omega = 0.4 fits R about twice the bulk estimate), a square
ring's non-uniform current, a round ring (`Lattice::with_disc`, not in the level format).

## 2026-11-02: the residual's cause (the owner: find it first)

Exact reference: `eddy2d.py` (finite differences, div(D grad B) = i w B, D = 1/(mu0 sigma(w))
with the Boltzmann sigma(w); the hole a huge D; B = 1 on the outer edge). Times B0^2 = 0.0225:
c = 100: 12.97 (0.4), 192 (1.6) against the game's 12.5, 184 (-3.4 %, -4 %); c = 7.3: 2.64,
3.25, 3.88, 4.93 at 0.4, 0.8, 1.6, 3.2 against the game's 2.64, 7.29, 17.4, 63.7.

Cause 1 (fixed): `field_work`/`induced_work` were booked at the push's velocity, before
`darwin_step` took back the induction's share of the impulse: (q^2/2m) E.dA_self/dt dt too
much per carrier and time, ~dt/2tau of the unscreened heat (density- and c-independent;
doubles with the mean free path halved: the walls-1 runs' 15 %). Now booked after
`darwin_step` (`Tube::book_work`). Rerun (`fixed`): 2.10, 4.21, 5.89, 17.4.

Caveat: the window t = 10-40 is not whole periods; the drive's work swings by
sqrt(1 + (wL/R)^2) of its mean, so each point carries up to ~17 % window error.
`--whole 1` now rounds the window to whole periods.

Remaining: ~1-2 % of the unscreened heat at 0.8-3.2. Under test: dt 0.006 (halfdt), the
self-field cell 0.125 (finecell; c/omega_p = 0.46 is under two cells of 0.25 in every
design so far), V-D8 (disc, c/omega_p 0.85).

V-D8 (disc, c/omega_p 0.85, dt 0.012): +0.9 %, +1.0 % (two seeds) against the cylinder.

Noise (ring c 7.3, 1.6, whole periods, dt 0.012, cell 0.25): seeds 21-24: 5.62, 4.88, 4.99,
4.59: mean 5.02 +- 0.22 (sd 0.43) against exact 3.88 (3.75 with the game's sigma).

Resolved design (walls 2, mfp 0.5, density 22, c 14.5: c/omega_p 0.87, L/R ~1.4, bend ~0.7;
B0 0.1, t from 5, whole periods, two seeds 21/22): 0.7: 8.99/10.41 (exact 9.32); 1.4:
13.52/15.22 (14.84); 2.8: 20.75/19.43 (18.36); c = 1e4 at 1.4: 70.7 (74.2). ~190 s per run
(7 at once). Residual at 2.8: +1.7 +- 0.5 = 0.6 % of the unscreened heat, as the first
design's at 1.6 (0.6 %): a common residual ~0.5 % of the unscreened heat, visible where the
screening leaves < 10 %. Under test: 4 seeds at dt 0.006 and 4 at cell 0.125 (first design).
