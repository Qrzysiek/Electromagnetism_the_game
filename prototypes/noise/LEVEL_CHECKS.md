# The thermal-noise level's two checks (2026-11-15)

DRUDE.md Phase 2 item 1 asked, before designing the level, for two measurements that could
sink it: the range of capacitances a player can reach, and how many `RC` the preview holds.
Both were measured on `proto2.py --thermal` (two strips 4 × 4 with a floating contact each,
A the reference, B the player's with a kicker above it; thermal-wall contacts, DEC-43), with
`generator tube-scan t2.json <placements> --trace 4` (`runs.out`: `t2_r1`, `t2_r2`,
`cscan_r1`). The runs shared the machine with the heavy tests' coverage pass (16 cores busy),
on 4 threads each, so their wall times say nothing about the preview's (below).

## The prototype with thermal-wall contacts

| run | spread ratio B/A | decay | pace | B's rms, t = 10–400 |
|---|---|---|---|---|
| refill contacts, preview (`n2_r1`, 2026-11-05) | 0.622 | 6.15 | 8.99 | |
| thermal walls, preview (`t2_r1`, weights 0.08) | 0.842 | 2.41 | 1.96 | 0.723 |
| thermal walls, verification (`t2_r2`, weights 0.04) | 0.746 | 2.20 | 1.78 | 0.484 |

- The pace over the decay is 0.81 (preview) and 0.81 (verification), where the refill
  contacts gave 1.46. D14 measures 0.990 ± 0.019 over many kicks; one kick's decay scatters
  by about 15 % (the table below).
- The rms halves in variance with the macroparticles' weight (0.484/0.723 = 0.67 against
  √½ = 0.71): kT is a macroparticle's, as DEC-42 says. The ratio cancels it: 0.842 against
  0.746, 12 % apart, about the ratio's expected statistical error (each spread `1/√(2N)`
  with `N = T/2τ ≈ 100` independent samples, about 7 %; the ratio of two about 10 %).
- `RC` is short: the pace is about 2 time units, so the window 10–400 holds about 200 of
  them, and the pace's own error `√(2τ/T)` is about 10 %.

## The capacitance's range: a grounded plate beside contact B

`cscan.json`: the empty board, then the player's grounded plate (4 × 0.5, vertical) beside
B's contact at gaps 1.75 (the nearest the plates' clearance allows on whole nodes), 2.75,
4.75 and 7.75.

| plate | spread ratio | decay | pace |
|---|---|---|---|
| none | 0.842 | 2.41 | 1.96 |
| gap 1.75 | 0.785 | 2.70 | 2.36 |
| gap 2.75 | 0.804 | 2.48 | 2.16 |
| gap 4.75 | 0.804 | 2.36 | 1.96 |
| gap 7.75 | 0.802 | 2.03 | 1.65 |

**The plate is not a usable capacitance control here.** Even at the nearest gap it raises
`C_B` by about 15 % against the empty board (`(0.842/0.785)²`), about 5 % against the far
plates (`(0.803/0.785)²`): the ratio moves by 7 % at most, less than its expected statistical
error (about 10 %, above). The decay and the pace, one kick and one window, scatter by about
±15 % between placements that can hardly differ (the far plates: 2.03–2.48 in the decay). The
runs share their seeds, so their scatter may understate the true error. `C_B` is dominated by
what the plate cannot change: the strip's metal at the contact's potential and the kicker
0.5 away (a parallel-plate estimate: the kicker adds ~0.3, the nearest plate ~0.18, per unit
length, to a total of ~1.2). In the z-invariant world a conductor's capacitance goes as one
over the logarithm of the distance to the nearest grounded one, so any plate raises it a
little and its gap matters less.

## What would make the level work (for the owner)

1. A capacitor as the player's control: a lumped capacitance on the floating contact
   (`circuits[].capacitance`, which the engine already has) chosen by the player from a
   list or a logarithmic slider, like a supply's. It spans a factor of ten easily, `kT/C`
   is literally a capacitor's reset noise (every camera pixel), and the two goals then pin
   `C` (the ratio) and `R` (a cut, through the pace) separately. A new control kind: the
   editor's tool, `check_editable`, the slider rules, a DEC entry.
2. A geometry where a plate dominates: a small contact far from the strip's end and the
   kicker, a longer plate parallel to its long face, plates allowed nearer to it. Perhaps a
   factor of two in `C`; to measure.
3. Longer windows to cut the scatter (the preview's 30 s: not measured, below).
4. The pace alone, through a cut: one control scaling one response, an introduction at
   most (the owner's rule).

Recommended: 1. Not decided here.

## The preview's time and the verdict's error bars

Measured at 07:40 with 14 of the 16 cores free (`generator check --no-cache t2.json`,
`t2_check.out`): **preview 24.4 s** to t = 420, inside the owner's 30 s; **verification
112.5 s**, within "a few minutes". (Under the coverage pass's load, on 4 threads, the same run
took about 450 s.)

The verification's error bars are larger than the estimate above: the spread ratio
0.746 ± 0.144 (19 %), the decay 2.20 ± 0.31 (14 %). A goal window has to span about two of
them each side for a comfortable verdict, so the control that sets `C` must move the ratio by
at least about 40 %: `C` by a factor of two or more. That rules out the plate (15 % at best)
and any modest geometry change; a capacitor chosen over a decade (option 1) moves the ratio by
a factor of about three. Longer windows (option 3) shrink the bars as one over the root of the
window: halving them needs four times the preview, past the 30 s.
