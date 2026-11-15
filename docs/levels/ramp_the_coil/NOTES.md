# Ramp the coil: notes

## The level's text on the induced field (for the owner; not changed here)

The description says the induced electric field "pushes the ion along a little". For every rate
that solves the level (all negative, `doc.tex` eq. 40.8) the induced field above the coil points
against the ion's motion: `E_x = rate · a_phi · (y − y_c)/rho < 0` with `rate < 0`, `a_phi > 0`
and the ion above the centre. Along the reference flight it does the work −0.0909 on the ion
(`T0 = 0.5`), which arrives at speed 0.905 instead of 1: it slows the ion by about a tenth of its
speed. "Along" can be read as "along its path", but a player will likely read "forward". A
wording such as "which here slows the ion a little" would match the physics.

Reproduce: `py -3.14 scripts/level_docs.py calc ramp_the_coil` (`workRef`, `speedEnd`).
