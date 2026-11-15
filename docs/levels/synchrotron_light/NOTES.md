# Synchrotron light: notes

## The reality remark on the loss per turn (for the owner; not changed here)

`scripts/reality.py`'s remark "A hugely radiating particle" says the particle "radiates 0.1 % of
its energy per turn at modest speed". In this level it radiates far more. At the launch the drag
power is `P = 2q⁴B²p²/(3m⁴c³) = 1.96e-3` (`B = 0.1612`, `p = 0.952`, `q = m = 1`, `c = 2`) and the
period 43.2, so the first turn takes 0.085: 19.7 % of the kinetic energy `T0 = 0.43`, or 1.9 % of
the total energy `γmc² = 4.43`. Over the reference flight (about twelve turns) it radiates 82 % of
`T0`. The "0.1 %" matches the fields.py design comment's general setting (`B ~ 0.1`, "spirals
inwards over ~30 turns"), not this level. A remark such as "radiates about a fifth of its kinetic
energy in its first turn" would match.

Reproduce: `py -3.14 scripts/level_docs.py calc synchrotron_light` (`powerZero`, `period`,
`lossTurn`, `radFrac`).
