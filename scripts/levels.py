"""All shipped levels, as a curriculum (single source of truth).

    python scripts/levels.py            # write levels/NN_*.json
    python scripts/levels.py --only 07  # rebuild one level
    python scripts/levels.py --text     # only the texts changed: rewrite them, check nothing

What each level represents and how far it is the real thing (its headline, remarks and the
labels of its electrodes and metals) is written in `scripts/reality.py`.

Order: by how hard the phenomenon is to understand and use. Every new element or
concept first gets an easy introduction level; within a chapter the possibility space
and the rarity of solutions grow on average (measured with `generator analyze`), and
so does the number of elements needed.

Each level is described in physical terms. Where a level has a designed reference setup
(`reference=`), detectors can be placed automatically around where the reference
flights actually land (`Auto` detectors): the script writes a probe version with wide
detector strips, runs `generator check` on it, and puts a detector of the requested size
around each landing point. Levels without a designed reference get one from the solver
(`generator solve --write`).

Units: k = 1, cell = 1. Particles are weakly charged (q ~ 1e-6) and electrodes strongly
charged (Q ~ 1e6) so that neglected radiation stays far below the numerical accuracy
(PHYSICS.md section 8). Magnetic strengths are in field units (PHYSICS.md 2.2): a coil of
radius a with kappa = mu0 I / 4 pi has B = 2 pi kappa / a at its centre, and a particle of
momentum p circles with radius p / (q B).
"""

from designs.build import main

if __name__ == "__main__":
    main()
