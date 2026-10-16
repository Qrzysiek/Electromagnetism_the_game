"""Companion of docs/levels/@SLUG@/doc.tex (Level @NUMBER@, @TITLE@).

Every number of the document, computed from the level file (doc.level) with the document's
formulas (their equation labels in comments), independently of the game: then compared with
the game's numbers (doc.flight, doc.result: game.json, written by
`python scripts/level_docs.py game @SLUG@`).

Rules (docs/levels/GROK.md): imports only from the standard library's math modules, numpy,
scipy, mpmath, sympy; no files, no processes; deterministic; under a minute. The game's
numbers enter only doc.compare (checked: the companion is run again with them changed, and
none of its own values may move). Run it alone with `python scripts/level_docs.py calc @SLUG@`.
"""
import mpmath as mp

# TEMPLATE. Tolerances, each set from an error estimate before comparing, with the reason:
# the game's verified flights are accurate to about 1e-10 (Foundations, section 6); an exact
# result here then gets 1e-9; an approximation gets its own error estimate.
TOL = 1e-9


def compute(doc):
    mp.mp.dps = 30
    lvl = doc.level
    shot = lvl["shots"][0]
    # TEMPLATE. The setup: every symbol of the Setup section, from the level file.
    c = mp.mpf(lvl["physics"]["c"])
    doc.value("c", c, 4)

    # TEMPLATE. The derivation and the steps by hand: doc.value for every number the text
    # prints (\val{name}), doc.table for tables, doc.data for figure data, doc.check for
    # internal consistency (e.g. a conserved quantity).

    # TEMPLATE. The comparison with the game's verified reference flight.
    g = doc.flight(0)
    doc.compare("Outcome", "outcome", "Arrived", g["outcome"])
