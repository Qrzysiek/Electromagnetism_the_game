"""The level designs (scripts/designs/) against the built levels: every design builds without
running anything, its file exists in the curriculum's order, and levels/curriculum.json is
the curriculum's."""
import json
import os
import sys

# Text-only mode, read once when the designs are imported: building a design then runs no
# generator (a level's reference is taken from its built file).
sys.argv = [sys.argv[0], "--text"]

from designs.common import ROOT  # noqa: E402
from designs.curriculum import ARCS, LEVELS  # noqa: E402
from designs.build import with_reality  # noqa: E402


def test_slugs_are_unique():
    slugs = [slug for slug, _ in LEVELS]
    assert len(slugs) == len(set(slugs))


def test_every_level_is_built_in_the_curriculums_order():
    built = sorted(f[:-5] for f in os.listdir(os.path.join(ROOT, "levels"))
                   if f.endswith(".json") and f[:2].isdigit())
    expected = sorted(f"{i:02d}_{slug}" for i, (slug, _) in enumerate(LEVELS, start=1))
    assert built == expected


def test_every_design_builds_with_its_reality():
    for slug, make in LEVELS:
        lvl = with_reality(slug, make())
        assert lvl["name"] and lvl["description"], slug
        assert lvl.get("reality", {}).get("represents"), slug


def test_curriculum_file_is_the_curriculums():
    with open(os.path.join(ROOT, "levels", "curriculum.json"), encoding="utf-8") as f:
        written = json.load(f)
    arcs = [{"name": name, "tiers": [{"name": t, "levels": [slug for slug, _ in lvls]}
                                     for t, lvls in tiers]}
            for name, tiers in ARCS]
    assert written == {"arcs": arcs}
