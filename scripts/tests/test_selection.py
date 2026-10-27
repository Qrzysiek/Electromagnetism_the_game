"""affected.py's selection rules on synthetic records (DEC-16): which tests each kind of
change selects, with coverage up to date and with coverage older than the verdict."""
import pytest

import affected

SNAP = "s1"
FILES = {"a.rs": ["a-code", "a-decl"], "b.rs": ["b-code", "b-decl"], "c.rs": ["c-code", "c-decl"]}


class Program:
    """The test program as the selection asks it: the tests and the levels' physics."""

    def tests(self):
        return ["golden::one", "golden::two"]

    def physics(self):
        return {"one": "p1", "two": "p2"}


def record(files, cov_snapshot=SNAP, level_file="f", physics="p1", golden="g"):
    return {"snapshot": SNAP, "cov_snapshot": cov_snapshot, "files": files, "result": "ok",
            "level_file": level_file, "physics": physics, "golden": golden}


@pytest.fixture
def tree(monkeypatch):
    """The working tree's state, as the test sets it."""
    state = {"universe": {k: list(v) for k, v in FILES.items()}, "global": "G",
             "golden": {"one": "g", "two": "g"}, "level_files": {"one": "f", "two": "f"}}
    monkeypatch.setattr(affected, "universe", lambda: state["universe"])
    monkeypatch.setattr(affected, "global_inputs", lambda: state["global"])
    monkeypatch.setattr(affected, "golden_entries", lambda: state["golden"])
    monkeypatch.setattr(affected, "sha", lambda p: state["level_files"].get(p.split("/")[-1][:-5]))
    return state


def records(one, two):
    return {"snapshots": {SNAP: {"universe": FILES, "global": "G"}},
            "tests": {"golden::one": one, "golden::two": two}}


def selected(recs):
    return sorted(n for n, _ in affected.select(recs, Program()))


def test_nothing_changed_selects_nothing(tree):
    assert selected(records(record(["a.rs"]), record(["b.rs"], physics="p2"))) == []


def test_a_body_selects_the_tests_that_executed_its_file(tree):
    tree["universe"]["a.rs"] = ["a-code-2", "a-decl"]
    assert selected(records(record(["a.rs"]), record(["b.rs"], physics="p2"))) == ["golden::one"]


def test_a_declaration_selects_every_test(tree):
    tree["universe"]["c.rs"] = ["c-code-2", "c-decl-2"]
    assert selected(records(record(["a.rs"]), record(["b.rs"], physics="p2"))) == ["golden::one", "golden::two"]


def test_a_new_file_selects_every_test(tree):
    tree["universe"]["d.rs"] = ["d-code", "d-decl"]
    assert selected(records(record(["a.rs"]), record(["b.rs"], physics="p2"))) == ["golden::one", "golden::two"]


def test_stale_coverage_is_selected_by_any_body(tree):
    # Test two's coverage predates its verdict: a body it did not execute then may matter now.
    tree["universe"]["a.rs"] = ["a-code-2", "a-decl"]
    recs = records(record(["a.rs"]), record(["b.rs"], cov_snapshot="s0", physics="p2"))
    assert selected(recs) == ["golden::one", "golden::two"]


def test_stale_coverage_without_a_change_is_not_selected(tree):
    recs = records(record(["a.rs"]), record(None, cov_snapshot=None, physics="p2"))
    assert selected(recs) == []


def test_a_reworded_level_is_not_selected(tree):
    tree["level_files"]["one"] = "f-reworded"  # the file changed, its physics did not
    assert selected(records(record(["a.rs"]), record(["b.rs"], physics="p2"))) == []


def test_a_levels_physics_selects_its_tests(tree, monkeypatch):
    tree["level_files"]["one"] = "f-2"
    monkeypatch.setattr(Program, "physics", lambda self: {"one": "p1-changed", "two": "p2"})
    assert selected(records(record(["a.rs"]), record(["b.rs"], physics="p2"))) == ["golden::one"]


def test_a_golden_entry_selects_its_golden_test(tree):
    tree["golden"]["two"] = "g-2"
    assert selected(records(record(["a.rs"]), record(["b.rs"], physics="p2"))) == ["golden::two"]


def test_a_global_input_selects_every_test(tree):
    tree["global"] = "G-2"
    assert selected(records(record(["a.rs"]), record(["b.rs"], physics="p2"))) == ["golden::one", "golden::two"]


def test_failed_and_unrecorded_tests_are_selected(tree):
    failed = dict(record(["a.rs"]), result="failed")
    recs = records(failed, record(["b.rs"], physics="p2"))
    del recs["tests"]["golden::two"]
    assert selected(recs) == ["golden::one", "golden::two"]
