"""The disk cleanup (DEC-37): what it deletes and what it must keep."""
import os

import cleanup


def touch(path, size=10, age_days=0.0, now=1e9):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "wb") as f:
        f.write(b"x" * size)
    t = now - age_days * cleanup.DAY
    os.utime(path, (t, t))


def names(items):
    return sorted(os.path.basename(it.path) for it in items)


def test_known_outputs_stay_and_the_rest_goes(tmp_path):
    t = tmp_path / "target"
    for known in ("release", "debug", "affected", "run-cache", "progress", "scratch"):
        touch(str(t / known / "f"))
    touch(str(t / "CACHEDIR.TAG"))
    touch(str(t / "claude-gen" / "CACHEDIR.TAG"))
    touch(str(t / "stray.log"))
    touch(str(t / "temp" / "build" / "x"))
    items = {os.path.basename(it.path): it for it in cleanup.leftovers(str(t))}
    assert sorted(items) == ["claude-gen", "stray.log", "temp"]
    assert items["claude-gen"].how == "target"
    assert items["stray.log"].how == "tree"
    assert items["temp"].size == 10


def test_builds_over_their_caps(tmp_path):
    t = tmp_path / "target"
    touch(str(t / "release" / "deps" / "a"), size=3000)
    touch(str(t / "debug" / "deps" / "a"), size=100)
    touch(str(t / "affected" / "build" / "a"), size=3000)
    caps = {"release": 2000 / cleanup.GB, "debug": 2000 / cleanup.GB,
            os.path.join("affected", "build"): 2000 / cleanup.GB}
    items = {os.path.relpath(it.path, t): it.how for it in cleanup.over_caps(str(t), caps)}
    assert items == {"release": "profile:release", os.path.join("affected", "build"): "target"}


def test_run_cache_records_older_than_the_generator_go(tmp_path):
    now = 1e9
    touch(str(tmp_path / "release" / "generator.exe"), age_days=1, now=now)
    touch(str(tmp_path / "run-cache" / "old.json"), age_days=2, now=now)
    touch(str(tmp_path / "run-cache" / "new.json"), age_days=0.5, now=now)
    items = cleanup.stale_runs(str(tmp_path / "run-cache"), str(tmp_path / "release" / "generator.exe"))
    assert names(items) == ["old.json"]
    # Without a generator nothing can be told stale.
    assert cleanup.stale_runs(str(tmp_path / "run-cache"), str(tmp_path / "none.exe")) == []


def test_scratch_old_and_uncited_goes_cited_stays(tmp_path):
    now = 1e9
    s = tmp_path / "scratch"
    touch(str(s / "old.py"), age_days=20, now=now)
    touch(str(s / "recent.py"), age_days=3, now=now)
    touch(str(s / "amp_proto.py"), age_days=30, now=now)
    touch(str(s / "scans" / "a.json"), age_days=30, now=now)
    touch(str(s / "busy" / "old.json"), age_days=30, now=now)
    touch(str(s / "busy" / "new.json"), age_days=1, now=now)
    os.utime(str(s / "scans"), (now - 30 * cleanup.DAY,) * 2)
    os.utime(str(s / "busy"), (now - 30 * cleanup.DAY,) * 2)
    items, kept = cleanup.stale_scratch(str(s), ["amp_*.py"], now)
    assert names(items) == ["old.py", "scans"]
    assert kept == ["amp_proto.py"]


def test_only_unlisted_worktree_folders_without_git_go(tmp_path):
    w = tmp_path / "worktrees"
    touch(str(w / "live" / ".git"))
    touch(str(w / "unlisted_live" / ".git"))
    touch(str(w / "leftover" / "scripts" / "x"))
    items = cleanup.orphan_worktrees(str(w), [str(w / "live")])
    assert names(items) == ["leftover"]


def test_remove_tree_deletes_read_only_files(tmp_path):
    touch(str(tmp_path / "d" / "ro"))
    os.chmod(str(tmp_path / "d" / "ro"), 0o444)
    assert cleanup.remove_tree(str(tmp_path / "d")) == []
    assert not os.path.exists(str(tmp_path / "d"))
