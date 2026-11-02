"""Disk cleanup: what takes space in the checkouts and what is left over (DEC-37; CLAUDE.md,
*Disk*).

    python scripts/cleanup.py [status]   each checkout's target/ by folder and what `run`
                                         would delete (read-only)
    python scripts/cleanup.py run        deletes it

Every checkout git lists (the main one and its worktrees) is cleaned alike:
- target/temp/: temporary target directories and outputs (the game's executable locked, an
  experiment): deleted whole.
- a folder or file under target/ not in KNOWN: left over (an extra CARGO_TARGET_DIR, a stray
  output): deleted. A tool that writes a new folder under target/ registers it in KNOWN.
- a build over its cap (CAPS): cargo never deletes a build it no longer uses (a dependency, a
  feature set or the toolchain changed), so the profiles only grow: cleaned (cargo clean), and
  in the main checkout the game and the generator rebuilt right after (the owner plays from
  it); elsewhere the next build rebuilds.
- target/run-cache/: records older than the checkout's generator (another engine wrote them:
  never hit again).
- target/scratch/: entries untouched for SCRATCH_DAYS days that no tracked file cites. A cited
  one belongs in prototypes/ (it is not versioned here, and goes with its worktree): reported.
- folders under .claude/worktrees/ that git no longer lists and that hold no `.git` (a removed
  worktree's leftovers): deleted. A worktree whose branch is merged into main with nothing
  uncommitted is reported: removing it (`git worktree remove`) is the owner's call.
"""
import argparse
import fnmatch
import os
import shutil
import stat
import subprocess
import sys
import time
from dataclasses import dataclass

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GB = float(1 << 30)
DAY = 86400.0
SCRATCH_DAYS = 14
# What may live under a checkout's target/, and who writes it.
KNOWN = {
    "release": "cargo's release profile: the game, the generator, the heavy tests",
    "debug": "cargo's dev profile: the fast tier, clippy",
    "tmp": "cargo's CARGO_TARGET_TMPDIR (integration tests)",
    "doc": "cargo doc",
    "affected": "scripts/affected.py: the shared records (main checkout), the instrumented build",
    "run-cache": "the generator's run cache (level::run_cache)",
    "test-passes": "the heavy tests' passes (level::run_cache)",
    "progress": "scripts/progress.py, runner.py: the jobs' logs and ledgers",
    "measured": "scripts/measured.py: the numbers the tests measured",
    "level-docs": "scripts/level_docs.py: the level documents' builds",
    "pairs": "scripts/pairs.py",
    "screenshots": "scripts/readme_shots.py: the captures before conversion",
    "scratch": "a session's scripts and outputs (pruned after SCRATCH_DAYS days)",
    "temp": "temporary target directories and outputs (deleted whole)",
    "CACHEDIR.TAG": "cargo",
    ".rustc_info.json": "cargo",
    ".future-incompat-report.json": "cargo",
    ".package-cache": "cargo",
}
# A build's cap in GB: about twice its fresh size (measured 2026-11-02, DEC-37: the release
# profile with the game, the generator and the heavy tests 2.4 GB, rebuilt in ~6.5 min; the
# dev one with the fast tier and clippy 4.4 GB, ~5.5 min; affected/build, the instrumented
# heavy tests, 0.2-0.4 GB).
CAPS = {"release": 5.0, "debug": 9.0, os.path.join("affected", "build"): 1.0}


@dataclass
class Item:
    """Something `run` deletes: its path, size in bytes, why, and how (`tree`: the folder or
    file; `profile:<name>`: cargo clean of that profile; `target`: cargo clean of that target
    directory)."""
    path: str
    size: int
    why: str
    how: str = "tree"


def tree_size(path: str) -> int:
    """Bytes under `path` (a file's own size), links not followed."""
    try:
        st = os.lstat(path)
    except OSError:
        return 0
    if not stat.S_ISDIR(st.st_mode):
        return st.st_size
    total = 0
    stack = [path]
    while stack:
        try:
            with os.scandir(stack.pop()) as it:
                for e in it:
                    try:
                        if e.is_dir(follow_symlinks=False):
                            stack.append(e.path)
                        else:
                            total += e.stat(follow_symlinks=False).st_size
                    except OSError:
                        pass
        except OSError:
            pass
    return total


def newest_mtime(path: str) -> float:
    """The latest modification under `path` (its own for a file)."""
    newest = os.lstat(path).st_mtime
    if os.path.isdir(path) and not os.path.islink(path):
        for d, dirs, files in os.walk(path):
            for n in dirs + files:
                try:
                    newest = max(newest, os.lstat(os.path.join(d, n)).st_mtime)
                except OSError:
                    pass
    return newest


def leftovers(target: str) -> list[Item]:
    """target/temp and every entry under target/ that KNOWN does not name."""
    items = []
    for name in sorted(os.listdir(target)):
        path = os.path.join(target, name)
        if name == "temp":
            items.append(Item(path, tree_size(path), "temporary (target/temp)"))
        elif name not in KNOWN:
            cargo = os.path.isfile(os.path.join(path, "CACHEDIR.TAG"))
            items.append(Item(path, tree_size(path),
                              "a leftover target directory" if cargo else "not a known output",
                              "target" if cargo else "tree"))
    return items


def over_caps(target: str, caps: dict[str, float]) -> list[Item]:
    """The builds over their caps."""
    items = []
    for rel, cap in caps.items():
        path = os.path.join(target, rel)
        size = tree_size(path)
        if size > cap * GB:
            how = f"profile:{rel}" if rel in ("release", "debug") else "target"
            items.append(Item(path, size, f"over its cap of {cap:g} GB (stale builds)", how))
    return items


def stale_runs(cache: str, generator: str) -> list[Item]:
    """Run-cache records older than the generator that would read them."""
    if not (os.path.isdir(cache) and os.path.isfile(generator)):
        return []
    built = os.stat(generator).st_mtime
    items = []
    for e in os.scandir(cache):
        if e.is_file() and e.stat().st_mtime < built:
            items.append(Item(e.path, e.stat().st_size, "written by an earlier generator"))
    return items


def cited_scratch(checkout: str) -> list[str]:
    """The scratch entries tracked files cite (`target/scratch/<pattern>`), as name patterns;
    prototypes/ not searched (its scripts name the paths they ran with, not citations)."""
    r = subprocess.run(["git", "grep", "-ohE", r"target/scratch/[A-Za-z0-9_.*-]+", "--",
                        ".", ":!prototypes"], cwd=checkout, capture_output=True, text=True)
    return sorted({c[len("target/scratch/"):].rstrip(".") for c in r.stdout.split()})


def stale_scratch(scratch: str, cited: list[str], now: float,
                  days: float = SCRATCH_DAYS) -> tuple[list[Item], list[str]]:
    """Scratch entries untouched for `days` that no citation matches, and the cited ones."""
    if not os.path.isdir(scratch):
        return [], []
    items, kept = [], []
    for name in sorted(os.listdir(scratch)):
        path = os.path.join(scratch, name)
        if any(fnmatch.fnmatch(name, c) for c in cited):
            kept.append(name)
        elif now - newest_mtime(path) > days * DAY:
            items.append(Item(path, tree_size(path), f"scratch untouched for {days:g} days"))
    return items, kept


def orphan_worktrees(folder: str, listed: list[str]) -> list[Item]:
    """Folders under .claude/worktrees that git does not list and that hold no `.git`."""
    if not os.path.isdir(folder):
        return []
    known = {os.path.normcase(os.path.abspath(p)) for p in listed}
    items = []
    for name in sorted(os.listdir(folder)):
        path = os.path.join(folder, name)
        if (os.path.isdir(path) and os.path.normcase(os.path.abspath(path)) not in known
                and not os.path.exists(os.path.join(path, ".git"))):
            items.append(Item(path, tree_size(path), "a removed worktree's leftovers"))
    return items


def git(args: list[str], cwd: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(["git", *args], cwd=cwd, capture_output=True, text=True)


def checkouts() -> list[tuple[str, str]]:
    """(path, branch) of every checkout git lists, the main one first."""
    r = git(["worktree", "list", "--porcelain"], ROOT)
    if r.returncode != 0:
        sys.exit("cleanup: git worktree list failed: " + r.stderr.strip())
    found: list[tuple[str, str]] = []
    path = ""
    for line in r.stdout.splitlines() + [""]:
        if line.startswith("worktree "):
            path = os.path.normpath(line[len("worktree "):])
        elif line.startswith("branch "):
            found.append((path, line[len("branch refs/heads/"):]))
            path = ""
        elif not line and path:
            found.append((path, ""))
            path = ""
    return found


def plan(path: str, main: bool, now: float) -> tuple[list[Item], list[str]]:
    """What to delete in one checkout, and remarks."""
    target = os.path.join(path, "target")
    items: list[Item] = []
    notes: list[str] = []
    if os.path.isdir(target):
        items += leftovers(target)
        items += over_caps(target, CAPS)
        exe = "generator.exe" if os.name == "nt" else "generator"
        items += stale_runs(os.path.join(target, "run-cache"), os.path.join(target, "release", exe))
        old, kept = stale_scratch(os.path.join(target, "scratch"), cited_scratch(path), now)
        items += old
        if kept:
            notes.append("tracked files cite target/scratch/ " + ", ".join(kept)
                         + ": move them to prototypes/ (scratch is not versioned)")
    if main:
        listed = [p for p, _ in checkouts()]
        items += orphan_worktrees(os.path.join(path, ".claude", "worktrees"), listed)
    return items, notes


def gb(n: float) -> str:
    return f"{n / GB:6.2f} GB"


def report(path: str, items: list[Item], notes: list[str]) -> None:
    target = os.path.join(path, "target")
    print(f"\n{path}")
    if os.path.isdir(target):
        sizes = sorted(((tree_size(os.path.join(target, n)), n) for n in os.listdir(target)),
                       reverse=True)
        print(f"  target/ {gb(sum(s for s, _ in sizes)).strip()}")
        for s, n in sizes:
            if s >= 0.01 * GB:
                print(f"  {gb(s)}  {n}  ({KNOWN.get(n, 'not known')})")
    groups: dict[tuple[str, str], list[Item]] = {}
    for it in items:
        groups.setdefault((os.path.dirname(it.path), it.why), []).append(it)
    for (folder, why), group in groups.items():
        size = sum(it.size for it in group)
        if len(group) > 3:
            print(f"  delete {gb(size)}  {os.path.relpath(folder, path)}: {len(group)} entries, {why}")
        else:
            for it in group:
                print(f"  delete {gb(it.size)}  {os.path.relpath(it.path, path)}: {why}")
    for n in notes:
        print("  note: " + n)


def remove_tree(path: str) -> list[str]:
    """Deletes a folder or file; the paths it could not delete (in use)."""
    failed: list[str] = []
    if os.path.isdir(path) and not os.path.islink(path):
        def onexc(_f: object, p: str, _e: BaseException) -> None:
            try:
                os.chmod(p, stat.S_IWRITE)
                os.remove(p)
            except OSError:
                failed.append(p)
        shutil.rmtree(path, onexc=onexc)
    else:
        try:
            os.remove(path)
        except OSError:
            failed.append(path)
    return failed


def delete(path: str, it: Item) -> bool:
    """Deletes one item; whether it all went."""
    if it.how == "tree":
        failed = remove_tree(it.path)
        for p in failed[:5]:
            print(f"    could not delete {p} (in use?)")
        return not failed
    if it.how == "target":
        cmd = ["cargo", "clean", "--target-dir", it.path]
    else:
        profile = it.how.split(":", 1)[1]
        cmd = ["cargo", "clean", "--release" if profile == "release" else "--profile=dev"]
    r = subprocess.run(cmd, cwd=path, capture_output=True, text=True)
    if r.returncode != 0:
        print("    " + (r.stderr.strip().splitlines() or ["cargo clean failed"])[-1])
    return r.returncode == 0


def main() -> None:
    ap = argparse.ArgumentParser(description="Disk cleanup of the checkouts (DEC-37).")
    ap.add_argument("action", nargs="?", default="status", choices=["status", "run"])
    args = ap.parse_args()
    now = time.time()
    total = 0
    all_co = checkouts()
    main_path = all_co[0][0]
    for k, (path, branch) in enumerate(all_co):
        items, notes = plan(path, k == 0, now)
        if k > 0 and branch:
            merged = git(["merge-base", "--is-ancestor", branch, "main"], path).returncode == 0
            clean = not git(["status", "--porcelain"], path).stdout.strip()
            if merged and clean:
                notes.append(f"{branch} is merged into main with nothing uncommitted: remove "
                             f"the worktree when done with it (git worktree remove {path})")
        report(path, items, notes)
        total += sum(it.size for it in items)
        if args.action == "run":
            rebuild = False
            failed = 0
            for it in items:
                if not delete(path, it):
                    failed += 1
                    print(f"  partly deleted {os.path.relpath(it.path, path)}")
                rebuild |= path == main_path and it.how == "profile:release"
            if items:
                print(f"  deleted {len(items) - failed} of {len(items)}")
            if rebuild:
                print("  rebuilding the game and the generator (cargo build --release)")
                subprocess.run(["cargo", "build", "--release", "-p", "game", "-p", "generator"],
                               cwd=path)
    verb = "freed" if args.action == "run" else "`run` would free"
    print(f"\n{verb} {gb(total).strip()}")


if __name__ == "__main__":
    main()
