"""The heavy tests a change can affect, selected by what they executed (DEC-16,
docs/TESTING.md).

    python scripts/affected.py status
        the heavy tests the working tree selects, and why (nothing is run)
    python scripts/affected.py run [--all] [--filter TEXT] [--part K/N] [--no-coverage]
                                   [--coverage-only] [--jobs N]
        runs the selected heavy tests (--all: every one; --filter: those whose name contains
        TEXT; --part: the Kth of N interleaved parts of all the tests): first their verdicts,
        with the plain test program at full speed, then their coverage, with the instrumented
        one (each test in its own process on one thread, as many at once as the machine has
        threads, the longest first). --no-coverage leaves the coverage for `refresh`;
        --coverage-only runs the instrumented program alone, for verdict and coverage at once
        (a first recording). One `[k/N] <test> ok|FAILED (<s> s)` line per run; exit 1 if
        one failed
    python scripts/affected.py refresh [--jobs N]
        the coverage of every test whose coverage is older than its verdict
    python scripts/affected.py note [--base REV]
        the change note for the commit message: the files changed against REV (default HEAD,
        with the untracked ones), what they affect, the heavy tests they select and the
        results recorded for them

The heavy tests are those of `crates/level/tests/per_level.rs` (`<check>::<level>`). A
test's record holds the state of every engine file when it last passed (a snapshot: each
file's code and its declarations, both hashed without comments or layout), the workspace
files its run executed (LLVM source-based coverage, `-C instrument-coverage`) and the
snapshot that coverage was taken at, the build's global inputs, its level's physics
(`Level::physics_json`, hashed) and, for `golden::`, the golden entry. A change to a
function's body selects the tests that executed its file; a change to the declarations
(constants, statics, types, signatures, macros, `const fn` bodies: what can act without
running, which coverage cannot see) or a new or removed file selects every test; a changed
global input, physics or golden entry its tests; a test that failed or was never recorded is
selected; and a test whose coverage is older than its last verdict (a code change may have
sent it into files it did not execute before) is selected by any change to the engine until
`refresh` brings its coverage up to date. The runs are deterministic: a test whose inputs are
bit for bit unchanged can only give its recorded result. Records live in `target/affected/`
(per checkout; a fresh one records everything on its first run).
"""

import argparse
import concurrent.futures
import datetime
import glob
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import threading
import time
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
AFF = os.path.join(ROOT, "target", "affected")
BUILD = os.path.join(AFF, "build")
PROF = os.path.join(AFF, "prof")
RECORDS = os.path.join(AFF, "records.json")
ENGINE = ["crates/physics/src", "crates/level/src"]
TEST_FILE = "crates/level/tests/per_level.rs"
GLOBAL_FILES = ["Cargo.lock", "Cargo.toml", "crates/physics/Cargo.toml",
                "crates/level/Cargo.toml", "rust-toolchain.toml"]
GOLDEN = "levels/golden_hashes.json"
DEPENDENTS = {"physics": ["physics", "level", "game", "generator"],
              "level": ["level", "game", "generator"], "game": ["game"],
              "generator": ["generator"]}


def rel(path):
    """A path relative to the repository, with forward slashes (None outside it)."""
    r = os.path.relpath(os.path.abspath(path), ROOT)
    return None if r.startswith("..") else r.replace("\\", "/")


def sha(path):
    try:
        with open(os.path.join(ROOT, path), "rb") as f:
            return hashlib.sha256(f.read()).hexdigest()
    except FileNotFoundError:
        return None


def test_only_files():
    """The engine's unit-test modules (`#[cfg(test)] mod name;`) and the files under them:
    not part of the heavy tests' program at all."""
    out = set()
    for d in ENGINE:
        for p in glob.glob(os.path.join(ROOT, d, "**", "*.rs"), recursive=True):
            lines = open(p, encoding="utf-8").read().split("\n")
            base = os.path.dirname(p)
            if os.path.basename(p) not in ("lib.rs", "mod.rs", "main.rs"):
                base = os.path.join(base, os.path.splitext(os.path.basename(p))[0])
            for i, line in enumerate(lines):
                if line.strip() != "#[cfg(test)]":
                    continue
                j = i + 1
                while j < len(lines) and lines[j].strip().startswith(("#[", "//")):
                    j += 1
                m = re.match(r"\s*(?:pub(?:\([^)]*\))?\s+)?mod\s+(\w+)\s*;", lines[j] if j < len(lines) else "")
                if m:
                    for c in (os.path.join(base, m.group(1) + ".rs"), os.path.join(base, m.group(1), "mod.rs")):
                        if os.path.exists(c):
                            out.add(rel(c))
                    sub = os.path.join(base, m.group(1))
                    out.update(rel(q) for q in glob.glob(os.path.join(sub, "**", "*.rs"), recursive=True))
    return out


def universe():
    """Every engine file (its unit tests' aside) and the test program's source: path ->
    (hash of the file, hash of its declarations)."""
    out = {}
    skip = test_only_files()
    paths = [rel(p) for d in ENGINE for p in glob.glob(os.path.join(ROOT, d, "**", "*.rs"), recursive=True)]
    for p in paths + [TEST_FILE]:
        if p in skip:
            continue
        text = open(os.path.join(ROOT, p), encoding="utf-8").read()
        code = " ".join(tokens(text))  # comments and layout aside
        out[p] = [hashlib.sha256(code.encode()).hexdigest(),
                  hashlib.sha256(declarations(text).encode()).hexdigest()]
    return out


TOKEN = re.compile(r"""
    (?P<ws>\s+)
  | (?P<line>//[^\n]*)
  | (?P<block>/\*)
  | (?P<raw>b?r(?P<hashes>\#*)")
  | (?P<str>b?"(?:\\.|[^"\\])*")
  | (?P<char>b?'(?:\\(?:x[0-9a-fA-F]{2}|u\{[0-9a-fA-F]+\}|.)|[^\\'])')
  | (?P<life>'[A-Za-z_][A-Za-z0-9_]*)
  | (?P<num>[0-9][0-9A-Za-z_.]*)
  | (?P<ident>[A-Za-z_][A-Za-z0-9_]*)
  | (?P<punct>.)
""", re.VERBOSE | re.DOTALL)


def tokens(text):
    """Rust's tokens without comments: strings, chars and lifetimes whole."""
    i, out = 0, []
    while i < len(text):
        m = TOKEN.match(text, i)
        kind = m.lastgroup if m.lastgroup != "hashes" else "raw"
        if m.group("raw"):
            end = text.index('"' + m.group("hashes"), m.end())
            out.append(text[i:end + 1 + len(m.group("hashes"))])
            i = end + 1 + len(m.group("hashes"))
            continue
        if kind == "block":
            depth, j = 1, m.end()
            while depth:
                if text.startswith("/*", j):
                    depth, j = depth + 1, j + 2
                elif text.startswith("*/", j):
                    depth, j = depth - 1, j + 2
                else:
                    j += 1
            i = j
            continue
        if kind not in ("ws", "line"):
            out.append(m.group(0))
        i = m.end()
    return out


def declarations(text):
    """A file's code without its functions' bodies (`const fn` bodies kept: they run at
    compile time, where coverage cannot see them). What is left (constants, statics,
    types, signatures, macros, attributes) can change a test without being executed, so a
    change to it selects every test; a change to a body only those that executed the file."""
    t = tokens(text)
    out, i = [], 0
    while i < len(t):
        tok = t[i]
        out.append(tok)
        if tok == "fn" and i + 1 < len(t) and re.match(r"[A-Za-z_]", t[i + 1]) \
                and not (i > 0 and t[i - 1] == "const"):
            # The signature up to the body's `{` (or `;` for a declaration), outside brackets.
            j, depth = i + 1, 0
            while j < len(t):
                if t[j] in "([":
                    depth += 1
                elif t[j] in ")]":
                    depth -= 1
                elif depth == 0 and t[j] in "{;":
                    break
                out.append(t[j])
                j += 1
            if j < len(t) and t[j] == "{":
                braces, k = 1, j + 1
                while braces:
                    braces += {"{": 1, "}": -1}.get(t[k], 0)
                    k += 1
                out.append("{}")
                i = k
                continue
            i = j
            continue
        i += 1
    return " ".join(out)


def global_inputs():
    h = hashlib.sha256()
    for f in GLOBAL_FILES:
        h.update(f.encode() + b"\0" + (sha(f) or "missing").encode() + b"\0")
    h.update(subprocess.run(["rustc", "-vV"], cwd=ROOT, capture_output=True, text=True).stdout.encode())
    return h.hexdigest()


def load_records():
    try:
        with open(RECORDS, encoding="utf-8") as f:
            records = json.load(f)
    except FileNotFoundError:
        return {"snapshots": {}, "tests": {}}
    for r in records["tests"].values():
        # Records from before verdicts and coverage were taken apart: their coverage was
        # taken with their verdict, by the instrumented program.
        if "cov_snapshot" not in r and r.get("files") is not None:
            r["cov_snapshot"] = r["snapshot"]
            r["seconds_cov"] = r.pop("seconds", None)
    return records


def save_records(records):
    os.makedirs(AFF, exist_ok=True)
    tmp = RECORDS + ".tmp"
    with open(tmp, "w", encoding="utf-8") as f:
        json.dump(records, f, indent=0, sort_keys=True)
    for _ in range(40):
        try:
            os.replace(tmp, RECORDS)
            return
        except PermissionError:  # Windows: a reader holds it
            time.sleep(0.05)
    raise SystemExit(f"could not write {RECORDS}")


def llvm_tool(name):
    sysroot = subprocess.run(["rustc", "--print", "sysroot"], capture_output=True, text=True).stdout.strip()
    host = next(l.split()[1] for l in subprocess.run(["rustc", "-vV"], capture_output=True,
                                                      text=True).stdout.splitlines()
                if l.startswith("host:"))
    exe = os.path.join(sysroot, "lib", "rustlib", host, "bin", name + (".exe" if os.name == "nt" else ""))
    if not os.path.exists(exe):
        raise SystemExit(f"{name} not found: rustup component add llvm-tools")
    return exe


def build(instrumented=True):
    """The test program: the instrumented one (its own target directory) or the plain one
    (cargo's usual target directory)."""
    os.makedirs(PROF, exist_ok=True)
    env = dict(os.environ)
    if instrumented:
        # Only the workspace's crates are instrumented (`coverage_rustc.py`, cargo's
        # RUSTC_WORKSPACE_WRAPPER): the dependencies' code, which `Cargo.lock` alone decides,
        # runs without counters. The build scripts and proc macros are instrumented too and
        # run during the build: their profiles go to the discarded ones (by default they land
        # in the working directory).
        wrapper = os.path.join(ROOT, "scripts", "coverage_rustc" + (".cmd" if os.name == "nt" else ".py"))
        env.update(RUSTC_WORKSPACE_WRAPPER=wrapper, CARGO_TARGET_DIR=BUILD,
                   LLVM_PROFILE_FILE=os.path.join(PROF, "build-%p.profraw"))
        env.pop("RUSTFLAGS", None)
    print(f"[affected] building the {'instrumented' if instrumented else 'plain'} test program", flush=True)
    r = subprocess.run(["cargo", "test", "--release", "-p", "level", "--test", "per_level", "--no-run",
                        "--message-format=json"], cwd=ROOT, env=env, capture_output=True, text=True,
                       encoding="utf-8", errors="replace")
    if r.returncode != 0:
        sys.stderr.write(r.stderr[-4000:])
        raise SystemExit("build failed")
    for line in r.stdout.splitlines():
        try:
            m = json.loads(line)
        except json.JSONDecodeError:
            continue
        if m.get("reason") == "compiler-artifact" and m.get("target", {}).get("name") == "per_level" \
                and m.get("executable"):
            return m["executable"]
    raise SystemExit("no test program in cargo's output")


class Program:
    """The instrumented test program, built when first needed. `fresh=False` (status, note)
    takes the last build if there is one: it only lists the tests and the levels' physics,
    and a change to the code that computes them selects the tests anyway."""

    def __init__(self, fresh=True):
        self.exe, self.fresh = None, fresh

    def get(self):
        if self.exe is None and not self.fresh:
            built = glob.glob(os.path.join(BUILD, "release", "deps", "per_level-*" + (".exe" if os.name == "nt" else "")))
            built = [b for b in built if not b.endswith((".d", ".pdb"))]
            if built:
                self.exe = max(built, key=os.path.getmtime)
        if self.exe is None:
            self.exe = build()
        return self.exe

    def quiet_env(self, **extra):
        os.makedirs(PROF, exist_ok=True)
        return dict(os.environ, LLVM_PROFILE_FILE=os.path.join(PROF, "discard-%p.profraw"), **extra)

    def tests(self):
        r = subprocess.run([self.get(), "--list", "--ignored"], capture_output=True, text=True,
                           env=self.quiet_env())
        return [l[:-len(": test")] for l in r.stdout.splitlines() if l.endswith(": test")]

    def physics(self):
        r = subprocess.run([self.get()], capture_output=True, text=True, env=self.quiet_env(EM_TEST_INPUTS="1"))
        return dict(l.split() for l in r.stdout.splitlines() if l.strip())


def level_of(test):
    return test.split("::", 1)[1]


def golden_entries():
    with open(os.path.join(ROOT, GOLDEN), encoding="utf-8") as f:
        return json.load(f)


def select(records, program, names=None):
    """(test, reason) for every heavy test the working tree selects, in name order."""
    now, glob_now = universe(), global_inputs()
    golden = golden_entries()
    names = names if names is not None else program.tests()
    physics = {}  # asked of the program only for levels whose file changed

    def physics_of(level):
        if not physics:
            physics.update(program.physics())
        return physics.get(level)

    out = []
    diffs = {}  # per snapshot: (files whose declarations changed, files whose bodies did)

    def diff(snap_id):
        if snap_id not in diffs:
            old = records["snapshots"][snap_id]["universe"]
            changed = {f for f in set(old) | set(now) if old.get(f) != now.get(f)}
            # A new or removed file counts as declarations: it may add an impl anywhere.
            decl = sorted(f for f in changed if (old.get(f) or [0, 0])[1] != (now.get(f) or [0, 0])[1])
            diffs[snap_id] = (decl, sorted(changed - set(decl)))
        return diffs[snap_id]

    for name in names:
        r = records["tests"].get(name)
        level = level_of(name)
        if r is None:
            out.append((name, "not recorded"))
            continue
        if r["result"] != "ok":
            out.append((name, "failed last time"))
            continue
        snap = records["snapshots"][r["snapshot"]]
        if snap["global"] != glob_now:
            out.append((name, "build inputs changed (Cargo files, toolchain)"))
            continue
        if sha(f"levels/{level}.json") != r["level_file"] and physics_of(level) != r["physics"]:
            out.append((name, "the level's physics changed"))
            continue
        if name.startswith("golden::") and golden.get(level) != r["golden"]:
            out.append((name, "its golden entry changed"))
            continue
        decl, bodies = diff(r["snapshot"])
        if decl:
            out.append((name, "declarations changed: " + listed(decl)))
            continue
        if r.get("files") is None or r.get("cov_snapshot") != r["snapshot"]:
            # Coverage older than the verdict: any change to the engine may concern it.
            if bodies:
                out.append((name, "coverage to refresh, and the engine changed: " + listed(bodies)))
            continue
        hit = [f for f in bodies if f in r["files"]]
        if hit:
            out.append((name, "executed code changed: " + listed(hit)))
    return out


def listed(files):
    more = f" and {len(files) - 3} more" if len(files) > 3 else ""
    return ", ".join(files[:3]) + more


def all_tests(records, program):
    """Every heavy test the program lists (records of tests no longer listed are dropped)."""
    names = program.tests()
    for gone in set(records["tests"]) - set(names):
        del records["tests"][gone]
    return names


# Checks whose coverage is another check's of the same level: `repeated::` computes
# `golden::`'s hash twice in one process, the same engine code on the same level (the runs
# are deterministic: the same path), so it executes exactly the files `golden::` does. Its
# verdict needs a run (the plain program's); its coverage needs none. The metal tubes'
# `repeated::` took 34 minutes and more on one instrumented thread.
DERIVED = {"repeated": "golden"}


def source_of(name):
    """The test whose coverage `name` takes, or None."""
    check, level = name.split("::", 1)
    return f"{DERIVED[check]}::{level}" if check in DERIVED else None


def derive_coverage(records, snap_id, names):
    """Each derived test among `names` takes its source's coverage where both were given at
    `snap_id`; the names left (their source's coverage is not current) are returned."""
    rest = []
    for n in names:
        r, src = records["tests"].get(n), records["tests"].get(source_of(n))
        if r and src and r.get("snapshot") == snap_id and src.get("cov_snapshot") == snap_id \
                and src.get("files") is not None:
            r.update({"files": list(src["files"]), "cov_snapshot": snap_id, "derived_from": source_of(n)})
        else:
            rest.append(n)
    save_records(records)
    return rest


def take_coverage(chosen, exe, records, snap_id, jobs):
    """The coverage of `chosen` [(test, why)]: the derived tests' from their sources, the
    others' by instrumented runs, the longest first. The failed tests' names."""
    direct = [(n, why) for n, why in chosen if source_of(n) is None]
    derived = [n for n, _ in chosen if source_of(n) is not None]
    direct.sort(key=lambda t: -expected_seconds(t[0], records["tests"].get(t[0])))
    print(f"[affected] coverage of {len(direct)} tests, the longest first"
          + (f" ({len(derived)} more take their source's)" if derived else ""), flush=True)
    # One rayon thread per test, as many tests as threads: the coverage counters are memory
    # shared by a process's threads, and threads contending for them ran a tube level slower
    # at 4 threads (23 s) than at 1 (14 s; plain build: 5.8 s at 1 thread, 1.6 s at 16).
    failed = run_tests(direct, exe, records, snap_id, coverage=True, jobs=jobs, threads=1)
    rest = derive_coverage(records, snap_id, derived)
    if rest:
        print(f"[affected] {len(rest)} whose source's coverage is not current: run instrumented", flush=True)
        failed += run_tests([(n, "coverage") for n in rest], exe, records, snap_id, coverage=True,
                            jobs=jobs, threads=1)
    return failed


def stale(records):
    """The tests whose coverage is older than their verdict (or missing)."""
    return sorted(n for n, r in records["tests"].items()
                  if r.get("files") is None or r.get("cov_snapshot") != r["snapshot"])


def expected_seconds(name, r):
    """A test's expected time in the instrumented program on one thread, for ordering the
    longest first: its last such run, else ten times its plain run (four threads, 2.5 times
    the instrumentation), else a guess by its level's kind and its check."""
    if r and r.get("seconds_cov"):
        return r["seconds_cov"]
    if r and r.get("seconds"):
        return 10 * r["seconds"]
    check, level = name.split("::", 1)
    try:
        with open(os.path.join(ROOT, "levels", f"{level}.json"), encoding="utf-8") as f:
            lvl = json.load(f)
    except OSError:
        return 0.0
    tube = lvl.get("tube") or {}
    kind = 600 if tube.get("metals") else 60 if tube else 30 if any(s.get("beam") for s in lvl.get("shots", [])) else 5
    return kind * {"repeated": 2, "radiation": 0.2}.get(check, 1)


def snapshot_now(records, before):
    snap_id = hashlib.sha256(json.dumps(before, sort_keys=True).encode()).hexdigest()[:16]
    records["snapshots"].setdefault(snap_id, {"universe": before, "global": global_inputs()})
    return snap_id


def run_tests(chosen, exe, records, snap_id, coverage, jobs, threads):
    """Runs `chosen` [(test, why)] in `exe`, `jobs` at a time, each its own process on
    `threads` rayon threads, and records each as it finishes: its verdict (and, with
    `coverage`, the files it executed). The failed tests' names."""
    golden = golden_entries()
    if coverage:
        profdata, cov = llvm_tool("llvm-profdata"), llvm_tool("llvm-cov")
    lock = threading.Lock()
    done, failed, total = [0], [], len(chosen)
    label = "coverage " if coverage else ""

    def one(item):
        name, why = item
        d = os.path.join(PROF, name.replace(":", "_"))
        shutil.rmtree(d, ignore_errors=True)
        os.makedirs(d)
        # One profile per process (`%p`; the merge pool `%m` writes nothing on Windows).
        env = dict(os.environ, LLVM_PROFILE_FILE=os.path.join(d, "%p.profraw"),
                   RAYON_NUM_THREADS=str(threads))
        if coverage:
            env["EM_TEST_FULL"] = "1"  # computed, not taken from a pass record: coverage sees it
        t0 = time.time()
        r = subprocess.run([exe, "--ignored", "--exact", name, "--test-threads", "1"], cwd=ROOT, env=env,
                           capture_output=True, text=True, encoding="utf-8", errors="replace")
        seconds = time.time() - t0
        ok = r.returncode == 0 and "1 passed" in r.stdout
        files = None
        if coverage:
            raws = glob.glob(os.path.join(d, "*.profraw"))
            if not raws:
                raise SystemExit(f"{name}: the instrumented run wrote no profile (exit {r.returncode})")
            merged = os.path.join(d, "t.profdata")
            subprocess.run([profdata, "merge", "-sparse", *raws, "-o", merged], check=True,
                           capture_output=True)
            export = subprocess.run([cov, "export", "-summary-only", f"-instr-profile={merged}", exe],
                                    capture_output=True, text=True, check=True).stdout
            files = sorted(p for f in json.loads(export)["data"][0]["files"]
                           if (p := rel(f["filename"])) is not None and f["summary"]["regions"]["covered"] > 0)
        shutil.rmtree(d, ignore_errors=True)
        level = level_of(name)
        with lock:
            rec = records["tests"].get(name) or {}
            was_ok = rec.get("result") == "ok" and rec.get("snapshot") == snap_id
            rec.update({"snapshot": snap_id, "result": "ok" if ok else "failed",
                        "level_file": sha(f"levels/{level}.json"), "physics": None,
                        "golden": golden.get(level) if name.startswith("golden::") else None,
                        "why": why, "when": datetime.datetime.now().isoformat(timespec="seconds")})
            rec.setdefault("files", None)
            if coverage:
                rec.update({"files": files, "cov_snapshot": snap_id, "seconds_cov": round(seconds, 1)})
            else:
                rec["seconds"] = round(seconds, 1)
            records["tests"][name] = rec
            save_records(records)
            done[0] += 1
            print(f"[{done[0]}/{total}] {label}{name} {'ok' if ok else 'FAILED'} ({seconds:.1f} s"
                  + (f", {len(files)} files)" if coverage else ")"), flush=True)
            if not ok:
                failed.append(name)
                if was_ok and coverage:
                    print(f"    {name} passed in the plain program and failed instrumented: "
                          "a result that depends on the build", flush=True)
                tail = [l for l in (r.stdout + r.stderr).splitlines()
                        if "panicked" in l or "assert" in l or "FAILED" in l][:8]
                print("\n".join("    " + l for l in tail), flush=True)

    with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as pool:
        for f in [pool.submit(one, item) for item in chosen]:
            f.result()
    return failed


def fill_physics(records, program):
    """The levels' physics hashes, asked of the program once (one run of it)."""
    physics = program.physics()
    for name, r in records["tests"].items():
        if r.get("physics") is None:
            r["physics"] = physics.get(level_of(name))
    save_records(records)


def run(args):
    records = load_records()
    before = universe()
    inst = Program()
    inst.get()
    if universe() != before:
        raise SystemExit("engine files changed during the build: run again")
    names = all_tests(records, inst)
    chosen = [(n, "--all") for n in names] if args.all else select(records, inst, names)
    if args.filter:
        chosen = [(n, why) for n, why in chosen if args.filter in n]
    if args.part:
        # Interleaved over the full sorted list (not the selection, which shrinks as a run
        # goes): an interrupted chunk, run again, takes the same tests.
        k, n = (int(x) for x in args.part.split("/"))
        index = {name: i for i, name in enumerate(sorted(names))}
        chosen = [(t, why) for t, why in chosen if index[t] % n == k - 1]
    print(f"[affected] {len(chosen)} of {len(names)} heavy tests selected", flush=True)
    snap_id = snapshot_now(records, before)
    failed = []
    # The verdicts with the plain program: all of them, or with --coverage-only those of the
    # tests whose coverage is another's (they need no instrumented run).
    plain_ones = [t for t in chosen if not args.coverage_only or source_of(t[0]) is not None]
    if plain_ones:
        plain = build(instrumented=False)
        if universe() != before:
            raise SystemExit("engine files changed during the build: run again")
        # The plain program as it runs these tests itself: four at a time on four threads each.
        cpus = os.cpu_count() or 4
        failed = run_tests(plain_ones, plain, records, snap_id, coverage=False, jobs=4,
                           threads=max(1, cpus // 4))
    if not args.no_coverage:
        todo = [(n, why) for n, why in chosen if n not in failed]
        # Coverage is still needed where it is older than the verdict just given.
        todo = [(n, why) for n, why in todo
                if records["tests"].get(n, {}).get("cov_snapshot") != snap_id]
        failed += take_coverage(todo, inst.get(), records, snap_id, args.jobs)
    fill_physics(records, inst)
    left = stale(records)
    print(f"[affected] {len(chosen) - len(set(failed))} ok, {len(set(failed))} failed"
          + (": " + ", ".join(sorted(set(failed))) if failed else "")
          + (f"; coverage to refresh: {len(left)} (python scripts/affected.py refresh)" if left else ""),
          flush=True)
    return 1 if failed else 0


def refresh(args):
    records = load_records()
    before = universe()
    inst = Program()
    inst.get()
    if universe() != before:
        raise SystemExit("engine files changed during the build: run again")
    all_tests(records, inst)
    snap_id = snapshot_now(records, before)
    # Coverage taken now stands for the verdict's state only if nothing changed since.
    todo = [n for n in stale(records) if records["tests"][n]["snapshot"] == snap_id]
    skipped = len(stale(records)) - len(todo)
    if skipped:
        print(f"[affected] {skipped} tests changed since their verdict: `run` them", flush=True)
    failed = take_coverage([(n, "refresh") for n in todo], inst.get(), records, snap_id, args.jobs)
    fill_physics(records, inst)
    return 1 if failed else 0


def changed_files(base):
    diff = subprocess.run(["git", "diff", "--name-only", base], cwd=ROOT, capture_output=True,
                          text=True).stdout.split()
    new = subprocess.run(["git", "ls-files", "--others", "--exclude-standard"], cwd=ROOT,
                         capture_output=True, text=True).stdout.split()
    return sorted(set(diff) | set(new))


def summarize(selected):
    """Selected tests grouped by check: `golden:: 109_ohms_law, 122_inverter`."""
    by = {}
    for name, _ in selected:
        check, level = name.split("::", 1)
        by.setdefault(check, []).append(level)
    return [f"  {c}:: {len(ls)}: " + ", ".join(ls[:12]) + (" ..." if len(ls) > 12 else "")
            for c, ls in sorted(by.items())]


def status(args):
    records = load_records()
    program = Program(fresh=False)
    names = program.tests()
    selected = select(records, program, names)
    print(f"{len(selected)} of {len(names)} heavy tests selected ({len(records['tests'])} recorded)")
    reasons = {}
    for name, why in selected:
        reasons.setdefault(why, []).append(name)
    for why, names in sorted(reasons.items(), key=lambda kv: -len(kv[1])):
        print(f"- {why}: {len(names)}")
    print("\n".join(summarize(selected)))
    left = stale(records)
    if left:
        print(f"coverage to refresh: {len(left)} (python scripts/affected.py refresh)")
    return 0


def note(args):
    files = changed_files(args.base)
    records = load_records()
    program = Program(fresh=False)
    names = program.tests()
    selected = select(records, program, names)
    crates = sorted({c for f in files if f.startswith("crates/")
                     for c in DEPENDENTS.get(f.split("/")[1], [f.split("/")[1]])})
    lines = [f"Affected (scripts/affected.py, against {args.base}):",
             f"- files changed: {len(files)}" + (": " + ", ".join(files[:8]) + (" ..." if len(files) > 8 else "")
                                                if files else ""),
             "- fast tier: " + (", ".join(crates) if crates else "no crate's code") + " (cargo test --workspace)"]
    if selected:
        lines.append(f"- heavy tests selected: {len(selected)} of {len(names)}")
        lines += summarize(selected)
    else:
        lines.append(f"- heavy tests selected: none of {len(names)} "
                     "(every recorded input unchanged, or run since)")
    ran = [r for r in records["tests"].values() if r.get("when")]
    if ran:
        last = max(r["when"] for r in ran)
        recent = [r for r in ran if r["when"][:10] == last[:10]]
        bad = sum(r["result"] != "ok" for r in recent)
        lines.append(f"- heavy tests run on {last[:10]}: {len(recent) - bad} ok, {bad} failed")
    print("\n".join(lines))
    return 0


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("status")
    s.add_argument("--base", default="HEAD")
    r = sub.add_parser("run")
    r.add_argument("--all", action="store_true")
    r.add_argument("--filter")
    r.add_argument("--part")
    r.add_argument("--no-coverage", action="store_true")
    r.add_argument("--coverage-only", action="store_true")
    f = sub.add_parser("refresh")
    f.add_argument("--jobs", type=int, default=os.cpu_count() or 4)
    r.add_argument("--jobs", type=int, default=os.cpu_count() or 4)
    n = sub.add_parser("note")
    n.add_argument("--base", default="HEAD")
    args = ap.parse_args()
    for f in glob.glob(os.path.join(PROF, "discard-*.profraw")) + glob.glob(os.path.join(PROF, "build-*.profraw")):
        try:
            os.remove(f)
        except OSError:
            pass
    sys.exit({"status": status, "run": run, "refresh": refresh, "note": note}[args.cmd](args))


if __name__ == "__main__":
    main()
