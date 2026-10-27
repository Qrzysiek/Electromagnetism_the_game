"""Long runs of independent items, controllable while they run (DEC-24, docs/PROGRESS.md).

    python scripts/runner.py run <job> --plan PLAN.json [--jobs N] [--idle]
                                 [--on-hang stop|continue] [--on-fail continue|stop]
        runs the plan's items, each its own process, N at a time; skips those that passed
        before with the same inputs (the job's ledger), and runs first those that stopped or
        failed it last time; one `[k/N] <item> <status> (<s> s)` line each, in
        target/progress/<job>.log too (`/progress` shows the job). A hung item stops the job
        (`--on-hang continue` goes on): the next items may well hang too.
    python scripts/runner.py ctl <job> pause | resume | stop | skip ITEM | jobs N
        while the job runs: pause (no new item starts; the running ones finish), resume, stop
        (the running items are killed and recorded as interrupted, the job ends), skip an item
        (killed if running, never started if pending; recorded as skipped), change how many
        run at once
    python scripts/runner.py state <job> [--save FILE] [--load FILE]
        the job's ledger (each item's last outcome) and its incidents (where it stopped, and
        why); --save copies the ledger, --load puts a saved one in its place (a later `run`
        then skips what it holds as passed, inputs unchanged)

When a job stops on a hang (or a failure, with `--on-fail stop`), what finished is in its
ledger, and its incident log (target/progress/<job>.incidents.jsonl) names the item, how long
it ran against its limit and the end of its output. Fix it and run the job again: the item
that stopped it runs first, every finished item is checked against its fingerprint (one whose
inputs the fix changed runs again), and the rest continue.

A plan (JSON): {"inputs": [the job's own inputs: a program's binary, ...],
                "timeout_factor": 4, "min_timeout": 60,
                "items": [{"id": "...", "cmd": [...], "inputs": [...], "estimate": seconds,
                           "env": {...}}, ...]}
Paths in `inputs` are files or glob patterns. An item's fingerprint is the hash of its command,
its environment additions and the contents of every input (the job's and its own): a rerun
after a change runs again exactly the items whose inputs changed, and the others count as
done. An item that runs longer than `timeout_factor` times its expected time (its last
passing run's, else its estimate; at least `min_timeout`) is killed as hung, recorded as a
timeout, and the job goes on. The ledger (target/progress/<job>.ledger.jsonl) takes each
outcome as it comes, so nothing finished is lost to a crash, a stop or a hang.

The library (`Item`, `Runner`) serves other scripts too (`affected.py`): an `on_done` callback
gets each finished item with its output.
"""
import argparse
import datetime
import glob
import hashlib
import json
import os
import shutil
import signal
import subprocess
import sys
import threading
import time
from dataclasses import dataclass, field

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DIR = os.path.join(ROOT, "target", "progress")
OK, FAILED, TIMEOUT, SKIPPED, INTERRUPTED = "ok", "failed", "timeout", "skipped", "interrupted"


@dataclass
class Item:
    id: str
    cmd: list
    inputs: list = field(default_factory=list)
    estimate: float | None = None
    env: dict = field(default_factory=dict)
    cwd: str = ROOT


def file_hashes(patterns):
    """`(path, sha256)` of every file the patterns name, sorted; a missing one hashes as such."""
    out = []
    for pattern in patterns:
        path = pattern if os.path.isabs(pattern) else os.path.join(ROOT, pattern)
        found = sorted(glob.glob(path)) if any(c in pattern for c in "*?[") else [path]
        if not found:
            out.append((pattern, "no match"))
        for p in found:
            try:
                with open(p, "rb") as f:
                    out.append((os.path.relpath(p, ROOT), hashlib.sha256(f.read()).hexdigest()))
            except OSError:
                out.append((os.path.relpath(p, ROOT), "missing"))
    return out


def fingerprint(item, job_inputs):
    h = hashlib.sha256()
    h.update(json.dumps([item.cmd, sorted(item.env.items())]).encode())
    for path, digest in file_hashes(list(job_inputs) + list(item.inputs)):
        h.update(f"{path}\0{digest}\0".encode())
    return h.hexdigest()


def ledger_path(job):
    return os.path.join(DIR, f"{job}.ledger.jsonl")


def incidents_path(job):
    return os.path.join(DIR, f"{job}.incidents.jsonl")


def read_incidents(job):
    try:
        with open(incidents_path(job), encoding="utf-8") as f:
            return [json.loads(line) for line in f if line.strip()]
    except FileNotFoundError:
        return []


def control_path(job):
    return os.path.join(DIR, f"{job}.ctl")


def read_ledger(job):
    """Each item's last entry."""
    out = {}
    try:
        with open(ledger_path(job), encoding="utf-8") as f:
            for line in f:
                if line.strip():
                    e = json.loads(line)
                    out[e["id"]] = e
    except FileNotFoundError:
        pass
    return out


def kill_tree(p):
    """The item's process and its children (a test program's, cargo's)."""
    if p.poll() is not None:
        return
    if os.name == "nt":
        subprocess.run(["taskkill", "/F", "/T", "/PID", str(p.pid)], capture_output=True)
    else:
        try:
            os.killpg(p.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
    try:
        p.wait(timeout=30)
    except subprocess.TimeoutExpired:
        pass


class Runner:
    """Runs items `jobs` at a time under the job's controls; see the module's documentation.
    `on_done(item, status, seconds, output)` is called for each finished item (in the runner's
    thread, one at a time); `ledger=False` leaves the outcomes to the caller (no done-before
    skips)."""

    def __init__(self, job, items, jobs=4, job_inputs=(), timeout_factor=4.0, min_timeout=60.0,
                 on_done=None, ledger=True, idle=False, quiet=False, on_hang="stop",
                 on_fail="continue"):
        self.job, self.items, self.jobs = job, list(items), jobs
        self.on_hang, self.on_fail = on_hang, on_fail
        self.job_inputs, self.timeout_factor, self.min_timeout = list(job_inputs), timeout_factor, min_timeout
        self.on_done, self.use_ledger, self.idle, self.quiet = on_done, ledger, idle, quiet
        self.results = {}
        os.makedirs(DIR, exist_ok=True)
        self.log_path = os.path.join(DIR, f"{job}.log")
        self.out_dir = os.path.join(DIR, f"{job}.items")

    def say(self, line):
        with open(self.log_path, "a", encoding="utf-8") as f:
            f.write(line + "\n")
        if not self.quiet:
            print(line, flush=True)

    def register(self, state, total, extra=None):
        """The job's record as `progress.py` keeps it (`/progress`)."""
        sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        import progress  # noqa: E402 (a sibling script)
        info = {"job": self.job, "command": ["runner", self.job], "log": self.log_path,
                "total": total, "each": self.each, "match": r"^\[\d+/\d+\] ",
                "started": self.started, "state": state, "chain": None, "pid": os.getpid()}
        info.update(extra or {})
        progress.save(self.job, info)

    def run(self):
        """Runs the items; returns {id: status} for those run now (done-before ones: "done")."""
        with open(self.log_path, "w", encoding="utf-8"):
            pass
        shutil.rmtree(self.out_dir, ignore_errors=True)
        os.makedirs(self.out_dir, exist_ok=True)
        ctl = control_path(self.job)
        with open(ctl, "w", encoding="utf-8"):
            pass  # a fresh control file: earlier commands do not carry over
        ledger = read_ledger(self.job) if self.use_ledger else {}
        prints = {it.id: fingerprint(it, self.job_inputs) for it in self.items} if self.use_ledger else {}
        pending = []
        for it in self.items:
            e = ledger.get(it.id)
            if e and e["status"] == OK and e["fingerprint"] == prints.get(it.id):
                self.results[it.id] = "done"
            else:
                pending.append(it)
        done_before = len(self.items) - len(pending)
        if done_before:
            self.say(f"[runner] {done_before} of {len(self.items)} passed before with the same inputs")
        # What stopped or failed the job last time runs first: whether the fix works shows at once.
        last = read_incidents(self.job)
        first = {e["item"] for e in last[-1:]} | {i for i, e in ledger.items()
                                                  if e["status"] in (FAILED, TIMEOUT)}
        pending.sort(key=lambda it: it.id not in first)
        if last and any(it.id in first for it in pending):
            e = last[-1]
            self.say(f"[runner] last stopped at {e['item']} ({e['status']} after {e['seconds']:.0f} s, "
                     f"{e['when']}): it runs first")
        total = len(pending)
        known = [self.expected(it, ledger) for it in pending if self.expected(it, ledger)]
        self.each = sum(known) / len(known) / max(1, self.jobs) if known else None
        self.started = time.time()
        self.register("running", total)
        running = {}  # id -> (item, process, start, output file handle)
        paused = stopping = False
        skip = set()
        offset = 0
        count = 0
        while pending or running:
            # The controls written since the last look.
            with open(ctl, encoding="utf-8") as f:
                f.seek(offset)
                commands = f.read()
                offset = f.tell()
            for line in commands.splitlines():
                words = line.split(maxsplit=1)
                if not words:
                    continue
                verb, arg = words[0], (words[1].strip() if len(words) > 1 else "")
                if verb == "pause":
                    paused = True
                    self.say("[runner] paused: no new item starts")
                    self.register("paused", total)
                elif verb == "resume":
                    paused = False
                    self.say("[runner] resumed")
                    self.register("running", total)
                elif verb == "stop":
                    stopping = True
                    self.say("[runner] stopping")
                elif verb == "skip":
                    skip.add(arg)
                elif verb == "jobs" and arg.isdigit():
                    self.jobs = max(1, int(arg))
                    self.say(f"[runner] {self.jobs} at a time")
            # Skipped before they start.
            for it in [it for it in pending if it.id in skip]:
                pending.remove(it)
                count += 1
                self.finish(it, SKIPPED, 0.0, "", count, total)
            # Running items: finished, skipped, hung, or stopped.
            for iid, (it, p, t0, out) in list(running.items()):
                seconds = time.time() - t0
                status = None
                if p.poll() is not None:
                    status = OK if p.returncode == 0 else FAILED
                elif stopping:
                    kill_tree(p)
                    status = INTERRUPTED
                elif iid in skip:
                    kill_tree(p)
                    status = SKIPPED
                else:
                    limit = self.limit(it, ledger)
                    if limit is not None and seconds > limit:
                        kill_tree(p)
                        status = TIMEOUT
                if status:
                    out.close()
                    with open(out.name, encoding="utf-8", errors="replace") as f:
                        output = f.read()
                    del running[iid]
                    count += 1
                    self.finish(it, status, seconds, output, count, total)
                    stop_here = (status == TIMEOUT and self.on_hang == "stop") or \
                        (status == FAILED and self.on_fail == "stop")
                    if stop_here and not stopping:
                        self.incident(it, status, seconds, self.limit(it, ledger), output)
                        stopping = True
            if stopping:
                for it in pending:
                    self.results[it.id] = INTERRUPTED
                pending = []
            # New items while there is room.
            while pending and not paused and not stopping and len(running) < self.jobs:
                it = pending.pop(0)
                out = open(os.path.join(self.out_dir, safe(it.id) + ".log"), "w", encoding="utf-8")
                env = dict(os.environ, **it.env)
                kwargs = {}
                if os.name == "nt":
                    kwargs["creationflags"] = subprocess.IDLE_PRIORITY_CLASS if self.idle else 0
                else:
                    kwargs["start_new_session"] = True
                cmd = list(it.cmd)
                if cmd[0] in ("python", "python3"):
                    cmd[0] = sys.executable
                p = subprocess.Popen(cmd, cwd=it.cwd, env=env, stdout=out, stderr=subprocess.STDOUT,
                                     **kwargs)
                running[it.id] = (it, p, time.time(), out)
            time.sleep(0.2)
        bad = [i for i, s in self.results.items() if s not in (OK, "done")]
        self.say(f"[runner] {sum(s == OK for s in self.results.values())} ok, "
                 f"{sum(s == 'done' for s in self.results.values())} done before, "
                 f"{len(bad)} not: " + (", ".join(f"{i} {self.results[i]}" for i in bad) if bad else "none"))
        self.register("done" if not bad else "failed", total,
                      {"finished": time.time(), "exit": 0 if not bad else 1})
        return self.results

    def incident(self, it, status, seconds, limit, output):
        """Where and why the job stops: the incident log, and a word on what to do."""
        entry = {"item": it.id, "status": status, "seconds": round(seconds, 1),
                 "limit": None if limit is None else round(limit, 1),
                 "when": datetime.datetime.now().isoformat(timespec="seconds"),
                 "cmd": it.cmd, "output_tail": output.splitlines()[-40:]}
        with open(incidents_path(self.job), "a", encoding="utf-8") as f:
            f.write(json.dumps(entry) + "\n")
        why = (f"no end after {seconds:.0f} s (its limit {limit:.0f} s)" if status == TIMEOUT
               else "it failed")
        self.say(f"[runner] stopping at {it.id}: {why}. What finished is kept; fix it and run "
                 f"the job again (it runs first; finished items whose inputs the fix changed run "
                 f"again). Its output's end: python scripts/runner.py state {self.job}")

    def expected(self, it, ledger):
        e = ledger.get(it.id)
        if e and e.get("status") == OK and e.get("seconds"):
            return e["seconds"]
        return it.estimate

    def limit(self, it, ledger):
        exp = self.expected(it, ledger)
        if exp is None:
            return None
        return max(self.min_timeout, self.timeout_factor * exp)

    def finish(self, it, status, seconds, output, count, total):
        self.results[it.id] = status
        if self.use_ledger:
            entry = {"id": it.id, "fingerprint": fingerprint(it, self.job_inputs), "status": status,
                     "seconds": round(seconds, 2),
                     "when": datetime.datetime.now().isoformat(timespec="seconds")}
            with open(ledger_path(self.job), "a", encoding="utf-8") as f:
                f.write(json.dumps(entry) + "\n")
        word = status if status != OK else "ok"
        extra = None
        if self.on_done:
            extra = self.on_done(it, status, seconds, output)
        self.say(f"[{count}/{total}] {it.id} {word} ({seconds:.1f} s{', ' + extra if extra else ''})")


def safe(item_id):
    return "".join(c if c.isalnum() or c in "-_." else "_" for c in item_id)


def control(job, words):
    """Appends a command to a running job's control file."""
    path = control_path(job)
    if not os.path.exists(path):
        sys.exit(f"{job}: no control file (is it running?)")
    with open(path, "a", encoding="utf-8") as f:
        f.write(" ".join(words) + "\n")
    print(f"{job}: {' '.join(words)}")


def state(job, save=None, load=None):
    if load:
        shutil.copy(load, ledger_path(job))
        print(f"{job}: ledger loaded from {load}")
    entries = read_ledger(job)
    counts = {}
    for e in entries.values():
        counts[e["status"]] = counts.get(e["status"], 0) + 1
    print(f"{job}: {len(entries)} items in the ledger: "
          + ", ".join(f"{n} {s}" for s, n in sorted(counts.items())))
    for e in entries.values():
        if e["status"] != OK:
            print(f"  {e['id']}: {e['status']} ({e['when']})")
    for e in read_incidents(job)[-3:]:
        print(f"incident {e['when']}: {e['item']} {e['status']} after {e['seconds']} s"
              + (f" (limit {e['limit']} s)" if e.get("limit") else "") + "; its output's end:")
        print("\n".join("    " + line for line in e["output_tail"][-12:]))
    if save:
        shutil.copy(ledger_path(job), save)
        print(f"{job}: ledger saved to {save}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("job")
    r.add_argument("--plan", required=True)
    r.add_argument("--jobs", type=int, default=4)
    r.add_argument("--idle", action="store_true")
    r.add_argument("--on-hang", choices=["stop", "continue"], default="stop")
    r.add_argument("--on-fail", choices=["continue", "stop"], default="continue")
    c = sub.add_parser("ctl")
    c.add_argument("job")
    c.add_argument("words", nargs="+")
    s = sub.add_parser("state")
    s.add_argument("job")
    s.add_argument("--save")
    s.add_argument("--load")
    args = ap.parse_args()
    if args.cmd == "ctl":
        control(args.job, args.words)
        return 0
    if args.cmd == "state":
        state(args.job, args.save, args.load)
        return 0
    with open(args.plan, encoding="utf-8") as f:
        plan = json.load(f)
    items = [Item(id=i["id"], cmd=i["cmd"], inputs=i.get("inputs", []), estimate=i.get("estimate"),
                  env=i.get("env", {})) for i in plan["items"]]
    results = Runner(args.job, items, jobs=args.jobs, job_inputs=plan.get("inputs", []),
                     timeout_factor=plan.get("timeout_factor", 4.0), min_timeout=plan.get("min_timeout", 60.0),
                     idle=args.idle, on_hang=args.on_hang, on_fail=args.on_fail).run()
    return 0 if all(s in (OK, "done") for s in results.values()) else 1


if __name__ == "__main__":
    sys.exit(main())
