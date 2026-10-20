"""Progress of long computations (the /progress skill: .claude/skills/progress/SKILL.md).

    python scripts/progress.py run <job> --total N [--each S] [--match REGEX] [--idle] -- CMD...
        runs CMD, its output (stdout and stderr) into target/progress/<job>.log; every output
        line matching REGEX counts as one finished item of N; --each is the estimated time
        per item, in seconds (before it starts: see `estimate` and the skill's timings)
    python scripts/progress.py chain <job> --plan PLAN.json [--idle]
        runs the plan's chunks one after another, each as `run` would (its log
        target/progress/<job>-<k>.log), recording each one that finishes; a failed chunk stops
        the chain; run again, it resumes at the first chunk not finished. The plan:
        {"chunks": [{"name": ..., "cmd": [...], "total": N, "each": S, "match": REGEX}, ...]}
    python scripts/progress.py watch <job> --log FILE --total N [--each S] [--match REGEX]
        tracks a computation started otherwise, from its log file
    python scripts/progress.py status [--all]
        each job and chain: done/total, elapsed, time per item measured against the
        estimate, ETA, its log's last line (read-only); --all: finished ones too
    python scripts/progress.py estimate <level.json> [--refine N]
        a metal tube level's run to its goals' end: carriers, steps, estimated seconds

The counts are read from the logs whenever the status is asked for, so a status is right
even if the wrapper stopped.
"""
import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DIR = os.path.join(ROOT, "target", "progress")
# A metal tube run's wall time per carrier and step on the owner's 16-thread machine at full
# speed (2026-10-20: Laminations' preview 6400 carriers x 3571 steps in 36 s, its
# verification 12 800 x 7143 in 124 s; Magnetoresistance 5200 x 2286 in 17 s, 10 400 x 4571
# in 76 s): 1.4-1.6e-6 s.
CARRIER_STEP_SECONDS = 1.5e-6


def job_path(job):
    return os.path.join(DIR, f"{job}.json")


def save(job, info):
    os.makedirs(DIR, exist_ok=True)
    tmp = job_path(job) + ".tmp"
    with open(tmp, "w", encoding="utf-8") as f:
        json.dump(info, f, indent=1)
    for _ in range(20):
        try:
            os.replace(tmp, job_path(job))
            return
        except PermissionError:  # Windows: the status may be reading it
            time.sleep(0.05)


def run_one(job, cmd, total, each, match, idle, chain=None):
    """Runs one job (its output in target/progress/<job>.log); its exit code."""
    os.makedirs(DIR, exist_ok=True)
    log = os.path.join(DIR, f"{job}.log")
    cmd = list(cmd)
    info = {"job": job, "command": cmd, "log": log, "total": total, "each": each,
            "match": match, "started": time.time(), "state": "running", "chain": chain}
    flags = subprocess.IDLE_PRIORITY_CLASS if idle and os.name == "nt" else 0
    # Windows: "python" may be a launcher alias CreateProcess cannot find; resolve the rest.
    if cmd[0] in ("python", "python3"):
        cmd[0] = sys.executable
    else:
        cmd[0] = shutil.which(cmd[0]) or cmd[0]
    with open(log, "w", encoding="utf-8") as out:
        p = subprocess.Popen(cmd, cwd=os.getcwd(), stdout=out, stderr=subprocess.STDOUT,
                             creationflags=flags)
        info["pid"] = p.pid
        save(job, info)
        code = p.wait()
    info.update(state="done" if code == 0 else "failed", exit=code, finished=time.time())
    save(job, info)
    print(f"{job}: {info['state']} (exit {code}) after {fmt(time.time() - info['started'])}",
          flush=True)
    return code


def run(args):
    sys.exit(run_one(args.job, args.cmd, args.total, args.each, args.match, args.idle))


def chain(args):
    """The plan's chunks in order, skipping those a previous run finished."""
    plan = json.load(open(args.plan, encoding="utf-8"))["chunks"]
    state_path = os.path.join(DIR, f"{args.job}.chain.json")
    state = {"plan": os.path.abspath(args.plan), "done": []}
    if os.path.exists(state_path):
        state["done"] = json.load(open(state_path, encoding="utf-8")).get("done", [])
    os.makedirs(DIR, exist_ok=True)
    state.update(pid=os.getpid(), chunks=len(plan), started=time.time(), state="running")
    write_chain(args.job, state)
    for k, c in enumerate(plan, start=1):
        name = c.get("name", str(k))
        if name in state["done"]:
            print(f"{args.job}: chunk {k} ({name}) already done", flush=True)
            continue
        code = run_one(f"{args.job}-{k}", c["cmd"], c.get("total"), c.get("each"),
                       c.get("match"), args.idle, chain=args.job)
        if code != 0:
            state.update(state=f"stopped at chunk {k} ({name}): exit {code}", finished=time.time())
            write_chain(args.job, state)
            sys.exit(code)
        state["done"].append(name)
        write_chain(args.job, state)
    state.update(state="done", finished=time.time())
    write_chain(args.job, state)
    print(f"{args.job}: all {len(plan)} chunks done", flush=True)


def write_chain(job, state):
    path = os.path.join(DIR, f"{job}.chain.json")
    tmp = path + ".tmp"
    with open(tmp, "w", encoding="utf-8") as f:
        json.dump(state, f, indent=1)
    os.replace(tmp, path)


def watch(args):
    info = {"job": args.job, "command": None, "log": os.path.abspath(args.log),
            "total": args.total, "each": args.each, "match": args.match,
            "started": os.path.getmtime(args.log) if args.since_log else time.time(),
            "state": "watching"}
    save(args.job, info)
    print(f"{args.job}: watching {info['log']}")


def fmt(s):
    if s is None:
        return "?"
    s = int(round(s))
    return f"{s // 3600}:{s % 3600 // 60:02d}:{s % 60:02d}" if s >= 3600 else f"{s // 60}:{s % 60:02d}"


def alive(pid):
    if pid is None:
        return False
    if os.name == "nt":
        out = subprocess.run(["tasklist", "/FI", f"PID eq {pid}", "/NH"], capture_output=True,
                             text=True).stdout
        return str(pid) in out
    try:
        os.kill(pid, 0)
        return True
    except OSError:
        return False


def measure(info):
    """Done items and their times (the log's modification time for the last one)."""
    log = info["log"]
    if not os.path.exists(log):
        return 0, None, ""
    pattern = re.compile(info["match"]) if info.get("match") else None
    done, last = 0, ""
    with open(log, encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.rstrip()
            if line:
                last = line
            if pattern and pattern.search(line):
                done += 1
    return done, os.path.getmtime(log), last


def status(args):
    if not os.path.isdir(DIR):
        print("no jobs")
        return
    now = time.time()
    rows = []
    for name in sorted(os.listdir(DIR)):
        if not name.endswith(".json") or name.endswith((".chain.json", ".plan.json")):
            continue
        try:
            info = json.load(open(os.path.join(DIR, name), encoding="utf-8"))
        except (json.JSONDecodeError, OSError):
            continue
        if not isinstance(info, dict) or "log" not in info:
            continue  # not a job's record
        state = info.get("state")
        if state == "running" and not alive(info.get("pid")):
            state = "lost (its process ended without a word)"
        if state in ("done", "failed") and not args.all and now - info.get("finished", now) > 3600:
            continue
        done, touched, last = measure(info)
        total, each = info.get("total"), info.get("each")
        end = info.get("finished", now)
        elapsed = end - info["started"]
        per = elapsed / done if done else None
        eta = None
        if state in ("running", "watching") and total:
            if per:
                eta = per * max(total - done, 0)
            elif each:
                eta = max(each * total - elapsed, 0)
        quiet = now - touched if touched else None
        rows.append((info["job"], state, done, total, elapsed, per, each, eta, quiet, last))
    lines = ["| job | state | done | elapsed | per item (estimate) | ETA | log quiet for | last line |",
             "|---|---|---|---|---|---|---|---|"]
    for job, state, done, total, elapsed, per, each, eta, quiet, last in rows:
        flag = ""
        if per and each and per > 1.5 * each:
            flag = " SLOW"
        lines.append(f"| {job} | {state} | {done}/{total or '?'} | {fmt(elapsed)} | "
                     f"{fmt(per)} ({fmt(each)}){flag} | {fmt(eta)} | {fmt(quiet)} | "
                     f"{last[:90].replace('|', '/')} |")
    # The chains: chunks done, the one running, and the rest's estimate.
    chains = []
    for name in sorted(os.listdir(DIR)):
        if not name.endswith(".chain.json"):
            continue
        st = json.load(open(os.path.join(DIR, name), encoding="utf-8"))
        job = name[: -len(".chain.json")]
        if st.get("state") == "running" and not alive(st.get("pid")):
            st["state"] = "lost (run the chain again to resume)"
        if st.get("state") == "done" and not args.all and now - st.get("finished", now) > 3600:
            continue
        try:
            plan = json.load(open(st["plan"], encoding="utf-8"))["chunks"]
        except (OSError, KeyError, json.JSONDecodeError):
            plan = []
        left = [c for c in plan if c.get("name") not in st.get("done", [])]
        rest = sum((c.get("each") or 0) * (c.get("total") or 0) for c in left[1:])
        current = next((r for r in rows if r[0] == f"{job}-{plan.index(left[0]) + 1}"), None) \
            if left else None
        eta = (current[7] or 0) + rest if current else rest
        chains.append(f"| {job} | {st.get('state')} | {len(st.get('done', []))}/{len(plan)} chunks | "
                      f"{left[0].get('name') if left else '-'} | {fmt(eta) if left else '-'} |")
    if chains:
        lines += ["", "| chain | state | done | current chunk | ETA (all chunks) |",
                  "|---|---|---|---|---|"] + chains
    print("\n".join(lines))


def estimate(args):
    lvl = json.load(open(args.level, encoding="utf-8"))
    t = lvl.get("tube")
    if not t:
        print("not a tube level")
        return
    metals = t.get("metals", [])
    if not metals:
        print("no metals: a beam's macroparticles depend on its emission; time one run")
        return
    ends = [g["end"] for key in ("goals", "potential_goals", "heat_goals") for g in t.get(key, [])]
    end = max(ends) if ends else lvl["physics"].get("t_max", 200.0)
    for refine in ([args.refine] if args.refine else [1, 2]):
        carriers = sum((m["corners"][1][0] - m["corners"][0][0])
                       * (m["corners"][1][1] - m["corners"][0][1])
                       * m["carriers_per_area"] for m in metals) * t.get("particles", 1.0) * refine
        steps = end / (t["step"] / refine)
        s = CARRIER_STEP_SECONDS * carriers * steps
        name = "preview" if refine == 1 else "verification" if refine == 2 else f"refine {refine}"
        print(f"{name}: {carriers:.0f} carriers x {steps:.0f} steps (to t = {end:g}) "
              f"~ {s:.0f} s at full speed")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="what", required=True)
    r = sub.add_parser("run")
    r.add_argument("job")
    r.add_argument("--total", type=int)
    r.add_argument("--each", type=float)
    r.add_argument("--match")
    r.add_argument("--idle", action="store_true")
    c = sub.add_parser("chain")
    c.add_argument("job")
    c.add_argument("--plan", required=True)
    c.add_argument("--idle", action="store_true")
    w = sub.add_parser("watch")
    w.add_argument("job")
    w.add_argument("--log", required=True)
    w.add_argument("--total", type=int)
    w.add_argument("--each", type=float)
    w.add_argument("--match")
    w.add_argument("--since-log", action="store_true",
                   help="time from the log's last change, not from now")
    s = sub.add_parser("status")
    s.add_argument("--all", action="store_true")
    e = sub.add_parser("estimate")
    e.add_argument("level")
    e.add_argument("--refine", type=int)
    # The command follows the first "--" (argparse's remainder took the options too).
    argv = sys.argv[1:]
    cmd = []
    if "--" in argv:
        k = argv.index("--")
        argv, cmd = argv[:k], argv[k + 1:]
    args = ap.parse_args(argv)
    if args.what == "run":
        if not cmd:
            ap.error("run: give the command after --")
        args.cmd = cmd
        run(args)
    elif args.what == "chain":
        chain(args)
    elif args.what == "watch":
        watch(args)
    elif args.what == "status":
        status(args)
    else:
        estimate(args)


if __name__ == "__main__":
    main()
