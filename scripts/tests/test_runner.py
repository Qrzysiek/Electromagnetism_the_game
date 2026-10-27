"""The runner (DEC-24): ledger and fingerprints, hung items, and the controls while it runs."""
import os
import sys
import threading
import time

import pytest

import runner


@pytest.fixture
def jobdir(tmp_path, monkeypatch):
    monkeypatch.setattr(runner, "DIR", str(tmp_path))
    monkeypatch.setattr(runner, "ROOT", str(tmp_path))
    monkeypatch.setattr(runner.Runner, "register", lambda self, *a, **k: None)
    return tmp_path


def py(code):
    return [sys.executable, "-c", code]


def item(iid, code, **kw):
    return runner.Item(id=iid, cmd=py(code), cwd=str(os.getcwd()), **kw)


def run(job, items, **kw):
    return runner.Runner(job, items, quiet=True, **kw).run()


def control_soon(job, words, after):
    """Writes a control command once the job's control file exists and `after` has passed."""
    def go():
        time.sleep(after)
        runner.control(job, words)
    t = threading.Thread(target=go)
    t.start()
    return t


def test_items_run_and_a_rerun_skips_them(jobdir):
    items = [item("a", "pass"), item("b", "pass")]
    assert run("j", items) == {"a": "ok", "b": "ok"}
    assert run("j", items) == {"a": "done", "b": "done"}


def test_a_changed_input_reruns_its_item_alone(jobdir):
    (jobdir / "a.txt").write_text("1")
    (jobdir / "b.txt").write_text("1")
    items = [item("a", "pass", inputs=["a.txt"]), item("b", "pass", inputs=["b.txt"])]
    run("j", items)
    (jobdir / "a.txt").write_text("2")
    assert run("j", items) == {"a": "ok", "b": "done"}


def test_a_failed_item_runs_again(jobdir):
    flag = jobdir / "flag"
    code = f"import sys, os; sys.exit(0 if os.path.exists({str(flag)!r}) else 1)"
    items = [item("f", code)]
    assert run("j", items) == {"f": "failed"}
    flag.write_text("")
    assert run("j", items) == {"f": "ok"}


def test_a_hang_stops_the_job_and_the_fix_continues_from_it(jobdir):
    # One at a time: "before" passes, "hang" hangs, "after" never starts.
    hang = "import time; time.sleep(60)"
    items = [item("before", "pass"), item("hang", hang, estimate=0.1), item("after", "pass")]
    t0 = time.time()
    results = run("j", items, jobs=1, timeout_factor=2.0, min_timeout=1.0)
    assert results == {"before": "ok", "hang": "timeout", "after": "interrupted"}
    assert time.time() - t0 < 30
    incident = runner.read_incidents("j")[-1]
    assert incident["item"] == "hang" and incident["status"] == "timeout"
    # Fixed (its command changed): it runs first, "before" is done, "after" follows.
    fixed = [item("before", "pass"), item("hang", "pass", estimate=0.1), item("after", "pass")]
    order = []
    r = runner.Runner("j", fixed, jobs=1, quiet=True,
                      on_done=lambda it, *a: order.append(it.id))
    assert r.run() == {"before": "done", "hang": "ok", "after": "ok"}
    assert order == ["hang", "after"]


def test_a_hang_may_be_passed_over(jobdir):
    items = [item("hang", "import time; time.sleep(60)", estimate=0.1), item("quick", "pass")]
    results = run("j", items, jobs=1, timeout_factor=2.0, min_timeout=1.0, on_hang="continue")
    assert results == {"hang": "timeout", "quick": "ok"}


def test_a_failure_may_stop_the_job(jobdir):
    items = [item("bad", "import sys; sys.exit(1)"), item("next", "pass")]
    assert run("j", items, jobs=1) == {"bad": "failed", "next": "ok"}
    assert run("j", items, jobs=1, on_fail="stop") == {"bad": "failed", "next": "done"}
    items2 = [item("bad2", "import sys; sys.exit(1)"), item("next2", "pass")]
    assert run("k", items2, jobs=1, on_fail="stop") == {"bad2": "failed", "next2": "interrupted"}


def test_skip_from_the_control_file(jobdir):
    items = [item("long", "import time; time.sleep(60)"), item("after", "pass")]
    t = control_soon("j", ["skip", "long"], 1.0)
    t0 = time.time()
    results = run("j", items, jobs=1)
    t.join()
    assert results == {"long": "skipped", "after": "ok"}
    assert time.time() - t0 < 30
    # Skipped is not passed: the next run tries it again.
    assert runner.read_ledger("j")["long"]["status"] == "skipped"


def test_stop_keeps_what_finished_and_the_rest_runs_next_time(jobdir):
    items = [item("first", "pass"), item("long", "import time; time.sleep(60)"), item("last", "pass")]
    t = control_soon("j", ["stop"], 1.5)
    results = run("j", items, jobs=1)
    t.join()
    assert results["first"] == "ok" and results["long"] == "interrupted"
    assert results["last"] == "interrupted"
    # Resumed: the finished one is done, the others run.
    fast = [item("first", "pass"), item("long", "pass"), item("last", "pass")]
    # (the same ids; "long" now quick, and its command changed, so its fingerprint did too)
    assert run("j", fast) == {"first": "done", "long": "ok", "last": "ok"}


def test_pause_holds_new_items_until_resumed(jobdir):
    stamp = jobdir / "second_started"
    items = [item("one", "import time; time.sleep(1.5)"),
             item("two", f"open({str(stamp)!r}, 'w').write(str(__import__('time').time()))")]
    t1 = control_soon("j", ["pause"], 0.5)
    t2 = control_soon("j", ["resume"], 3.5)
    t0 = time.time()
    results = run("j", items, jobs=1)
    t1.join()
    t2.join()
    assert results == {"one": "ok", "two": "ok"}
    # "one" ended at ~1.5 s; "two" waited for the resume at ~3.5 s.
    assert float(stamp.read_text()) - t0 > 3.0


def test_a_saved_ledger_loaded_skips_what_it_holds(jobdir):
    items = [item("a", "pass")]
    run("j", items)
    saved = jobdir / "saved.jsonl"
    runner.state("j", save=str(saved))
    os.remove(runner.ledger_path("j"))
    assert run("j", items) == {"a": "ok"}
    os.remove(runner.ledger_path("j"))
    runner.state("j", load=str(saved))
    assert run("j", items) == {"a": "done"}
