"""The plan of Laminations' layout scan: previews, then verifications (a chain)."""
import json
import os

ROOT = r"C:\Users\qrzys\Documents\Projects\Electromagnetism_the_game\.claude\worktrees\tubes-next"
S = os.path.join(ROOT, "target", "scratch")
GEN = os.path.join(ROOT, r"target\claude-gen\release\generator.exe")
level = os.path.join(S, "draft_laminations.json")
placements = os.path.join(S, "lam20_placements.json")
plan = {"chunks": [
    {"name": "previews", "cmd": [GEN, "tube-scan", level, placements, "--refine", "1"],
     "total": 9, "each": 60, "match": r"^\d+: "},
    {"name": "verifications", "cmd": [GEN, "tube-scan", level, placements, "--refine", "2"],
     "total": 9, "each": 180, "match": r"^\d+: "},
]}
path = os.path.join(ROOT, "target", "progress", "lam20.plan.json")
os.makedirs(os.path.dirname(path), exist_ok=True)
json.dump(plan, open(path, "w", encoding="utf-8"), indent=1)
print(path)
