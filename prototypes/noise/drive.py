"""The prototype with its floating contact driven instead (tunable, no circuit): the strip's
resistance from the current at a supply."""
import json, sys
import proto
name, length, width, end, seed = sys.argv[1], float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]), int(sys.argv[5])
lvl = proto.build(length, width, 0.0, end, None, seed, 25.0)
lvl["electrodes"][1]["tunable"] = True
lvl["tube"]["circuits"] = []
json.dump(lvl, open(name + ".json", "w", encoding="utf-8"), indent=1, ensure_ascii=False)
node = lvl["electrodes"][1]["center"]
json.dump([[{"node": node, "kind": "supply", "value": v}] for v in (4.0, -4.0)], open(name + "_pl.json", "w"))
