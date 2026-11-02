"""A scratch copy of Eddy currents with copper's regime (design only; the level comes from
levels.py), and swing placements to scan."""
import json
import math

l = json.load(open("levels/115_eddy_currents.json", encoding="utf-8"))
l["name"] = "Eddy currents in copper (scratch)"
m = l["tube"]["metals"][0]
m["mean_free_path"] = 0.3
m["density"] = 2.7
m["core_radius"] = 0.0025
l["tube"]["solenoid_drive"]["omega"] = 0.13 * math.pi
l["tube"]["step"] = 0.033
l["limits"]["solenoid_fields"] = [-1.5, 0.43, 1.5]
for e in l["reference_solution"]:
    e["value"] = 0.43
json.dump(l, open("target/scratch/copper.json", "w", encoding="utf-8"), indent=1)
node = l["reference_solution"][0]["node"]
json.dump([[{"node": node, "kind": "supply", "value": b}] for b in [0.1, 0.2, 0.3, 0.43, 0.6, 0.8, 1.2]],
          open("target/scratch/copper_b.json", "w"))
print("ok", node)
