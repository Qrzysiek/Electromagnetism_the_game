"""scripts/levels.py: Laminations 20 x 12, the goal at most 0.65 (measured layouts)."""
p = "scripts/levels.py"
s = open(p, encoding="utf-8").read()


def sub(a, b):
    global s
    assert s.count(a) == 1, (s.count(a), a[:70])
    s = s.replace(a, b)


a = s.index("    # Designed with `generator tube-scan` (target/scratch/lam_scan16.py, lam_seeds.py): a")
b = s.index("    return level(\n        \"Laminations\",", a)
s = s[:a] + '''    # Designed with `generator tube-scan` (target/scratch/lam20_*.py): a transformer's core
    # 20 x 12 on the solenoid's axis (30 carriers per area: 14 400 in the verification), the
    # field fixed (B0 = 0.68 at omega 0.125), three cuts of any length and angle (the owner:
    # a bigger plate, cuts not only across the whole metal). Verification (local Ohm's law):
    # uncut 2.833 (3.060); three full cuts along, evenly 0.516 (0.277), at y = 7, 11, 14
    # 0.597, at 6, 10, 14 0.751; two along and one across 0.407; two along 0.755 (0.475);
    # one along and one across 0.791 (0.765); one along 1.260 (0.991); corner to corner 1.435;
    # one across 1.633 (1.692); two half-length cuts along 1.745, two staggered 16 long 1.737.
    # The long cuts' faces carry current across as capacitors (PHYSICS.md 2.11): along them
    # the gain falls short of Ohm's law. The goal, at most 0.65, takes all three cuts across
    # the whole core, along it, about evenly (the best preview of three along read 0.609).
''' + s[b:]
sub('''angle; you have three cuts. Make the core take at most 0.6 of heat per unit time "
        "(averaged over $t = 25$ to $125$).",''', '''angle; you have three cuts. Make the core take at most 0.65 of heat per unit time "
        "(averaged over $t = 25$ to $125$).",''')
sub('''                  heat_goals=[heat_goal(0, 0.0, 0.6, 25.0, 125.0)],''',
    '''                  heat_goals=[heat_goal(0, 0.0, 0.65, 25.0, 125.0)],''')
sub('''        reference=[cut(2, 8, 22, 8), cut(2, 12, 22, 12)])''',
    '''        reference=[cut(2, 8, 22, 8), cut(2, 12, 22, 12), cut(12, 4, 12, 16)])''')
open(p, "w", encoding="utf-8", newline="\n").write(s)
print("ok")
