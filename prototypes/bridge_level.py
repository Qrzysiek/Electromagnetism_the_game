"""Adds the Wheatstone bridge (metals arc, Intermediate, after Magnetoresistance) to
scripts/levels.py and its reality entry to scripts/reality.py. Run once."""

LEVEL = '''

def wheatstone_bridge():
    # Designed with `generator tube-scan` (target/scratch/bridge_*.py): Magnetoresistance's two
    # strips on one fixed supply at their right ends, each left end to ground through an
    # ideal resistor (R0 the long strip's, R2 the short wide one's); the output, a voltmeter
    # between the two left ends, is zero when X / Y = R2 / R0 (X, Y the short and the long
    # strip's resistances). From the ends' potentials (V 10, R0 10, R2 2.8, preview,
    # t = 10-80): X / Y 0.22, 0.28, 0.30, 0.33, 0.44 at B = 8, 10, 12, 14, 16 (X 2.6 to 6.5,
    # Y 11.8 to 14.8: the long strip's own magnetoresistance, smaller), the output -2.0 at
    # B = 0, 0.00 at 10, +1.04 at 16. The output moves ~0.17 per unit of B near the balance,
    # and the preview's and verification's differed by up to 0.13 (the long strip's end
    # noisier, ~0.1, through R0 = 10 than the short one's, ~0.05, through 2.8); halving the
    # resistors (5 and 1.5, the same ratio), the supply 12 and the window to t = 150 made it
    # monotonic (preview -0.47, -0.36, -0.07, 0.03, 0.16, 0.36, 0.46 at B = 9..15).
    # Verification at 10.5, 12, 13.5: -0.270 +- 0.006, 0.061 +- 0.049, 0.238 +- 0.146 (the
    # two ends' noise largely common: they moved together by ~0.1). Goal: the voltmeter
    # within +-0.3, B about 10.3-13 (either sign: the bridge reads |B|).
    return level(
        "Wheatstone bridge",
        "The classic way to measure a small change of resistance, here a field sensor. "
        "Magnetoresistance's two strips share a supply of 12 at their right ends: the long "
        "one, whose resistance barely changes in a field, and the short wide one, whose "
        "resistance grows with the field's strength. Each strip's left end goes to ground "
        "through a resistor, 5 for the long strip and 1.5 for the short one, and a voltmeter "
        "reads the difference between those two ends. It reads zero, the bridge balanced, "
        "when the strips' resistances stand in the resistors' ratio: the short one 1.5/5 of "
        "the long one. Set your solenoid so that the voltmeter reads between -0.3 and 0.3 "
        "(averaged over $t = 10$ to $150$). Real sensors work this way: magnetoresistors in a "
        "bridge, its small output read by an amplifier.",
        grid=(24, 20), c=None, t_max=170.0,
        # 0, 1: the long strip's left and right ends; 2, 3: the short wide strip's. The right
        # ends on the supply; the left ends to ground through 5 and 1.5.
        electrodes=[contact(3, 15, 4), contact(21, 15, 4),
                    contact(9, 6, 10), contact(15, 6, 10)],
        solenoid_fields=slider(16.0, 12.0, 6.0),
        tube=tube(0, None, 0, 0, 0, 0, step=METAL_STEP, solenoid=(23, 19),
                  circuits=[tube_circuit(1, dc(12.0), resistance=0.0, linked=[3]),
                            tube_circuit(0, dc(0.0), resistance=5.0),
                            tube_circuit(2, dc(0.0), resistance=1.5)],
                  potentials=[potential_goal(0, "mean", -0.3, 0.3, 10.0, 150.0, reference=2)],
                  metals=[metal((4, 13), (20, 17), [0, 1], seed=11, **METAL),
                          metal((10, 1), (14, 11), [2, 3], seed=12, **METAL)]),
        reference=[supply(23, 19, 12.0)])
'''

REALITY = '''
REALITY["wheatstone_bridge"] = R(
    "A Wheatstone bridge reading a magnetoresistor: a field sensor",
    kept("A bridge balances on a ratio",
         "It reads zero when its arms' resistances stand in the same ratio, so a small change "
         "shows as a voltage from zero: how strain gauges and field sensors are read."),
    kept("Shorted Hall field, longer paths",
         "The short wide strip's contacts short out its Hall field, its carriers' paths bend, "
         "and it resists more in either field direction: InSb field plates in sensor bridges."),
    simpler("Two strips and two ideal resistors",
            "A full sensor bridge has four matched elements whose drifts cancel; here two arms "
            "are ideal lumped resistors, so the long strip's own small magnetoresistance shifts "
            "the balance a little."),
    METAL_CLASSICAL, METAL_Z, TUBE_CIRCUITS,
    electrodes=["long strip, left end (to ground through 5; the voltmeter's +)",
                "long strip, right end (the supply, 12)",
                "wide strip, left contact (to ground through 1.5; the voltmeter's −)",
                "wide strip, right contact (the supply, 12)"],
    metals=["long strip (barely magnetoresistive)", "short wide strip (magnetoresistive)"])
'''


def sub(s, old, new):
    assert s.count(old) == 1, old
    return s.replace(old, new)


def main():
    p = "scripts/levels.py"
    s = open(p, encoding="utf-8").read()
    assert "def wheatstone_bridge" not in s
    a = s.index("def magnetoresistance():")
    b = s.index("\n\n\ndef ", a)
    s = s[:b] + LEVEL.rstrip("\n") + "\n" + s[b:]
    s = sub(s, '''        ("magnetoresistance", magnetoresistance),
''', '''        ("magnetoresistance", magnetoresistance),
        ("wheatstone_bridge", wheatstone_bridge),
''')
    open(p, "w", encoding="utf-8", newline="\n").write(s)
    p = "scripts/reality.py"
    s = open(p, encoding="utf-8").read()
    assert 'REALITY["wheatstone_bridge"]' not in s
    a = s.index('REALITY["magnetoresistance"] = R(')
    b = s.index("\n\nREALITY[", a)
    s = s[:b] + "\n" + REALITY.rstrip("\n") + s[b:]
    open(p, "w", encoding="utf-8", newline="\n").write(s)
    print("wheatstone bridge added")


if __name__ == "__main__":
    main()
