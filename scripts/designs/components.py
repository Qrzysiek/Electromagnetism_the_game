"""Inside the components: vacuum tubes, metals, transistors; AC and induction."""

import math
from designs.common import (
    METAL,
    METAL_STEP,
    contact,
    current_goal,
    cut,
    dc,
    heat_goal,
    level,
    metal,
    plate,
    player_plate,
    potential_goal,
    projectile,
    sine,
    slider,
    supply,
    tube,
    tube_circuit,
)
from designs.matter import pulse

# Inside the components: vacuum tubes (docs/TUBES.md, PHYSICS.md 2.11). z-invariant:
# the electrodes are long prisms, the electrons long lines of charge, the currents per
# unit length; Newtonian (tube voltages give v ~ 0.03 c in SI; here c is infinite).

def vacuum_diode():
    # The planar diode of crates/level/tests/tubes.rs: plates 10 long, 1 thick, faces 7
    # apart, inside the tube's grounded envelope (PHYSICS.md 2.11). Preview currents
    # (`generator tube-scan`): 1.00 at 22, 1.15 at 24, 1.29 at 26, 1.45 at 28; the
    # verification's are ~7 % lower (I ~ V^1.5) with an error of ~10 % (the envelope made
    # the convergence first order, PHYSICS.md 2.11). The goal 0.9 to 1.5 is verified for
    # about 23.6 to 28 (1.0 to 1.4 left only 25.3 to 26.9); the reference 26: 1.203 +-
    # 0.124.
    return level(
        "Vacuum diode",
        "Inside a vacuum tube. The left plate is a hot cathode: its right face emits "
        "electrons (blue dots), as many as its field pulls off. The right plate, the "
        "anode, is on your power supply. Set its voltage so that a current between 0.9 "
        "and 1.5 reaches the anode (averaged over $t = 10$ to $20$). The electrons in "
        "flight repel the ones behind them: this space charge chokes the current, which "
        "grows as $I \\propto V^{3/2}$ (Child and Langmuir), faster than in proportion "
        "to the voltage. The potential map shows the space charge as a flat region in "
        "front of the cathode, where the field is pulled down to zero.",
        grid=(24, 16), c=None, t_max=200.0,
        electrodes=[plate(8, 8, 10, thickness=1.0, angle_deg=90.0),
                    plate(16, 8, 10, thickness=1.0, angle_deg=90.0, tunable=True)],
        supplies=slider(70.0, 13.0, 26.0, 52.0),
        tube=tube(0, 1, 0.9, 1.5, 10.0, 20.0, emit_toward_deg=0.0),
        reference=[supply(16, 8, 26.0)])



def sideways():
    # Designed with `generator tube-scan` (PHYSICS.md 2.11, V8): a diode (30 across a gap of
    # 9: Hull's planar cutoff sqrt(2V)/d = 0.86) in the player's solenoid field, and a
    # collector at 10 above the top end of the gap; inside the grounded envelope. Preview
    # collector current: 0.47 without field (the anode 0.96), nothing with B negative (the
    # drift goes down); rising as the electrons curve upwards, 0.83 at B = 0.8, 1.28 at
    # 0.95, 1.34 at 1.1, 1.27 at 1.2, 1.15 at 1.3, while the anode's falls to 0 from 1.1
    # (the cutoff; below it, 0.13 at 0.95, comes in at the gap's ends, where the
    # translation symmetry behind Hull's theorem is broken); above ~1.35 the electrons
    # curl back onto the cathode, which takes them (0.87 at 1.4, 0.50 at 1.5). The goal
    # >= 1.1 asks for B in ~0.95-1.3: it solved on 0.925-1.225, 3.8 % of the solenoid's +-4
    # (the showcase rule's review, 2026-10-20). At least 0.9 (still nearly twice the 0.47
    # without field) asks for ~0.82-1.38, between Hull's cutoff and the electrons' curling
    # back onto the cathode.
    return level(
        "Sideways",
        "Around this tube is a solenoid: switch on its supply (the yellow coil at the top "
        "right) to fill the tube with a magnetic field along its length, out of the board "
        "when positive. Electrons crossing it curve, and above a critical field (Hull's "
        "cutoff) they cannot reach the anode at all: they drift along the cathode instead, "
        "sideways, perpendicular to both the electric and the magnetic field. Send the "
        "electrons into the collector above the gap: at least 0.9 must land there "
        "(averaged over $t = 20$ to $40$).",
        grid=(24, 22), c=None, t_max=60.0,
        electrodes=[plate(6, 11, 10, thickness=1.0, angle_deg=90.0),
                    plate(16, 11, 10, thickness=1.0, angle_deg=90.0, kind="potential", value=30.0),
                    plate(9, 19, 6, thickness=1.0, kind="potential", value=10.0)],
        solenoid_fields=slider(4.0, -1.1, 0.6, 1.1),
        tube=tube(0, 2, 0.9, 10.0, 20.0, 40.0, emit_toward_deg=0.0, solenoid=(23, 21)),
        reference=[supply(23, 21, 1.1)])



def bipolar_diode():
    # Designed with `generator tube-scan` (PHYSICS.md 2.11): Vacuum diode's tube with the
    # anode's face emitting ions under the space-charge limit. The ions' q/m = 1/16 (16
    # times heavier per charge than the electrons, not 1836 or more: in a static field the
    # mass changes neither the paths nor the electrons' gain, only the time scale and the
    # ion current, sqrt(m_e/m_i) of the electrons' in a planar diode; heavier ions would
    # take longer to settle and longer runs). Preview over t = 20-30 (both currents settled:
    # at q/m = 1/25 the electrons' was 2.540 there and 2.553 over 30-50): electrons into the
    # anode 0.613 at 10 V, 1.723 at 20, 2.558 at 26, 3.496 at 32 (V^1.5; Vacuum diode's
    # 1.29 at 26 without ions: the ions' charge doubles it, 1.86 in Langmuir's infinite
    # planar diode), ions into the cathode 0.148, 0.410, 0.619, 0.830 (0.24 of the
    # electrons'; the planar law: 0.25). The goals, electrons >= 2.0 and ions <= 0.75, ask
    # for about 22 to 29.5 V in the preview; without the ions the electron goal would need
    # 41 V, far beyond the wear limit.
    return level(
        "Bipolar diode",
        "The same diode, but now the anode's face is coated with a salt that, hot, emits "
        "positive ions (orange), as many as its field pulls off, just as the cathode emits "
        "electrons. The ions fly the other way, and their positive charge cancels part of "
        "the electrons' cloud: more electrons get through than in the vacuum diode at the "
        "same voltage (1.86 times as many in an infinitely wide diode: Langmuir, 1929). But "
        "the ions strike the cathode and wear it away. Set the anode's voltage so that at "
        "least 2.0 of electrons reach the anode while at most 0.75 of ions strike the "
        "cathode (both averaged over $t = 20$ to $30$). Both currents grow as "
        "$I \\propto V^{3/2}$.",
        grid=(24, 16), c=None, t_max=200.0,
        electrodes=[plate(8, 8, 10, thickness=1.0, angle_deg=90.0),
                    plate(16, 8, 10, thickness=1.0, angle_deg=90.0, tunable=True)],
        supplies=slider(40.0, 13.0, 26.0, 33.0),
        tube=tube(0, 1, 2.0, 10.0, 20.0, 30.0, emit_toward_deg=0.0,
                  ions=(1, 0.0625, 180.0),
                  goals=[current_goal(0, 0.0, 0.75, 20.0, 30.0, carriers="ions")]),
        reference=[supply(16, 8, 26.0)])



def through_the_cloud():
    # Designed with `generator tube-scan` and verified with `check` (PHYSICS.md 2.11),
    # inside the tube's grounded envelope. The owner (2026-10-12) found the flow at 30 V
    # too weak for the ion to cross a cloud: the anode is now at 60 (a current of 2.40
    # without plates, 2.5 times the former, its space charge twice), and the ion (q/m =
    # 0.02, launched at t = 10 when the flow is established) at 1.8: slower, at 60-70 V it
    # ends on the cathode whatever single plate helps. Without a plate it leaves the top at
    # x = 4.8, left of the detector (8 to 14). Of 90 single plates, 24 bring it in:
    # positive ones below the gap, left of its path, kick it rightwards (the original
    # idea) but draw electrons (2.7-3.7 at +30 to +75), and one above the cathode's top
    # end (6, 18) at +65/+75 pushes it right as it leaves the gap (2.85). Verified with
    # the anode at 60: (9, 2) vertical at +30, current 2.898 +- 0.036, the ion arrives at
    # t = 21.23 and 21.20 (the reference); (8, 2) horizontal at +45, 2.956 +- 0.018, t =
    # 21.33 and 21.32; (6, 18) vertical at +65, 2.830 +- 0.073, t = 23.80 and 23.85. Too
    # strong a plate below draws too much current (3.2 and more at +55 to +75).
    # The showcase rule's review (2026-10-20): with the reference plate the anode's current
    # is 1.75, 2.26, 2.58, 2.90, 3.23, 3.59, 4.14, 4.88, 5.56 at 40, 50, 55, 60, 65, 70, 80,
    # 90, 100 V, the ion leaving the top at x = 9.66, 8.42, 7.78, 7.15, 6.53, 5.91, 4.59, 3.41
    # (the detector 8 to 14) and at 100 turned back onto the cathode: the supply's +-100
    # shows it all. With the current goal 2.5-3.1 the supply solved on 54-61 (3.5 %); at
    # 2.3-3.3 the ion bounds it above (~62) and the current below (~51).
    # Measured: the supply 51-61 (5.0 %), the plate 30-49, next to +-100's half's middle (50):
    # the plates span +-120 (its middle 60; the interval 7.9 %).
    return level(
        "Through the cloud",
        "An ion (a long line of positive charge, the large dot) is fired up through the gap "
        "of a working diode, after its electrons have started to flow. Keep the diode "
        "running, with a current between 2.3 and 3.3 into the anode (averaged over "
        "$t = 15$ to $30$), and bring the ion into the detector at the top. The anode's "
        "field pushes positive charges toward the cathode, and the electrons' cloud pulls "
        "on the ion as it passes: at the voltage the current needs, the ion leaves the "
        "tube far left of the detector. Your two plates can push it back on course, but "
        "they act on the electrons as well, and the ion's own charge disturbs the flow it "
        "crosses.",
        grid=(24, 20), c=None, t_max=40.0,
        electrodes=[plate(6, 10, 10, thickness=1.0, angle_deg=90.0),
                    plate(16, 10, 10, thickness=1.0, angle_deg=90.0, tunable=True)],
        supplies=slider(100.0, 40.0, 60.0, 80.0),
        max_plates=2, plate_voltages=slider(120.0, -40.0, 20.0, 30.0, 45.0, 65.0),
        plate_size=(3.0, 0.5, 4.0),
        tube=tube(0, 1, 2.3, 3.3, 15.0, 30.0, emit_toward_deg=0.0,
                  projectiles=[projectile(11, 0, 90.0, 1.8, 0.2, 10.0, launch=10.0,
                                          detector=((8, 18), (14, 20)))]),
        reference=[supply(16, 10, 60.0), player_plate(9, 2, 30.0, angle_deg=90.0)])



def electron_gun():
    # Designed with `generator tube-scan`, verified with `check` (PHYSICS.md 2.11). The
    # cathode is temperature limited (0.2 per unit length of its 8-long face: at most 1.6),
    # so pulling harder gains nothing. Inside the tube's grounded envelope the drift
    # region's potential peaks on the axis, which focuses the beam a little: with the
    # former 0.1 it reached the collector unaided (0.113 +- 0.042 for the goal 0.08), so
    # the cathode was made hotter, its beam's space charge spreading it harder. Without
    # plates now 0.155 +- 0.026 (the goal 0.24: missed). Two horizontal plates beside the
    # beam after the slit focus it (preview 0.25-0.31 at +2 to +8, x = 17 to 23).
    # Verified: (21, 13/7) at +6 0.294 +- 0.008 (the reference); (17, 13/7) at +6 0.288
    # +- 0.023; (19, 13/7) at +4 0.294 +- 0.019. At 0 the lens is too strong (0.192 +-
    # 0.025) and at +12 too weak (0.161 +- 0.012): both missed.
    # The plates' range (the showcase rule's review, 2026-10-20): the collector's current
    # (preview) with both plates at V: 0.180, 0.253, 0.281, 0.293, 0.273, 0.196, 0.155 at
    # V = 0, 2, 4, 6, 8, 10, 14; nothing at -10 to -40 (the beam turned back), 0.125, 0.101,
    # 0.098 at +20, +30, +40 (the plates, above the anode's 30, collect it). One plate alone
    # steers the beam off: each slider solves on ~3-7.5 with the other at 6, 2.8 % of +-80.
    # Every behaviour lies within +-40 (beyond it only more of the same), so the plates span
    # +-40 (that interval 5.6 %).
    return level(
        "Electron gun",
        "An electron gun: the hot cathode (left) gives at most 0.2 per unit length of its "
        "face (its temperature limits it: a stronger field pulls no more), the anode's slit "
        "lets a thin beam through, and the small collector on the right should catch it. "
        "But the beam's electrons repel one another: it spreads and misses. Place two plates "
        "beside the beam to focus it, so that at least 0.24 reaches the collector (averaged "
        "over $t = 20$ to $40$). A plate at a lower potential than the beam's surroundings "
        "pushes the electrons back toward the axis; too low, and it throws them back.",
        grid=(32, 20), c=None, t_max=60.0,
        electrodes=[plate(4, 10, 8, thickness=1.0, angle_deg=90.0),
                    plate(12, 14, 6, thickness=1.0, angle_deg=90.0, kind="potential", value=30.0),
                    plate(12, 6, 6, thickness=1.0, angle_deg=90.0, kind="potential", value=30.0),
                    plate(28, 10, 2, thickness=1.0, angle_deg=90.0, kind="potential", value=30.0)],
        max_plates=2, plate_voltages=slider(40.0, -10.0, 3.0, 6.0, 12.0, 20.0),
        plate_size=(3.0, 0.5, 4.0),
        tube=tube(0, 3, 0.24, 5.0, 20.0, 40.0, emit_toward_deg=0.0, saturation=0.2),
        reference=[player_plate(21, 13, 6.0), player_plate(21, 7, 6.0)])



def ion_diode():
    # Designed with `generator tube-scan` (PHYSICS.md 2.11): an inverted coaxial diode, the
    # cathode the ring (mid-diameter 10, wall 1: inner radius b = 4.5), the anode the rod
    # (radius a = 1) emitting ions (q/m = 1/16, as in Bipolar diode) on every face, in the
    # player's solenoid field. Hull's cutoff from r0 to r1 for a charge released at rest:
    # B = sqrt(2 m V/q) 2 r1/|r0^2 - r1^2|: for the electrons (b to a) 0.80 at 30 V, for the
    # ions (a to b) 4 * 2b/2a = 18 times that: the ions cross, the electrons are insulated
    # in a sheath along the ring, whose charge draws more ions (the magnetically insulated
    # ion diode). The sheath keeps spreading towards the rod (the impedance collapse of real
    # ion diodes), so the goal's window is a fixed time after switch-on. A first design (with
    # the game's usual 150 macroparticles per emitting face) was shelved: the sheath's ion
    # current depended mostly on the number of macroparticles (x1 to x8: 8.21, 7.42, 6.90,
    # 6.66 at 28 V, B = 1.8) and the verdicts' bars were 15-25 %. With 8 times as many (the
    # tree code makes them affordable: PHYSICS.md 2.11) it converges in their number (x16:
    # 8.30 against 8.19 at 30 V, B = 1.7), and the rest of the preview's error was the
    # emitters' release layer, which the sheath magnifies: emitters meshed twice as finely
    # bring the preview within 7.3 % of the verification at B = 1.7 (4.8 % at 2.4): bars of
    # ~10 %, a preview of ~3.4 s and a verification of ~13 s. (V, B) grid in the preview
    # (ions onto the ring / electrons leaking into the rod, t = 20-30):
    #   V \ B  1.5         1.7        1.9        2.2        2.6        3.0
    #   24     5.55/0      4.82/0     4.39/0     3.97/0     3.64/0     3.46/0
    #   27     7.09/0      6.12/0     5.52/0     4.92/0     4.47/0     4.21/0
    #   30     8.87/0.006  7.78/0     6.80/0     6.01/0     5.40/0     5.05/0
    #   33    11.84/0.094  9.88/0.002 8.31/0     7.20/0     6.42/0     5.95/0
    #   36    14.81/0.171 11.32/0    10.06/0     8.55/0     7.53/0     6.95/0
    # The boost is largest just above the cutoff (1.76 times the strong field's at 30 V).
    # The showcase rule's review (2026-10-20): with ions 5.5 to 8 the supply solved on 27-32
    # (3.1 % of +-80; the +-80 shows the leak above ~50). Widened to 5.0-8.5 (not 9: the ion
    # current climbs as B falls towards the electrons' cutoff, 8.87 at B = 1.5, and the
    # field's interval would then take in 1.5, its half's middle). Measured: the supply
    # 26-33.5 (4.7 %), the field 1.55-3.0, up to +-3's end (5.05 at B = 3): the field spans
    # +-8 (its half's middle 4; the ions bent harder towards their own cutoff, ~14).
    # The goals (ions 5.5 to 8, leak <= 0.05) leave a diagonal strip, B roughly in proportion
    # to sqrt(V): at the reference (30, 1.9) V ~ 27-32 and B ~ 1.66-2.55; full voltage and
    # field (~8.1) overshoot it, the sliders' middles fall short. The leak the description
    # promises is within both sliders (the owner, 2026-10-13: whatever a level says can
    # happen must be possible): at 30 V the electrons pour into the rod with the field off
    # (17.4 at B = 0, 0.76 at 1.2, 0.006 at 1.5), and at B = 1.9 above ~50 V (0.066 at 50,
    # 0.25 at 60, 3.6 at 80): hence the supply's +-80. A grounded envelope (a thin ring
    # around the cathode) closes the tube: without it the open 2D boundary left a field
    # outside the cathode, which drew electrons off its outer face into an arc drifting
    # round it (E x B) off the board.
    return level(
        "Ion diode",
        "An ion diode turned inside out: the ring is the cathode, emitting electrons from "
        "its inner face, and the rod in its middle, on your supply, emits ions. Without a "
        "magnetic field the electrons rush to the rod and waste the supply's current. Your "
        "solenoid (the yellow coil, top right) curls them: above Hull's critical field they "
        "cannot reach the rod, and stay in a sheath along the ring. The heavy ions still "
        "cross, and the sheath's negative charge pulls more of them across. Deliver an ion "
        "current between 5 and 8.5 to the ring while at most 0.05 of electrons leak to the "
        "rod (averaged over $t = 20$ to $30$). The critical field grows with the voltage, "
        "as $B \\propto \\sqrt{V}$: too weak a field for your voltage lets the electrons "
        "leak; the sheath helps the ions most just above it. The thin outer ring is the "
        "tube's grounded envelope.",
        grid=(24, 16), c=None, t_max=100.0,
        electrodes=[plate(12, 8, 10, thickness=1.0, shape="ring"),
                    plate(12, 8, 2, thickness=2.0, shape="rod", tunable=True),
                    plate(12, 8, 14.5, thickness=0.5, shape="ring")],
        supplies=slider(80.0, 27.0, 30.0, 33.0),
        solenoid_fields=slider(8.0, 1.7, 1.9, 2.2),
        tube=tube(0, 1, 0.0, 0.05, 20.0, 30.0, solenoid=(23, 15), step=0.05,
                  ions=(1, 0.0625), particles=8.0, emitter_mesh=2.0,
                  goals=[current_goal(0, 5.0, 8.5, 20.0, 30.0, carriers="ions")]),
        reference=[supply(12, 8, 30.0), supply(23, 15, 1.9)])



def rectifier():
    # Designed with `generator tube-scan`, verified with `check` (PHYSICS.md 2.11): a
    # cathode between two anodes 6 away, each on an AC supply (period 20) whose setting is
    # its amplitude, the right one's half a period behind; the cathode is the output,
    # through R = 30 and C = 1 to ground, and charges positive as it loses electrons.
    # Output over t = 60-100 (preview): both anodes at 20/30/40: mean 6.0/9.9/13.9, swing
    # 0.34/0.55/0.77 (full wave); one at 20/30/40: 3.7/6.2/9.0, swing 1.7/2.6/3.5 (half
    # wave); both, but of opposite signs (in phase): 4.5/7.4/10.5, swing 3.3/5.0/6.7. The
    # goal (mean 8.5 to 11.5, ripple below 1.2, over t = 40-80, after the output has
    # settled) needs both anodes in opposite phase: each supply alone, the other at the
    # reference, solves on 24 to 36 (preview; `generator window`), the verification's
    # error (+-0.7 on the mean) narrowing it (at 30/30: 9.26 +- 0.66, at the edge). The
    # reference 32 sits in the middle.
    omega = 2.0 * math.pi / 20.0
    return level(
        "Rectifier",
        "A rectifier tube: a hot cathode in the middle, an anode on each side. Each anode is "
        "on its own AC power supply, a sine of period 20 whose amplitude you set; the right "
        "one's runs half a period behind the left one's. Electrons flow from the cathode to "
        "an anode only while that anode is positive: the tube conducts one way. The cathode "
        "is connected to ground through a load resistor (30) and a capacitor (1): losing "
        "electrons, it charges positive, and its potential is the rectifier's output, DC "
        "from AC. Make a smooth supply: the cathode's mean potential between 8.5 and 11.5, "
        "and its ripple (half its peak-to-peak) below 1.2, over $t = 40$ to $80$.",
        grid=(24, 16), c=None, t_max=200.0,
        electrodes=[plate(12, 8, 10, thickness=1.0, angle_deg=90.0),
                    plate(5, 8, 10, thickness=1.0, angle_deg=90.0, tunable=True),
                    plate(19, 8, 10, thickness=1.0, angle_deg=90.0, tunable=True)],
        supplies=slider(80.0, 16.0, 32.0, 48.0),
        tube=tube(0, None, 0.0, 0.0, 0.0, 0.0,
                  circuits=[tube_circuit(1, sine(0.0, omega), resistance=0.0),
                            tube_circuit(2, sine(0.0, omega, phase=math.pi), resistance=0.0),
                            tube_circuit(0, dc(0.0), resistance=30.0, capacitance=1.0)],
                  potentials=[potential_goal(0, "mean", 8.5, 11.5, 40.0, 80.0),
                              potential_goal(0, "swing", 0.0, 1.2, 40.0, 80.0)]),
        reference=[supply(5, 8, 32.0), supply(19, 8, 32.0)])



def triode():
    # Designed with `generator tube-scan`, verified with `check` (PHYSICS.md 2.11). A
    # cathode 6 long, five grid wires (0.6 x 0.5, 2 apart) 1.25 before its face on one
    # sine source (amplitude 1, period 40) whose offset, the bias, the player sets, and an
    # anode 6 long through R = 60 to B+ = 100. The anode's own capacitance (it faces the
    # grounded envelope and the grid) shunts the load: g_m/(omega C) capped the gain at ~1.3
    # with a 14-long anode at period 20, so the anode is short and the signal slow. Output
    # amplitude (verification, over one period t = 40-80, the particles softened at an
    # eighth of a segment): 1.96 at bias -1.2, 2.37 at -1, 2.76 at -0.8, 3.12 at -0.6,
    # 3.49 at -0.4, 3.79 at -0.2, 4.09 at 0, the preview ~7 % above (bars +-0.3 to
    # +-0.36): the goal 2.2 to 3.8 is verified for a bias of about -0.9 to -0.4 (at 0 and
    # -0.2 the gain is too high, at -1 too low); cut off below about -2 (only the grid's
    # capacitive feedthrough, 0.32, is left).
    omega = 2.0 * math.pi / 40.0
    wires = [plate(7, y, 0.6, thickness=0.5, angle_deg=90.0, tunable=(y == 8))
             for y in (8, 4, 6, 10, 12)]
    return level(
        "Triode",
        "A triode: between the hot cathode (left) and the anode (right) a grid of five thin "
        "wires, all on one supply: a small signal, a sine of amplitude 1 and period 40, "
        "around a steady potential, the bias, which you set. Close to the cathode, where the "
        "electrons are slow, the grid controls the current far more strongly than the "
        "anode does. The anode takes its current through a resistor (60) from a fixed 100, "
        "so as the current swings, the anode's potential swings too, and much more than "
        "the grid's: an amplifier. Bias the grid so that the anode's potential oscillates "
        "with an amplitude between 2.2 and 3.8 (over $t = 40$ to $80$). Too negative a "
        "bias cuts the current off.",
        grid=(24, 16), c=None, t_max=200.0,
        electrodes=[plate(5, 8, 6, thickness=1.0, angle_deg=90.0),
                    plate(12, 8, 6, thickness=1.0, angle_deg=90.0)] + wires,
        supplies=slider(2.6, -1.4, -0.7, 0.7),
        tube=tube(0, None, 0.0, 0.0, 0.0, 0.0, emit_toward_deg=0.0,
                  circuits=[tube_circuit(2, sine(1.0, omega), resistance=0.0,
                                         supply_offset=True, linked=(3, 4, 5, 6)),
                            tube_circuit(1, dc(100.0), resistance=60.0)],
                  potentials=[potential_goal(1, "amplitude", 2.2, 3.8, 40.0, 80.0,
                                             omega=omega)]),
        reference=[supply(7, 8, -0.7)])



def klystron():
    # Designed with `generator tube-scan` (PHYSICS.md 2.11). A weak beam (a temperature-
    # limited cathode, 0.02 per unit length of its 4-long face, 400 macroparticles per unit
    # time: 1/400 of a period apart) from a Wehnelt gun (the cathode between two plates at
    # its potential, the accelerating gap into the grounded drift tube's mouth a
    # converging lens: most of the beam reaches the collector) flies past a buncher
    # (plates 14-20 on a sine of amplitude V1, period 1: its entrance and exit gaps
    # modulate the electrons' speeds) and a catcher 21 further (one plate, a tank: L = 0.051
    # in series with R = 0.04 to ground, resonant at omega = 2 pi with the plate's own
    # capacitance, 0.497 from the boundary elements; Q ~ 8, rung up in 2.5). The player
    # sets the cathode's potential -V0 (the beam's energy) and V1. The beam's space charge
    # is negligible: the bunching is ballistic. Output amplitude against the drive at
    # V0 = 50 (preview): 0.073 at V1 = 8, 0.111 at 14, 0.121 at 20, 0.109 at 26
    # (overbunching); with the buncher at a fixed offset the J1 shape was measured to 2 %
    # (J1(X)/J1max with X = 0.092 V1 fitted at the peak: 0.306, 0.711, 0.956 at V1 = 4,
    # 10, 24 for the measured 0.300, 0.694, 0.959). The output's phase moves with the
    # drift angle omega L/v0, 7.5 degrees per volt near 50 (theta0/(2 V0) = 7.6: phase
    # pushing): 122.7 at V0 = 40, 169.6 at 45, -151.4 at 50, -114.4 at 55, -82.7 at 60
    # (V1 = 20); the amplitude hardly depends on V0. At V0 = 54: the amplitude 0.072 at
    # V1 = 8, 0.088 at 10, 0.095 at 11, 0.102 at 12, 0.111 at 14 (phase -125); past the
    # optimum 0.103 at 28, 0.096 at 30, 0.089 at 32, 0.079 at 34 (phase -113 to -103).
    # Refine 1 to 2: phases within 2.6 degrees, amplitudes 5 %. A flat optimum and those
    # bars leave no narrow window at saturation, so the goal is a klystron backed off below
    # it (as one feeding an accelerator is, so that its output follows its drive):
    # amplitude 0.085-0.105 and phase -135 to -105, which ask for V0 ~ 52-56 and V1 ~
    # 9.7-12.5, or overbunched ~27.5-33; the next island of phase lies beyond the slider.
    omega = 2.0 * math.pi
    P = plate
    return level(
        "Klystron",
        "A klystron amplifier. The gun on the left shoots a weak electron beam into a "
        "grounded tube; you set the cathode's voltage, and with it the beam's speed. On its "
        "way the beam passes the buncher, plates on a small sine (period 1) whose amplitude "
        "you set: crossing its gaps speeds some electrons up and slows others down. Further "
        "on the fast ones catch up with the slow ones and the beam gathers into bunches, one "
        "per period, which induce a current in the catcher plate, a tuned circuit ringing at "
        "the same frequency. Run it as a linear amplifier, below saturation: the catcher "
        "must oscillate with an amplitude between 0.085 and 0.105, and with a phase, "
        "relative to the buncher's signal, between $-135°$ and $-105°$ (over $t = 10$ to "
        "$20$). A stronger signal bunches the beam tighter, up to a point: too strong a one "
        "overbunches it, the fast electrons overtaking the slow ones before the catcher. "
        "The phase depends on how long the bunches take to arrive.",
        grid=(48, 20), c=None, t_max=40.0,
        electrodes=[
            P(3, 10, 4, thickness=1.0, angle_deg=90.0, tunable=True),
            P(3, 13, 4, thickness=1.0), P(3, 7, 4, thickness=1.0),
            P(10, 13, 6, thickness=1.0), P(10, 7, 6, thickness=1.0),
            P(17, 13, 6, thickness=1.0, tunable=True), P(17, 7, 6, thickness=1.0),
            P(28, 13, 14, thickness=1.0), P(28, 7, 14, thickness=1.0),
            P(38, 13, 4, thickness=1.0),
            P(38, 7, 4, thickness=1.0),
            P(43, 13, 4, thickness=1.0), P(43, 7, 4, thickness=1.0),
            P(47, 10, 4, thickness=1.0, angle_deg=90.0)],
        supplies=slider(80.0, -30.0, -54.0, -60.0, 11.0, 30.0),
        tube=tube(0, None, 0.0, 0.0, 0.0, 0.0, emit_toward_deg=0.0, saturation=0.02,
                  particle_rate=400.0,
                  circuits=[tube_circuit(0, dc(0.0), resistance=0.0, linked=(1, 2)),
                            tube_circuit(5, sine(0.0, omega), resistance=0.0, linked=(6,)),
                            tube_circuit(10, dc(0.0), resistance=0.04, inductance=0.051)],
                  potentials=[potential_goal(10, "amplitude", 0.085, 0.105, 10.0, 20.0,
                                             omega=omega),
                              potential_goal(10, "phase", -135.0, -105.0, 10.0, 20.0,
                                             omega=omega)]),
        reference=[supply(3, 10, -54.0), supply(17, 13, 11.0)])




def magnetron():
    # Designed with `generator tube-scan --trace` (PHYSICS.md 2.11; docs/TUBES.md). A strapped
    # eight-vane magnetron in the z-invariant model: the cathode a rod (radius a = 2) on the
    # player's supply, eight vanes (tips at b = 4.5; centres on grid nodes, so the diagonal
    # ones are shorter) and the anode block, a ring (inner radius 8.5). The vanes at 0, 90,
    # 180 and 270 degrees and the block are grounded; the diagonal ones are strapped onto one
    # node and a tank to ground (L = 0.036 in series with R = 0.002, and 1.34 to ground: with
    # the vanes' own 1.337 a cold resonance at 3.22, Q ~ 55; the added capacitance halves the
    # electrons' pull on the frequency, which without it hopped between branches near 3.0 and
    # 3.6). Its oscillation is the pi mode: four crests around the vanes, turning at omega/4.
    # The oscillation starts from the switch-on transient (the cathode's step rings the tank
    # to ~4) and the particles' noise, so its start-up time depends on the resolution (t ~ 20
    # to 80); its saturated swing does not (45 V, B = 6: preview 14.6, verification 15.6, x4
    # particles 16.0): the goal's window is late, after every resolution has started. An
    # injected priming signal did not make the start reliable (below threshold the hub
    # absorbs it) and is not used. Swing over t = 90-120 (preview; frequency from the
    # up-crossings 3.03-3.11 where it oscillates):
    #   V \ B   5.5    6      6.5
    #   37.5    0.57   0.55   0.32
    #   40      8.50   0.68   0.36
    #   42.5    3.13   9.03   0.41
    #   45      0.80  14.86   1.01
    #   47.5    5.16* 14.55  15.99
    #   50      4.20* 3.47   21.51
    #   52.5    1.53   1.61  22.84
    #   55      0.82   1.50   5.47   (* the other branch, near 3.6)
    # The oscillating band follows the Buneman-Hartree condition, the electrons' drift in
    # step with the wave, V = Omega B (b'^2 - a^2)/2 - Omega^2 b'^2/2 with Omega = omega/4 =
    # 0.77 and an effective anode radius b' ~ 5 (between the vanes' tips and their roots):
    # 37, 41 and 45 at B = 5.5, 6 and 6.5, the band's lower edges. Far above it the electrons
    # outrun the wave and the oscillation stops (B = 6: none at 50, below Hull's cutoff 58.7);
    # below Hull's critical field they cross to the vanes (at 45 V: 34 into the vanes at
    # B = 4.4; the cutoff is ~4.5, Hull's with b' = 5). Just above it a second solution:
    # at 45 V, B = 4.58-4.67 the tank swings 13.6 (omega 3.29) while 8.4 flows into the
    # vanes, seven times the main band's current (B = 4.8-5.25: the other branch, 5.8-7.8,
    # short of the goal). `generator window` at the reference: the supply 43.3-48.0 (2.9 %
    # of the range), the solenoid 5.83-6.25 and that pocket (and their mirrors: the field's
    # sign only turns the hub the other way); every middle verified.
    # Its output (the owner, 2026-10-13: the radiation on the map, from all the charges): the
    # tube moved to the middle of a 96 x 96 board (about a wavelength, 2 pi c/omega ~ 82 at
    # c = 40, across) and an antenna, a rod of diameter 2 twenty cells from the axis, is
    # strapped to the tank's vanes; the field solver (level::tube_em) shows its waves. Its
    # radiation resistance, measured with `generator tube-field`, is the tank's load: at R =
    # 0.002 the board radiated 0.91 of the load's R<I^2>, at R = 0.0018 the two agree to ~3 %
    # (test magnetron_antenna_radiates_its_load). Redesigned windows: the supply 42.0-47.3,
    # the solenoid 5.83-6.33 (the near-cutoff pocket shrank to a point at 4.58).
    return level(
        "Magnetron",
        "A magnetron, the heart of a microwave oven and of wartime radar. Electrons leave "
        "the cathode in the middle (your supply sets its voltage below the vanes) and, in "
        "your solenoid's field (the yellow coil, top right), curl around it instead of "
        "crossing to the vanes. The vanes are strapped in two sets: one grounded with the "
        "block, the other through a tuned circuit ringing near $\\omega = 3$. Its "
        "oscillation is a wave of four crests running around the vanes at $\\omega/4$; "
        "electrons drifting in step with it give it energy and gather into four spokes "
        "that turn with it. Make the tuned circuit oscillate with a swing of at least 10 "
        "(half its peak-to-peak, over $t = 90$ to $120$). The electrons drift at $E/B$: "
        "to keep pace with the wave the voltage must rise with the field. Too weak a field "
        "lets them straight onto the vanes; too fast a drift outruns the wave. The rod "
        "outside the block is the antenna, strapped to the tuned set of vanes: it radiates "
        "the microwaves: the magnetic or energy flow map shows them spreading over the "
        "board once the oscillation has grown out of the noise, after about $t = 50$ (the "
        "mouse wheel zooms in on the spokes, a drag pans).",
        grid=(96, 96), c=None, t_max=150.0,
        electrodes=[plate(48, 48, 4, thickness=4.0, shape="rod", tunable=True),
                    plate(54, 48, 3.0, thickness=0.6, angle_deg=0.0),
                    plate(52, 52, 2.313708, thickness=0.6, angle_deg=45.0),
                    plate(48, 54, 3.0, thickness=0.6, angle_deg=90.0),
                    plate(44, 52, 2.313708, thickness=0.6, angle_deg=135.0),
                    plate(42, 48, 3.0, thickness=0.6, angle_deg=180.0),
                    plate(44, 44, 2.313708, thickness=0.6, angle_deg=225.0),
                    plate(48, 42, 3.0, thickness=0.6, angle_deg=270.0),
                    plate(52, 44, 2.313708, thickness=0.6, angle_deg=315.0),
                    plate(48, 48, 18, thickness=1.0, shape="ring"),
                    plate(68, 48, 2, thickness=2.0, shape="rod")],
        supplies=slider(80.0, -38.0, -45.0, -52.0),
        solenoid_fields=slider(10.0, 5.5, 6.0, 6.5),
        tube=tube(0, None, 0.0, 0.0, 0.0, 0.0, solenoid=(95, 95), step=0.05,
                  circuits=[tube_circuit(1, dc(0.0), resistance=0.0, linked=(3, 5, 7, 9)),
                            tube_circuit(2, dc(0.0), resistance=0.0018, inductance=0.036,
                                         capacitance=1.34, linked=(4, 6, 8, 10))],
                  potentials=[potential_goal(2, "swing", 10.0, 1000.0, 90.0, 120.0)],
                  output=(1, 40.0)),
        reference=[supply(48, 48, -45.0), supply(95, 95, 6.0)])


def ohms_law():
    # Designed with `generator tube-scan` (PHYSICS.md 2.11, the Drude metal): a strip
    # 16 x 4 between two contacts. Preview currents into the right contact over t = 5-35
    # (the contacts emitting as thermal walls): 0.91 at 8, 1.72 at 16, 2.43 at 24, 3.17 at
    # 32, 3.84 at 40, -2.61 at -24: in proportion to within a slow fall of I/V (0.1135 to
    # 0.096: the carriers warm in the field); the bulk theory gives 0.1175 V. The goal
    # 2.6-3.2: the supply solves on 26-32 (`generator window`).
    return level(
        "Ohm's law",
        "Inside a metal. The rectangle is a conductor seen end-on: its electrons (blue; a "
        "few traced with their paths) fly at their thermal speed and bounce off the "
        "metal's ions (the dots), a zigzag. Your supply, on the contact at the right end, "
        "adds a slow drift along the field to that random motion: a current, in "
        "proportion to the voltage, $I = V/R$ (Ohm's law). Make a current between 2.6 and "
        "3.2 flow into the right contact (averaged over $t = 10$ to $40$). At every "
        "collision the electrons hand the ions the energy the field gave them: the Joule "
        "heat, shown below. A classical model of a metal (its real electrons are a "
        "quantum gas): see the physics notes.",
        grid=(24, 16), c=None, t_max=200.0,
        electrodes=[contact(3, 8, 4), contact(21, 8, 4, tunable=True)],
        supplies=slider(40.0, 15.0, 30.0, 35.0),
        tube=tube(0, None, 0, 0, 0, 0, step=METAL_STEP,
                  goals=[current_goal(1, 2.6, 3.2, 10.0, 40.0, carriers="net")],
                  metals=[metal((4, 6), (20, 10), [0, 1], seed=1, **METAL)]),
        reference=[supply(21, 8, 30.0)])



def long_and_thin():
    # Designed with `generator tube-scan`: two strips of the same metal in series, the upper
    # one 16 x 4 (L/W = 4), the lower one 8 x 8 (L/W = 1). Supply 1 on the upper strip's
    # left end, supply 2 on the junction (the upper strip's right end, strapped to the lower
    # strip's left end), the lower strip's right end grounded. Each strip alone (the
    # contacts emitting as thermal walls): the upper 2.63 at 26.5, the lower 2.32 at 5.5:
    # R = 10.1 and 2.37, a ratio of 4.3 (L/W gives 4; the thin strip's stronger field warms
    # its carriers). The goals 2.0-2.8 each: at the reference (32, 5.5) supply 1 solves on
    # 28-38, the junction on about 4.7-6.6 (the textbook ratio, V2 = V1/5 = 6.4, passes).
    return level(
        "Long and thin",
        "Two strips of the same metal: the upper one long and thin, the lower one short and "
        "wide. A strip's resistance grows with its length and falls with its width, "
        "$R = L/(\\sigma W)$ (per unit length along it): the upper one has about four times "
        "the lower one's. They are in series: your first supply holds the upper strip's "
        "left end, your second the junction (the upper strip's right end, wired to the lower "
        "strip's left end), and the lower strip's right end is grounded. Make a current "
        "between 2.0 and 2.8 flow through each strip (into its left end, averaged over "
        "$t = 10$ to $40$).",
        grid=(24, 16), c=None, t_max=200.0,
        electrodes=[contact(17, 5, 8), contact(3, 12, 4, tunable=True),
                    contact(21, 12, 4, tunable=True), contact(7, 5, 8)],
        supplies=slider(40.0, 5.5, 15.0, 32.0),
        tube=tube(0, None, 0, 0, 0, 0, step=METAL_STEP,
                  circuits=[tube_circuit(2, dc(0.0), resistance=0.0, linked=[3])],
                  goals=[current_goal(1, 2.0, 2.8, 10.0, 40.0, carriers="net"),
                         current_goal(3, 2.0, 2.8, 10.0, 40.0, carriers="net")],
                  metals=[metal((4, 10), (20, 14), [1, 2], seed=4, **METAL),
                          metal((8, 1), (16, 9), [3, 0], seed=5, **METAL)]),
        reference=[supply(3, 12, 32.0), supply(21, 12, 5.5)])



def potentiometer():
    # Designed with `generator tube-scan`: a strip 18 x 4, its left end grounded, a side
    # contact on its top a third of the way along (x 8-10). Preview with the right end at
    # 30: the side contact's current (positive in) -0.41 at 8, -0.16 at 9, 0.17 at 10,
    # 0.45 at 11, 0.65 at 12: none at about 9.5 (a uniform field gives 10; the ohmic
    # contacts' layers take a little). The goals: 2.4-3.0 along the strip (the right end
    # at about 27-33) and -0.35..0.35 through the side contact (about 8.3-10.6 at 30;
    # +-0.25 left the verification uncertain: the side contact's current carries the
    # strip's thermal noise, an error bar of ~0.18 with 50 macroparticles per area over
    # t = 10-40; now 100 over 10-80). With the contacts emitting as a thermal wall (they
    # had been heat sinks, PHYSICS.md 2.11) the side contact's current at 30 is -0.22 at 8,
    # -0.02 at 9, 0.31 at 10: none at about 9.05 (V1/3 = 10: the contacts take ~1 V).
    # The side contact's current verified with a bar of ~0.3 (+-0.35 left it uncertain):
    # the window +-0.5, about 7.5-10.6 at 30. The showcase rule's review (2026-10-20): it solved
    # on 6.75-10.25, 4.4 % of +-40: the window +-0.7 (the side current changes by ~0.27 per
    # volt, still a null to find).
    return level(
        "Potentiometer",
        "A longer strip of metal with a third contact on its side, a third of the way "
        "along, on your second supply. Inside a conductor carrying a steady current the "
        "field is uniform: the potential falls evenly from one end to the other (the "
        "contacts take a little). Make a "
        "current between 2.4 and 3.0 flow along the strip into its right end, and set the "
        "side contact so that no current flows through it (between $-0.7$ and $0.7$, "
        "averaged over $t = 10$ to $80$): it must sit at the potential the strip has "
        "there. A potentiometer measures a voltage this way, drawing no current from it.",
        grid=(24, 16), c=None, t_max=200.0,
        electrodes=[contact(2, 8, 4), contact(22, 8, 4, tunable=True),
                    contact(9, 11, 2, vertical=False, tunable=True)],
        supplies=slider(40.0, 10.0, 15.0, 30.0),
        tube=tube(0, None, 0, 0, 0, 0, step=METAL_STEP,
                  goals=[current_goal(1, 2.4, 3.0, 10.0, 80.0, carriers="net"),
                         current_goal(2, -0.7, 0.7, 10.0, 80.0, carriers="net")],
                  metals=[metal((3, 6), (21, 10), [0, 1, 2], seed=2,
                                **dict(METAL, carriers_per_area=100.0))]),
        reference=[supply(22, 8, 30.0), supply(9, 11, 9.0)])



def hall_probe():
    # Designed with `generator tube-scan` (100 macroparticles per area: the floating
    # contacts' voltage carries the resistor's thermal noise, ~0.1 over the window): a
    # strip 16 x 4, the right end on the supply, two floating contacts on its sides at
    # mid-length (C = 1 to ground). At 30 (I ~3.1) the top's potential against the
    # bottom's: -0.58 at B = 0.4, -1.15 at 0.7, -1.44 at 1.0, -2.1 at 1.2, -2.33 at 1.6,
    # +0.75 at -0.7 (electrons pushed to the top; about -0.15 at B = 0). Below the weak-field
    # R_H I B = (4/pi) I B/rho as omega_c tau grows (0.4 B: the Hall factor falls towards 1)
    # and the contacts short the strip's ends. The goals: 2.7-3.5 into the right end (about
    # 26-34) and the top 1.0-1.9 below the bottom (B about 0.6-1.15 at 30). With the
    # contacts emitting as a thermal wall: -0.90 at B = 0.9, -1.44 at 1.2, -1.68 at 1.5,
    # +0.22 at 0. `generator window` at the reference (30, 1.3): the supply solves on
    # 26-34, the solenoid on 1.0-1.4 (its slider +-4, the middles at +-2), over t = 10-80
    # (the floating contacts' noise).
    # The solenoid's full range (the owner, 2026-10-20: every behaviour the model can
    # show): +-16, from straight paths between collisions to circles 0.12 across
    # (omega_c tau = 8; the step turns the carriers by at most 0.56 rad).
    # Over it at 30 (target/scratch/b_range_scan.py; preview): current 3.14, 3.12, 3.06,
    # 3.00, 2.91, 2.73, 2.55, 2.21, 1.96 and the top against the bottom -0.45, -1.05, -2.43,
    # -3.68, -4.67, -6.61, -8.12, -10.57, -12.52 at B = 0.5, 1, 2, 3, 4, 6, 8, 12, 16 (the
    # Hall voltage ever more sublinear). The goal moved to -6.0 to -3.5: B about 2.9-5.5.
    hall = dict(METAL, carriers_per_area=100.0)
    return level(
        "Hall probe",
        "A strip of metal carrying a current, in your solenoid's field (B along the "
        "strip's length, out of the board when positive; the yellow coil, top right). The "
        "field pushes the moving electrons sideways (the Lorentz force) until the charge "
        "they pile on one side holds the rest back: a voltage across the strip, the Hall "
        "voltage, $V_H \\approx I B/(n q)$ (per unit length). The two contacts on the "
        "strip's sides read it, floating: no current flows through them. Try the field's "
        "whole range: weak, the electrons' paths barely bend; strong, they curl into small "
        "circles between collisions. Make a current between 2.7 and 3.5 flow into the right "
        "end and the top contact sit between 3.5 and 6 below the bottom one (averaged over "
        "$t = 10$ to $80$). Which way must the field point? Electrons are negative.",
        grid=(24, 16), c=None, t_max=200.0,
        electrodes=[contact(3, 8, 4), contact(21, 8, 4, tunable=True),
                    contact(12, 11, 2, vertical=False), contact(12, 5, 2, vertical=False)],
        supplies=slider(40.0, 15.0, 30.0, 35.0),
        solenoid_fields=slider(16.0, -4.2, 3.0, 4.2, 5.0),
        tube=tube(0, None, 0, 0, 0, 0, step=METAL_STEP, solenoid=(23, 15),
                  circuits=[tube_circuit(2, dc(0.0), capacitance=1.0),
                            tube_circuit(3, dc(0.0), capacitance=1.0)],
                  goals=[current_goal(1, 2.7, 3.5, 10.0, 80.0, carriers="net")],
                  potentials=[potential_goal(2, "mean", -6.0, -3.5, 10.0, 80.0,
                                             reference=3)],
                  metals=[metal((4, 6), (20, 10), [0, 1, 2, 3], seed=3, **hall)]),
        reference=[supply(21, 8, 30.0), supply(23, 15, 4.2)])



def positive_carriers():
    # Designed with `generator tube-scan`: two strips 16 x 4 in the player's solenoid, the
    # upper of the levels' metal (electrons), the lower the same but its carriers positive
    # (q/m = +1), each with its left end grounded, its right end on a supply and floating
    # side contacts at mid-length (C = 1). The Hall voltage's sign follows the carriers' for
    # a given current: with the carriers moving the same way the two strips' Hall voltages
    # agree in sign, so the tops both sit above the bottoms only with opposite supplies.
    # Preview at (30, B = -1.3, -15) over t = 10-80: upper current 3.10, top above bottom
    # 1.52; lower current 1.66, top above bottom 1.21 (both 0.3 lower at B = 0 for the upper
    # strip, 0 for the lower; the lower's potential shifts the upper's reading by ~0.2). The
    # floating contacts' noise over the window: 0.11 (upper) and 0.15 (lower) between seeds
    # and resolutions (50 macroparticles per area: two strips at 100 would pass the
    # verification's carrier limit). At B = -1.3 the lower Hall voltage at -12 (0.65-0.84)
    # did not verify against its noise; at B = -1.8: upper 2.23 +- 0.10 at 30, lower 0.89
    # at -10 and 1.23 at -14 (two seeds, both resolutions). `generator window`: the upper
    # supply solves on 28-34, the lower on -16 to -8, the solenoid on -2.5 to -1.5: its
    # slider is +-7 (at +-5 the half's middle, -2.5, solved).
    # The solenoid's full range (the owner, 2026-10-20: every behaviour the model can
    # show): +-16, from straight paths between collisions to circles 0.12 across
    # (omega_c tau = 8; the step turns the carriers by at most 0.56 rad).
    # Over it at (30, B, -12) (preview): upper current 3.11, 3.06, 3.06, 2.98, 2.85, 2.68,
    # 2.50, 2.22, 1.95; upper Hall 0.40, 1.12, 2.31, 3.57, 4.69, 6.12, 7.63, 10.47, 11.99;
    # lower 0.57, 0.73, 1.14, 1.47, 2.05, 2.80, 3.15, 4.44, 4.92 at B = -0.5, -1, -2, -3,
    # -4, -6, -8, -12, -16. The goals moved: current 2.6-3.5, upper 3.4-6.0, lower 1.3-2.6:
    # B about -5.6 to -2.9.
    material = dict(METAL, carriers_per_area=50.0)
    return level(
        "Positive carriers",
        "Two strips of different metals in your solenoid's field, each with two floating "
        "contacts on its sides that read its Hall voltage. The upper one conducts by "
        "electrons; the lower one's carriers are positive (orange). Some real metals, such "
        "as zinc and beryllium, behave so: band theory explains it, and in this classical "
        "model their carriers simply carry a positive charge. The field reaches strengths "
        "that curl the carriers into small circles. Make a current between 2.6 and 3.5 flow "
        "into the upper strip's right end, and make both strips' top contacts sit above "
        "their bottom ones: the upper by 3.4 to 6, the lower by 1.3 to 2.6 (averaged over "
        "$t = 10$ to $80$). Which way must each current flow?",
        grid=(24, 20), c=None, t_max=200.0,
        electrodes=[contact(3, 14, 4), contact(21, 14, 4, tunable=True),
                    contact(12, 17, 2, vertical=False), contact(12, 11, 2, vertical=False),
                    contact(3, 5, 4), contact(21, 5, 4, tunable=True),
                    contact(12, 8, 2, vertical=False), contact(12, 2, 2, vertical=False)],
        supplies=slider(40.0, -12.0, 15.0, 30.0),
        solenoid_fields=slider(16.0, -4.0, 3.0, 4.0),
        tube=tube(0, None, 0, 0, 0, 0, step=METAL_STEP, solenoid=(23, 19),
                  circuits=[tube_circuit(2, dc(0.0), capacitance=1.0),
                            tube_circuit(3, dc(0.0), capacitance=1.0),
                            tube_circuit(6, dc(0.0), capacitance=1.0),
                            tube_circuit(7, dc(0.0), capacitance=1.0),
                            tube_circuit(4, dc(0.0), resistance=0.0)],
                  goals=[current_goal(1, 2.6, 3.5, 10.0, 80.0, carriers="net")],
                  potentials=[potential_goal(2, "mean", 3.4, 6.0, 10.0, 80.0, reference=3),
                              potential_goal(6, "mean", 1.3, 2.6, 10.0, 80.0, reference=7)],
                  metals=[metal((4, 12), (20, 16), [0, 1, 2, 3], seed=6, **material),
                          metal((4, 3), (20, 7), [4, 5, 6, 7], seed=7, q_over_m=1.0,
                                **material)]),
        reference=[supply(21, 14, 30.0), supply(23, 19, -4.0), supply(21, 5, -12.0)])



# Why Hall probes are semiconductors (the owner, 2026-10-20): two strips of the same size and
# the same conductance in one field, the upper the levels' metal (few carriers that fly far: it
# resembles n-type germanium at 1e16 per cm^3 in these fields, game scales.rs), the lower four
# times as many carriers colliding four times as often (density 8, mean free path 0.25: sigma =
# rho tau the same; n-type silicon at 1e17 per cm^3, within x1.21). At the same current the
# lower's carriers drift four times slower, the magnetic push on each is four times weaker, and
# so is its Hall voltage (V_H = IB/(nq)); at the same field they turn a quarter as far between
# collisions (omega_c tau). Copper is the trend a million times further. The step follows the
# dense strip's plasma oscillation (omega_p = 10.0: 0.019, omega_p dt = 0.19).
DENSE = dict(density=8.0, thermal_speed=2.0, mean_free_path=0.25, core_radius=0.002,
             carriers_per_area=50.0)

DENSE_STEP = 0.019



def few_or_many():
    # Designed with `generator tube-scan` (target/scratch/few_scratch.py): at (30, B, 30)
    # (preview, t = 10-80) over B = 0, 1, 2, 3, 4.2, 6, 9, 12, 16: upper current 3.18, 3.16,
    # 3.14, 3.03, 2.93, 2.77, 2.53, 2.29, 2.04; lower 3.25, 3.40, 3.33, 3.29, 3.25, 3.10, 3.13,
    # 2.94, 2.68; upper Hall (top less bottom) -0.19, -1.58, -3.13, -3.91, -5.22, -6.73, -9.27,
    # -10.7, -13.3; lower +0.76, +0.02, -0.81, -0.85, -1.45, -1.80, -2.85, -4.00, -4.97 (its
    # floating contacts scatter ~0.4: its macroparticles are four times heavier). Per unit IB
    # the Hall voltages are 0.42 and 0.11: four times apart, as the densities. Verified at
    # (30, 3.8, 30): currents 3.00 +- 0.05 and 3.32 +- 0.16, Hall -4.69 +- 0.08 and -1.37 +-
    # 0.19. The Hall readings' error bars grow with B (upper 0.08, 0.51, 0.74 and lower 0.29,
    # 0.36, 0.70 at B = 2.5, 5, 5.5): with goals -6 to -3.5 and -2 to 0 the solenoid solved on
    # 3.0-4.5 only (4.7 %), so -7 to -3 and -2.5 to 0. `generator window`: the solenoid
    # solves on 2.0-6.0 (12.5 % of +-16; 1.5 and 6.5-8.5 fail), the upper supply on 28-35
    # (8.8 % of +-40), the lower on 25-31 (7.5 %; both measured with the narrower Hall goals,
    # their edges set by the current goals). The solenoid's +-16 spans straight paths to
    # circles (omega_c tau 6.4 in the upper strip, 1.6 in the lower).
    material = dict(METAL, carriers_per_area=50.0)
    return level(
        "Few carriers or many",
        "Two strips of the same size and the same resistance in your solenoid's field, each "
        "with two floating contacts on its sides that read its Hall voltage. The upper one's "
        "carriers are few and fly far between collisions (like germanium's); the lower one has "
        "four times as many, colliding four times as often (like silicon's). Make a current "
        "between 2.7 and 3.5 flow into each strip's right end, and the upper strip's top "
        "contact sit between 3 and 7 below its bottom one, while the lower strip's sits at "
        "most 2.5 below (averaged over $t = 10$ to $80$). Carrying the same current, the many "
        "carriers drift slowly, the field pushes each one gently, and the Hall voltage is "
        "small: real Hall probes are semiconductors with few, fast carriers (InSb, GaAs), "
        "while copper's carriers are about a million times as many, and its Hall voltage a "
        "million times smaller: microvolts.",
        grid=(24, 20), c=None, t_max=200.0,
        electrodes=[contact(3, 14, 4), contact(21, 14, 4, tunable=True),
                    contact(12, 17, 2, vertical=False), contact(12, 11, 2, vertical=False),
                    contact(3, 5, 4), contact(21, 5, 4, tunable=True),
                    contact(12, 8, 2, vertical=False), contact(12, 2, 2, vertical=False)],
        supplies=slider(40.0, 15.0, 30.0, 35.0),
        solenoid_fields=slider(16.0, -3.8, 3.0, 3.8, 5.0),
        tube=tube(0, None, 0, 0, 0, 0, step=DENSE_STEP, solenoid=(23, 19),
                  circuits=[tube_circuit(2, dc(0.0), capacitance=1.0),
                            tube_circuit(3, dc(0.0), capacitance=1.0),
                            tube_circuit(6, dc(0.0), capacitance=1.0),
                            tube_circuit(7, dc(0.0), capacitance=1.0),
                            tube_circuit(4, dc(0.0), resistance=0.0)],
                  goals=[current_goal(1, 2.7, 3.5, 10.0, 80.0, carriers="net"),
                         current_goal(5, 2.7, 3.5, 10.0, 80.0, carriers="net")],
                  potentials=[potential_goal(2, "mean", -7.0, -3.0, 10.0, 80.0, reference=3),
                              potential_goal(6, "mean", -2.5, 0.0, 10.0, 80.0, reference=7)],
                  metals=[metal((4, 12), (20, 16), [0, 1, 2, 3], seed=6, **material),
                          metal((4, 3), (20, 7), [4, 5, 6, 7], seed=7, **DENSE)]),
        reference=[supply(21, 14, 30.0), supply(23, 19, 3.8), supply(21, 5, 30.0)])



def magnetoresistance():
    # Designed with `generator tube-scan` (target/scratch/design_mr_intro.py): a strip 16 x 4
    # (current along it) and one 4 x 10 (short along the current, wide across it), both on
    # one supply, in the player's solenoid. In the long strip the Hall field holds the
    # carriers back from the sides and they flow straight; the short wide strip's contacts
    # short that field: the carriers' paths bend and its resistance grows with |B|
    # (geometric magnetoresistance; with L/W = 0.4 short of the Corbino limit 1 + (0.4 B)^2).
    # Preview over t = 10-80, the currents into the two right contacts: at 7 V the long
    # strip 0.73, 0.69, 0.68, 0.63, 0.63 at B = 3..7 (its own, intrinsic magnetoresistance:
    # tau ~ 1/v), the short one 5.43, 4.62, 4.08, 3.64, 3.04 (about 7.2 at B = 0); at 6 V
    # and 8 V in proportion. A run takes ~16 s, at B = 6 as at 0. `generator window`: B solves on
    # 4.5-5.5 and -5.5 to -4.5 (its mirror's middle, -5, not verified: near the goals' edges);
    # with the long strip's goal 0.62-0.76 the supply solved on 7.0-7.5 only (1.2 % of +-20):
    # widened to 0.58-0.78, the slider +-10.
    # The solenoid's full range (the owner, 2026-10-20: every behaviour the model can
    # show): +-16, from straight paths between collisions to circles 0.12 across
    # (omega_c tau = 8; the step turns the carriers by at most 0.56 rad).
    # Over it at 7 V (preview): long strip 0.795, 0.772, 0.738, 0.683, 0.679, 0.626, 0.573,
    # 0.547, 0.522, 0.452; short 6.96, 6.09, 5.41, 4.63, 4.06, 3.61, 2.68, 2.09, 1.56, 1.18 at
    # B = 1, 2, 3, 4, 5, 6, 8, 10, 13, 16. At B = 5 both currents go as the supply (long
    # 0.096 V, short 0.58 V), so its window is the overlap of the two goals' ratios: with
    # 0.58-0.78 and 3.2-4.6 it solved on 6.5-7.5 only (5 %; 6.0 gave 0.550, 8.0 gave 4.68).
    # Widened to 0.56-0.80 and 3.2-4.8: the supply about 5.8-8.2, B about 3.8-6.6 (the half's
    # middle, 8, clear).
    return level(
        "Magnetoresistance",
        "Two strips of the same metal on one supply, in your solenoid's field: a long one and "
        "a short, wide one. In the long strip the charge the field piles on its sides (the "
        "Hall voltage) pushes the carriers back, and they flow straight on. The short strip's "
        "contacts short that charge out: its carriers' paths bend and lengthen, and it "
        "resists more, whichever way the field points; at full strength its carriers circle "
        "on the spot between collisions. Make between 0.56 and 0.8 flow into the long "
        "strip's right end and between 3.2 and 4.8 into the short strip's (averaged over "
        "$t = 10$ to $80$).",
        grid=(24, 20), c=None, t_max=200.0,
        # 0, 1: the long strip's contacts; 2, 3: the short wide strip's.
        electrodes=[contact(3, 15, 4), contact(21, 15, 4, tunable=True),
                    contact(9, 6, 10), contact(15, 6, 10)],
        supplies=slider(10.0, 7.0, 4.0),
        solenoid_fields=slider(16.0, -5.0, 3.0, 5.0),
        tube=tube(0, None, 0, 0, 0, 0, step=METAL_STEP, solenoid=(23, 19),
                  circuits=[tube_circuit(1, dc(0.0), resistance=0.0, linked=[3])],
                  goals=[current_goal(1, 0.56, 0.80, 10.0, 80.0, carriers="net"),
                         current_goal(3, 3.2, 4.8, 10.0, 80.0, carriers="net")],
                  metals=[metal((4, 13), (20, 17), [0, 1], seed=11, **METAL),
                          metal((10, 1), (14, 11), [2, 3], seed=12, **METAL)]),
        reference=[supply(21, 15, 7.0), supply(23, 19, 5.0)])


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
    # within +-0.3: B solves on 11-13 (6.3 %; 10.5 too, 10.75 uncertain at the edge) and
    # -13.5 to -10.5 (9.4 %; `generator window`, every 0.5; the bridge reads |B|).
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




def strength_and_sign():
    # The metals' master (task 6): a magnetic-field sensor of two parts in series on one
    # supply, joined at a floating junction. The short wide strip of Magnetoresistance (4 x
    # 10, electrons) resists more as |B| grows, whichever way B points; the Hall bar 16 x 4
    # carries positive carriers (Positive carriers' lower strip), its Hall voltage odd in B.
    # The current pins the supply, the junction's potential (the wide strip's share of it)
    # |B|, and the Hall voltage's sign B's sign.
    # Designed with `generator tube-scan` (target/scratch/level_scratch.py): at 30 over B =
    # -16, -8, -5, -3, 0, 3, 5, 8, 12, 16 (preview, t = 10-80) the current 1.20, 2.02, 2.48,
    # 2.72, 2.98, 2.75, 2.43, 2.05, 1.58, 1.16, the junction 7.44, 5.10, 4.16, 3.37, 2.86, 3.34,
    # 4.06, 5.00, 5.94, 6.48 (both even in B: the two strips' magnetoresistance) and the Hall
    # voltage -10.9, -7.25, -5.43, -3.69, -0.11, 3.67, 5.36, 6.97, 9.34, 10.9 (odd: the top
    # above the bottom for B out of the board, the current then flowing leftwards through
    # positive carriers). Then mirrored left to right (the wide strip under the Hall bar's
    # right end, the wire between them short and clear of the board's legend): the current
    # flows rightwards, and the top sits above the bottom for B into the board; at B = 5,
    # before the mirror, 2.40, 3.98, 5.15. Verified at (30, -5): 2.47 +- 0.06, 4.065 +-
    # 0.005, 5.43 +- 0.19. `generator window`: the supply solves on 27-32 (6.2 % of +-40; 26
    # and 33 fail on the junction), the solenoid on -6.0 to -3.5 (7.8 % of +-16; -6.5 to -9
    # and -3 to -2 fail). The junction over the current is the wide strip's resistance, whatever the
    # supply: 0.96, 1.22, 1.67, 2.44, 3.75, 5.59 at |B| = 0, 3, 5, 8, 12, 16. Goals: the
    # current 2.1-2.9, the junction 3.6-4.5 (their ratio 1.24-2.14: |B| about 3-7, clear of
    # the half's middle, 8), the Hall voltage 3-8 (the sign).
    return level(
        "Strength and sign",
        "A magnetic field sensor in two parts, in series on your supply. The current enters "
        "the Hall bar's left end (its carriers are positive), leaves its right end into the "
        "junction, and crosses a short, wide strip of the electron metal to ground; that "
        "strip resists more as the field grows, whichever way the field points. Make a "
        "current between 2.1 and 2.9 flow through the wide strip, the junction sit between "
        "3.6 and 4.5, and the Hall bar's top contact sit 3 to 8 above its bottom one "
        "(averaged over $t = 10$ to $80$). The wide strip reads the field's strength but not "
        "its sign; the Hall bar reads its sign.",
        grid=(24, 20), c=None, t_max=200.0,
        # 0: the wide strip's bottom contact (grounded); 1: the Hall bar's left end (the
        # supply); 2, 3: its side contacts; 4, 5: the junction (the Hall bar's right end
        # wired to the wide strip's top contact, straight below it).
        electrodes=[contact(17, 2, 10, vertical=False), contact(3, 15, 4, tunable=True),
                    contact(12, 18, 2, vertical=False), contact(12, 12, 2, vertical=False),
                    contact(21, 15, 4), contact(17, 8, 10, vertical=False)],
        supplies=slider(40.0, 15.0, 30.0, 35.0),
        solenoid_fields=slider(16.0, -5.0, 3.0, 5.0),
        tube=tube(0, None, 0, 0, 0, 0, step=METAL_STEP, solenoid=(23, 19),
                  circuits=[tube_circuit(2, dc(0.0), capacitance=1.0),
                            tube_circuit(3, dc(0.0), capacitance=1.0),
                            tube_circuit(4, dc(0.0), capacitance=1.0, linked=[5])],
                  goals=[current_goal(5, 2.1, 2.9, 10.0, 80.0, carriers="net")],
                  potentials=[potential_goal(4, "mean", 3.6, 4.5, 10.0, 80.0),
                              potential_goal(2, "mean", 3.0, 8.0, 10.0, 80.0, reference=3)],
                  metals=[metal((4, 13), (20, 17), [1, 2, 3, 4], seed=13, q_over_m=1.0,
                                **METAL),
                          metal((12, 3), (22, 7), [5, 0], seed=14, **METAL)]),
        reference=[supply(3, 15, 30.0), supply(23, 19, -5.0)])



# Inside the components: transistors (docs/DRUDE.md, the owner approved the plan
# 2026-10-20). The field effect: a gate beside a thin channel of few carriers pushes them out
# or draws more in, and the channel's current follows. A metal cannot be gated (it screens
# within an atom); the model's lightly doped channel can (its Debye length 0.56, a quarter of
# its width). The channel's material: density 1, mean free path 0.5 (n-type silicon at 1e16
# per cm^3 by its own ratios, within x1.1: game scales.rs, `resemblance`).
CHANNEL = dict(density=1.0, thermal_speed=2.0, mean_free_path=0.5, core_radius=0.005,
               carriers_per_area=400.0)

CHANNEL_STEP = 0.04



def gate():
    # Designed with `generator tube-scan` (target/scratch/gate_design*.py): a channel 8 x 2
    # between a grounded source and a drain held at 10, gates 1.5 above and below it on the
    # player's supply. The channel charges through its own resistance (R C ~ L^2: a channel 16
    # long had not settled by t = 30; 8 long settles within ~10). Convergence: with the gates
    # 0.5 away and 100 macroparticles per area the preview's current was 0.342 at gate 0
    # against the verification's 0.220 (and 0.243 at refine 3); the macroparticles were the
    # cause (twice as many: 0.262; four times: 0.266; the verification with twice as many:
    # 0.256; the softening and the contacts' segments hardly mattered), and the gates moved
    # 1.5 away (three preview segments). Depletion needs more: the current at -10, -15, -25,
    # -35 was 0.176, 0.163, 0.150, 0.159 at 200 per area (a plateau, an artefact), 0.182, 0.136,
    # 0.071, 0.050 at 400, 0.181, 0.139, 0.069, 0.027 at 800, 0.171, 0.134, 0.065, 0.023 at the
    # verification with 800: below threshold the current crosses the gate's barrier
    # thermally, and coarse macroparticles' grainy field rocks the barrier (at -46: 0.078,
    # 0.052, 0.007 at 200, 400, 800). So 400 per area: converged within 10 % to about -25;
    # pinch-off beyond -40 is shown as it is, its current not converged at the preview. At
    # 200 per area over the slider: 0.957, 0.598, 0.262 at +100, 40, 0, the channel's carriers
    # 13575, 6917, 2969 (doubled at 400: at +60 the verification's ~36000 pass the particle
    # cap and pause emission, where nothing solves). Goal 0.08-0.16: the gate about -13 to
    # -24.
    return level(
        "Gate",
        "The principle of the transistor. A thin strip of lightly doped silicon carries a "
        "current from the drain, held at 10, to the grounded source. Two gate electrodes beside "
        "it, on your supply, never touch it: a negative gate pushes the strip's carriers out "
        "and the current fades, until none flows; a positive one draws more in from the "
        "contacts. Set the gate so that between 0.08 and 0.16 flows into the drain, a third to "
        "two thirds of what flows with the gate at 0 (averaged over $t = 20$ to $60$). A metal screens any field within an atom, so no gate reaches "
        "inside copper; billions of these gates, each over a thin oxide on a silicon channel, "
        "switch in every chip.",
        grid=(20, 16), c=None, t_max=100.0,
        # 0: the source (grounded); 1: the drain (held at 10); 2, 3: the gates (your supply).
        electrodes=[contact(5, 8, 2), contact(15, 8, 2),
                    plate(10, 11, 6, thickness=1.0, tunable=True),
                    plate(10, 5, 6, thickness=1.0)],
        supplies=slider(60.0, -18.0, 0.0, 20.0),
        tube=tube(0, None, 0, 0, 0, 0, step=CHANNEL_STEP,
                  circuits=[tube_circuit(1, dc(10.0), resistance=0.0),
                            tube_circuit(2, dc(0.0), resistance=0.0, linked=[3])],
                  goals=[current_goal(1, 0.08, 0.16, 20.0, 60.0, carriers="net")],
                  metals=[metal((6, 7), (14, 9), [0, 1], seed=3, **CHANNEL)]),
        reference=[supply(10, 11, -18.0)])



def inverter():
    # Designed with `generator tube-scan --trace` (target/scratch/fet_inverter*.py): Gate's
    # channel with its drain fed from 25 through a load of 30, the gate at -60
    # (off) until a pulse (rise 2) brings it to the player's supply at t = 100 for 100.
    # The drain's potential (preview, window means): at V_DD 25 off 24.08 (the channel leaks
    # 0.03 at -60; at -40 it leaked 0.09), on 14.35, 17.14, 20.33 with the gate at 0, -10,
    # -20: the low level about 0.3 per volt of gate. The output settles within ~40 of each
    # switch (the channel's and the load's R C); at the switch the gate's step jumps the drain
    # up by ~4 for a few units (their capacitance carries charge across: the gate-drain
    # feedthrough). With load 20 the swing was 4.9, with 40 7.1 (and slower). With V_DD on a
    # second slider (goals 22-26.5 and 16-19) it solved on 24-26 only (3.75 %: the low level
    # follows V_DD one for one), the gate on -14 to -6: so V_DD is the level's, and this an
    # introduction (the load is its new piece). Goals: the output 22.5-25.5 while off (t =
    # 60-100), 15-19 while on (t = 160-200): the gate's on level about -17 to -3.
    return level(
        "Inverter",
        "The simplest logic gate. The last level's channel, its drain now fed from 25 through "
        "a load resistor: while the channel is off, almost no current flows and the drain "
        "sits near 25 (a logical 1); while it conducts, the current through the load pulls "
        "the drain down (a 0). The gate starts at -60, off; at $t = 100$ a pulse brings it to "
        "your supply's level, for 100. Make the drain sit between 22.5 and 25.5 while the "
        "gate is off (averaged over $t = 60$ to $100$) and between 15 and 19 while it "
        "conducts ($t = 160$ to $200$). Watch the drain jump as the gate switches: the gate "
        "pushes charge onto it through the capacitance between them.",
        grid=(20, 16), c=None, t_max=260.0,
        # 0: the source (grounded); 1: the drain, fed from 25 through the load; 2, 3: the
        # gates (the pulse, its top on your supply).
        electrodes=[contact(5, 8, 2), contact(15, 8, 2),
                    plate(10, 11, 6, thickness=1.0, tunable=True),
                    plate(10, 5, 6, thickness=1.0)],
        supplies=slider(40.0, -10.0, 10.0),
        tube=tube(0, None, 0, 0, 0, 0, step=CHANNEL_STEP,
                  circuits=[tube_circuit(1, dc(25.0), resistance=30.0),
                            tube_circuit(2, pulse(-60.0, 2.0, delay=100.0, width=100.0,
                                                  fall=2.0, low=-60.0),
                                         resistance=0.0, linked=[3])],
                  potentials=[potential_goal(1, "mean", 22.5, 25.5, 60.0, 100.0),
                              potential_goal(1, "mean", 15.0, 19.0, 160.0, 200.0)],
                  metals=[metal((6, 7), (14, 9), [0, 1], seed=3, **CHANNEL)]),
        reference=[supply(10, 11, -10.0)])



def amplifier():
    # Designed with `generator tube-scan` (target/scratch/amp_*.py): Gate's channel, a sine
    # of amplitude 3 on both gates around the player's supply (the bias), the drain fed from
    # V_DD through R_L. The DC characteristics (fet_dc.py, preview, the drain through 20):
    # at V_D ~ 39 g_m ~ 0.018 between gates -20 and -10 and the output conductance ~0.0026
    # (intrinsic gain ~7); at V_D ~ 12 g_m ~ 0.009 against 0.012 (a resistor): gain needs
    # the drain high, so a load of 200 from V_DD 100 put it at ~41 at the bias -15. The drain
    # starts at V_DD and settles with tau ~ 40 (R_par C: C_out ~ 0.3, the drain contact's
    # 2D capacitance and the channel's coupling to the gates), so the gain-bandwidth
    # product g_m / C ~ 0.06: at the bias -15 the gain (window 2 periods after 100) was 1.14,
    # 1.41 at periods 100, 200, and ~1.9 at 400 (traced). Shorter channels were faster but
    # gained less (L 6, load 150: 1.11, 1.37; L 4, load 100: 0.67, 0.72), a load of 300
    # (V_DD 130) no more at period 300 (1.73). At period 300, measured over one period after
    # 120, the gain against the bias: -30 0.59, -25 0.55 (the preview's unconverged pinch-off
    # plateau, the drain stuck near 59), -20 1.62, -15 1.69, -10 1.43, -5 1.02, 0 0.78,
    # 8 0.42. The lag 60-70 degrees on top of the inversion (the output pole and the gate-drain
    # feedthrough, a right-half-plane zero). Preview against verification (refine 2): the
    # gain at -22 0.92 / 1.61 (the preview's pinch-off not converged: the drain stuck at 56
    # against 59.5), -15 1.69 / 1.69, -9 1.36 / 1.42. The verdict's error is the two runs'
    # difference, so with the goal 4.5-7 (gain 1.5) it solved on -20 to -13 only (`generator
    # window`, every 1: 6 %; -22 uncertain, the preview's plateau). V_DD 120 keeps the drain
    # in saturation to higher biases, where g_m is larger: the gain (preview) -22 0.85, -18
    # 1.88, -14 1.93, -10 1.75, -6 1.59, -2 1.21. It solves on -20 to -5 (every 1; 12.5 %
    # of +-60); the reference -13 verified 5.95 +- 0.91 (preview 5.35).
    period = 300.0
    omega = 2.0 * math.pi / period
    return level(
        "Amplifier",
        "A small signal made larger: the heart of every radio and audio amplifier. The gates "
        "now carry a slow sine of amplitude 3 (period 300) around your supply's level, the "
        "bias, and the drain is fed from 120 through a load of 200. As the gate rises the "
        "channel conducts more and the drain falls: the drain swings with the signal, "
        "inverted, and at a good bias by more. Bias it so that the drain's oscillation has an "
        "amplitude between 4.5 and 7 (a gain of 1.5 or more), measured over one period "
        "from $t = 150$, once the drain has settled from 120. Too negative, and the pinched "
        "channel's current barely follows the gate; too positive, and the drain falls so low "
        "that the channel is a plain resistor. The drain also lags: its capacitance charges "
        "through the load.",
        grid=(20, 16), c=None, t_max=460.0,
        # 0: the source (grounded); 1: the drain (the output), fed from 120 through the load;
        # 2, 3: the gates (the signal around your supply).
        electrodes=[contact(5, 8, 2), contact(15, 8, 2),
                    plate(10, 11, 6, thickness=1.0, tunable=True),
                    plate(10, 5, 6, thickness=1.0)],
        supplies=slider(60.0, -13.0, 0.0, 20.0),
        tube=tube(0, None, 0, 0, 0, 0, step=CHANNEL_STEP,
                  circuits=[tube_circuit(1, dc(120.0), resistance=200.0),
                            tube_circuit(2, sine(3.0, omega), resistance=0.0,
                                         supply_offset=True, linked=[3])],
                  potentials=[potential_goal(1, "amplitude", 4.5, 7.0, 150.0, 450.0,
                                             omega=omega)],
                  metals=[metal((6, 7), (14, 9), [0, 1], seed=3, **CHANNEL)]),
        reference=[supply(10, 11, -13.0)])



def ringing_tank():
    # Designed with `generator tube-scan` (target/scratch/rlc_proto.py): Magnetoresistance's
    # long strip (16 x 4, R ~ 8.8 at DC) as the resistor of a parallel tank: its left contact
    # fed through a coil L = 7.5 (no series resistance) by the player's step at t = 5, a
    # capacitor C = 1.5 from it to ground, its right contact grounded. Ideal: omega0 =
    # 1/sqrt(LC) = 0.298, Q = R sqrt(C/L) = 3.9. Traced (step 10): the end swings to 16.8,
    # 4.8, 13.5, 7.3 at t = 16.8, 29.6, 40.0, 52.5 about 10, a period of ~23 (the contact and
    # the strip add ~0.3 to C: 1/sqrt(7.5 x 1.8) gives 23.1) and Q ~ 6 (each period 0.59 of
    # the last: the strip resists more than its DC 8.8 at these frequencies). With a series
    # resistance 0.5 Q was ~2.4. The amplitude at the ringing's omega 2 pi / 23 over three
    # periods (t = 5-74): per volt 0.460, 0.469, 0.473 at 11, 14, 16.5 (the strip warms a
    # little at larger steps); preview against verification 4.998 / 5.056, 6.373 / 6.569,
    # 7.683 / 7.798. Goal 4.5-9 on a +-50 slider: the step solves on 10.5-18.5 (`generator
    # window`, every 0.5; 9.5-10 too, 10.25 uncertain at the goal's edge, an isolated 19.5)
    # and -18 to -14 at least (-10 not): 8-9 %.
    period = 23.0
    return level(
        "Ringing tank",
        "A tuned circuit, the heart of every radio: a coil (L = 7.5) and a capacitor "
        "(C = 1.5) meet at a metal strip's left end, and the strip's right end is grounded. At "
        "$t = 5$ your supply steps the coil's far end from 0 to its level. The coil's current "
        "overshoots, and the strip's end rings about the supply's level with a period of "
        "about 23 (the strip and its contact add a little to the capacitor). The strip is "
        "the circuit's resistor: its carriers' collisions turn each swing's energy into "
        "heat, and the ringing dies away. Make the ringing's amplitude between 4.5 and 9 "
        "(its oscillation at the period 23, over three periods from $t = 5$ to $74$).",
        grid=(24, 20), c=None, t_max=80.0,
        # 0: the strip's left end, the tank's node (the coil from your supply, the capacitor
        # to ground); 1: its right end (grounded).
        electrodes=[contact(3, 10, 4, tunable=True), contact(21, 10, 4)],
        supplies=slider(50.0, 14.0, 6.0),
        tube=tube(0, None, 0, 0, 0, 0, step=METAL_STEP,
                  circuits=[tube_circuit(0, pulse(0.0, 1.0, delay=5.0), resistance=0.0,
                                         inductance=7.5, capacitance=1.5)],
                  potentials=[potential_goal(0, "amplitude", 4.5, 9.0, 5.0, 74.0,
                                             omega=2.0 * math.pi / period)],
                  metals=[metal((4, 8), (20, 12), [0, 1], seed=11, **METAL)]),
        reference=[supply(3, 10, 14.0)])



def field_tuned_damping():
    # Designed with `generator tube-scan` (target/scratch/ftd_proto.py): Ringing tank's
    # parallel tank (no series resistance) with Magnetoresistance's short wide strip (4 along
    # the current, 10 across, R ~ 1 at B = 0 to ~6 at 16) as its resistor, in the player's
    # solenoid: Q = R sqrt(C/L) from ~0.5 to ~3 with L 6, C 1.5. The tank rings with a period
    # of ~20 (ideal 18.8: the strip and its contact add to C). The amplitude at 2 pi / 20 over
    # three periods (t = 5-65), step 6, 50 carriers per area (preview): 0.47, 0.75, 1.11,
    # 1.44, 1.46, 1.86, 2.06 at B = 0, 4, 8, 10, 12, 14, 16. Preview against verification, step
    # 8: 1.78 / 1.84, 2.08 / 2.20, 2.78 / 2.46 at 9.5, 11.5, 13.5: the error grows with Q (the
    # strip's thermal noise rings the tank at its own frequency, the more the higher its Q).
    # With 100 per area and the step 10: 2.19 / 2.21, 2.51 / 2.63, 3.15 / 2.95 (errors 0.04,
    # 0.18, 0.29). Goal 2.15-3.0: its lower bound keeps the window clear of the half's middle
    # (8 +- 0.8). B solves on 9.5-13 (every 0.5; 9 and 13.5 not: 10.9 %) and -13 to -11 (-9
    # not).
    period = 20.0
    return level(
        "Field-tuned damping",
        "Ringing tank's circuit, its resistor now Magnetoresistance's short wide strip, in "
        "your solenoid. With no field the strip barely resists, and across the tank it "
        "shorts the ringing out: the end hardly swings after the step to 10 at $t = 5$. In a "
        "field its carriers' paths bend and it resists more, and the tank rings longer: the "
        "field sets the tank's Q, whichever way it points. Make the ringing's amplitude "
        "between 2.15 and 3 (at the period 20, over three periods from $t = 5$ to $65$). The "
        "strip's own thermal noise rings the tank too, more the longer it rings.",
        grid=(24, 20), c=None, t_max=70.0,
        # 0: the strip's left contact, the tank's node (the coil from the step, the capacitor
        # to ground); 1: its right contact (grounded).
        electrodes=[contact(9, 10, 10), contact(15, 10, 10)],
        solenoid_fields=slider(16.0, 11.0, 4.0),
        tube=tube(0, None, 0, 0, 0, 0, step=METAL_STEP, solenoid=(23, 19),
                  circuits=[tube_circuit(0, pulse(10.0, 1.0, delay=5.0), resistance=0.0,
                                         inductance=6.0, capacitance=1.5)],
                  potentials=[potential_goal(0, "amplitude", 2.15, 3.0, 5.0, 65.0,
                                             omega=2.0 * math.pi / period)],
                  metals=[metal((10, 5), (14, 15), [0, 1], seed=12,
                                **dict(METAL, carriers_per_area=100.0))]),
        reference=[supply(23, 19, 11.0)])



def lock_in_hall():
    # Designed with `generator tube-scan` (target/scratch/lockin_proto.py): Hall probe's strip
    # (16 x 4 at the board's centre, 100 carriers per area, floating side contacts C = 1) with
    # a steady current from the player's supply on its right end, in the player's solenoid
    # swinging as B0 sin(EDDY_OMEGA t) (period 50.3). The top contact against the bottom one,
    # its amplitude and phase at the field's frequency over two periods from t = 10
    # (preview): (V, B0) = (30, 4) 4.89 at 160 deg, (-30, 4) 4.89 at -15 deg, (30, 8) 7.85 at
    # 160; at -30: 2.85, 3.80, 5.97, 6.45 at B0 = 2, 3, 5, 6 (sublinear above ~5); at B0 = 4:
    # 3.24, 6.26 at -20, -38. The phase lags the ideal 0 / 180 by ~17 deg for either sign of
    # the current (the contacts' C charging through the strip, ~2.5 time units; an inductive
    # pickup would have shifted the two signs opposite ways): the induced field circles the
    # board's centre and the contacts sit on its symmetry line. Preview against
    # verification: (-30, 3.5) 4.40 / 4.43, (-30, 2.7) 3.64 / 3.47, (-30, 5) 5.97 / 5.81,
    # (-24, 3.5) 3.78 / 3.54; the phases within 3 deg. The amplitude goal's ratio sets both
    # sliders' windows (the voltage goes as I B0): 4.9-7 puts the supply's between the
    # guessable 22 and 38 and the field's below 7.2 (sublinear near 7). The phase goal
    # -60..30 (in step with the field) takes the current's direction. Over two periods with
    # the goal 4.9-7 B0 solved on 4.5-6.0 only (4.7 %; the supply -36 to -26, 12.5 %): three
    # periods (the error bars ~0.8 times) and the goal from 4.7: the reference verified
    # 5.77 +- 0.26 at -17.5 +- 1.4 deg; B0 solves on 4.0-6.5 (every 0.5; 3.5 and 7 not:
    # 7.8 %), the supply on -36 to -24 (-37 and -22 not: 15 %); the opposite signs fail the
    # phase (supply +30 or B0 -5).
    period = 2.0 * math.pi / EDDY_OMEGA
    start, end = 10.0, 10.0 + 3.0 * period
    hall = dict(METAL, carriers_per_area=100.0)
    return level(
        "Lock-in Hall",
        "Hall probe's strip with a steady current (your supply on its right end), in your "
        "solenoid swinging as $B_0 \\sin(0.125\\,t)$, a period of about 50 (your second "
        "supply sets $B_0$). The Hall voltage across the strip swings with the field, in "
        "step with it or against it as the current's direction says. A lock-in amplifier "
        "picks out a signal's swing at one frequency, its size and its phase, from the noise "
        "and from anything else the strip does: here the top contact against the bottom one "
        "at the field's frequency, over three periods from $t = 10$. Make its amplitude "
        "between 4.7 and 7, in step with the field (its phase between -60 and 30 degrees: it "
        "lags a little, as the contacts charge through the strip). Which way must the "
        "current flow?",
        grid=(24, 16), c=None, t_max=end + 5.0,
        # 0: the strip's left end (grounded); 1: its right end (your supply); 2, 3: the top
        # and bottom side contacts (floating voltmeters).
        electrodes=[contact(3, 8, 4), contact(21, 8, 4, tunable=True),
                    contact(12, 11, 2, vertical=False), contact(12, 5, 2, vertical=False)],
        supplies=slider(40.0, -30.0, 30.0),
        solenoid_fields=slider(16.0, 5.0, 3.0),
        tube=tube(0, None, 0, 0, 0, 0, step=METAL_STEP, solenoid=(23, 15),
                  solenoid_drive=sine(1.0, EDDY_OMEGA),
                  circuits=[tube_circuit(2, dc(0.0), capacitance=1.0),
                            tube_circuit(3, dc(0.0), capacitance=1.0)],
                  potentials=[potential_goal(2, "amplitude", 4.7, 7.0, start, end,
                                             omega=EDDY_OMEGA, reference=3),
                              potential_goal(2, "phase", -60.0, 30.0, start, end,
                                             omega=EDDY_OMEGA, reference=3)],
                  metals=[metal((4, 6), (20, 10), [0, 1, 2, 3], seed=11, **hall)]),
        reference=[supply(21, 8, -30.0), supply(23, 15, 5.0)])



# The eddy-current levels' drive: the solenoid's field swings as B0 sin(wt), period ~50: long
# against the metal's settling (a few tenths) and its slow polarizations (a few units), and
# slow enough that a cut's two faces, a Debye length apart, carry little current across as
# a capacitor (at period 25 a core cut along its length took 61 % more heat than Ohm's law
# gives: the cut's 1/(wC) was comparable to a strip's resistance; at period 50, 14 %).
EDDY_OMEGA = 0.125



def cold_electrode(x, y):
    """A small floating electrode (zero charge) that emits nothing (saturation 0): the cathode
    a tube level with only isolated metals needs. Floating, it cannot draw charge from this
    open world's "infinity", which made a metal beside a grounded electrode run away
    (docs/OPEN_PROBLEMS.md, OP-1)."""
    return plate(x, y, 0.5, thickness=0.5, kind="charge", value=0.0)



def eddy_currents():
    # Designed with `generator tube-scan` (target/scratch/eddy_window.py): an isolated plate
    # 12 x 12 on the solenoid's axis, its field B0 sin(0.125 t). Local Ohm's law (PHYSICS.md
    # 2.11, D8): heat sigma <Bdot^2> int u dA = 0.47 x 728.75 x (0.125 B0)^2 / 2 = 2.68 B0^2;
    # the carriers' linear response holds up to Bdot ~ 0.1 (qEl/kT = 0.1 at the edges'
    # middles): B0 up to ~0.8. Measured (preview / verification): 0.744 / 0.753 at B0 = 0.55,
    # 1.153 / 1.181 at 0.68, 1.504 / 1.595 at 0.8 (5-7 % below local Ohm's law: the walls,
    # PHYSICS.md 2.11): about 2.55 B0^2. With the goal 1.0-1.4 the slider solved on 0.66-0.72
    # only (`generator window`, every 0.0225; the verdict's error straddles the goal's edges
    # beyond): widened to 0.9-1.5. Then the swing's full range (the owner, 2026-10-20): +-3;
    # heat over B0^2 2.31, 2.51, 2.47, 2.48, 2.39 at B0 = 0.5, 1, 1.5, 2, 3 (previews; the
    # carriers barely warm), so the goal 1.2-4.0 solves on B0 ~ 0.69-1.26.
    return level(
        "Eddy currents",
        "A metal plate alone in your solenoid, whose field now swings back and forth, "
        "$B = B_0 \\sin \\omega t$. A changing magnetic field induces an electric field "
        "that circulates about the solenoid's axis (Faraday's law), and in the metal it drives "
        "the carriers round in loops: eddy currents, which heat it, the faster the field "
        "changes the more (as $\\dot B^2$). Set the swing $B_0$ so that the plate takes "
        "between 1.2 and 4 of heat per unit time (averaged over $t = 25$ to $125$).",
        grid=(24, 20), c=None, t_max=200.0,
        electrodes=[cold_electrode(1, 1)],
        solenoid_fields=slider(3.0, 0.95),
        tube=tube(0, None, 0, 0, 0, 0, step=METAL_STEP, saturation=0.0, solenoid=(23, 19),
                  solenoid_drive=sine(1.0, EDDY_OMEGA),
                  heat_goals=[heat_goal(0, 1.2, 4.0, 25.0, 125.0)],
                  metals=[metal((6, 4), (18, 16), [], seed=21, **METAL)]),
        reference=[supply(23, 19, 0.95)])



# Eddy currents in Ohm's regime (the owner, 2026-10-20: Eddy currents' plate reads as a doped
# semiconductor, omega_c tau 0.38 at its answer; make one like copper's). The same plate with
# carriers that collide far more often, in a field that swings faster, so that every ratio
# deciding its eddy currents is small, as in copper: omega tau = 0.049 (the swing slow against
# the collisions), omega_c tau <= 0.18 over the whole slider (0.052 at the answer; copper's
# 0.0044 per tesla), omega/(4 pi sigma) = 0.10 (the edge charge keeps up with the swing), the
# plate 40 mean free paths across (the other's 12). The metal's own ratios, though, are a lightly
# doped semiconductor's: omega_p tau 0.70 and a screening length of 1.14 mean free paths match
# n-type silicon at 1e16 per cm^3 within x1.08 (game scales.rs, `resemblance`), where copper's are
# 400 and 0.0014 (real metals' omega_p tau 50-500: the step must follow the plasma oscillation,
# so the model affords about 1). Named for silicon (the owner's choice). The mean free path 0.3
# (tau = 0.12), cores of radius 0.0025 (a/2l = 0.004, the model's limit 0.005; 96 000 of them),
# density 2.7 (sigma = 0.32), the step 0.033 (omega_p dt = 0.19).
SILICON = dict(density=2.7, thermal_speed=2.0, mean_free_path=0.3, core_radius=0.0025,
               carriers_per_area=50.0)

SILICON_STEP = 0.033

# 13 whole cycles of the heat (it goes as Bdot^2, period pi/omega) in the goal's window.
SILICON_OMEGA = 0.13 * math.pi



def eddy_silicon():
    # Designed with `generator tube-scan` (target/scratch/copper_scratch.py): heat (preview)
    # 0.114, 0.471, 1.02, 1.11, 1.19, 1.67, 2.06, 4.42, 7.35, 16.4 at B0 = 0.1, 0.2, 0.3, 0.32,
    # 0.35, 0.38, 0.43, 0.6, 0.8, 1.2: 11.5 B0^2, local Ohm's law's 11.1 (sigma <Bdot^2> int u
    # dA; 40 mean free paths across, its walls matter less than the other plate's, 5-7 %
    # below), the preview scattering ~10 % at small swings. With the goal 1.2-4.0 B0 solved on
    # 0.40-0.55 (5.0 % of +-1.5): from 1.0. Verified 2.067 +- 0.012 at 0.43 (the other plate's
    # +- 0.28: the shorter collision time averages the thermal noise out). The swing's +-1.5
    # spans heat 0 to ~26 with omega_c tau <= 0.18.
    return level(
        "Eddy currents in silicon",
        "The last plate's carriers flew far between collisions, so that the field turned their "
        "paths by a good part of a radian on the way ($\\omega_c \\tau \\approx 0.4$; the Scales "
        "section). This plate is lightly doped silicon's: its carriers collide far more often, "
        "in a field that swings faster. The field barely turns them between collisions "
        "($\\omega_c \\tau \\approx 0.05$), so the eddy currents follow Ohm's law alone, as they "
        "do in copper at any laboratory field, and their heat grows as $B_0^2$. Set the swing "
        "$B_0$ so that the plate takes between 1 and 4 of heat per unit time (averaged over "
        "$t = 25$ to $125$).",
        grid=(24, 20), c=None, t_max=200.0,
        electrodes=[cold_electrode(1, 1)],
        solenoid_fields=slider(1.5, 0.43),
        tube=tube(0, None, 0, 0, 0, 0, step=SILICON_STEP, saturation=0.0, solenoid=(23, 19),
                  solenoid_drive=sine(1.0, SILICON_OMEGA),
                  heat_goals=[heat_goal(0, 1.0, 4.0, 25.0, 125.0)],
                  metals=[metal((6, 4), (18, 16), [], seed=21, **SILICON)]),
        reference=[supply(23, 19, 0.43)])



def laminations():
    # Designed with `generator tube-scan` (target/scratch/lam20_*.py): a transformer's core
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
    return level(
        "Laminations",
        "A transformer's core in a field that swings back and forth: the eddy currents waste "
        "its energy as heat. They run in loops across the whole core; cut it, and each loop "
        "must close within a narrower piece, so the currents are weaker and the heat falls, "
        "the narrower the pieces, the more. Drag across the core to cut it, any length and "
        "angle; you have three cuts. Make the core take at most 0.65 of heat per unit time "
        "(averaged over $t = 25$ to $125$).",
        grid=(24, 20), c=None, t_max=200.0,
        electrodes=[cold_electrode(1, 1)],
        max_cuts=3,
        tube=tube(0, None, 0, 0, 0, 0, step=METAL_STEP, saturation=0.0, b_z=0.68,
                  solenoid_drive=sine(1.0, EDDY_OMEGA),
                  heat_goals=[heat_goal(0, 0.0, 0.65, 25.0, 125.0)],
                  metals=[metal((2, 4), (22, 16), [], seed=122,
                                **dict(METAL, carriers_per_area=30.0))]),
        reference=[cut(2, 8, 22, 8), cut(2, 12, 22, 12), cut(12, 4, 12, 16)])
