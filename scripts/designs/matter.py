"""Arcs: metals, delivering beams (direction and energy goals), electrodes and circuits."""

from designs.common import (
    Auto,
    M,
    box,
    charge,
    circle_coil,
    dc,
    drive,
    level,
    magnet,
    metal_sphere,
    plate,
    player_plate,
    shot,
    slider,
    supply,
)

# =======================================================================================
# Metals (PHYSICS.md 2.6). A metal sphere is an equipotential: held at a potential V it
# acts from outside like a charge V a at its centre, plus the charges that nearby charges
# induce on it (their images). Values: V = 4e5 on a sphere of radius 2.5 is 1e6 of charge.

def high_voltage_dome():
    return level(
        "High-voltage dome",
        "A metal sphere held at high voltage, like the dome of a Van de Graaff generator. "
        "From outside it pushes like a big charge at its centre, but it is a conductor: a "
        "charge you put near it pulls the dome's charge towards itself. The dome lifts the "
        "beam; add a charge to bring it into the detector.",
        shots=[shot(1e-6, 1.0, (0, 10), 0.0, 0.5, box((27, 16, 30, 19)))],
        conductors=[metal_sphere(15, 4, 2.5, "potential", 3e5)],
        max_charges=3, magnitudes=[m * M for m in (0.25, 0.5, 1, 2)])



def image_charge():
    return level(
        "Image charge",
        "A large grounded sphere. Whatever charge you place near it, the metal answers with "
        "an opposite image charge inside: near the surface your charge is almost cancelled, "
        "far away the pair acts like a dipole. Steer the beam around the sphere.",
        shots=[shot(1e-6, 1.0, (0, 6), 0.0, 0.5, box((27, 13, 30, 17)))],
        conductors=[metal_sphere(14, 10, 3.5, "grounded")],
        max_charges=3, magnitudes=[m * M for m in (0.5, 1, 2, 4)])



def polarized_sphere():
    return level(
        "Polarised sphere",
        "An isolated, uncharged metal sphere next to a strong charge. It stays neutral, but "
        "its charges separate: the side facing the charge takes the opposite sign, the far "
        "side the same sign. The sphere becomes a dipole that moves with every charge you "
        "add. Bring the beam to the detector.",
        shots=[shot(1e-6, 1.0, (0, 4), 0.0, 0.5, box((27, 15, 30, 18)))],
        elements=[charge(8, 16, 4 * M)],
        conductors=[metal_sphere(16, 10, 2.5, "charge", 0.0)],
        max_charges=3, magnitudes=[m * M for m in (0.5, 1, 2, 4)], region=(4, 2, 24, 18))



# =======================================================================================
# Delivering beams: detectors that also require a direction or an energy on arrival, as
# the entrance of a next stage does (PHYSICS.md 6.1).

def injection():
    return level(
        "Injection",
        "The detector is the entrance of the next accelerator stage: the beam must not only "
        "hit it, it must enter moving along the axis (+x, within ±8°), or the next stage "
        "loses it. The beam already reaches the entrance, but at 20°: straighten it. The "
        "cone shows what is accepted.",
        shots=[shot(1e-6, 1.0, (0, 5), 20.0, 0.5, box((28, 12, 30, 17), direction=(0.0, 8.0)))],
        max_charges=3, magnitudes=[m * M for m in (0.25, 0.5, 1, 2)])



def soft_landing():
    # T0 = 0.5; the ions must arrive with 0.02 <= T <= 0.08: a potential hill of
    # 0.42-0.48 T0 at the detector, without reflecting them.
    return level(
        "Soft landing",
        "Ions for surface science must land gently, or they would smash the surface: they "
        "may arrive with at most a quarter of their launch energy. The two electrodes at "
        "the target already brake them; add charges so that they land slowly, on the "
        "target, without turning back.",
        # A retarding electrode pair at the target (as deceleration optics have); the
        # player sets the final braking and keeps the beam on the target.
        shots=[shot(1e-6, 1.0, (0, 10), 0.0, 0.5, box((29, 9, 30, 11), kinetic=(0.02, 0.12)))],
        elements=[charge(28, 6, 0.6 * M), charge(28, 14, 0.6 * M)],
        max_charges=3, magnitudes=[m * M for m in (0.1, 0.2, 0.3, 0.5, 0.75, 1)],
        region=(12, 2, 26, 18))



def collimator():
    return level(
        "Collimator",
        "Three rays leave the source at −10°, 0° and +10°. All three must reach the "
        "detector travelling parallel to the axis (within ±3°): a parallel beam, as optics "
        "downstream expect. Only charges of the particle's sign are available.",
        shots=[shot(1e-6, 1.0, (0, 10), a, 0.5, box((28, 6, 30, 14), direction=(0.0, 3.0)))
               for a in (-10.0, 0.0, 10.0)],
        max_charges=4, magnitudes=[m * M for m in (0.25, 0.5, 1, 2)], signs=(True, False),
        region=(4, 1, 20, 19))



# =======================================================================================
# Electrodes (PHYSICS.md 2.7): real plates and apertures with their fringe fields,
# computed by the boundary element method. V is the electrode potential; for the
# particle q = 1e-6, a potential difference dV costs q dV = 1e-6 dV of kinetic energy
# (T0 = 0.5 for the usual beams).

def deflection_plates():
    return level(
        "Deflection plates",
        "Real deflection plates, as in a cathode-ray tube: two metal plates at +V and −V. "
        "Between them the field is nearly uniform, but it bulges out at the ends (the "
        "fringe field) and the plates' charge rearranges when you bring charges near. Add "
        "a charge to bring the beam into the detector.",
        shots=[shot(1e-6, 1.0, (0, 10), 0.0, 0.5, box((27, 14, 30, 18)))],
        electrodes=[plate(11, 13, 8, kind="potential", value=-3e4),
                    plate(11, 7, 8, kind="potential", value=3e4)],
        max_charges=3, magnitudes=[m * M for m in (0.25, 0.5, 1, 2)], region=(16, 1, 26, 19))



def power_supply():
    # The top plate's supply: -120k, -60k, 60k or 120k lands the beam about 4, 2, -2 or
    # -4 cells from the axis; off (grounded) it flies straight.
    return level(
        "Power supply",
        "Instead of placing charges, turn a knob: the top deflection plate is connected to "
        "a power supply (yellow frame). Click the plate to switch it on, and set its "
        "potential with the slider in the panel. A positive plate pushes positive ions away, a negative one pulls "
        "them. Bring the beam into the detector.",
        shots=[shot(1e-6, 1.0, (0, 10), 0.0, 0.5, Auto((27, 0, 30, 20), "x", 2))],
        electrodes=[plate(11, 13, 8, tunable=True), plate(11, 7, 8)],
        supplies=slider(3.4e5, -1.2e5, -6e4, 6e4, 1.2e5), reference=[supply(11, 13, -1.2e5)])



# Circuits (PHYSICS.md 2.9-2.10): a source through a resistor (with an inductor or a
# switch) drives a plate; the plates of "Power supply" (the top one's capacitance against
# the grounded lower plate is C = 2.63, from the boundary elements). The plate levels are
# Newtonian (c = infinity), where the circuits' quasi-static model is exact; the ions
# (T0 = 0.5, m = 1) fly at 1 cell per time unit and cross the plates (x = 7 to 15) at t
# = 7 to 15 after their launch.

def charging_a_plate():
    # Through R = 3 the top plate charges with tau = RC = 7.9: as the ions cross it holds
    # 59 to 85 % of its supply's voltage. -120k lands the beam at y = 12.49, -60k at 11.22
    # (in "Power supply", charged at once, -60k gave 12). The supply solves from -140k
    # to -85k (generator window); the slider's ends (+-310k) and their middles lie
    # outside.
    return level(
        "Charging a plate",
        "The top deflection plate is connected to its power supply through a resistor R. "
        "A plate is a capacitor: its charge has to flow in through R, so its potential "
        "rises as $V(t) = V_s(1 - e^{-t/RC})$ towards the supply's $V_s$, with the time "
        "constant $\\tau = RC$, here about 8 (C comes from the plates' shape). The ions "
        "cross the plates while it is still charging: the Circuits panel shows V(t), with "
        "the launch marked. Set the supply so that the beam lands in the detector.",
        shots=[shot(1e-6, 1.0, (0, 10), 0.0, 0.5, box((27, 12, 30, 13)))],
        electrodes=[plate(11, 13, 8, tunable=True, drive=drive(dc(), 3.0)), plate(11, 7, 8)],
        supplies=slider(3.1e5, -1.2e5, -6e4, 6e4, 1.2e5), reference=[supply(11, 13, -1.2e5)],
        c=None)



def chopper():
    # The top plate carries no charge until a switch connects it to its supply through
    # R = 0.4 at t = 30 (tau = 1.05): the ion launched at t = 0 crosses the plates at
    # t = 7 to 15 and flies straight on; the one launched at t = 40 meets a fully
    # charged plate: 120k lands it at y = 5.54, 60k at 8.23 (none: 10). The supply
    # solves from 90k to 165k; the slider's ends (+-360k) and their middles lie outside.
    return level(
        "Chopper",
        "A switch connects the top plate to its power supply at t = 30, through a small "
        "resistor: the plate charges within about one time unit, so the switch turns the "
        "deflector on. Two ions are launched, at t = 0 and at t = 40. Set the supply so "
        "that each reaches its own detector: the first flies past the plates before the "
        "switch closes, the second after.",
        shots=[shot(1e-6, 1.0, (0, 10), 0.0, 0.5, box((27, 9, 30, 11))),
               shot(1e-6, 1.0, (0, 10), 0.0, 0.5, box((27, 5, 30, 7)), time=40.0)],
        electrodes=[plate(11, 13, 8, tunable=True, drive=drive(dc(), 0.4, switch=(False, 30.0))),
                    plate(11, 7, 8)],
        supplies=slider(3.6e5, -1.2e5, -6e4, 6e4, 1.2e5), reference=[supply(11, 13, 1.2e5)],
        c=None)



def ringing_plate():
    # With L = 4.66 in series (and R = 0.2) the top plate (C = 2.63) rings: omega =
    # 1/sqrt(LC) = 0.286, the first maximum at t = pi/omega = 11, overshooting the
    # supply's voltage by 79 %. The ions cross the plates at t = 7 to 15, around it:
    # 120k lands the beam at y = 3.98 (6 cells down; "Power supply"'s static plate gave
    # 4), 60k at 7.60. The supply solves from 107.5k to 170k (generator window); the
    # slider's ends (+-205k) and their middles (+-102.5k) lie outside.
    return level(
        "Ringing plate",
        "Now an inductor L sits between the supply and the top plate. With the plate's "
        "capacitance it makes an LC circuit: the charge rushes in, overshoots and swings "
        "around the supply's voltage at $\\omega = 1/\\sqrt{LC}$, slowly damped by the "
        "small resistance. The first swing peaks just as the ions pass, at almost twice "
        "the supply's voltage (see the Circuits panel). Set the supply so that the beam "
        "lands in the detector.",
        shots=[shot(1e-6, 1.0, (0, 10), 0.0, 0.5, box((27, 3, 30, 5)))],
        electrodes=[plate(11, 13, 8, tunable=True, drive=drive(dc(), 0.2, l=4.66)),
                    plate(11, 7, 8)],
        supplies=slider(2.05e5, -1.4e5, -7e4, 7e4, 1.4e5), reference=[supply(11, 13, 1.4e5)],
        c=None)



def induction_kick():
    # Jackson's Faraday level (5.15) with the coil on a circuit: its current rises
    # through R from 0, kappa(t) = kappa_inf (1 - e^{-t/tau}), tau = L/R = 20, with L =
    # N/c^2 its self-inductance (N = 448.5 for radius 8 and wire 0.1). A driven coil
    # needs a finite c: c = 1e4, where the circuit's quasi-static parameter (light time
    # across the arena / tau) is 1.5e-4, its corrections of order its square (the
    # first-order retardation of a closed current vanishes), and the ion (v ~ 1)
    # Newtonian to 1e-8. The
    # canonical angular momentum holds p_phi = -q kappa(t) A_phi(rho) for the ion
    # starting at rest: the kick follows the flux, not its rate. kappa_inf = 1e6 (V/R =
    # kappa_inf c^2 = 1e14, R = L/tau = 2.24e-7, V = 2.24e7): at 4 cells from the axis,
    # where A_phi ~ 1.6 per unit kappa, a kick p ~ 1.6. As in Faraday's level the
    # detector asks for T >= 0.3, which only the induced field supplies: three charges
    # of at most 0.2M at least 4 cells away give at most 0.15.
    tau, n_coil, c = 20.0, 448.5, 1e4
    r = n_coil / c ** 2 / tau
    coil = circle_coil(15, 10, 8.0, 0.0)
    coil["drive"] = drive(dc(1e6 * c ** 2 * r), r)
    return level(
        "Induction kick",
        "The coil's current now comes from a circuit: switched on at t = 0 through a "
        "resistor, it rises as $I(t) = (V/R)(1 - e^{-t/\\tau})$ with $\\tau = L/R$, L the "
        "coil's own inductance. While it rises, the changing flux induces an electric "
        "field around the axis that drives the particle at rest; once the current "
        "settles, the induction stops. Faraday's law in integral form: the kick depends "
        "only on how much the flux changed, not on how fast. Bring the particle into the "
        "detector with at least 0.3 of kinetic energy.",
        shots=[shot(1e-6, 1.0, (15, 6), 0.0, 1e-6, box((18, 11, 22, 15), kinetic=(0.3, 10.0)))],
        coils=[coil], c=c, t_max=200.0,
        max_charges=3, magnitudes=[m * M for m in (0.05, 0.1, 0.2)], region=(8, 10, 22, 17))



def pulse(high, rise, delay=0.0, width=1e4, fall=0.0, low=0.0):
    """A trapezoid source (PHYSICS.md 2.9): `low` until `delay`, up to `high` over `rise`,
    held for `width`, down over `fall` (its top the player's power supply sets)."""
    return {"kind": "pulse", "low": low, "high": high, "delay": delay, "rise": rise,
            "width": width, "fall": fall}



def betatron(core=1.34e8, guide=2.3e3, kinetic=(8.6, 11.3), core_slider=3.1e8,
             guide_slider=8.3e3, rise=300.0, detector=None):
    # Wideroe's betatron: an ion at rest at R = 8 between a core coil (radius 4) and a
    # guide coil (radius 12), concentric. Its canonical angular momentum is 0, so
    # p = q A_phi(R, t); a circle at R needs p = q R B_z(R, t). With both coils' strengths
    # rising in proportion that holds at every t iff A_phi(R) = R B_z(R) for the mixture:
    # the core must rise 2.057 times as fast as the guide (exact loop integrals; field
    # index n = -3.8 < 1, radially stable). The guide ramps at the player's rate g; the
    # core follows a trapezoid source through a fast circuit (L/R = 2, c = 1e4): its
    # current rises in step with the source for `rise`, then holds. The rise needs
    # kappa_top / rise = 2.057 g, i.e. a source top V = 2.057 g rise N / tau (N = 189.5,
    # its self-flux). When the core stops, the flux through the orbit stops growing while
    # the guide's field rises on: the orbit shrinks onto the detector between the coils,
    # with the energy the betatron gave it. A wrong ratio moves the orbit during the rise:
    # inward onto the detector too slow, or outward into the guide. The core's supply is
    # its source's top voltage (`supply_voltages`), the guide's its ramp rate. Measured
    # (generator trace): the orbit holds at r = 8.00 +- 0.04 through the rise while T
    # climbs to 4.6, then shrinks onto the detector (r = 6 to 7) at t = 567, T ~ 10.
    # Solving intervals (generator window): the core 130M to 142.5M (the energy is its
    # flux's), the guide 1.48k to 3.39k (and islands up to 3.59k: a wrong ratio settles
    # the orbit at another radius, which the extraction still reaches); the sliders
    # (+-310M, +-8.3k) keep 0, their ends and middles outside.
    c = 1e4
    n_core = 189.5
    tau = 2.0
    r = n_core / c ** 2 / tau
    core_coil = circle_coil(15, 15, 4.0, 0.0)
    core_coil["drive"] = drive(pulse(0.0, rise), r)
    core_coil["tunable"] = True
    guide_coil = circle_coil(15, 15, 12.0, 0.0)
    guide_coil["tunable"] = True
    detector = detector or box((8, 14, 9, 16), kinetic=kinetic)
    return level(
        "Betatron",
        "Wideroe's betatron accelerates without any electrode: a changing flux induces an "
        "electric field around it, and a magnetic field holds the particle on its circle. "
        "The particle starts at rest between two coils. The inner one (the core, on a "
        "circuit whose source rises and then holds) supplies most of the flux, the outer "
        "one (the guide, ramped) the field at the orbit. The orbit keeps its radius only if "
        "the flux inside it grows twice as fast as the guide field times the orbit's area, "
        "$\\Phi(R) = 2\\pi R^2 B(R)$: then the particle gains energy turn after turn on "
        "the same circle. When the core's source stops rising, the flux stops growing but "
        "the guide's field rises on, and the orbit shrinks onto the detector. Set the two "
        "supplies so that the particle arrives there with its energy in the window.",
        grid=(30, 30), c=c, t_max=900.0,
        shots=[shot(1e-6, 1.0, (23, 15), 90.0, 1e-9, detector)],
        coils=[core_coil, guide_coil],
        supplies=slider(core_slider, core), coil_rates=slider(guide_slider, guide),
        reference=[supply(19, 15, core), supply(27, 15, guide)])



def pulse_sorter(refs=(1.46e5, 9.6e4, -3.36e5), detectors=None):
    # Arc 5 finale: three ions, launched at t = 0, 30 and 60, pass three deflection plate
    # pairs in a row (x = 5-9, 16-20, 27-31; gaps of 8 cells), each top plate on its own
    # circuit and its own supply slider. A charges through R = 20 (tau ~ 25 with these
    # small plates' C ~ 1.25), B is switched on at t = 40, C rings (L = 130: period ~ 80).
    # Measured deflections per volt (generator check, linear at 30 kV), each ion relative
    # to the largest: A 0.28, 0.82, 1; B 0, 0.96, 1; C 1, 0.89, 0.13. Each ion weighs the
    # three circuits differently, so the three landing points fix the three supplies and
    # no slider is redundant; the player must read the circuits' plots. The reference
    # (+146k, +96k, -336k: 0.32, 0.21 and 0.73 of the slider's +-460k, away from its
    # guessable points) lands the ions at y = 12.3, 9.3 and 6.5, top to bottom in launch
    # order, every flight clear of the plates; the detectors are placed around them.
    # Small, low plates (4 x 2 cells, as in the microscope column): six full-size ones
    # took 5.5 s to set up, over the sandbox budget.
    def pair(x, d):
        return [plate(x, 14, 4, height=2.0, tunable=True, drive=d), plate(x, 6, 4, height=2.0)]
    electrodes = (pair(7, drive(dc(), 20.0)) + pair(18, drive(dc(), 0.4, switch=(False, 40.0)))
                  + pair(29, drive(dc(), 0.5, l=130.0)))
    detectors = detectors or [Auto((38, 0, 40, 20), "x", 1, own=True) for _ in range(3)]
    shots = [shot(1e-6, 1.0, (0, 10), 0.0, 0.5, d, time=t)
             for t, d in zip((0.0, 30.0, 60.0), detectors)]
    return level(
        "Pulse sorter",
        "Finale of the circuits: three ions, launched at t = 0, 30 and 60, fly through "
        "three pairs of deflection plates, and each top plate hangs on its own circuit "
        "with its own power supply. The first charges through a resistor, the second is "
        "switched on at t = 40, the third rings with an inductor. Each ion meets the "
        "three plates at a different stage of their circuits (the Circuits panel shows "
        "each plate's V(t), with the launches marked). Set the three supplies so that "
        "every ion lands in its own detector.",
        grid=(40, 20), shots=shots, electrodes=electrodes, c=None, t_max=160.0,
        supplies=slider(4.6e5, -3e5, -2e5, -1e5, 1e5, 2e5, 3e5),
        reference=[supply(x, 14, v) for x, v in zip((7, 18, 29), refs) if v])



def tune_the_lens():
    # The middle aperture's halves have a supply each; only both at 300k focus the three
    # rays into the detector (100k and 200k focus too weakly, 400k too strongly).
    apertures = []
    for x, tunable in ((10, False), (14, True), (18, False)):
        apertures += [plate(x, 15, 6, height=3.0, angle_deg=90.0, tunable=tunable),
                      plate(x, 5, 6, height=3.0, angle_deg=90.0, tunable=tunable)]
    return level(
        "Tune the lens",
        "A real Einzel lens: three metal apertures, the outer ones grounded. The two halves "
        "of the middle aperture have a power supply each. At the right voltage the lens "
        "focuses the three diverging rays into the small detector; too little and they "
        "stay apart, too much and they cross before it.",
        shots=[shot(1e-6, 1.0, (0, 10), a, 0.5, Auto((28, 0, 30, 20), "x", 1))
               for a in (-6.0, 0.0, 6.0)],
        electrodes=apertures, supplies=slider(4.9e5, 1e5, 2e5, 3e5, 4e5),
        reference=[supply(14, 15, 3e5), supply(14, 5, 3e5)])



def build_a_deflector():
    return level(
        "Build a deflector",
        "Now you build the electrodes. Place a metal plate at a potential of your choice "
        "near the beam (R turns it, Q/E changes its potential). A plate at a positive "
        "potential pushes positive ions away. Plates keep a cell away from other metal "
        "and from the other elements.",
        shots=[shot(1e-6, 1.0, (0, 10), 0.0, 0.5, Auto((27, 0, 30, 20), "x", 1))],
        max_plates=3, plate_voltages=slider(2.3e5, -1e5, -7.5e4, -5e4, 5e4, 7.5e4, 1e5),
        region=(5, 3, 22, 17), reference=[player_plate(10, 12, 7.5e4)])



def shielding():
    return level(
        "Shielding",
        "A strong charge next to the beam line throws the beam out of the arena. A "
        "grounded metal plate between them screens its field: the charge it induces on the "
        "plate cancels much of the field on the far side. Place a grounded plate.",
        shots=[shot(1e-6, 1.0, (0, 10), 0.0, 0.5, Auto((27, 0, 30, 20), "x", 2))],
        elements=[charge(15, 14, 2 * M)],
        max_plates=3, plate_voltages=[0.0], region=(5, 3, 25, 17),
        reference=[player_plate(15, 12, 0.0)])



def glass(x, y, length, thickness, eps, height=4.0, angle_deg=0.0):
    """A dielectric box of relative permittivity eps (PHYSICS.md 2.7)."""
    return {"center": [x, y, 0], "length": length, "thickness": thickness, "height": height,
            "angle_deg": angle_deg, "permittivity": eps}



def ferrite(x, y, length, thickness, mu, height=4.0, angle_deg=0.0):
    """A ferrite box (high-mu insulator) of relative permeability mu (PHYSICS.md 2.7)."""
    return {"center": [x, y, 0], "length": length, "thickness": thickness, "height": height,
            "angle_deg": angle_deg, "permeability": mu}



def ferrite_shield():
    # Ferrites' introduction: a strong magnet (2M at (15, 15)) beside the beam line,
    # behind a ferrite plate (mu = 1000, 1.5 cells thick, 2 high, from x = 10 to 20).
    # Measured (generator check): without the plate the magnet bends the beam 2.27 cells
    # up at the right edge; with it 0.98 (a taller plate shields more: 0.50 at height 4;
    # mu = 5 still 1.50). The detector (one cell) is placed around the reference's landing
    # with the plate; without it the beam lands ~1.3 cells higher, which the build checks
    # (NEEDS_BODY: the reference fails without the ferrite).
    return level(
        "Ferrite shield",
        "A ferrite is a magnetically soft insulator: it carries no current, but its "
        "magnetic domains line up with a field and draw the field lines into it, so "
        "behind a ferrite plate a magnet's field is much weaker (a magnetic shield; "
        "electric fields pass, as through any insulator). The strong magnet behind the "
        "plate still nudges the beam a little. Place magnets to bring it into the "
        "detector.",
        shots=[shot(1e-6, 1.0, (0, 10), 0.0, 0.5, Auto((28, 0, 30, 20), "x", 1))],
        elements=[magnet(15, 15, 2 * M)],
        ferrites=[ferrite(15, 12, 10, 1.5, 1000.0, height=2.0)],
        max_magnets=3, strengths=[m * M for m in (0.5, 1, 2)], region=(4, 1, 26, 9),
        reference=[magnet(17, 7, -1 * M)])



def glass_screen():
    # Shielding's charge (2M at (15, 14)) behind a glass slab (eps = 4, 1.5 cells thick,
    # 2 high, from x = 10 to 20) instead of a metal plate. The glass's bound charge
    # weakens the field beyond it but does not cancel it, as metal does: the beam is still
    # thrown off, less, and the player's charges must make up the rest. The build checks
    # that the reference fails without the glass (NEEDS_BODY): the glass matters
    # (measured: the solver's one charge, -1M at (22, 2), brings the beam in with it; the
    # beam crashes into that charge without it). Small and low (setup 0.8 s; a 12 x 2 x 4
    # slab took 3.7 s).
    return level(
        "Glass screen",
        "Shielding with glass instead of metal. A dielectric is an insulator: its charges "
        "cannot flow, but its molecules polarize, and the bound charge on its surfaces "
        "weakens the field that crosses it, here by roughly the permittivity ε = 4, never "
        "to zero as metal does. The strong charge behind the glass still pushes the beam "
        "off; place charges to bring it into the detector.",
        shots=[shot(1e-6, 1.0, (0, 10), 0.0, 0.5, box((27, 8, 30, 12)))],
        elements=[charge(15, 14, 2 * M)],
        dielectrics=[glass(15, 12, 10, 1.5, 4.0, height=2.0)],
        max_charges=3, magnitudes=[m * M for m in (0.25, 0.5, 1)], region=(4, 1, 26, 9))



def real_einzel_lens():
    # Three apertures (plate pairs with a gap of 4 cells around y = 10) at x = 10, 14,
    # 18; the outer ones grounded, the middle one at V. For positive ions a positive V
    # decelerates and focuses.
    apertures = []
    for x, kind, value in ((10, "grounded", None), (14, "potential", 2.5e5),
                           (18, "grounded", None)):
        apertures += [plate(x, 15, 6, height=3.0, angle_deg=90.0, kind=kind, value=value),
                      plate(x, 5, 6, height=3.0, angle_deg=90.0, kind=kind, value=value)]
    return level(
        "Real Einzel lens",
        "The Einzel lens of arc 1, as it is really built: three metal apertures, the "
        "middle one at high voltage. Its focus is not quite on the detector, and the rays "
        "leave at ±6° now. Add charges to bring all three rays into the small "
        "detector.",
        shots=[shot(1e-6, 1.0, (0, 10), a, 0.5, box((28, 9, 30, 11))) for a in (-6.0, 0.0, 6.0)],
        electrodes=apertures,
        max_charges=4, magnitudes=[m * M for m in (0.1, 0.2, 0.3, 0.5)], region=(20, 2, 27, 18))
