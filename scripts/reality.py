"""What each level represents, and how far it is the real thing (the owner, 2026-10-21:
"leaving it out gives wrong impression on what it represents").

Each level's entry (`REALITY[slug]`, attached by `levels.py` to the level file as
`reality`, and shown folded under the description, its headline in view) has

- `represents`: the real thing the level shows, in a few words (at most 80 characters);
- remarks, each a short line (at most 70 characters) and its explanation:
  `kept` what carries over from the real thing, `simpler` what is idealized, scaled or
  approximated, `omitted` what the model leaves out (at least three, at least one `kept`
  and one of the others; test `every_level_states_its_reality`);
- `electrodes`, `metals`: each one's function in the device, in the level's order (shown
  with its number wherever it is named).

The numerical model's own approximations (boundary elements, macroparticles, time steps)
are the game's "Physics model and its limits" (`level::model`), not repeated here. A text
change rebuilds nothing: `python scripts/levels.py --text`.
"""


def kept(line, detail):
    return ("kept", line, detail)


def simpler(line, detail):
    return ("simplified", line, detail)


def omitted(line, detail):
    return ("omitted", line, detail)


def R(represents, *remarks, electrodes=(), metals=()):
    return {"represents": represents, "remarks": list(remarks),
            "electrodes": list(electrodes), "metals": list(metals)}


REALITY = {}

# =======================================================================================
# The particle levels: the plane z = 0 of a 3D world.

COULOMB = kept(
    "Exact Coulomb fields",
    "The charges' fields are the exact 3D inverse-square fields, and the particle follows the "
    "exact equation of motion in them (relativistic where the level has a speed of light).")
FIXED = simpler(
    "Your charges are held fixed",
    "The particle pulls on them as hard as they push on it, but they stay put, as if "
    "infinitely heavy. Real instruments shape their fields with electrodes held at set "
    "voltages, not with point charges.")
PLANE = simpler(
    "Motion in one plane",
    "Everything is symmetric about the plane z = 0, so the particle never leaves it. A real "
    "beam also has a width across the plane, where the same fields act differently "
    "(focusing in the plane can defocus across it).")
VACUUM = omitted(
    "Residual gas",
    "The arena is a perfect vacuum. In a real tube the beam hits gas molecules, which scatter "
    "and slow it and, once ionized, screen the fields.")
ONE_AT_A_TIME = simpler(
    "One particle at a time",
    "The shots fly separately and never meet, as in a weak beam; in an intense beam the "
    "particles repel each other (the Beams arc shows it).")
SPHERE_MAGNETS = simpler(
    "Magnets are magnetized spheres",
    "A uniformly magnetized sphere has exactly a dipole's field outside it, which keeps the "
    "field simple and exact; real magnets are blocks, rings and pole pieces, and an "
    "instrument's magnet makes a field as uniform as it can.")
LOOP_COIL = simpler(
    "A coil is one thin loop",
    "The coil is a single loop of thin wire carrying a steady current, its field the exact "
    "Biot–Savart field (perpendicular to the plane in it, and far from uniform); a real "
    "magnet has many turns, iron and pole pieces shaped for a uniform field.")
POINT_ANTENNA = simpler(
    "Antennas are point dipoles",
    "Each antenna is an oscillating point dipole with its exact retarded fields (near, "
    "induction and radiation zones); real RF deflectors are resonant cavities or plates "
    "driven by a generator.")
SLOW_LIGHT = simpler(
    "Light is slow here",
    "c is a few cells per time unit, so that the board spans less than a wavelength of the "
    "RF and the particle crosses it in a few periods; the physics depends only on such "
    "ratios (v/c, the wavelength against the board), which real devices share.")

REALITY["first_bend"] = R(
    "Steering a charged particle with Coulomb's force",
    COULOMB,
    kept("Like repels, unlike attracts",
         "The force falls as the inverse square of the distance, so a charge bends the path "
         "strongly only where the particle passes close by."),
    FIXED, PLANE)

REALITY["slingshot"] = R(
    "A charge turned back by an attracting one",
    COULOMB,
    kept("Orbits are conic sections",
         "Past a single attracting charge the path is a hyperbola, as for a comet passing the "
         "Sun: it turns back only if it passes very close or slowly."),
    FIXED, PLANE)

REALITY["geiger_marsden"] = R(
    "Rutherford scattering: alpha particles off a gold nucleus",
    kept("Rutherford's law holds exactly",
         "The deflection follows tan(θ/2) = qQ/(2T₀b). For the Coulomb force quantum mechanics "
         "gives the same scattering as classical mechanics, which is why Rutherford's "
         "classical analysis was right."),
    simpler("One alpha aimed at one nucleus",
            "Geiger and Marsden fired many alphas through a gold foil thousands of atoms "
            "thick; each met nuclei at random distances and about 1 in 8000 bounced back. "
            "Here you aim one alpha at one nucleus."),
    simpler("The nucleus does not recoil",
            "It is held fixed. A gold nucleus, 49 times heavier than an alpha, recoils a "
            "little: an alpha bouncing straight back loses about 8 % of its energy to it."),
    omitted("The atom's electrons",
            "They screen the nucleus only far from it and are far too light to deflect the "
            "alpha; they are left out."),
    simpler("Real distances are femtometres",
            "A 5 MeV alpha comes within about 45 fm of a gold nucleus head-on; the cell-sized "
            "board is the same physics scaled up (the Scales section)."))

REALITY["around_the_wall"] = R(
    "Beam transport: steering a beam into the next stage",
    COULOMB,
    simpler("The wall is a row of charges",
            "A row of fixed point charges stands in for a charged obstacle; a real wall, a "
            "conductor or an insulator, would also carry charges induced by its "
            "surroundings."),
    FIXED, PLANE)

REALITY["thomson_crt"] = R(
    "Thomson's cathode-ray tube (1897): electrons deflected electrically",
    kept("The deflection falls as 1/T",
         "Through a given field a faster particle spends less time and is kicked less: the "
         "spot of an oscilloscope moves in proportion to the deflecting voltage over the "
         "beam's energy."),
    simpler("Point charges, not plates",
            "Thomson deflected the rays between two parallel metal plates at a voltage, a "
            "nearly uniform field; here you build the field from point charges."),
    omitted("The magnetic half of the experiment",
            "Thomson balanced the electric deflection against a magnetic one to measure the "
            "electrons' charge-to-mass ratio; only the electric part is here."),
    omitted("Residual gas",
            "Hertz had failed to see the electric deflection: ionized gas in his tube "
            "screened the plates. Thomson's better vacuum made it visible; here the vacuum "
            "is perfect."),
    PLANE)

REALITY["twin_beams"] = R(
    "One field setup serving two beams",
    COULOMB, ONE_AT_A_TIME, FIXED)

REALITY["hemispherical_analyzer"] = R(
    "A hemispherical electron energy analyzer (as in XPS)",
    kept("The field of a point charge",
         "Between two concentric hemispheres the field is exactly that of a point charge at "
         "their centre, which is what your charge makes: the orbits are the same Kepler "
         "ellipses, and the energies separate as in the instrument."),
    simpler("No hemispheres",
            "A real analyzer's field exists only between its two hemispherical electrodes, "
            "entered and left through slits; here the point charge's field fills the board."),
    omitted("The input lens",
            "Real analyzers slow the electrons in a lens before the hemispheres (a fixed pass "
            "energy) and image a band of energies onto a position-sensitive detector."),
    PLANE)

REALITY["reflectron"] = R(
    "The reflectron: the ion mirror of time-of-flight mass spectrometers",
    kept("Ions turn round in a retarding field",
         "Faster ions go deeper into the mirror before turning, so ions of different energies "
         "come back along different paths and at different times."),
    simpler("A mirror of point charges",
            "A real reflectron is a stack of ring electrodes making a uniform retarding field; "
            "here you build one from point charges."),
    omitted("The timing",
            "The instrument tells masses apart by their arrival times, and the mirror's point "
            "is that faster ions catch up with slower ones; here only where the ions arrive "
            "counts."),
    PLANE)

REALITY["einzel_lens"] = R(
    "An einzel lens of electron and ion optics",
    kept("A decelerating lens focuses",
         "Rays far from the axis are pushed back towards it more strongly than those near it, "
         "so the bundle converges, as in the electrostatic lenses of microscopes and ion "
         "sources."),
    simpler("Point charges, not cylinders",
            "A real einzel lens is three coaxial cylinders or apertures, the middle one at "
            "another voltage: round, so it focuses in every direction. Charges in one plane "
            "focus only in that plane."),
    simpler("Three rays",
            "A real beam is a continuum of rays, and a lens focuses its outer rays shorter "
            "(spherical aberration); the detector's size here tolerates some."),
    COULOMB)

REALITY["fast_lane"] = R(
    "Relativistic momentum: why fast particles are hard to bend",
    kept("Exact relativistic motion",
         "The momentum grows as γmv and the force changes the momentum, exactly at any speed: "
         "near c a particle gains momentum but hardly speed."),
    SLOW_LIGHT,
    kept("Radiation negligible",
         "The particle's charge is small so that it radiates far less than the accuracy, as "
         "a real beam particle does in fields like these."),
    FIXED)

REALITY["beta_spectrometer"] = R(
    "A beta-ray spectrometer for relativistic electrons",
    kept("The relativistic orbit condition",
         "A circular orbit needs |qQ| = γmv²R: at these speeds the Newtonian estimate would "
         "be wrong by tens of per cent."),
    simpler("Electric, not magnetic",
            "Real beta spectrometers bend the electrons in magnetic fields (semicircular "
            "180° focusing or lens types); a central electric field bends them here."),
    simpler("Two energies",
            "Beta decay gives electrons a continuous spread of energies (the neutrino takes "
            "the rest); two chosen energies stand for it."),
    PLANE)

REALITY["first_coil"] = R(
    "Magnetic bending: a charge circling in a coil's field",
    kept("The magnetic force does no work",
         "q v × B is perpendicular to the velocity: it turns the particle but never changes "
         "its speed (watch the energy bars)."),
    LOOP_COIL, FIXED)

REALITY["dempster"] = R(
    "Dempster's magnetic mass spectrometer (1918)",
    kept("180° focusing",
         "Ions leaving at slightly different angles meet again after half a turn in a "
         "magnetic field, the principle of Dempster's instrument."),
    simpler("A loop's field, not a uniform one",
            "Dempster bent the ions between the flat poles of an electromagnet, a uniform "
            "field; the large loop's field falls off from its centre, so r = √(2mT)/(qB) "
            "holds only roughly."),
    simpler("Masses 1, 2 and 4",
            "Real isotopes differ in mass by a few per cent; masses differing by factors of "
            "two separate widely."),
    simpler("Accelerated by point charges",
            "Dempster accelerated the ions with a voltage across slits; your charges play "
            "that part."),
    PLANE)

REALITY["wien_filter"] = R(
    "A Wien filter (velocity selector)",
    kept("Crossed fields pass v = E/B",
         "Where the electric push qE balances the magnetic qvB the particle flies straight; "
         "slower ones are pushed one way, faster ones the other."),
    simpler("Fields that only roughly balance",
            "A real Wien filter has parallel plates (a uniform E) inside a magnet's uniform B, "
            "over the same length; here B is a loop's field and E comes from point charges, "
            "so they balance only approximately along the path."),
    PLANE)

REALITY["first_magnet"] = R(
    "A permanent magnet's field bending a beam",
    kept("A dipole's field",
         "Outside a uniformly magnetized sphere the field is exactly a dipole's: in the plane "
         "it points across the plane and falls as 1/r³, so only close passes bend much."),
    SPHERE_MAGNETS,
    simpler("Magnets never change",
            "A magnet keeps its magnetization whatever field it is in (an ideally hard "
            "magnet), and is held in place."))

REALITY["calutron"] = R(
    "Lawrence's calutron: separating uranium isotopes magnetically",
    kept("Heavier ions bend less",
         "At equal energy the radius r = √(2mT)/(qB) grows as the square root of the mass, "
         "which is how the calutrons separated uranium-235 from uranium-238."),
    simpler("Masses 1, 2 and 4",
            "U-235 and U-238 differ by 1.3 % in mass: after a semicircle of 1.2 m radius "
            "their beams landed about 1.5 cm apart."),
    simpler("Magnets you place",
            "The calutrons bent the ions in the uniform field of a huge electromagnet; here "
            "dipole magnets make a field far from uniform."),
    omitted("Space charge",
            "The real beams were intense and pushed themselves apart; here the ions fly one "
            "at a time (Calutron with space charge shows it)."),
    PLANE)

REALITY["build_wien_filter"] = R(
    "A velocity selector built from charges and magnets",
    kept("Crossed fields pass v = E/B",
         "Where the electric and magnetic pushes balance, one speed passes straight; slower "
         "and faster particles are pushed to opposite sides."),
    simpler("Point sources only",
            "Uniform crossed fields balance everywhere; the fields of a few charges and "
            "dipole magnets balance only approximately along the path, and the exit's "
            "tolerance absorbs the rest."),
    SPHERE_MAGNETS, PLANE)

REALITY["stern_gerlach"] = R(
    "The Stern–Gerlach experiment (1922): spin splits a beam in two",
    kept("The force of a field gradient",
         "A neutral particle with a magnetic moment feels no Lorentz force, only the pull of "
         "the field's gradient, m ∇B_z, as the silver atoms did."),
    simpler("Two classical spin states",
            "Each state flies as its own shot, its moment fixed along z. In quantum mechanics "
            "each atom is in a superposition, and the field picks one of two outcomes: the "
            "splitting into exactly two beams is the quantum part, put in by hand here."),
    simpler("The magnet",
            "Stern and Gerlach used an electromagnet with a knife-edge pole over a grooved "
            "one for a strong gradient; the beams split by a fraction of a millimetre. A "
            "spherical magnet's dipole field here."),
    omitted("The oven's speed spread",
            "The atoms came from an oven with a broad spread of speeds, which smeared each "
            "beam; here each state has one speed."))

REALITY["stray_field"] = R(
    "Keeping a beam on target against a stray field",
    kept("A uniform field acts exactly",
         "The stray field pushes the particle as a real one would; the motion is Newtonian, "
         "as at these speeds."),
    simpler("A uniform field, on or off",
            "Stray fields in a lab come from cables and supplies and vary in space and time; "
            "here a uniform field is either on or off."),
    omitted("Shielding",
            "A real beam line would sit in a grounded metal pipe that keeps outside electric "
            "fields out; here the setup must work despite the field."),
    FIXED)

REALITY["mains_hum"] = R(
    "Mains interference on a beam, and an imaging lens",
    kept("Imaging cancels the launch angle",
         "A lens that images the source onto the detector brings every ray from the source "
         "there, whatever its angle."),
    kept("The same field everywhere at once",
         "The field changes in step over the whole board (c = ∞), as mains fields do over a "
         "lab: their wavelength is thousands of kilometres."),
    simpler("The hum is fast here",
            "Real hum (50 or 60 Hz) is frozen during a flight of nanoseconds and shifts each "
            "particle's path as a slowly changing static field; here it changes during the "
            "flight, which acts as a random launch angle."),
    PLANE)

REALITY["earths_field"] = R(
    "A cathode-ray tube's sensitivity to the Earth's magnetic field",
    kept("Bending grows with the flight time",
         "A weak field bends the path as the square of the time spent in it, so slow "
         "particles over long paths are the most sensitive."),
    simpler("Only the field across the plane",
            "The Earth's field (about 50 µT, which moves a picture tube's beam by "
            "millimetres) also has parts along the tube, which would push the beam out of "
            "the plane; only the part across the plane is here."),
    simpler("Corrected with charges",
            "Real tubes were corrected with small magnets and coils around the neck, or "
            "shielded; here charges do it."),
    PLANE)

REALITY["rf_kick"] = R(
    "Deflection by a radio-frequency field",
    kept("The kick depends on the phase",
         "The antenna's field reverses every half period, so what a passing particle gets "
         "depends on when it passes."),
    POINT_ANTENNA, SLOW_LIGHT)

REALITY["rf_separator"] = R(
    "An RF separator, as at CERN for kaons and pions",
    kept("Time-dependent fields tell bunches apart",
         "Particles that reach an oscillating field at different phases are kicked "
         "differently, which no static field can do."),
    simpler("Identical particles, delayed",
            "Real separators sort particles of equal momentum but different mass: after a "
            "long drift their different speeds bring them to the deflecting cavities at "
            "different phases. Here two identical particles leave half a period apart."),
    POINT_ANTENNA, SLOW_LIGHT)

REALITY["streak_camera"] = R(
    "A streak camera: time turned into position",
    kept("A sweeping field maps time to place",
         "What arrives early is deflected one way and what arrives late the other, so the "
         "screen shows the arrival times side by side."),
    simpler("Real streak cameras",
            "Light hits a photocathode and its photoelectrons are swept across a phosphor "
            "screen by a voltage ramp on deflection plates, resolving picoseconds; here an "
            "antenna's sinusoidal field sweeps bunches a third of a period apart."),
    POINT_ANTENNA)

REALITY["synchrotron_light"] = R(
    "Synchrotron radiation damping: a circling charge spirals in",
    kept("Radiation reaction",
         "The particle radiates as Larmor's formula says and loses that energy (the "
         "Landau–Lifshitz force), exact in classical electrodynamics while the loss per turn "
         "is small."),
    simpler("A hugely radiating particle",
            "Its charge to mass is chosen so that it radiates 0.1 % of its energy per turn "
            "at modest speed. A real electron this slow would radiate almost nothing; real "
            "synchrotron light comes from electrons with γ in the thousands, beamed forward."),
    omitted("Re-acceleration",
            "Electron rings replace the lost energy with RF cavities every turn; here nothing "
            "does, so the orbit shrinks."),
    LOOP_COIL)

REALITY["tune_the_rf"] = R(
    "Choosing an RF frequency to separate bunches",
    kept("The phase difference is ωΔt",
         "Two bunches Δt apart meet the field at phases ωΔt apart; they are kicked oppositely "
         "near half a period."),
    POINT_ANTENNA, SLOW_LIGHT)

# The metal spheres, electrodes, dielectrics and ferrites.

IDEAL_SUPPLY = simpler(
    "Ideal supplies",
    "A plate on a supply sits at its voltage whatever charge that takes; a real supply has "
    "an internal resistance and leads, and charging a plate takes a current for a while "
    "(the Circuits arc).")
PLATE_FIELDS = kept(
    "Real plates and their fringe fields",
    "The electrodes are metal boxes standing on the plane, their surface charge solved "
    "exactly (boundary elements): uniform between plates, bulging out past their ends, "
    "and rearranging when charges come near, as in a real instrument.")
SLITS = simpler(
    "Slits, not round apertures",
    "Each aperture is a pair of plates a few cells high, a slit; a real einzel lens is "
    "round (apertures or tubes) and focuses in every direction at once, a slit lens only "
    "across the slit.")
LUMPED = simpler(
    "Ideal lumped circuits",
    "Resistors, capacitors, inductors and switches are ideal; the leads have no "
    "resistance, inductance or delay, and the circuit acts on the plate as if the field "
    "spread instantly (quasi-static, here exact: the motion is Newtonian).")

REALITY["high_voltage_dome"] = R(
    "The high-voltage dome of a Van de Graaff generator",
    kept("A conductor is an equipotential",
         "The dome's charge arranges itself so that its whole surface has one potential: "
         "from outside it acts like a charge at its centre, plus the images of the charges "
         "near it."),
    simpler("Held at its voltage",
            "An ideal supply holds the dome at its potential; a real dome is charged by a "
            "belt and loses charge into the air (corona) until the two balance."),
    FIXED, PLANE)

REALITY["image_charge"] = R(
    "Image charges: a grounded metal sphere",
    kept("The images are exact",
         "For a sphere the method of images is exact: a charge q at distance d from a "
         "grounded sphere of radius a acts with the image −qa/d at a²/d from its centre."),
    kept("The metal answers at once",
         "A metal rearranges its charge within femtoseconds, instantly for a beam."),
    omitted("The particle's own image",
            "The particle's charge also induces an image that pulls it towards the metal; "
            "for its tiny charge this matters only within a hair of the surface, and it is "
            "left out."),
    PLANE)

REALITY["polarized_sphere"] = R(
    "Electrostatic induction: a neutral metal sphere polarized",
    kept("Neutral, but polarized",
         "An isolated conductor keeps its net charge (here zero) but its charges separate: "
         "the near side takes the opposite sign, the far side the same, exactly as images "
         "give."),
    kept("Every charge moves the induced charge",
         "The induced charge answers every charge you add, so the sphere's field changes "
         "with your setup, as real conductors do."),
    FIXED, PLANE)

REALITY["injection"] = R(
    "Injecting a beam into the next accelerator stage",
    kept("Position and angle both count",
         "Each stage of an accelerator accepts particles only within a range of positions "
         "and angles, its acceptance; a beam outside it is lost."),
    simpler("A box and a cone",
            "A real acceptance is a region of phase space (positions and angles together, "
            "an ellipse); here it is a box with an angle cone, for one particle."),
    FIXED, PLANE)

REALITY["soft_landing"] = R(
    "Soft landing: slowing ions gently onto a surface",
    kept("A potential hill takes the energy",
         "An ion arrives with its launch energy less q times the potential it climbed: the "
         "braking electrodes take the energy back without stopping it short."),
    simpler("Brakes made of charges",
            "Real soft-landing setups slow the ions with electrodes biased close to the ions' "
            "energy, landing them at a few eV per charge so that fragile molecules survive; "
            "point charges brake them here."),
    omitted("The surface",
            "The target is a detector: what the ion does on arrival (sticking, breaking "
            "apart, sputtering) is not modelled."),
    PLANE)

REALITY["collimator"] = R(
    "Collimating a diverging beam into a parallel one",
    kept("A lens with the source at its focus",
         "Rays from a point at a lens's focal point leave parallel: the same optics as for "
         "light."),
    simpler("A lens of point charges",
            "Real collimating lenses are electrodes shaped for small aberrations; a few "
            "point charges make a lens with large ones, which the angle tolerance absorbs."),
    FIXED, PLANE)

REALITY["deflection_plates"] = R(
    "The deflection plates of a cathode-ray tube",
    PLATE_FIELDS,
    simpler("Plates at fixed voltages",
            "The plates are held at ±V by ideal supplies; their charge follows what you add "
            "nearby, their potential does not."),
    FIXED, PLANE,
    electrodes=["upper deflection plate (−V)", "lower deflection plate (+V)"])

REALITY["power_supply"] = R(
    "Deflection set by a power supply's knob",
    kept("The voltage sets the deflection",
         "The sideways push grows in proportion to the plate's voltage, its sign with the "
         "voltage's: an oscilloscope's deflection works so."),
    PLATE_FIELDS, IDEAL_SUPPLY, PLANE,
    electrodes=["deflection plate (on your supply)", "deflection plate (grounded)"])

REALITY["charging_a_plate"] = R(
    "An RC circuit charging a deflection plate",
    kept("The plate is a capacitor",
         "Its capacitance comes from its shape (the boundary elements); its voltage rises as "
         "V_s(1 − e^(−t/RC)), so a beam passing early sees part of it. In an oscilloscope "
         "the plates charge in nanoseconds, which limits how fast a signal it can follow."),
    LUMPED, PLATE_FIELDS, PLANE,
    electrodes=["deflection plate (charged through R)", "deflection plate (grounded)"])

REALITY["chopper"] = R(
    "A beam chopper: a deflector switched on at a set time",
    kept("A switch and a small RC",
         "Once the switch closes, the plate charges within RC, about one time unit: what "
         "passes before is undeflected, what passes after is swept aside, as in the pulsed "
         "deflectors that chop accelerator beams into bunches."),
    simpler("An ideal switch",
            "It closes in no time, without bounce or arc; real fast choppers use "
            "semiconductor switches and still take nanoseconds."),
    LUMPED, PLANE,
    electrodes=["deflection plate (switched on at t = 30)", "deflection plate (grounded)"])

REALITY["ringing_plate"] = R(
    "An LC circuit ringing on a deflection plate",
    kept("LC resonance and overshoot",
         "The charge rushes in through L, overshoots and swings at ω = 1/√(LC), damped "
         "slowly by R; fast pulsed deflectors ring the same way through their leads' stray "
         "inductance."),
    LUMPED, PLATE_FIELDS,
    electrodes=["deflection plate (on an LC circuit)", "deflection plate (grounded)"])

REALITY["induction_kick"] = R(
    "Faraday's law: an ion kicked by a coil's rising flux",
    kept("The kick follows the flux change",
         "While the current rises, the changing flux induces an electric field around the "
         "axis; the momentum it gives depends only on how much the flux changed (the "
         "canonical momentum), not on how fast."),
    kept("The coil's own inductance",
         "L comes from the coil's own flux, and the current rises with τ = L/R."),
    simpler("Light is very fast here",
            "c = 10⁴ makes the coil's field follow its current at once (quasi-static), as "
            "in any coil much smaller than a wavelength."),
    LOOP_COIL)

REALITY["betatron"] = R(
    "The betatron: accelerating by induction (Wideröe's 2:1 rule)",
    kept("The 2:1 condition",
         "The orbit keeps its radius only if the flux through it grows twice as fast as the "
         "guide field times its area; then the induced field accelerates the particle turn "
         "after turn on the same circle."),
    simpler("Two loops, not a shaped magnet",
            "Real betatrons (Kerst's first in 1940, for electrons) use an iron magnet whose "
            "poles shape both the central flux and the guide field; two current loops make "
            "them here."),
    simpler("Stable only in the plane",
            "Here the guide field grows outward (field index −3.8), which holds the orbit's "
            "radius but would throw a real particle out of the plane: real betatrons need "
            "an index between 0 and 1. The plane's symmetry hides that."),
    PLANE)

REALITY["pulse_sorter"] = R(
    "Timed deflectors: sorting ions by their circuits' timing",
    kept("Each circuit has its own V(t)",
         "A plate charging through R, one switched on, one ringing through L: each ion "
         "meets each plate at a different stage, as fast deflectors in a beam line do."),
    LUMPED, ONE_AT_A_TIME, PLANE,
    electrodes=["plate A (charged through R)", "plate A (grounded)",
                "plate B (switched on at t = 40)", "plate B (grounded)",
                "plate C (ringing through L)", "plate C (grounded)"])

REALITY["tune_the_lens"] = R(
    "An einzel lens: three apertures, the middle one at a voltage",
    kept("Weak, right and strong focusing",
         "The focal length shortens as the middle electrode's voltage rises: too little and "
         "the rays stay apart, too much and they cross before the detector."),
    PLATE_FIELDS, SLITS,
    electrodes=["entrance aperture, upper half", "entrance aperture, lower half",
                "middle aperture, upper half (your supply)",
                "middle aperture, lower half (your supply)",
                "exit aperture, upper half", "exit aperture, lower half"])

REALITY["build_a_deflector"] = R(
    "Electrostatic deflection by plates you place",
    PLATE_FIELDS,
    simpler("Plates at ideal fixed voltages",
            "Each plate you place sits at its chosen potential whatever charge that takes, "
            "and stays where you put it."),
    PLANE)

REALITY["shielding"] = R(
    "Electrostatic shielding by a grounded plate",
    kept("Induced charge cancels the field",
         "The charge the field induces on the grounded plate cancels much of it beyond the "
         "plate."),
    simpler("A plate, not an enclosure",
            "A plate screens only partly: the field reaches round its edges. A closed "
            "grounded box (a Faraday cage) would screen completely."),
    PLATE_FIELDS, PLANE)

REALITY["ferrite_shield"] = R(
    "Magnetic shielding by a ferrite plate",
    kept("A soft magnet draws the field in",
         "A material of high permeability (μ = 1000) carries the field lines through itself, "
         "so beyond it the magnet's field is much weaker; electric fields pass, as through "
         "any insulator."),
    simpler("An ideal ferrite",
            "Its permeability is constant: no saturation (real ferrites saturate near "
            "0.3–0.5 T), no hysteresis, no dependence on frequency; it is magnetized by the "
            "magnets and coils only."),
    SPHERE_MAGNETS, PLANE)

REALITY["glass_screen"] = R(
    "Partial screening by a dielectric (glass)",
    kept("Bound charge weakens the field",
         "The glass's molecules polarize, and the bound charge on its surfaces weakens the "
         "field that crosses it, by roughly its permittivity, never to zero as metal does."),
    simpler("Static, linear glass",
            "ε = 4 at every field and frequency, with no breakdown and no losses."),
    omitted("Charging up",
            "Real insulators near beams collect stray charge on their surfaces, which "
            "deflects the beam (why beam lines avoid exposed insulators); this glass "
            "never charges."),
    FIXED, PLANE)

REALITY["real_einzel_lens"] = R(
    "An einzel lens as it is built: three metal apertures",
    kept("A real lens's field",
         "The apertures' surface charge is solved exactly, so the lens has its true focal "
         "length and aberrations: its focus is not quite on the detector, as a real lens's "
         "rarely is."),
    SLITS, FIXED,
    electrodes=["entrance aperture, upper half (grounded)",
                "entrance aperture, lower half (grounded)",
                "middle aperture, upper half (high voltage)",
                "middle aperture, lower half (high voltage)",
                "exit aperture, upper half (grounded)", "exit aperture, lower half (grounded)"])

REALITY["two_stages"] = R(
    "An instrument's stages, passed in order",
    COULOMB,
    simpler("A stage is a box",
            "Real stages are lenses, filters and slits that act on the beam; here a stage "
            "is a region the particle must pass, acting on nothing."),
    FIXED, PLANE)

# Beams: many particles at once.

SPACE_CHARGE = kept(
    "Space charge, exactly",
    "The particles repel each other with the exact Coulomb force (with c = ∞ that is the "
    "whole interaction), so an intense beam spreads as a real one does.")
BEAM_SAMPLE = simpler(
    "A beam of a few particles",
    "A dozen or so particles stand for a beam of billions: each carries the charge of many "
    "(its mass and energy scaled alike, so that its own path is unchanged), so their "
    "repulsion matches an intense beam's. The spreads are a fixed sample of the beam's "
    "distribution, the same every time.")
NEUTRALIZATION = omitted(
    "Neutralization by the residual gas",
    "In a real beam line the beam ionizes the remaining gas, and the electrons it frees "
    "cancel much of a positive beam's space charge; in this perfect vacuum nothing does.")
LENS_OF_CHARGES = simpler(
    "Lenses made of point charges",
    "Real ion optics use shaped electrodes with small aberrations; a few point charges make "
    "lenses and deflectors with large ones, which the detectors' tolerances absorb.")

REALITY["space_charge"] = R(
    "Space charge: an intense beam spreading under its own repulsion",
    SPACE_CHARGE, BEAM_SAMPLE, NEUTRALIZATION, PLANE)

REALITY["stern_gerlach_beam"] = R(
    "The Stern–Gerlach experiment with a real beam",
    kept("The force of a field gradient",
         "Each atom's moment is pulled along the field's gradient; the state pulled towards "
         "the magnet passes closer and fans out more, as in the real splitting's uneven "
         "lines."),
    kept("Neutral atoms barely interact",
         "Their magnetic dipole forces on each other are tiny; in this Newtonian level they "
         "are none."),
    simpler("Two classical spin states",
            "Each state flies as its own beam, its moment fixed along z; the quantum "
            "superposition and the measurement's choice of one of two outcomes are put in "
            "by hand."),
    simpler("A narrow speed spread",
            "2 % in energy; the oven's atoms had a broad thermal spread, which slits "
            "narrowed only in direction."))

REALITY["collimated_beam"] = R(
    "Collimating an intense, diverging beam",
    SPACE_CHARGE, BEAM_SAMPLE, LENS_OF_CHARGES, PLANE)

REALITY["velocity_selector"] = R(
    "A Wien filter for an intense beam",
    kept("Crossed fields pass v = E/B",
         "The middle speed flies straight where the electric and magnetic pushes balance."),
    SPACE_CHARGE, BEAM_SAMPLE, LOOP_COIL)

REALITY["beam_preparation"] = R(
    "Preparing a beam: collimation, then delivery",
    SPACE_CHARGE, BEAM_SAMPLE,
    simpler("A stage is a box",
            "The gate is a region the beam must cross parallel; real stages are lenses and "
            "slits that also remove what is off axis."),
    PLANE)

REALITY["relativistic_beam"] = R(
    "A relativistic beam: magnetic attraction holds it together",
    kept("Net repulsion 1/γ²",
         "Charges moving side by side attract magnetically, which cancels all but 1/γ² of "
         "their electric repulsion: at 0.8 c, 36 %. Fast beams hold together for this "
         "reason."),
    kept("Radiation reaction",
         "Each particle feels its own radiation's reaction (Landau–Lifshitz): with these "
         "scaled charges it radiates too much to neglect."),
    simpler("Fields from where the others are now",
            "Each feels the others' fields computed from their present motion, exact in the "
            "velocities and checked against the exact retarded fields; real fields come from "
            "where the others were a light-time earlier."),
    BEAM_SAMPLE)

REALITY["chromatic_aberration"] = R(
    "Chromatic aberration of an electrostatic lens",
    kept("Faster ions focus further away",
         "A lens bends a faster particle less, so a spread of energies smears the focus "
         "along the axis, as in every charged-particle lens."),
    simpler("An exaggerated energy spread",
            "8 % in energy, to make the aberration large; real sources have spreads from a "
            "fraction of an eV (field emitters) to tens of eV (plasma sources) on beams of "
            "kilovolts."),
    LENS_OF_CHARGES,
    simpler("The ions do not interact",
            "Each ion flies alone: the aberration shows without space charge."))

REALITY["real_analyzer"] = R(
    "A hemispherical analyzer with a real source's angular spread",
    kept("Focusing to first order only",
         "Half a turn in a central field brings rays of slightly different angles together "
         "to first order; the image still broadens as the square of the angle, which limits "
         "a real analyzer's resolution."),
    simpler("No hemispheres",
            "The point charge's field fills the board; a real analyzer's exists only between "
            "its electrodes, entered through a slit."),
    PLANE)

REALITY["crt_earth_field"] = R(
    "A cathode-ray tube adjusted for the Earth's field",
    kept("A weak field bends a slow beam",
         "Over the flight the field turns the electrons a little, the more the longer they "
         "take: enough to move a spot."),
    simpler("Only the field across the plane",
            "The real field (about 50 µT) also has parts along the tube; only the part "
            "across the plane is here."),
    simpler("Corrected with charges",
            "Real tubes were corrected with small magnets and coils around the neck."),
    PLANE)

REALITY["calutron_space_charge"] = R(
    "The calutron at production currents",
    kept("Space charge spreads the spots",
         "The ions repel each other, so each isotope's spot grows into its neighbour's."),
    BEAM_SAMPLE, NEUTRALIZATION,
    kept("Radiation included",
         "With the scaled charges the ions would radiate more than may be neglected, so "
         "their radiation's reaction is included."))

REALITY["beam_pipe"] = R(
    "A beam pipe screening a steering charge",
    kept("The wall's induced charge screens",
         "The grounded wall takes up charge that cancels much of your charge's field inside "
         "the pipe."),
    simpler("One wall, not a pipe",
            "The pipe is a single grounded plate along one side; a real pipe surrounds the "
            "beam, and steering is done with magnets or with electrodes inside it."),
    PLATE_FIELDS, PLANE,
    electrodes=["beam pipe wall (grounded)"])

REALITY["soft_landing_current"] = R(
    "Soft landing with an intense ion beam",
    kept("Space charge grows as the beam slows",
         "The same current packed into a slower beam is a denser charge: braked ions crowd "
         "and blow up just before the target."),
    BEAM_SAMPLE,
    omitted("The surface",
            "The target is a detector: what the ions do on arrival is not modelled."),
    PLANE)

REALITY["sorting_station"] = R(
    "An ion-optical bench: lens, energy sorter, lens",
    COULOMB, LENS_OF_CHARGES, FIXED, PLANE)

REALITY["microscope_column"] = R(
    "The column of a charged-particle microscope",
    kept("Condenser, crossover, deflector",
         "A lens images the source through a small crossover, and a deflector steers the "
         "beam behind it, as in a microscope's column."),
    simpler("Electrostatic, with positive particles",
            "Most electron microscopes focus with magnetic lenses; electrostatic columns "
            "like this serve ion beams (as in focused-ion-beam instruments). The particles "
            "here are positive."),
    simpler("The ion pump is a fixed charge",
            "Its high voltage acts on the beam as a point charge; a real pump is shielded "
            "from the column."),
    SLITS,
    electrodes=["condenser entrance aperture, upper", "condenser entrance aperture, lower",
                "condenser middle aperture, upper (your supply)",
                "condenser middle aperture, lower (your supply)",
                "condenser exit aperture, upper", "condenser exit aperture, lower",
                "deflection plate (your supply)"])

REALITY["mass_spectrometer"] = R(
    "A mass spectrometer: an electrostatic lens and a magnetic sector",
    kept("Electric optics ignore the mass",
         "At equal energy every mass follows the same path through electric fields (only "
         "T/q matters), so one lens focuses all of them."),
    kept("Magnets sort by momentum",
         "In a magnetic field the momentum √(2mT) decides the radius, so the masses "
         "separate."),
    SPHERE_MAGNETS, FIXED, PLANE)

REALITY["rf_beam_line"] = R(
    "A beam line: steering robust to a stray field, then an RF separator",
    kept("Only a changing field separates identical bunches",
         "Two identical bunches half a period apart differ only in when they arrive; an "
         "oscillating field tells them apart."),
    POINT_ANTENNA,
    simpler("A uniform stray field, on or off",
            "The building's stray field is uniform and either on or off."),
    PLANE)

REALITY["isotope_separator"] = R(
    "An isotope separator: collimation, then a magnetic sector",
    kept("The lighter isotope bends more",
         "At equal energy the lighter ion has less momentum, so a magnetic field bends it "
         "more: both bend the same way."),
    SPACE_CHARGE, BEAM_SAMPLE,
    simpler("Masses 1 and 2",
            "Real isotopes differ by a few per cent in mass; a factor of two separates "
            "widely."),
    SPHERE_MAGNETS)

# The Jackson series.

UNIFORM = kept(
    "Uniform static fields, as in the book",
    "The fields across the arena are exactly uniform and static, as Jackson assumes, and "
    "the motion in them is exact.")

REALITY["jackson_exb_drift"] = R(
    "E×B drift in crossed uniform fields (Jackson §12.3)",
    kept("Drift E×B/B² for any charge",
         "A particle starting at rest rolls along a cycloid and drifts at E×B/B² whatever "
         "its charge and mass: the drift that moves plasma in the ionosphere, magnetrons and "
         "Hall thrusters."),
    UNIFORM,
    simpler("Guiding centres follow equipotentials",
            "Only where your charge's field changes little over a gyration's size; the "
            "motion itself is computed exactly."),
    FIXED)

DIPOLE_EARTH = kept(
    "The Earth as a dipole",
    "Near the Earth its field is close to a dipole's, which a magnetized sphere at the "
    "centre gives exactly.")
EQUATOR_ONLY = omitted(
    "Bouncing between the poles",
    "The particles stay in the equatorial plane. Real trapped particles also spiral along "
    "the field lines and bounce between mirror points in the north and the south.")
TEST_PARTICLES = simpler(
    "Test particles in given fields",
    "The particles move in fields that they do not change; a real magnetosphere is shaped "
    "by its plasma's own currents (the magnetopause, the ring current, the tail).")

REALITY["jackson_van_allen"] = R(
    "Drift around a dipole Earth: the Van Allen belts (Jackson Pr. 12.9)",
    kept("Protons and electrons drift opposite ways",
         "The field grows towards the Earth, so a gyrating particle's orbit is tighter on "
         "the near side and it drifts around the Earth, the two signs in opposite senses: "
         "the ring current."),
    DIPOLE_EARTH,
    simpler("Large gyrations",
            "Here a gyration is 0.4 of the distance to the Earth; a real belt particle "
            "gyrates over kilometres thousands of kilometres out, and drifts round in "
            "minutes to hours."),
    EQUATOR_ONLY)

REALITY["jackson_gradient_drift"] = R(
    "Gradient drift in a nonuniform magnetic field (Jackson §12.4)",
    kept("Drift along lines of equal |B|",
         "A gyrating orbit is tighter where the field is stronger, so the particle drifts "
         "sideways, the two signs in opposite directions, as in the textbook's guiding-centre "
         "picture."),
    UNIFORM,
    simpler("A row of magnets makes the gradient",
            "The textbook's gradient is gentle against a gyration; this one changes over a "
            "few gyrations, and the exact motion shows the drift and its wobble."),
    SPHERE_MAGNETS)

REALITY["jackson_runaway"] = R(
    "Crossed fields with |E| > c|B|: no frame to drift in (Jackson Pr. 12.5)",
    kept("No frame without an electric field",
         "When |E| > c|B| no moving frame cancels E (E² − c²B² is the same in every frame), "
         "and the particle runs away along E, exactly as relativity says."),
    SLOW_LIGHT,
    kept("Radiation negligible",
         "The particle's charge is small and every field large (only qE and qB matter), so "
         "that its radiation stays below the accuracy."),
    SPHERE_MAGNETS)

REALITY["throw_a_charge"] = R(
    "Momentum passed between two free charges by their fields",
    kept("Equal and opposite pushes",
         "Each charge's field pushes the other, and the momentum one loses the other gains, "
         "without contact."),
    simpler("Instantaneous forces",
            "With c = ∞ each feels the other's push at once and Newton's third law holds at "
            "every moment. Really the push travels at the speed of light and the field holds "
            "momentum on the way, which matters only near c."),
    simpler("A slingshot launch",
            "You give the free charge its speed and direction by hand."),
    PLANE)

REALITY["jackson_recoil"] = R(
    "Equal masses recoil at right angles (Jackson §13.1)",
    kept("90° for equal masses, any force",
         "Momentum and energy conservation alone send equal masses off at right angles after "
         "an elastic collision with one at rest; cloud-chamber photographs of alphas hitting "
         "helium nuclei show the same fork."),
    simpler("Newtonian",
            "At relativistic speeds the angle closes below 90°; here c = ∞."),
    simpler("A distant push, not a contact",
            "The two never touch: the Coulomb force does the whole collision."),
    PLANE)

REALITY["ramp_the_coil"] = R(
    "A coil whose current its supply ramps",
    kept("The field follows the current",
         "The field grows in proportion to the current, and a rising field induces an "
         "electric field around the axis (Faraday's law), both acting on the ion."),
    simpler("A prescribed ramp",
            "The current rises at the set rate whatever the coil's inductance; a real supply "
            "must push it against L dI/dt (Induction kick has that circuit)."),
    LOOP_COIL)

REALITY["jackson_faraday"] = R(
    "Faraday's law: an electric field induced by a rising current (Jackson §5.15)",
    kept("∮E·dl = −dΦ/dt",
         "The changing flux drives the particle around the axis with no charge nearby, and "
         "the growing field tightens its orbit (the flux through it is an adiabatic "
         "invariant, §12.5)."),
    simpler("Quasi-static",
            "The field follows the current at once (c = ∞ for the motion); a real coil is "
            "quasi-static too as long as it is much smaller than a wavelength."),
    LOOP_COIL)

REALITY["jackson_knock_on"] = R(
    "Energy transfer in a Coulomb collision (Jackson Pr. 13.1)",
    kept("T(b) = T_max/(1 + (b/b_min)²)",
         "The closer the pass, the harder the kick; close passes throw the target forward, "
         "distant ones sideways."),
    simpler("A mass ratio of 40",
            "A fast ion in matter knocks out electrons thousands of times lighter than "
            "itself: the knock-on electrons (δ rays) along particle tracks."),
    simpler("One target",
            "In matter the ion passes many electrons at every distance, and the sum of their "
            "kicks is its stopping power (Bethe's formula)."),
    PLANE)

REALITY["jackson_stormer"] = R(
    "Störmer's forbidden region: cosmic rays turned away by the Earth (Jackson §12.1)",
    kept("Canonical angular momentum is conserved",
         "The dipole field is symmetric about its axis, so r p_φ + q r A_φ cannot change and "
         "a slow particle cannot reach the equator, whatever its aim."),
    DIPOLE_EARTH,
    simpler("Only the equatorial plane",
            "In 3D the forbidden region is a ring around the equator: cosmic rays reach the "
            "Earth more easily near the poles, why their cutoff depends on latitude."),
    FIXED)

REALITY["jackson_magnetosphere"] = R(
    "Solar-wind particles drifting into a dipole magnetosphere",
    kept("E×B and gradient drifts together",
         "The crossed fields carry both signs the same way; near the Earth the gradient "
         "drift splits them, protons one way, electrons the other."),
    simpler("The solar wind as uniform crossed fields",
            "In the Earth's frame the flowing solar wind carries an electric field −v×B; here "
            "it is uniform everywhere, even next to the Earth, where the real magnetosphere "
            "shields it."),
    TEST_PARTICLES, EQUATOR_ONLY)

STRONG_PARTICLE = simpler(
    "A strongly charged particle",
    "The particle carries a charge comparable to the sphere's, so that the charge it induces "
    "on the metal pulls it noticeably; for a real ion that pull matters only within "
    "nanometres of a surface. Newtonian (c = ∞), so nothing radiates.")
IMAGES = kept(
    "Images, exactly",
    "The metal's induced charge is the exact image series plus a fitted surface charge, "
    "its boundary error below 1e-10: the conductor is an exact equipotential.")

REALITY["jackson_own_image"] = R(
    "A charge pulled by its own image in a grounded sphere (Jackson §2.2)",
    IMAGES,
    kept("The pull grows steeply near the surface",
         "The image −qR/d at R²/d pulls with q²Rd/(d² − R²)², exactly."),
    STRONG_PARTICLE)

REALITY["jackson_two_spheres"] = R(
    "Two conducting spheres: images of images (Jackson Pr. 2.6)",
    IMAGES,
    kept("The neutral sphere becomes a dipole",
         "Each sphere's charge induces an image in the other, and so on: the neutral one "
         "polarizes towards the charged one, and both answer your charges and the "
         "particle."),
    FIXED, PLANE)

REALITY["jackson_slot"] = R(
    "A grounded wall with an opening (Jackson §3.13)",
    kept("The field leaks through and fades",
         "Through the opening the screened charge's field leaks out and falls off within a "
         "few widths of it."),
    simpler("A finite wall with a slot",
            "Jackson's wall is an infinite plane with a round hole; this one is two plates "
            "8 long and 4 high with a slot between them, so the field also reaches round "
            "their ends and over their top."),
    PLATE_FIELDS, FIXED,
    electrodes=["wall, lower half (grounded)", "wall, upper half (grounded)"])

REALITY["jackson_golden_ratio"] = R(
    "Repulsion that turns to attraction near a charged sphere (Jackson Pr. 2.4)",
    kept("The force vanishes at 1.618 radii",
         "For a charge equal to the sphere's, the image wins inside the golden ratio of the "
         "radius from the centre: far away it repels, close by it attracts."),
    IMAGES, STRONG_PARTICLE)

REALITY["jackson_multipoles"] = R(
    "The multipole expansion: far away only the lowest moment counts (Jackson §4.1)",
    kept("Q/b, p/b², qd²/b³",
         "A beam passing at distance b is turned by a net charge as 1/b, by a dipole as 1/b², "
         "by a quadrupole as 1/b³, exactly as the expansion says."),
    COULOMB, FIXED)

REALITY["jackson_sphere_slalom"] = R(
    "Charged conducting spheres: repulsion far, image attraction near",
    kept("Each sphere repels, then captures",
         "Carrying the particle's own charge, a sphere repels it from afar and pulls it in "
         "inside 1.618 radii, and each polarizes the others."),
    IMAGES, STRONG_PARTICLE)

QUANTUM_ATOM = omitted(
    "Quantum mechanics",
    "Real atoms do not collapse: quantum mechanics keeps the electron in stationary states, "
    "and it radiates only when it jumps between them. This classical collapse is the failure "
    "that made quantum theory necessary.")
FAST_COLLAPSE = simpler(
    "A collapse within a few turns",
    "The atom is scaled so that it radiates away its orbit in a few turns (c = 2 cells per "
    "time unit); a classical hydrogen atom would circle about a hundred thousand times "
    "before collapsing, in about 10⁻¹¹ s. The nucleus is held fixed.")
LL_REACTION = kept(
    "Radiation reaction",
    "The electron radiates as Larmor's formula says and loses that energy through the "
    "Landau–Lifshitz force, exact in classical electrodynamics while the reaction is small "
    "against the Coulomb force.")

REALITY["jackson_classical_atom"] = R(
    "The classical atom's radiative collapse (Jackson Pr. 16.2)",
    LL_REACTION,
    kept("Close passes radiate most",
         "The radiated power grows as the acceleration squared, so an orbit that dips close "
         "to the nucleus loses its energy faster."),
    FAST_COLLAPSE, QUANTUM_ATOM)

REALITY["jackson_circularization"] = R(
    "Elliptic orbits round off as they radiate (Jackson Pr. 16.3)",
    LL_REACTION,
    kept("Ellipses circularize",
         "An elliptic orbit radiates most near the nucleus, which takes away more energy "
         "than angular momentum there: the orbit becomes rounder as it decays."),
    FAST_COLLAPSE, QUANTUM_ATOM)

REALITY["jackson_three_orbits"] = R(
    "Three classical orbits radiating away (Jackson Ch. 16)",
    LL_REACTION,
    kept("Shape decides the decay",
         "Small and eccentric orbits radiate fastest; a wide circle alone would take five "
         "times too long."),
    FAST_COLLAPSE, QUANTUM_ATOM)

THOMSON_ATOM = simpler(
    "Thomson's atom",
    "The positive charge is a uniform sphere the electron moves through (J. J. Thomson's "
    "model, 1904): inside, the pull grows linearly with the distance, a spring. Rutherford "
    "found the positive charge in a tiny nucleus instead, but a bound electron's response to "
    "weak fields is still described by such a spring (the Lorentz oscillator of optics). "
    "The sphere is rigid and fixed.")

REALITY["jackson_bound_charge"] = R(
    "An electron bound by a spring: Thomson's atom (Jackson §16.7)",
    kept("A harmonic pull inside a uniform charge",
         "Inside a uniformly charged sphere the field grows linearly from the centre, so "
         "the electron circles with one period at any radius; outside the pull is "
         "Coulomb's 1/r², no spring."),
    THOMSON_ATOM,
    omitted("Ionization energies",
            "A real atom's electron leaves when it gains its binding energy in a quantum "
            "jump; here it leaves when it climbs out of the classical well."))

REALITY["jackson_resonance"] = R(
    "Resonant driving of a bound electron (Jackson §16.8)",
    kept("Only the resonant frequency builds up",
         "A drive at ω₀ pushes the electron in step every cycle and its swing grows; at "
         "other frequencies the pushes cancel out, as in the absorption of light by atoms."),
    THOMSON_ATOM,
    omitted("Damping",
            "Without radiation reaction here the swing grows until the electron escapes; a "
            "real oscillator's damping (its radiation, collisions) limits it and gives the "
            "resonance its width."),
    POINT_ANTENNA)

REALITY["jackson_bound_knock"] = R(
    "A passing ion kicks a bound electron (Jackson Pr. 13.2)",
    kept("Quick passes kick, slow ones give back",
         "A pass shorter than the electron's period is a sharp kick; a slow one lets the "
         "electron follow and hand the energy back (an adiabatic pass)."),
    THOMSON_ATOM,
    simpler("A mass ratio of 40",
            "A real ion is thousands of times heavier than an electron."),
    PLANE)

REALITY["jackson_spectroscopy"] = R(
    "Spectroscopy: atoms told apart by their resonances",
    kept("Each atom answers its own frequency",
         "A drive at one atom's natural frequency builds up its electron's swing and barely "
         "moves the other's: how spectral lines identify atoms."),
    THOMSON_ATOM,
    omitted("Damping",
            "No radiation reaction here, so the resonances have no width: a real line's "
            "width comes from the atom's radiation and its collisions."),
    POINT_ANTENNA)

RADIATION_EXACT = kept(
    "The radiation is exact for the motion",
    "What the receiver collects is computed from the particle's computed path with the exact "
    "far-zone fields (Liénard; the spectrum from Jackson's radiation integral), and the "
    "particle feels its radiation's reaction (Landau–Lifshitz).")
IN_PLANE_RECEIVER = simpler(
    "Directions in the plane only",
    "The receiver averages over directions in the plane of the motion; the real radiation "
    "also goes out of the plane, where its pattern differs.")
MILD_GAMMA = simpler(
    "γ = 3 and slow light",
    "c is 2 cells per time unit and the particle's γ is 3, so that its cone (1/γ, about 19°) "
    "and its flashes fit the board. Real synchrotron sources run electrons at γ in the "
    "thousands, beamed into fractions of a milliradian, out to X-rays.")

REALITY["jackson_beaming"] = R(
    "Relativistic beaming: a fast charge radiates forward (Jackson §14.3)",
    kept("A cone of half-angle about 1/γ",
         "An accelerated relativistic charge sends almost all its radiation forward, into a "
         "narrow cone around its velocity: none while it flies straight."),
    RADIATION_EXACT, MILD_GAMMA, IN_PLANE_RECEIVER)

REALITY["jackson_critical_frequency"] = R(
    "The critical frequency of synchrotron radiation (Jackson §14.6)",
    kept("Tighter bends, higher frequencies",
         "A bend lights the receiver for a time of about R/(γ³c), so its spectrum reaches "
         "about (3/2)γ³c/R: why synchrotron light from tight bends of fast electrons reaches "
         "X-rays."),
    RADIATION_EXACT, MILD_GAMMA, IN_PLANE_RECEIVER)

REALITY["jackson_quiet_turn"] = R(
    "Radiation follows the velocity: turning quietly (Jackson §14.2)",
    kept("Gentle turns are quiet",
         "The power grows as the acceleration squared while the time spent heading for the "
         "receiver grows only with the radius: wide turns send less there."),
    RADIATION_EXACT, MILD_GAMMA, IN_PLANE_RECEIVER)

REALITY["jackson_thomson"] = R(
    "Thomson scattering of light by a fast charge: inverse Compton (Jackson §14.8)",
    kept("Doppler up to about 4γ²ω",
         "A charge flying into a light wave meets it faster and sends the scattered light "
         "forward squeezed: how laser light becomes X-rays in inverse Compton sources."),
    RADIATION_EXACT,
    omitted("The photon's recoil",
            "Classical scattering leaves out the Compton recoil, valid while the light's "
            "photon energy, seen from the charge, is far below its rest energy; at higher "
            "energies the recoil lowers the scattered frequency."),
    IN_PLANE_RECEIVER)

REALITY["jackson_cross_section"] = R(
    "Why the sky is blue: light scattered by a bound electron (Jackson §16.8)",
    kept("Rayleigh, resonance, Thomson",
         "The cross-section σ_T ω⁴/((ω₀² − ω²)² + (ωΓ)²) grows as ω⁴ far below the "
         "resonance (Rayleigh), peaks at it, and tends to a free electron's above it "
         "(Thomson), with Γ from the electron's own radiation."),
    THOMSON_ATOM,
    simpler("One atom",
            "The sky's colour is the scattering of countless air molecules, whose resonances "
            "lie in the ultraviolet, adding up incoherently; here one atom scatters."),
    IN_PLANE_RECEIVER)

REALITY["jackson_braking"] = R(
    "Bremsstrahlung: the flash of a sudden stop (Jackson §15.2)",
    kept("Zero straight ahead, peaked at cos θ = β",
         "A stopped charge radiates as sin²θ/(1 − β cos θ)², every frequency alike, as "
         "Jackson derives."),
    simpler("An instant stop",
            "Real bremsstrahlung comes from deflections by nuclei, over a short but finite "
            "time, so the spectrum is flat only below the inverse of that time; and in "
            "quantum theory no photon carries more than the particle's kinetic energy."),
    RADIATION_EXACT, IN_PLANE_RECEIVER)

REALITY["jackson_undulator"] = R(
    "An undulator: a wiggled beam's radiation adds up (Jackson §14.7)",
    kept("ω₁ = 2γ²ω_u/(1 + K²/2)",
         "Seen from straight ahead the wiggles' radiation adds up at one frequency, Doppler-"
         "squeezed, and the line grows with the number of periods."),
    simpler("Magnets are spheres",
            "A real undulator is two rows of magnet blocks with poles, periods of centimetres, "
            "for electrons with γ in the thousands, making X-ray lines; here a row of "
            "dipole spheres."),
    RADIATION_EXACT, MILD_GAMMA)

RING_FIELDS = simpler(
    "Uniform crossed fields",
    "The charges circle in uniform crossed fields given from outside; their own fields add "
    "at the receiver, and they act on each other.")

REALITY["jackson_in_step"] = R(
    "Coherent radiation: two charges in step (Jackson Pr. 14.23)",
    kept("Fields add, energies square",
         "Two charges circling in step send twice the field, so four times the energy: why "
         "bunched beams radiate coherently (free-electron lasers)."),
    RADIATION_EXACT, RING_FIELDS, IN_PLANE_RECEIVER)

REALITY["jackson_quiet_ring"] = R(
    "Two charges on one circle cancel their radiation (Jackson Pr. 14.23)",
    kept("A dipole that never turns",
         "Charges on opposite sides of one circle, circling together, have a dipole moment "
         "that stays zero: their radiation at ω₀ cancels."),
    kept("Relativistic velocity addition",
         "The partner must move at the drift added relativistically to its circling "
         "velocity; with Galilean addition the cancellation fails by 30 times."),
    RADIATION_EXACT, RING_FIELDS)

REALITY["jackson_ring_of_four"] = R(
    "A ring of charges: why a steady current does not radiate (Jackson Pr. 14.23)",
    kept("N charges radiate at multiples of Nω₀",
         "Equal charges spaced evenly round a circle cancel every harmonic below Nω₀; a ring "
         "of many, a steady current, radiates almost nothing, unlike one orbiting charge."),
    RADIATION_EXACT, RING_FIELDS)

# =======================================================================================
# Inside the components: the z-invariant world (tubes, metals, transistors).

TUBE_SPACE_CHARGE = kept(
    "Space charge from the electrons themselves",
    "Every electron in flight repels the others, and the charge it induces on the electrodes "
    "is solved at every step: the space-charge limit (Child and Langmuir) comes out of the "
    "simulation, not from a formula.")
TUBE_Z = simpler(
    "Infinitely long electrodes",
    "The tube is the same all along its length: its electrodes are long prisms, its electrons "
    "long lines of charge, its currents per unit length. A real tube is a few centimetres "
    "long and has end effects; the textbooks' planar and coaxial laws assume this symmetry "
    "too.")
TUBE_NEWTON = kept(
    "Newtonian electrons",
    "At tube voltages electrons reach a few per cent of c, where Newtonian motion is accurate "
    "to a fraction of a per cent; their magnetic forces on each other (of order v²/c²) are "
    "left out.")
COLD_EMISSION = simpler(
    "Electrons released at rest",
    "The cathode gives off as many electrons as its field pulls off (the space-charge limit), "
    "released at rest. A real hot cathode emits them with thermal energies of about 0.1 eV, "
    "which lets a little current flow against a small retarding field and moves the "
    "potential's minimum slightly out from its face.")
SECONDARIES = omitted(
    "Secondary electrons and gas",
    "An electron that hits an electrode stops there. Real electrodes also knock out secondary "
    "electrons (the tetrode's trouble), get hot, and the residual gas ionizes; all left out.")
LIGHT_IONS = simpler(
    "Light ions",
    "The ions are 16 times heavier per charge than the electrons, not 1836 times or more: in "
    "a static field the mass changes neither the paths nor the electrons' gain, only the time "
    "scale and the ion current (√(m_e/m_i) of the electrons' in a planar diode).")
TUBE_CIRCUITS = simpler(
    "Lumped circuits per unit length",
    "Supplies, resistors, capacitors and coils are ideal lumped parts, per unit length of the "
    "tube; they drive the electrodes and feel the current the electrons induce (Ramo).")

REALITY["vacuum_diode"] = R(
    "A planar vacuum diode: the Child–Langmuir law",
    TUBE_SPACE_CHARGE,
    kept("I ∝ V^(3/2)",
         "The current grows faster than the voltage: the space-charge law of every "
         "space-charge-limited diode."),
    TUBE_Z, COLD_EMISSION,
    omitted("The cathode's temperature limit",
            "This cathode gives whatever the field pulls; a real one's temperature caps the "
            "current (Richardson's law: the Electron gun level)."),
    electrodes=["cathode (hot, emitting)", "anode (your supply)"])

REALITY["sideways"] = R(
    "Hull's magnetic cutoff: the magnetron's principle",
    kept("Hull's critical field",
         "Above a field of about √(2V m/q)/d the electrons curve back before the anode and "
         "drift along the cathode, at E/B across both fields."),
    kept("Broken symmetry at the ends",
         "Below the cutoff a little current still reaches the anode at the gap's ends, where "
         "the translation symmetry behind Hull's theorem is broken."),
    TUBE_Z,
    simpler("A uniform field from the solenoid",
            "The solenoid round the tube makes an exactly uniform field along it; a real "
            "magnet's field is only nearly uniform."),
    electrodes=["cathode (hot)", "anode (at 30)", "collector (at 10)"])

REALITY["bipolar_diode"] = R(
    "A bipolar diode: ions cancelling the electrons' space charge (Langmuir, 1929)",
    kept("More electrons pass",
         "The ions' positive charge cancels part of the electrons' cloud, so more electrons "
         "pass at the same voltage (1.86 times in an infinite planar diode)."),
    LIGHT_IONS, TUBE_Z,
    omitted("The cathode's wear",
            "The ions that strike the cathode wear it away in a real tube (sputtering); here "
            "the goal only limits their current."),
    electrodes=["cathode (hot, emits electrons)", "anode (your supply; its coating emits ions)"])

REALITY["through_the_cloud"] = R(
    "An ion crossing a working diode's electron cloud",
    kept("The ion and the flow act on each other",
         "The cloud pulls the ion, and the ion's charge disturbs the flow it crosses; your "
         "plates act on both."),
    simpler("The ion is a line of charge",
            "A real ion is a point, its field falling as 1/r² instead of 1/r: it would disturb "
            "the flow far less, and be pulled less."),
    TUBE_SPACE_CHARGE, TUBE_Z,
    electrodes=["cathode (hot)", "anode (your supply)"])

REALITY["electron_gun"] = R(
    "An electron gun with a temperature-limited cathode",
    kept("Temperature-limited emission",
         "The cathode gives at most 0.2 per unit length of its face, as a real hot cathode "
         "gives at most its Richardson current: a stronger field pulls no more."),
    kept("A beam spread by its own charge",
         "The beam's electrons repel one another and it spreads; two plates at a lower "
         "potential make a lens that refocuses it."),
    simpler("A sheet beam through a slit",
            "Infinitely long electrodes make a flat ribbon of a beam; a real gun makes a round "
            "beam through a round aperture."),
    TUBE_Z,
    electrodes=["cathode (hot, temperature-limited)", "anode, upper jaw of the slit (at 30)",
                "anode, lower jaw of the slit (at 30)", "collector (at 30)"])

REALITY["ion_diode"] = R(
    "A magnetically insulated ion diode, as in pulsed-power ion sources",
    kept("Electrons insulated, ions not",
         "Hull's cutoff for the heavy ions is 18 times the electrons', so a field that holds "
         "the electrons in a sheath lets the ions cross, and the sheath's charge draws more "
         "of them."),
    kept("The sheath spreads",
         "The electron sheath keeps spreading towards the anode, as in real ion diodes (their "
         "impedance collapse): the goal's window is a fixed time after switch-on."),
    simpler("A steady supply",
            "Real magnetically insulated diodes run on pulsed power, megavolts for tens of "
            "nanoseconds; here a steady supply."),
    LIGHT_IONS, TUBE_Z,
    electrodes=["cathode (ring; emits electrons inwards)",
                "anode (rod, your supply; emits ions)", "envelope (grounded)"])

REALITY["rectifier"] = R(
    "A full-wave rectifier tube with a smoothing capacitor",
    kept("It conducts one way",
         "Electrons flow to an anode only while it is above the cathode, so two anodes on "
         "opposite phases charge the output on both halves of the cycle: full-wave "
         "rectification, as in the twin-anode rectifiers of old power supplies."),
    kept("The ripple of an RC filter",
         "The capacitor holds the output between the peaks and the load drains it: the "
         "ripple falls as RC grows."),
    TUBE_CIRCUITS, TUBE_Z, COLD_EMISSION,
    electrodes=["cathode (hot; the output, through R and C to ground)",
                "left anode (AC supply)", "right anode (AC supply, half a period behind)"])

REALITY["triode"] = R(
    "A triode amplifier biased by its grid",
    kept("The grid controls the current",
         "Close to the cathode, where the electrons are slow, the grid's potential governs the "
         "current far more than the anode's: the amplification factor."),
    kept("Gain capped by the anode's capacitance",
         "The anode's own capacitance shunts its load, which is why the signal is slow and the "
         "anode short (as real tubes are limited at high frequency)."),
    simpler("Five wires as the grid",
            "A real grid is a fine helix or mesh round the cathode; here five parallel wires."),
    TUBE_CIRCUITS, TUBE_Z,
    electrodes=["cathode (hot)", "anode (through R = 60 from 100)",
                "grid wire (your bias and the signal)", "grid wire", "grid wire",
                "grid wire", "grid wire"])

REALITY["klystron"] = R(
    "A two-cavity klystron amplifier",
    kept("Velocity modulation and bunching",
         "The buncher speeds some electrons up and slows others; downstream the fast catch "
         "the slow, the bunches' current follows the Bessel function J₁, and too strong a "
         "drive overbunches."),
    kept("Phase pushing",
         "The output's phase moves with the beam's transit time, about 7.5° per volt here, as "
         "the theory's θ₀/(2V₀)."),
    simpler("Plates and a tank, not cavities",
            "A real klystron's buncher and catcher are resonant cavities with narrow gaps the "
            "beam crosses; here plates, the catcher on a lumped tuned circuit."),
    simpler("A weak beam",
            "Its space charge is negligible and the bunching purely ballistic; in high-power "
            "klystrons the bunches' own repulsion (plasma oscillations) limits it."),
    TUBE_Z,
    electrodes=["cathode (your supply: −V₀)", "Wehnelt plate, upper (at the cathode's potential)",
                "Wehnelt plate, lower (at the cathode's potential)",
                "drift tube mouth, upper (grounded)", "drift tube mouth, lower (grounded)",
                "buncher, upper (your sine)", "buncher, lower (your sine)",
                "drift tube, upper (grounded)", "drift tube, lower (grounded)",
                "catcher, upper (grounded)", "catcher (the tuned circuit)",
                "drift tube end, upper (grounded)", "drift tube end, lower (grounded)",
                "collector (grounded)"])

REALITY["magnetron"] = R(
    "A strapped eight-vane magnetron, as in microwave ovens and radar",
    kept("Hull's cutoff and the Buneman–Hartree condition",
         "The electrons drift at E/B; it oscillates where that drift keeps pace with the wave "
         "running round the vanes, the band the Buneman–Hartree condition predicts."),
    kept("Spokes and start from noise",
         "Electrons in step with the π mode give it energy and gather into four turning "
         "spokes; the oscillation grows out of the noise, as in a real magnetron."),
    simpler("Vanes on a lumped tank",
            "A real magnetron's resonators are cavities cut into the anode block between the "
            "vanes; here the vanes are strapped onto one lumped tuned circuit."),
    simpler("The radiation shown, not felt",
            "The antenna's microwaves are drawn from the run's currents with light at 40 cells "
            "per time unit (one-way: the electrons stay Newtonian); real magnetrons feed a "
            "waveguide at a few gigahertz."),
    TUBE_Z,
    electrodes=["cathode (rod, your supply)", "vane at 0° (grounded set)",
                "vane at 45° (tuned set)", "vane at 90° (grounded set)",
                "vane at 135° (tuned set)", "vane at 180° (grounded set)",
                "vane at 225° (tuned set)", "vane at 270° (grounded set)",
                "vane at 315° (tuned set)", "anode block (grounded)",
                "antenna (strapped to the tuned set)"])

# The metals (the Drude model).

METAL_EMERGES = kept(
    "Ohm's law emerges",
    "The current, the Joule heat and the Hall voltage are not formulas: they come out of the "
    "carriers' motion and collisions, σ = nq²τ/m verified against the Boltzmann equation.")
METAL_CLASSICAL = simpler(
    "A classical metal, a semiconductor by its ratios",
    "The carriers are classical particles at the lattice's temperature bouncing off fixed ion "
    "cores (the Drude model). A real metal's electrons are a degenerate quantum gas whose own "
    "ratios (plasma oscillations per collision, screening per free path) the model cannot "
    "afford; its carriers have a lightly doped semiconductor's ratios (the Scales section "
    "names which), which conducts by the same laws.")
METAL_Z = simpler(
    "A conductor seen end-on",
    "The strip is infinitely long along z, its carriers lines of charge moving as one: "
    "currents are per unit length, and the noise is larger than a real conductor's, whose "
    "carriers move independently along its length.")
METAL_SCATTERING = simpler(
    "Hard-disc collisions, a fixed-temperature lattice",
    "Carriers bounce off randomly placed hard discs and leave with a fresh thermal speed; the "
    "lattice takes the heat and never warms. Real carriers scatter off lattice vibrations and "
    "impurities at rates that depend on energy and temperature.")
METAL_NO_SELF_FIELD = omitted(
    "The currents' own magnetic field",
    "The carriers' currents make no magnetic field of their own: no self-inductance and no "
    "skin effect, valid while the skin depth is far larger than the conductor. The only "
    "magnetic field is the solenoid's.")
CONTACTS = simpler(
    "Ideal ohmic contacts",
    "A contact takes in the carriers that reach it and gives off as many as keep the layer in "
    "front of it neutral, at the lattice's temperature: a perfect reservoir. Real contacts "
    "have a resistance, and between different materials a barrier.")
VOLTMETER = simpler(
    "Ideal voltmeters",
    "The floating side contacts read the potential as an ideal voltmeter does, through a "
    "capacitance of 1 to ground, drawing no current.")
IDLE_ELECTRODE = "inert electrode (floating, emits nothing; the model needs one)"

REALITY["ohms_law"] = R(
    "Ohm's law inside a conductor (the Drude model)",
    METAL_EMERGES,
    kept("The Joule heat",
         "At every collision the carriers hand the lattice the energy the field gave them: "
         "the heat is the field's work on the current."),
    METAL_CLASSICAL, METAL_Z, CONTACTS,
    electrodes=["left contact (grounded)", "right contact (your supply)"],
    metals=["resistor"])

REALITY["long_and_thin"] = R(
    "Resistance from shape: R = L/(σW)",
    kept("Length raises, width lowers",
         "Two strips of one metal: the long thin one resists about four times as much (4.3: "
         "its stronger field warms its carriers a little)."),
    METAL_EMERGES, METAL_CLASSICAL,
    simpler("An ideal wire",
            "The junction's two contacts are joined by a wire without resistance."),
    electrodes=["short wide strip, right end (grounded)", "long strip, left end (supply 1)",
                "long strip, right end: the junction (supply 2)",
                "short wide strip, left end: the junction (wired)"],
    metals=["long thin resistor", "short wide resistor"])

REALITY["potentiometer"] = R(
    "A potentiometer: a voltage measured without drawing current",
    kept("The potential falls evenly",
         "Inside a conductor carrying a steady current the field is uniform, so the potential "
         "falls evenly along it; the contacts take a little."),
    kept("The null method",
         "A contact at the strip's own potential there draws no current: a voltage measured "
         "without disturbing it."),
    METAL_CLASSICAL, CONTACTS,
    electrodes=["strip's left end (grounded)", "strip's right end (supply 1)",
                "side contact (supply 2)"],
    metals=["resistive strip"])

REALITY["hall_probe"] = R(
    "A Hall probe: the voltage across a current in a magnetic field",
    kept("V_H ≈ IB/(nq)",
         "The field pushes the moving carriers sideways until the charge they pile up holds "
         "the rest back; at strong fields the voltage grows ever more slowly (the Hall factor "
         "falls towards 1, and the end contacts short it)."),
    kept("A semiconductor, as real probes are",
         "By its ratios the strip is n-type germanium (the Scales section); real Hall probes "
         "are semiconductors for the reason Few carriers or many shows."),
    VOLTMETER, METAL_Z,
    electrodes=["strip's left end (grounded)", "strip's right end (your supply)",
                "Hall contact, top (floating)", "Hall contact, bottom (floating)"],
    metals=["Hall plate"])

REALITY["positive_carriers"] = R(
    "The Hall effect of positive carriers (holes)",
    kept("The Hall sign follows the carriers'",
         "For the same current, carriers of opposite sign move opposite ways, are pushed to "
         "the same side, and give Hall voltages of opposite sign."),
    simpler("Positive classical carriers",
            "In zinc, beryllium and p-type semiconductors the positive sign comes from holes, "
            "empty states in a nearly full band that move as positive charges; here the "
            "carriers simply carry a positive charge over a negative background."),
    VOLTMETER, METAL_CLASSICAL,
    electrodes=["electron strip, left end (grounded)", "electron strip, right end (supply 1)",
                "electron strip, top Hall contact", "electron strip, bottom Hall contact",
                "positive strip, left end (grounded)", "positive strip, right end (supply 2)",
                "positive strip, top Hall contact", "positive strip, bottom Hall contact"],
    metals=["electron metal (negative carriers)", "positive-carrier metal"])

REALITY["few_or_many"] = R(
    "Why Hall probes are semiconductors: few fast carriers or many slow ones",
    kept("Fewer carriers, a larger Hall voltage",
         "The same conductance σ = nq²τ/m from few carriers flying far or many colliding "
         "often; at the same current the few drift faster, are pushed harder, and give four "
         "times the Hall voltage."),
    simpler("Ratios, not materials",
            "The strips have germanium-like and silicon-like ratios; copper's carriers, a "
            "million times as many, lie beyond what the model can afford."),
    VOLTMETER, METAL_CLASSICAL,
    electrodes=["few-carrier strip, left end (grounded)",
                "few-carrier strip, right end (supply 1)",
                "few-carrier strip, top Hall contact", "few-carrier strip, bottom Hall contact",
                "many-carrier strip, left end (grounded)",
                "many-carrier strip, right end (supply 2)",
                "many-carrier strip, top Hall contact",
                "many-carrier strip, bottom Hall contact"],
    metals=["few fast carriers (germanium-like)", "many slow carriers (silicon-like)"])

REALITY["magnetoresistance"] = R(
    "Geometric magnetoresistance",
    kept("Shorted Hall field, longer paths",
         "In a long strip the Hall field holds the carriers on course; a short wide strip's "
         "contacts short it out, the carriers' paths bend and lengthen, and it resists more "
         "whichever way the field points: the principle of InSb magnetoresistors."),
    kept("The model's own magnetoresistance",
         "Even the long strip resists a little more in a field: with τ = ℓ/v the faster "
         "carriers are deflected differently from the slow ones."),
    METAL_CLASSICAL, METAL_SCATTERING,
    electrodes=["long strip, left end (grounded)", "long strip, right end (your supply)",
                "wide strip, left contact (grounded)", "wide strip, right contact (your supply)"],
    metals=["long strip", "short wide strip"])

REALITY["strength_and_sign"] = R(
    "A two-part magnetic field sensor: strength and sign",
    kept("Even and odd in the field",
         "A magnetoresistor reads |B| whatever its sign; a Hall bar's voltage is odd in B and "
         "reads its sign, as in real field sensors."),
    VOLTMETER, METAL_CLASSICAL, CONTACTS,
    electrodes=["wide strip, bottom contact (grounded)", "Hall bar, left end (your supply)",
                "Hall bar, top Hall contact", "Hall bar, bottom Hall contact",
                "Hall bar, right end (the junction)", "wide strip, top contact (the junction)"],
    metals=["Hall bar (positive carriers)", "magnetoresistor (electrons)"])

REALITY["eddy_currents"] = R(
    "Eddy currents heating a conducting plate",
    kept("Faraday's law, then Ohm's",
         "The changing field induces a circulating electric field, which drives the carriers "
         "round in loops; the heat grows as Ḃ², as the torsion-function law gives."),
    simpler("Carriers that turn between collisions",
            "Here the field turns the carriers' paths by a good part of a radian between "
            "collisions (ω_c τ ≈ 0.4): a semiconductor's regime; in copper it barely turns "
            "them (Eddy currents in silicon)."),
    METAL_NO_SELF_FIELD, METAL_CLASSICAL,
    electrodes=[IDLE_ELECTRODE],
    metals=["plate"])

REALITY["eddy_silicon"] = R(
    "Eddy currents in Ohm's regime, as in copper",
    kept("Ohm's law alone",
         "The field barely turns the carriers between collisions (ω_c τ ≈ 0.05), the swing is "
         "slow against them, and the heat grows as B₀², as in copper at any laboratory "
         "field."),
    simpler("Silicon's ratios, not copper's",
            "The plate's own ratios are lightly doped silicon's; copper's (400 plasma "
            "oscillations per collision) are beyond the model, and only its regime is "
            "reproduced."),
    METAL_NO_SELF_FIELD, METAL_SCATTERING,
    electrodes=[IDLE_ELECTRODE],
    metals=["plate (silicon-like)"])

REALITY["laminations"] = R(
    "A laminated transformer core",
    kept("Narrow pieces, small loops",
         "Cut into strips, the core's eddy loops close within each strip and the heat falls "
         "as the strips narrow (N laminations of one metal take about 1/N² of it)."),
    simpler("Cuts are ideal insulating slits",
            "Real laminations are sheets coated with varnish; a cut here is a wall whose faces, "
            "a screening length apart, still carry a little current across as a capacitor."),
    omitted("A magnetic core",
            "A real core is iron, whose permeability concentrates the flux and whose "
            "hysteresis wastes energy too; this one is a non-magnetic conductor in a given "
            "field."),
    METAL_NO_SELF_FIELD,
    electrodes=[IDLE_ELECTRODE],
    metals=["transformer core"])

# The transistors: a gated channel of the Drude metal.

FET_EMERGES = kept(
    "No transistor law is programmed",
    "Depletion, accumulation, pinch-off, the current's saturation at high drain voltage, the "
    "current over the gate's barrier below threshold and the channel's charging all come out "
    "of the carriers' motion in the fields.")
FET_NO_HOLES = omitted(
    "Holes, and so inversion",
    "Only electrons: no holes, no generation or recombination. A MOSFET's inversion layer and "
    "the pn junctions at its source and drain cannot form, so this is a uniformly doped, "
    "normally-on channel that the gate depletes: a junctionless (depletion-mode) transistor, "
    "not an enhancement-mode MOSFET.")
FET_NO_OXIDE = omitted(
    "The oxide",
    "The gate faces the channel across empty space, and the permittivity is 1 everywhere "
    "(silicon's is 11.7, its oxide's 3.9). No work-function difference (gate 0 is flat band), "
    "no oxide or interface charge, no traps, no gate leakage.")
FET_SCATTERING = simpler(
    "Simple collisions, a lattice that never warms",
    "The carriers bounce off random hard discs, not phonons and impurities: there is no "
    "saturation velocity (at high fields the drift keeps growing, about as √E), and the "
    "device never heats itself.")
FET_DOPANTS = simpler(
    "Dopants as a smooth background",
    "The donors are a uniform positive background, always fully ionized: they scatter "
    "nothing, and there are no random dopant fluctuations.")
FET_QUANTUM = omitted(
    "Quantum effects",
    "No band structure, no quantum confinement in the thin channel, no tunnelling. The "
    "classical Boltzmann statistics is right at this light doping.")
FET_SLICE = simpler(
    "A slice of an infinitely wide device",
    "The device is infinitely wide along z (the usual cut in device simulation), but each "
    "carrier is a line of charge moving as one: its noise and grainy field are far larger "
    "than a real device's, why the current deep below threshold is not converged.")
FET_SCALE = simpler(
    "Voltages in thermal voltages",
    "kT/q is 4 here, so the gate at −18 is 4.5 kT/q, about 0.12 V at room temperature; a "
    "real gate swings about 1 V, 40 kT/q. The channel matches n-type silicon at about 10¹⁶ "
    "per cm³ by its own ratios.")
FET_CONTACTS = simpler(
    "Ideal contacts, no n+ regions",
    "The source and drain are perfect reservoirs of the channel's own doping; real ones are "
    "heavily doped regions with a resistance of their own.")

REALITY["gate"] = R(
    "A field-effect transistor: a junctionless, depletion-mode FET",
    kept("The field effect",
         "A gate beside a channel a few screening lengths thin empties it (depletion) or fills "
         "it from the contacts (accumulation), and the current follows: what a metal, "
         "screening within an atom, can never do."),
    FET_EMERGES, FET_NO_HOLES, FET_NO_OXIDE, FET_SCATTERING, FET_DOPANTS, FET_QUANTUM,
    FET_SLICE, FET_SCALE, FET_CONTACTS,
    electrodes=["source (grounded)", "drain (held at 10)",
                "gate, above the channel (your supply)", "gate, below the channel (your supply)"],
    metals=["channel (lightly doped, n-type)"])

REALITY["inverter"] = R(
    "A resistor-load inverter: the simplest logic gate",
    kept("Logic levels from a load line",
         "Off, the channel passes almost nothing and the drain sits near the supply (a 1); on, "
         "its current through the load pulls the drain down (a 0)."),
    kept("Gate–drain feedthrough",
         "At each switch the gate's step jumps the drain for a moment through the capacitance "
         "between them, as in real logic."),
    simpler("Off needs a negative gate",
            "A depletion transistor conducts at 0, so the gate must go to −60 to switch it off; "
            "real logic uses enhancement-mode transistors, off at 0, in complementary pairs "
            "(CMOS) that draw no current when idle."),
    FET_NO_HOLES, FET_NO_OXIDE, FET_SLICE,
    electrodes=["source (grounded)", "drain (fed from 25 through the load)",
                "gate, above the channel (the pulse; its top your supply)",
                "gate, below the channel (the pulse)"],
    metals=["channel (lightly doped, n-type)"])

REALITY["amplifier"] = R(
    "A common-source amplifier with a resistor load",
    kept("Voltage gain from a load line",
         "The gate's small swing changes the channel's current, and the load turns that "
         "change into a larger swing of the drain: g_m times the load and the channel's own "
         "resistance in parallel."),
    kept("The bias point",
         "Too negative, the pinched channel's current barely follows the gate; too positive, "
         "the drain falls so low that the channel acts as a plain resistor: the gain peaks "
         "in between, as in every transistor stage."),
    kept("Inverted, and lagging",
         "The drain falls as the gate rises. It lags because its capacitance charges through "
         "the load, so a faster signal would be amplified less: the gain–bandwidth product, "
         "here g_m/C about 0.06."),
    simpler("A gain under 2",
            "A real stage gains 10 to 100: a long channel in saturation hardly feels its "
            "drain. This short, thick one does (its output resistance about 400 against the "
            "load's 200), and the drain's capacitance, slow to charge, takes some of the rest."),
    FET_NO_HOLES, FET_NO_OXIDE, FET_SLICE, TUBE_CIRCUITS,
    electrodes=["source (grounded)", "drain, the output (fed from 120 through the load)",
                "gate, above the channel (the signal around your bias)",
                "gate, below the channel (the signal)"],
    metals=["channel (lightly doped, n-type)"])
