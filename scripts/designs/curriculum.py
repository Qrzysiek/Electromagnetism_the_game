"""The curriculum: the arcs in order (`ARCS`), every level (`LEVELS`), and `levels/curriculum.json`."""

import json
import os
from designs.common import ROOT
from designs.charges import (
    around_the_wall,
    beta_spectrometer,
    einzel_lens,
    fast_lane,
    first_bend,
    geiger_marsden,
    hemispherical_analyzer,
    reflectron,
    slingshot,
    thomson_crt,
    twin_beams,
)
from designs.fields import (
    build_wien_filter,
    calutron,
    dempster,
    earths_field,
    first_coil,
    first_magnet,
    mains_hum,
    rf_kick,
    rf_separator,
    stern_gerlach,
    stray_field,
    streak_camera,
    synchrotron_light,
    tune_the_rf,
    wien_filter,
)
from designs.matter import (
    betatron,
    build_a_deflector,
    charging_a_plate,
    chopper,
    collimator,
    deflection_plates,
    ferrite_shield,
    glass_screen,
    high_voltage_dome,
    image_charge,
    induction_kick,
    injection,
    polarized_sphere,
    power_supply,
    pulse_sorter,
    real_einzel_lens,
    ringing_plate,
    shielding,
    soft_landing,
    tune_the_lens,
)
from designs.instruments import (
    beam_pipe,
    beam_preparation,
    calutron_space_charge,
    chromatic_aberration,
    collimated_beam,
    crt_earth_field,
    isotope_separator,
    mass_spectrometer,
    microscope_column,
    real_analyzer,
    relativistic_beam,
    rf_beam_line,
    soft_landing_current,
    sorting_station,
    space_charge,
    stern_gerlach_beam,
    two_stages,
    velocity_selector,
)
from designs.jackson import (
    jackson_beaming,
    jackson_bound_charge,
    jackson_bound_knock,
    jackson_braking,
    jackson_circularization,
    jackson_classical_atom,
    jackson_critical_frequency,
    jackson_cross_section,
    jackson_exb_drift,
    jackson_faraday,
    jackson_golden_ratio,
    jackson_gradient_drift,
    jackson_in_step,
    jackson_knock_on,
    jackson_magnetosphere,
    jackson_multipoles,
    jackson_own_image,
    jackson_quiet_ring,
    jackson_quiet_turn,
    jackson_recoil,
    jackson_resonance,
    jackson_ring_of_four,
    jackson_runaway,
    jackson_slot,
    jackson_spectroscopy,
    jackson_sphere_slalom,
    jackson_stormer,
    jackson_thomson,
    jackson_three_orbits,
    jackson_two_spheres,
    jackson_undulator,
    jackson_van_allen,
    ramp_the_coil,
    throw_a_charge,
)
from designs.components import (
    amplifier,
    bipolar_diode,
    eddy_currents,
    eddy_silicon,
    electron_gun,
    few_or_many,
    field_tuned_damping,
    gate,
    hall_probe,
    inverter,
    ion_diode,
    klystron,
    laminations,
    lock_in_hall,
    long_and_thin,
    magnetoresistance,
    magnetron,
    ohms_law,
    positive_carriers,
    potentiometer,
    rectifier,
    ringing_tank,
    sideways,
    strength_and_sign,
    through_the_cloud,
    triode,
    vacuum_diode,
    wheatstone_bridge,
)

# Arcs: (name, [(tier, [(slug, function)])]). The level files are numbered in this order;
# `levels/curriculum.json` tells the game each level's arc and tier.
ARCS = [
    ("Charges: steering and optics", [
        ("Introduction", [
            ("first_bend", first_bend),
            ("slingshot", slingshot),
            ("geiger_marsden", geiger_marsden),
            ("twin_beams", twin_beams),
            ("two_stages", two_stages),
            ("injection", injection),
        ]),
        ("Intermediate", [
            ("thomson_crt", thomson_crt),
            ("around_the_wall", around_the_wall),
            ("einzel_lens", einzel_lens),
            ("reflectron", reflectron),
            ("collimator", collimator),
            ("hemispherical_analyzer", hemispherical_analyzer),
            ("soft_landing", soft_landing),
        ]),
        ("Master", [
            ("sorting_station", sorting_station),
        ]),
    ]),
    ("Metal and electrodes", [
        ("Introduction", [
            ("high_voltage_dome", high_voltage_dome),
            ("polarized_sphere", polarized_sphere),
            ("image_charge", image_charge),
            ("deflection_plates", deflection_plates),
            ("power_supply", power_supply),
            ("build_a_deflector", build_a_deflector),
            ("shielding", shielding),
            ("glass_screen", glass_screen),
        ]),
        ("Intermediate", [
            ("tune_the_lens", tune_the_lens),
            ("real_einzel_lens", real_einzel_lens),
            ("beam_pipe", beam_pipe),
        ]),
        ("Master", [
            ("microscope_column", microscope_column),
        ]),
    ]),
    ("Relativity and magnetism", [
        ("Introduction", [
            ("fast_lane", fast_lane),
            ("first_coil", first_coil),
            ("first_magnet", first_magnet),
            ("ferrite_shield", ferrite_shield),
            ("stern_gerlach", stern_gerlach),
        ]),
        ("Intermediate", [
            ("beta_spectrometer", beta_spectrometer),
            ("dempster", dempster),
            ("wien_filter", wien_filter),
            ("calutron", calutron),
            ("build_wien_filter", build_wien_filter),
        ]),
        ("Master", [
            ("mass_spectrometer", mass_spectrometer),
        ]),
    ]),
    ("Time: noise, radio frequency and radiation", [
        ("Introduction", [
            ("stray_field", stray_field),
            ("rf_kick", rf_kick),
            ("ramp_the_coil", ramp_the_coil),
            ("synchrotron_light", synchrotron_light),
        ]),
        ("Intermediate", [
            ("mains_hum", mains_hum),
            ("earths_field", earths_field),
            ("crt_earth_field", crt_earth_field),
            ("rf_separator", rf_separator),
            ("tune_the_rf", tune_the_rf),
            ("streak_camera", streak_camera),
        ]),
        ("Master", [
            ("rf_beam_line", rf_beam_line),
        ]),
    ]),
    ("Circuits", [
        ("Introduction", [
            ("charging_a_plate", charging_a_plate),
            ("chopper", chopper),
            ("ringing_plate", ringing_plate),
        ]),
        ("Intermediate", [
            ("induction_kick", induction_kick),
            ("betatron", betatron),
        ]),
        ("Master", [
            ("pulse_sorter", pulse_sorter),
        ]),
    ]),
    ("Beams", [
        ("Introduction", [
            ("space_charge", space_charge),
            ("stern_gerlach_beam", stern_gerlach_beam),
            ("relativistic_beam", relativistic_beam),
        ]),
        ("Intermediate", [
            ("collimated_beam", collimated_beam),
            ("velocity_selector", velocity_selector),
            ("beam_preparation", beam_preparation),
            ("chromatic_aberration", chromatic_aberration),
            ("real_analyzer", real_analyzer),
            ("calutron_space_charge", calutron_space_charge),
            ("soft_landing_current", soft_landing_current),
        ]),
        ("Master", [
            ("isotope_separator", isotope_separator),
        ]),
    ]),
    ("Jackson: charges in fields and collisions", [
        ("Introduction", [
            ("jackson_exb_drift", jackson_exb_drift),
            ("jackson_van_allen", jackson_van_allen),
            ("throw_a_charge", throw_a_charge),
            ("jackson_faraday", jackson_faraday),
        ]),
        ("Intermediate", [
            ("jackson_knock_on", jackson_knock_on),
            ("jackson_gradient_drift", jackson_gradient_drift),
            ("jackson_runaway", jackson_runaway),
            ("jackson_recoil", jackson_recoil),
            ("jackson_stormer", jackson_stormer),
        ]),
        ("Master", [
            ("jackson_magnetosphere", jackson_magnetosphere),
        ]),
    ]),
    ("Jackson: electrostatics", [
        ("Introduction", [
            ("jackson_own_image", jackson_own_image),
            ("jackson_two_spheres", jackson_two_spheres),
        ]),
        ("Intermediate", [
            ("jackson_slot", jackson_slot),
            ("jackson_golden_ratio", jackson_golden_ratio),
            ("jackson_multipoles", jackson_multipoles),
        ]),
        ("Master", [
            ("jackson_sphere_slalom", jackson_sphere_slalom),
        ]),
    ]),
    ("Jackson: radiation damping", [
        ("Introduction", [
            ("jackson_classical_atom", jackson_classical_atom),
        ]),
        ("Intermediate", [
            ("jackson_circularization", jackson_circularization),
        ]),
        ("Master", [
            ("jackson_three_orbits", jackson_three_orbits),
        ]),
    ]),
    ("Jackson: bound charges", [
        ("Introduction", [
            ("jackson_bound_charge", jackson_bound_charge),
            ("jackson_resonance", jackson_resonance),
        ]),
        ("Intermediate", [
            ("jackson_bound_knock", jackson_bound_knock),
        ]),
        ("Master", [
            ("jackson_spectroscopy", jackson_spectroscopy),
        ]),
    ]),
    ("Jackson: radiation", [
        ("Introduction", [
            ("jackson_beaming", jackson_beaming),
            ("jackson_critical_frequency", jackson_critical_frequency),
            ("jackson_in_step", jackson_in_step),
        ]),
        ("Intermediate", [
            ("jackson_quiet_turn", jackson_quiet_turn),
            ("jackson_thomson", jackson_thomson),
            ("jackson_cross_section", jackson_cross_section),
            ("jackson_braking", jackson_braking),
            ("jackson_quiet_ring", jackson_quiet_ring),
            ("jackson_ring_of_four", jackson_ring_of_four),
        ]),
        ("Master", [
            ("jackson_undulator", jackson_undulator),
        ]),
    ]),
]



ARCS.append(("Inside the components: vacuum tubes", [
    ("Introduction", [
        ("vacuum_diode", vacuum_diode),
        ("sideways", sideways),
        ("bipolar_diode", bipolar_diode),
    ]),
    ("Intermediate", [
        ("through_the_cloud", through_the_cloud),
        ("electron_gun", electron_gun),
        ("ion_diode", ion_diode),
        ("rectifier", rectifier),
        ("triode", triode),
        ("klystron", klystron),
    ]),
    ("Master", [
        ("magnetron", magnetron),
    ]),
]))


ARCS.append(("Inside the components: metals", [
    ("Introduction", [
        ("ohms_law", ohms_law),
        ("long_and_thin", long_and_thin),
        ("potentiometer", potentiometer),
    ]),
    ("Intermediate", [
        ("hall_probe", hall_probe),
        ("positive_carriers", positive_carriers),
        ("few_or_many", few_or_many),
        ("magnetoresistance", magnetoresistance),
        ("wheatstone_bridge", wheatstone_bridge),
        ("eddy_currents", eddy_currents),
        ("eddy_silicon", eddy_silicon),
        ("laminations", laminations),
    ]),
    ("Master", [
        ("strength_and_sign", strength_and_sign),
    ]),
]))


ARCS.append(("Inside the components: transistors", [
    ("Introduction", [
        ("gate", gate),
        ("inverter", inverter),
    ]),
    ("Intermediate", [
        ("amplifier", amplifier),
    ]),
]))


ARCS.append(("AC and induction", [
    ("Introduction", [
        ("ringing_tank", ringing_tank),
    ]),
    ("Intermediate", [
        ("field_tuned_damping", field_tuned_damping),
        ("lock_in_hall", lock_in_hall),
    ]),
]))


LEVELS = [lv for _, tiers in ARCS for _, lvls in tiers for lv in lvls]



def write_curriculum():
    """`levels/curriculum.json`: every level's arc and tier, by slug, in order."""
    arcs = [{"name": name, "tiers": [{"name": t, "levels": [slug for slug, _ in lvls]}
                                     for t, lvls in tiers]}
            for name, tiers in ARCS]
    path = os.path.join(ROOT, "levels", "curriculum.json")
    with open(path, "w", newline="\n", encoding="utf-8") as f:
        f.write(json.dumps({"arcs": arcs}, indent=2, ensure_ascii=False) + "\n")
