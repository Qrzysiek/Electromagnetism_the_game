# Prototypes

The scripts, scans and draft tests that the tracked documents cite as the record of how a
choice was made: level designs' `generator tube-scan` runs (`scripts/designs/components.py`),
the prototypes behind decisions (`docs/DECISIONS.md`), exact references and explorations
(`docs/HANDOFF.md`, `docs/DRUDE.md`, `docs/DARWIN.md`, `docs/OPEN_PROBLEMS.md`, `PHYSICS.md`).

They were written in a session's `target/scratch/` (not versioned, and gone with its
worktree) and moved here when a tracked file cited them (DEC-37; CLAUDE.md, *Disk*). They
are not maintained: each ran from `target/scratch/` against the engine and levels of its date
(paths are relative to it, the `.rs` drafts belong in a crate's `tests/`). To rerun one,
copy it back into `target/scratch/` and expect to adapt it.

- `inductance/`: the exact eddy currents of a square ring (`eddy2d.py`) and of nested rings
  (`transformer2d.py`), the self-field's scans (`NOTES.md`), Self-inductance's designs.
- `joule/`: the heat at rest's runs (PHYSICS.md §2.12); `low_voltage_work.rs`, the field's
  work against `I V` at a low supply (DEC-33's history; a test for `crates/level/tests/`);
  `contact_bias.rs`, the contacts' layers: the field's work at rest against the resolution's
  parts, its profile by columns (with `Lattice::account_size` temporarily ten cells a unit),
  and the macroscopic heat against `I V` (DEC-41, PHYSICS.md §2.11).
- `noise/`: the thermal-noise level's prototypes (docs/HANDOFF.md, *Next*; DRUDE.md Phase 2
  item 1): `runs.out` (each run's `generator tube-scan` result line, its trace left out:
  `n2_r1` and `n2_r2` the two-goal prototype at the preview's and the verification's
  resolution, `sq*` a square-wave kick), `proto.py` (one strip, a floating contact with a
  lumped capacitance, optionally a grounded plate beside it), `proto2.py` (two strips, a
  reference contact and the player's with a kicker: the spread-ratio and decay goals; it
  needs the measures of the same date, not yet committed), `drive.py` (the strip's
  resistance from a driven contact), `acf.py` and `fold.py` (a `tube-scan --trace`'s spread,
  correlation time and folded response to a square-wave kick). `EQUIPARTITION.md`
  (2026-11-15): a floating contact's noise against `kT/C`, 18–43 % high with the refill
  contacts, its two causes, and the thermal wall that fixed it (DEC-43); `equipartition.rs`
  (its runs, a test to copy into `crates/physics/tests/`) and `reservoir_contacts.patch`
  (the reservoirs it tried, regular, Poisson and Poisson with the release's Boltzmann
  factor, exact or to first order: an experiment on `emission.rs` from before
  `emit_thermal`, which `git apply` puts back on the commit that precedes it).
  `LEVEL_CHECKS.md` (2026-11-15): the level's two checks with the thermal-wall contacts
  (`proto2.py --thermal`; `runs.out` `t2_r1`, `t2_r2`, `cscan_r1`; `cscan.json` the plate
  placements): `RC` short enough, the plate too weak a capacitance control.
- `timing/`: the tube levels' timing audit (docs/HANDOFF.md, the owner's timing rules):
  `timing_audit.py` (`generator check` with a fresh run cache: preview and verification
  times), `field_audit.py` (`generator tube-field`: the field maps' rate).
