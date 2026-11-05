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
- `timing/`: the tube levels' timing audit (docs/HANDOFF.md, the owner's timing rules):
  `timing_audit.py` (`generator check` with a fresh run cache: preview and verification
  times), `field_audit.py` (`generator tube-field`: the field maps' rate).
