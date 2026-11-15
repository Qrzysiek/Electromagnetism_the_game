# Mains hum: notes

## The reference is not an imaging lens (for the owner; not changed here)

The level's reality entry says it represents "an imaging lens", with the remark "Imaging cancels
the launch angle" as *kept* (and LEVELS.md's row says "imaging"). The reference does not image
the source onto the detector. On the line x = 28, without the hum (`doc.tex` §4.4, eq. 42.7):

- the arrival height's sensitivity to the start's transverse velocity is 12.5 (28 on the empty
  board): the charges halve the angle's effect; an imaging lens would make it 0;
- its sensitivity to the start's height is −0.006: the detector sits at the charges' focus for
  rays that start parallel, so the hum's shift of the start (phases 0° and 180°) hardly matters;
- tilted flights are kept on target by a fold, not by imaging: a flight aimed down passes close
  to the lower charge (22, 5) and is thrown back up; the arrival height is lowest, 8.795, at a
  transverse velocity of −0.045 and rises on both sides (10.62 at −0.1, passing 1.54 from the
  charge; 16.29 at −0.12, passing 1.20). Arrivals slightly below the detector still enter it
  through its bottom face as they climb.

The four listed phases arrive (the closest 0.02 inside the top edge), and so do all 72 phases
of a 5° sweep, the lowest exactly on the bottom face. Either the remark could describe what the
reference does (the start's height focused out, the downward angles folded back by the lower
charge), or the level could get a reference that images (∂y_e/∂θ = 0 at the detector), if one
exists within four charges in the region; this document does not search for one.

Reproduce: `py -3.14 scripts/level_docs.py calc mains_hum` (`dydv`, `dydy`, `foldV`, `foldY`,
`yAtMinus`, `closestMinus`, `margin`, `nSweepOk`).
