use super::*;

/// Every shipped (and custom) level can be rebuilt by hand in the sandbox: all its
/// values are within reach of the editors. New fields are covered at compile time by
/// the exhaustive destructuring above.
#[test]
fn all_levels_are_reproducible_in_the_sandbox() {
    let dir = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../levels");
    let mut checked = 0;
    for d in [dir.clone(), dir.join("custom")] {
        let Ok(entries) = std::fs::read_dir(&d) else {
            continue;
        };
        for e in entries {
            let p = e.unwrap().path();
            if !level::is_level_file(&p) {
                continue;
            }
            let level = Level::from_json(&std::fs::read_to_string(&p).unwrap()).unwrap();
            assert_eq!(check_editable(&level), Ok(()), "{}", p.display());
            checked += 1;
        }
    }
    assert!(checked > 0);
}

/// A ring's hole (`MetalSpec::hole`) is editable: its corners on the grid strictly inside the
/// metal pass, a corner on the metal's edge or off the grid does not.
#[test]
fn ring_holes_are_editable() {
    // By its slug: the curriculum renumbers the files.
    let base = level::shipped("eddy_currents");
    let with = |hole: [level::Node; 2]| {
        let mut l = base.clone();
        l.tube.as_mut().expect("a tube level").metals[0].hole = Some(hole);
        check_editable(&l)
    };
    assert_eq!(with([[10, 8, 0], [14, 12, 0]]), Ok(()));
    assert!(
        with([[6, 8, 0], [14, 12, 0]]).is_err(),
        "on the metal's edge"
    );
    assert!(with([[10, 8, 0], [40, 12, 0]]).is_err(), "off the grid");
}
