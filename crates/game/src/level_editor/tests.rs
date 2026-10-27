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
