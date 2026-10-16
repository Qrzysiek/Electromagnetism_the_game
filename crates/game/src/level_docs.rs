//! The level documents (docs/levels/README.md): each shipped level's physics derived from
//! first principles, typeset by `scripts/level_docs.py build` into `docs/levels/pdf/`
//! (generated, not committed). The level panel opens them with the system's PDF viewer.

use std::path::{Path, PathBuf};

/// The Foundations every level document builds on.
pub const FOUNDATIONS: &str = "foundations";

/// A level's document: typeset, written but not typeset, or none yet.
pub enum Doc {
    Ready(PathBuf),
    NotBuilt(String),
    Missing,
}

/// The level's name in `levels/curriculum.json`: its file name without the number
/// (`01_first_bend.json` -> `first_bend`); documents are filed under it, so renumbering the
/// levels does not move them.
pub fn slug(level_path: &Path) -> Option<String> {
    let stem = level_path.file_stem()?.to_str()?;
    let (number, rest) = stem.split_once('_')?;
    number
        .chars()
        .all(|c| c.is_ascii_digit())
        .then(|| rest.to_owned())
}

/// `docs/levels/`, next to the levels' folder.
fn docs_dir() -> Option<PathBuf> {
    Some(crate::levels_dir().parent()?.join("docs").join("levels"))
}

/// The document of `slug`.
pub fn find(slug: &str) -> Doc {
    let Some(dir) = docs_dir() else {
        return Doc::Missing;
    };
    let pdf = dir.join("pdf").join(format!("{slug}.pdf"));
    if pdf.is_file() {
        Doc::Ready(pdf)
    } else if dir.join(slug).join("doc.tex").is_file() {
        Doc::NotBuilt(format!("python scripts/level_docs.py build {slug}"))
    } else {
        Doc::Missing
    }
}

/// The document of a level file (none for custom levels: `levels/custom/`).
pub fn for_level(level_path: &Path) -> Doc {
    let shipped = level_path
        .parent()
        .and_then(Path::file_name)
        .is_some_and(|d| d == "levels");
    match slug(level_path) {
        Some(s) if shipped => find(&s),
        _ => Doc::Missing,
    }
}

/// Opens a PDF with the system's default viewer (its links to the other documents work there).
pub fn open(pdf: &Path) {
    if let Err(e) = opener::open(pdf) {
        bevy::log::warn!("could not open {}: {e}", pdf.display());
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn slug_drops_the_number() {
        assert_eq!(
            slug(Path::new("levels/01_first_bend.json")).as_deref(),
            Some("first_bend")
        );
        assert_eq!(
            slug(Path::new("levels/111_potentiometer.json")).as_deref(),
            Some("potentiometer")
        );
        assert_eq!(slug(Path::new("levels/custom/my_level.json")), None);
    }
}
