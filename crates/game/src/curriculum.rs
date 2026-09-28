//! The curriculum: each shipped level's arc and tier (introduction, intermediate, master),
//! from `levels/curriculum.json`, which `scripts/levels.py` writes next to the levels
//! (docs/CURRICULUM.md). Custom levels have no place in it.

use std::path::{Path, PathBuf};

/// Where a level stands in the curriculum.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Place {
    /// Arc number, from 1.
    pub arc: usize,
    pub arc_name: String,
    pub tier: String,
}

impl Place {
    /// "Arc 1 · Charges: steering and optics".
    pub fn arc_title(&self) -> String {
        format!("Arc {} · {}", self.arc, self.arc_name)
    }
}

/// The place of each level file (`NN_<slug>.json`), looked up in the `curriculum.json`
/// of its directory.
pub fn places(paths: &[PathBuf]) -> Vec<Option<Place>> {
    let mut read: Vec<(PathBuf, Option<serde_json::Value>)> = Vec::new();
    paths
        .iter()
        .map(|p| {
            let dir = p.parent()?.to_path_buf();
            if !read.iter().any(|(d, _)| *d == dir) {
                let json = std::fs::read_to_string(dir.join("curriculum.json"))
                    .ok()
                    .and_then(|t| serde_json::from_str(&t).ok());
                read.push((dir.clone(), json));
            }
            let json = read.iter().find(|(d, _)| *d == dir)?.1.as_ref()?;
            place_of(json, &slug(p)?)
        })
        .collect()
}

/// `first_bend` from `.../01_first_bend.json`.
fn slug(path: &Path) -> Option<String> {
    let stem = path.file_stem()?.to_str()?;
    let (number, rest) = stem.split_once('_')?;
    number
        .chars()
        .all(|c| c.is_ascii_digit())
        .then(|| rest.to_string())
}

fn place_of(json: &serde_json::Value, slug: &str) -> Option<Place> {
    for (a, arc) in json.get("arcs")?.as_array()?.iter().enumerate() {
        for tier in arc.get("tiers")?.as_array()? {
            let listed = tier.get("levels")?.as_array()?;
            if listed.iter().any(|l| l.as_str() == Some(slug)) {
                return Some(Place {
                    arc: a + 1,
                    arc_name: arc.get("name")?.as_str()?.to_string(),
                    tier: tier.get("name")?.as_str()?.to_string(),
                });
            }
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Every shipped level has a place, and the files are numbered in curriculum order:
    /// arcs ascend, and within an arc the tiers come in their listed order.
    #[test]
    fn every_shipped_level_has_a_place_in_order() {
        let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../levels");
        let mut paths: Vec<PathBuf> = std::fs::read_dir(&dir)
            .unwrap()
            .filter_map(|e| e.ok().map(|e| e.path()))
            .filter(|p| level::is_level_file(p))
            .collect();
        paths.sort();
        let places = places(&paths);
        let mut last = (0, String::new());
        let mut seen_tiers: Vec<(usize, String)> = Vec::new();
        for (p, place) in paths.iter().zip(&places) {
            let place = place
                .as_ref()
                .unwrap_or_else(|| panic!("{} is not in curriculum.json", p.display()));
            assert!(place.arc >= last.0, "{}: arcs out of order", p.display());
            let key = (place.arc, place.tier.clone());
            if key != last {
                assert!(
                    !seen_tiers.contains(&key),
                    "{}: tier {} of arc {} is split",
                    p.display(),
                    place.tier,
                    place.arc
                );
                seen_tiers.push(key.clone());
            }
            last = key;
        }
    }
}
