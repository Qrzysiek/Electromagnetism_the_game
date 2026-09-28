//! Math in level texts: `$...$` is LaTeX. The formulas are rendered once, at build time,
//! by MathJax (`scripts/math/render.mjs`, run by `scripts/levels.py`) to SVGs in
//! `levels/math/`, with a manifest of their sizes; here they are placed inline with the
//! text, coloured like it. A formula without an SVG (e.g. typed into a sandbox level
//! that has not been rebuilt) is shown as its LaTeX source.

use std::collections::HashMap;
use std::sync::{Mutex, OnceLock};

use bevy_egui::egui;

/// A rendered formula: its SVG file and its size in ex (the height of an x).
struct Formula {
    file: String,
    width_ex: f32,
    height_ex: f32,
}

/// The formulas rendered at build time, by their LaTeX source.
fn manifest() -> &'static HashMap<String, Formula> {
    static MANIFEST: OnceLock<HashMap<String, Formula>> = OnceLock::new();
    MANIFEST.get_or_init(|| {
        let path = crate::levels_dir().join("math").join("manifest.json");
        let Ok(text) = std::fs::read_to_string(path) else {
            return HashMap::new();
        };
        let Ok(serde_json::Value::Object(map)) = serde_json::from_str(&text) else {
            return HashMap::new();
        };
        map.into_iter()
            .filter_map(|(tex, v)| {
                #[allow(clippy::cast_possible_truncation)]
                let num = |k: &str| v.get(k)?.as_f64().map(|x| x as f32);
                Some((
                    tex,
                    Formula {
                        file: v.get("file")?.as_str()?.to_string(),
                        width_ex: num("width_ex")?,
                        height_ex: num("height_ex")?,
                    },
                ))
            })
            .collect()
    })
}

/// The texture of formula `file` in `color`, rasterized by resvg at `size` physical pixels
/// (cached; the SVG's `currentColor` becomes the text's colour). Rasterizing at the exact
/// size keeps the thin strokes crisp (a scaled-down large raster looked faint and soft).
fn texture(
    ctx: &egui::Context,
    file: &str,
    color: egui::Color32,
    size: [usize; 2],
) -> Option<egui::TextureHandle> {
    type Key = (String, [u8; 4], [usize; 2]);
    static CACHE: OnceLock<Mutex<HashMap<Key, egui::TextureHandle>>> = OnceLock::new();
    let key = (file.to_string(), color.to_array(), size);
    let mut cache = CACHE.get_or_init(Default::default).lock().ok()?;
    if let Some(t) = cache.get(&key) {
        return Some(t.clone());
    }
    let path = crate::levels_dir().join("math").join(file);
    let [r, g, b, _] = color.to_array();
    let svg = std::fs::read_to_string(path)
        .ok()?
        .replace("currentColor", &format!("#{r:02x}{g:02x}{b:02x}"));
    let tree =
        resvg::usvg::Tree::from_data(svg.as_bytes(), &resvg::usvg::Options::default()).ok()?;
    let [w, h] = size;
    let mut pixmap = resvg::tiny_skia::Pixmap::new(
        u32::try_from(w.max(1)).ok()?,
        u32::try_from(h.max(1)).ok()?,
    )?;
    let natural = tree.size();
    #[allow(clippy::cast_precision_loss)]
    let transform = resvg::tiny_skia::Transform::from_scale(
        w as f32 / natural.width(),
        h as f32 / natural.height(),
    );
    resvg::render(&tree, transform, &mut pixmap.as_mut());
    let image = egui::ColorImage::from_rgba_premultiplied([w.max(1), h.max(1)], pixmap.data());
    let handle = ctx.load_texture(
        format!("math/{}", key.0),
        image,
        egui::TextureOptions::LINEAR,
    );
    cache.insert(key, handle.clone());
    Some(handle)
}

/// Pieces of a text: plain text and `$...$` formulas (an unmatched `$` is text).
#[derive(Debug, PartialEq, Eq)]
pub enum Piece<'a> {
    Text(&'a str),
    Math(&'a str),
}

pub fn pieces(text: &str) -> Vec<Piece<'_>> {
    let mut out = Vec::new();
    let mut rest = text;
    while let Some(open) = rest.find('$') {
        let Some(len) = rest[open + 1..].find('$') else {
            break;
        };
        if open > 0 {
            out.push(Piece::Text(&rest[..open]));
        }
        out.push(Piece::Math(&rest[open + 1..open + 1 + len]));
        rest = &rest[open + 2 + len..];
    }
    if !rest.is_empty() {
        out.push(Piece::Text(rest));
    }
    out
}

/// A wrapped paragraph of `text` with its formulas typeset, in the style of `style`
/// applied to each word (e.g. italics).
pub fn paragraph(ui: &mut egui::Ui, text: &str, style: impl Fn(egui::RichText) -> egui::RichText) {
    // One ex of the formulas (MathJax's unit, its x-height) matched to the x-height of the
    // interface font, about 0.52 em of the body text.
    let em = ui
        .style()
        .text_styles
        .get(&egui::TextStyle::Body)
        .map_or(14.0, |f| f.size);
    let ex = 0.52 * em;
    let ppp = ui.ctx().pixels_per_point();
    let color = ui.visuals().text_color();
    ui.horizontal_wrapped(|ui| {
        ui.spacing_mut().item_spacing.x = 0.0;
        for piece in pieces(text) {
            match piece {
                Piece::Text(t) => {
                    // Word by word, so that the row wraps between words; each keeps the
                    // space that follows it.
                    for word in t.split_inclusive(' ') {
                        ui.label(style(egui::RichText::new(word)));
                    }
                }
                Piece::Math(tex) => {
                    let rendered = manifest().get(tex).and_then(|f| {
                        let size = egui::vec2(f.width_ex * ex, f.height_ex * ex);
                        #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
                        let px = [
                            (size.x * ppp).round() as usize,
                            (size.y * ppp).round() as usize,
                        ];
                        Some((texture(ui.ctx(), &f.file, color, px)?, size))
                    });
                    if let Some((tex_handle, size)) = rendered {
                        ui.add(egui::Image::new((tex_handle.id(), size)))
                            .on_hover_text(tex);
                    } else {
                        ui.label(egui::RichText::new(tex).monospace());
                    }
                }
            }
        }
    });
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pieces_split_text_and_formulas() {
        assert_eq!(
            pieces("a $x^2$ b $y$"),
            vec![
                Piece::Text("a "),
                Piece::Math("x^2"),
                Piece::Text(" b "),
                Piece::Math("y"),
            ]
        );
        assert_eq!(pieces("costs $5"), vec![Piece::Text("costs $5")]);
        assert_eq!(pieces("$x$"), vec![Piece::Math("x")]);
    }
}
