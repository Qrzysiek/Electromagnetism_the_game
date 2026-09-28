// Renders LaTeX formulas to standalone SVG with MathJax (run by scripts/levels.py).
//
//   node scripts/math/render.mjs <formulas.json> <out-dir>
//
// <formulas.json> is a JSON array of TeX strings (the $...$ of level descriptions). For
// each it writes <out-dir>/<key>.svg (key: the first 16 hex digits of the TeX's SHA-256)
// and <out-dir>/manifest.json, mapping each TeX string to its file and its size in ex
// (MathJax's unit: the height of an x) with the baseline offset, so that the game can
// place it inline with text. SVGs no longer in the list are removed. The paths are
// inline (no font cache) and coloured `currentColor`, which the game replaces with the
// colour of the surrounding text.

import { readFileSync, writeFileSync, readdirSync, unlinkSync, mkdirSync } from 'node:fs';
import { createHash } from 'node:crypto';
import { join } from 'node:path';
import { mathjax } from 'mathjax-full/js/mathjax.js';
import { TeX } from 'mathjax-full/js/input/tex.js';
import { SVG } from 'mathjax-full/js/output/svg.js';
import { liteAdaptor } from 'mathjax-full/js/adaptors/liteAdaptor.js';
import { RegisterHTMLHandler } from 'mathjax-full/js/handlers/html.js';
import { AllPackages } from 'mathjax-full/js/input/tex/AllPackages.js';

const [list, outDir] = process.argv.slice(2);
const formulas = [...new Set(JSON.parse(readFileSync(list, 'utf8')))].sort();
mkdirSync(outDir, { recursive: true });

const adaptor = liteAdaptor();
RegisterHTMLHandler(adaptor);
const doc = mathjax.document('', {
  InputJax: new TeX({ packages: AllPackages }),
  OutputJax: new SVG({ fontCache: 'none' }),
});

const key = (tex) => createHash('sha256').update(tex).digest('hex').slice(0, 16);
const manifest = {};
for (const tex of formulas) {
  const node = doc.convert(tex, { display: false });
  const svg = adaptor.firstChild(node);
  let text = adaptor.outerHTML(svg);
  if (text.includes('data-mjx-error') || text.includes('merror')) {
    console.error(`MathJax error in: ${tex}`);
    process.exit(1);
  }
  const ex = (name) => parseFloat(adaptor.getAttribute(svg, name));
  const width = ex('width');
  const height = ex('height');
  const style = adaptor.getAttribute(svg, 'style') || '';
  const valign = parseFloat((style.match(/vertical-align:\s*(-?[\d.]+)ex/) || [0, '0'])[1]);
  // Absolute size for rasterizers (1 ex = 16 px; the game scales to its font).
  text = text
    .replace(/ width="[\d.]+ex"/, ` width="${(width * 16).toFixed(3)}px"`)
    .replace(/ height="[\d.]+ex"/, ` height="${(height * 16).toFixed(3)}px"`)
    .replace(/ style="[^"]*"/, '');
  const k = key(tex);
  writeFileSync(join(outDir, `${k}.svg`), text + '\n');
  manifest[tex] = { file: `${k}.svg`, width_ex: width, height_ex: height, baseline_ex: -valign };
}
writeFileSync(join(outDir, 'manifest.json'), JSON.stringify(manifest, null, 1) + '\n');
const keep = new Set(Object.values(manifest).map((m) => m.file));
for (const f of readdirSync(outDir)) {
  if (f.endsWith('.svg') && !keep.has(f)) unlinkSync(join(outDir, f));
}
console.log(`${formulas.length} formulas rendered to ${outDir}`);
