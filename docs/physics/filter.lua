-- Pandoc filter for PHYSICS.md -> PHYSICS.tex (see main.tex and scripts/physics_pdf.py).
-- Backtick spans are either code (paths, identifiers: monospace, breakable) or formulas
-- written in Unicode (upright, in the text font); fenced blocks are displayed formulas.
-- Tables get column widths from their contents and a smaller font. Headings keep their
-- plain text (code included) for the PDF bookmarks.

local function escape(s)
  s = s:gsub("\\", "\0")
  s = s:gsub("([{}$&#_%%])", "\\%1")
  s = s:gsub("~", "\\textasciitilde{}")
  s = s:gsub("%^", "\\textasciicircum{}")
  s = s:gsub("%z", "\\textbackslash{}")
  return s
end

local ROOTS = { "crates/", "scripts/", "levels/", "docs/", "src/", "tests/", "target/" }
local EXTS = { ".rs", ".py", ".wgsl", ".wls", ".json", ".md", ".toml", ".tex", ".lua" }

local function is_code(s)
  if s:find("::", 1, true) then return true end
  for _, e in ipairs(EXTS) do
    if s:find(e, 1, true) then return true end
  end
  for _, r in ipairs(ROOTS) do
    if s:find(r, 1, true) then return true end
  end
  -- snake_case identifiers (not subscripts such as B_z, F_RR, T_edge).
  if s:find("%l%l_%l") then return true end
  if s:find("^cargo ") or s:find("^python ") or s:find("^EM_") or s:find("^%-%-") then
    return true
  end
  -- A single identifier (sample, Outcome, retarded_time, run()).
  if #s >= 4 and s:find("^[%a_][%w_]*$") and s:find("%l%l") then return true end
  return s:find("^[%a_][%w_]*%(%)$") ~= nil
end

local function code(el)
  local t = el.text
  if is_code(t) then
    -- Allow line breaks after separators in long paths and names.
    local e = escape(t):gsub("([/:%.])", "%1\\allowbreak{}"):gsub("(\\_)", "%1\\allowbreak{}")
    return pandoc.RawInline("latex", "\\code{" .. e .. "}")
  end
  return pandoc.RawInline("latex", "\\formula{" .. escape(t) .. "}")
end

local function code_block(el)
  local lines = {}
  for line in (el.text .. "\n"):gmatch("(.-)\n") do
    -- Keep indentation and runs of spaces (continuation lines are aligned by them).
    local lead = line:match("^( *)")
    local body = escape(line:sub(#lead + 1)):gsub("  +", function(sp)
      return string.rep("\\ ", #sp)
    end)
    table.insert(lines, (#lead > 0 and ("\\hspace*{" .. (#lead * 0.5) .. "em}") or "") .. body)
  end
  return pandoc.RawBlock("latex",
    "\\begin{formulablock}\n" .. table.concat(lines, "\n") .. "\n\\end{formulablock}")
end

-- Column widths proportional to (mean cell length)^0.7, at least 6 % each.
local function table_widths(tbl)
  local n = #tbl.colspecs
  local sums, counts = {}, {}
  for i = 1, n do sums[i], counts[i] = 0, 0 end
  local function add(rows)
    for _, row in ipairs(rows) do
      for i, cell in ipairs(row.cells) do
        if i <= n then
          sums[i] = sums[i] + #pandoc.utils.stringify(cell.contents)
          counts[i] = counts[i] + 1
        end
      end
    end
  end
  add(tbl.head.rows)
  for _, body in ipairs(tbl.bodies) do add(body.body) end
  local w, total = {}, 0
  for i = 1, n do
    w[i] = math.max((sums[i] / math.max(counts[i], 1)) ^ 0.7, 1)
    total = total + w[i]
  end
  local s = 0
  for i = 1, n do
    w[i] = math.max(w[i] / total, 0.06)
    s = s + w[i]
  end
  for i = 1, n do tbl.colspecs[i] = { tbl.colspecs[i][1], w[i] / s * 0.98 } end
  return {
    pandoc.RawBlock("latex", "\\begingroup\\small"),
    tbl,
    pandoc.RawBlock("latex", "\\endgroup"),
  }
end

local SECTIONS = { "section", "subsection", "subsubsection", "paragraph", "subparagraph" }

return {
  -- First, the headings' plain text (with their code) for the bookmarks.
  { Header = function(h)
      h.attributes["pdf"] = pandoc.utils.stringify(h.content)
      return h
    end },
  { Code = code, CodeBlock = code_block, Table = table_widths },
  { Header = function(h)
      local body = pandoc.write(pandoc.Pandoc({ pandoc.Plain(h.content) }), "latex")
      return pandoc.RawBlock("latex", string.format("\\%s{\\texorpdfstring{%s}{%s}}\\label{%s}",
        SECTIONS[h.level] or "subparagraph", body:gsub("%s+$", ""),
        escape(h.attributes["pdf"]), h.identifier))
    end },
}
