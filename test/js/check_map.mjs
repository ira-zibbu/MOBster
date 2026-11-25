// Load a generated MOBster map in jsdom (no network: jsdom fetches no external
// resources unless told to), hover every glyph, and print a JSON report on stdout.
//   node test/js/check_map.mjs path/to/mobster_map.html
import { readFileSync } from "node:fs";
import { JSDOM, VirtualConsole } from "jsdom";

const html = readFileSync(process.argv[2], "utf8");
const errors = [];
const vc = new VirtualConsole();
vc.on("jsdomError", (e) => errors.push(String(e.message || e)));
vc.on("error", (e) => errors.push(String(e)));

const dom = new JSDOM(html, { runScripts: "dangerously", pretendToBeVisual: true, virtualConsole: vc });
const { window } = dom;
const doc = window.document;

const out = { errors, rendered: !!(window.MOBSTER_MAP && window.MOBSTER_MAP.rendered) };
if (out.rendered) {
  const data = window.MOBSTER_MAP.data;
  const pts = [...doc.querySelectorAll("path.point")];
  out.n_points = pts.length;
  out.ring_counts = {};
  out.points = [];
  for (const p of pts) {
    const ring = p.getAttribute("data-ring");
    out.ring_counts[ring] = (out.ring_counts[ring] || 0) + 1;
    // Glyph position is relative to its record: CircosJS wraps each record's glyphs in
    // <g class="block" transform="rotate(start°)">. Compose both to get page coordinates.
    const m = /translate\(\s*([-\d.e]+)\s*,\s*([-\d.e]+)\s*\)/.exec(p.getAttribute("transform"));
    const blk = p.closest("g.block");
    const rm = blk && /rotate\(\s*([-\d.e]+)/.exec(blk.getAttribute("transform") || "");
    const rot = rm ? Number(rm[1]) * Math.PI / 180 : 0;
    const lx = m ? Number(m[1]) : NaN, ly = m ? Number(m[2]) : NaN;
    const gx = lx * Math.cos(rot) - ly * Math.sin(rot), gy = lx * Math.sin(rot) + ly * Math.cos(rot);
    p.dispatchEvent(new window.MouseEvent("mouseover", { bubbles: true }));
    const tip = doc.querySelector("div.circos-tooltip");
    const cells = tip ? [...tip.querySelectorAll("td")].map((td) => td.textContent) : [];
    const fields = {};
    for (let k = 0; k + 1 < cells.length; k += 2) fields[cells[k]] = cells[k + 1];
    p.dispatchEvent(new window.MouseEvent("mouseout", { bubbles: true }));
    out.points.push({
      i: Number(p.getAttribute("data-i")), ring,
      x: gx, y: gy,
      fill: p.getAttribute("fill"), d: p.getAttribute("d"), tooltip: fields,
    });
  }
  // record arcs: CircosJS gives each block a <g> whose class is the block id (r0, r1, ...)
  out.blocks = [...doc.querySelectorAll(".cs-layout g")].map((g) => g.getAttribute("class"))
    .filter((c) => /^[rfq]\d+$/.test(c || ""));
  // composite chords (paired map): geometry hooks + tooltip on hover
  out.chords = [...doc.querySelectorAll("path.chord")].map((c) => {
    c.dispatchEvent(new window.MouseEvent("mouseover", { bubbles: true }));
    const tip = doc.querySelector("div.chord-tooltip");
    const cells = tip ? [...tip.querySelectorAll("td")].map((td) => td.textContent) : [];
    const fields = {};
    for (let k = 0; k + 1 < cells.length; k += 2) fields[cells[k]] = cells[k + 1];
    const shown = tip ? tip.style.display : "";
    const note = tip && tip.querySelector(".moved") ? tip.querySelector(".moved").textContent : "";
    c.dispatchEvent(new window.MouseEvent("mouseout", { bubbles: true }));
    return {
      i: Number(c.getAttribute("data-i")), inverted: c.getAttribute("data-inverted") === "1",
      widened: c.getAttribute("data-widened") === "1",
      note,
      s: c.getAttribute("data-s").split(",").map(Number), t: c.getAttribute("data-t").split(",").map(Number),
      d: c.getAttribute("d"), fill: c.getAttribute("fill"), tooltip: fields, shown,
    };
  });
  out.external_refs = [...doc.querySelectorAll("[src],[href]")].map((e) => e.getAttribute("src") || e.getAttribute("href"))
    .filter((u) => /^(https?:)?\/\//.test(u));
  out.meta = data.meta;
}
process.stdout.write(JSON.stringify(out));
