"""Interactive circular map of MOBster calls (docs/viz_instructions.md).

Writes ONE self-contained HTML file: call data inlined as JSON, D3 and CircosJS
inlined from vendor/, no network access needed to render it. The page template is
a package data file (template.html), filled by plain token substitution.

Collision nudging (§6) is done here in Python rather than in the browser, so it is
deterministic and unit-testable; the geometry constants below are passed to the
page so both sides agree on radii and arc sizes.
"""

from __future__ import annotations

import argparse
import csv
import html
import json
import logging
import math
import sys
from importlib import resources
from pathlib import Path

from .. import __version__
from ..inputs import InputError, load_fasta

logger = logging.getLogger("MOBster")

# --------------------------------------------------------------------------- geometry
SIZE = 820                       # SVG width = height, px
LAYOUT_INNER, LAYOUT_OUTER = 330, 340
GAP_RAD = 0.04                   # angular gap after each record (CircosJS `gap`)

# One ring per call, outermost inward (§3). `hover_px` is the glyph's hover radius:
# glyph centres closer than 2*hover_px on the same ring are nudged apart (§6).
RINGS = {
    "NEW": dict(radius=312, color="#d1495b", shape="circle", size=56, hover_px=4.5),
    # Composite-unit copies (README deviation 11): also novel, so next to NEW, in a
    # second accent colour and a different shape so it reads apart from NEW.
    "COMPOSITE": dict(radius=290, color="#7b5ea7", shape="diamond", size=56, hover_px=4.5),
    "ANC": dict(radius=268, color="#6c7a89", shape="circle", size=36, hover_px=3.8),
    "AMB": dict(radius=246, color="#c4c4c4", shape="tick", size=0, hover_px=3.5),
}
RING_LABELS = {"NEW": "NEW: novel transposition",
               "COMPOSITE": "COMPOSITE: copy in a new composite unit", "ANC": "ANC: present in reference",
               "AMB": "AMB: not yet resolved"}

# Tooltip fields (§5), in display order: (CSV column, label).
TOOLTIP_FIELDS = [
    ("element_name", "Element"), ("element_family", "Family"), ("element_cluster", "Cluster"),
    ("query_record", "Record"), ("element_start", "Start"), ("element_end", "End"),
    ("call", "Call"), ("status", "Status"),
    ("midpoint_gap", "midpoint_gap"), ("anc_residual", "anc_residual"),
    ("tsd_estimate", "tsd_estimate"),
    ("composite_partner", "Partner"), ("source_cargo", "Cargo source"),
]


def is_partial(element_type: str) -> bool:
    return element_type.strip().lower() in ("p", "partial")


# --------------------------------------------------------------------------- nudging

def _isotonic(y: list[float]) -> list[float]:
    """Least-squares non-decreasing fit (pool adjacent violators)."""
    blocks: list[list[float]] = []          # [sum, count]
    for v in y:
        blocks.append([v, 1])
        while len(blocks) > 1 and blocks[-2][0] / blocks[-2][1] > blocks[-1][0] / blocks[-1][1]:
            s, c = blocks.pop()
            blocks[-1][0] += s
            blocks[-1][1] += c
    out = []
    for s, c in blocks:
        out += [s / c] * c
    return out


def nudge(positions: list[float], min_sep: float, lo: float, hi: float) -> list[float]:
    """Move points the least (in squared distance) so neighbours are >= min_sep apart.

    Solved exactly: with y_i = x_i - i*min_sep the spacing constraint becomes
    "y non-decreasing", i.e. isotonic regression, and the record bounds become a
    clip on y. Returns positions in the input order. Points already far apart do
    not move.
    """
    n = len(positions)
    if n < 2:
        return list(positions)
    order = sorted(range(n), key=lambda i: positions[i])
    y = [positions[i] - k * min_sep for k, i in enumerate(order)]
    fit = _isotonic(y)
    y_hi = hi - (n - 1) * min_sep
    if y_hi < lo:                       # more glyphs than fit on the record: best effort
        logger.warning("map: %d glyphs cannot be separated on a record of this size", n)
        y_hi = lo
    fit = [min(max(v, lo), y_hi) for v in fit]
    out = [0.0] * n
    for k, i in enumerate(order):
        out[i] = fit[k] + k * min_sep
    return out


def min_sep_bp(total_len: int, n_records: int, radius: float, sep_px: float) -> float:
    """Genomic distance at which two glyphs on a ring of `radius` are `sep_px` apart on screen.

    Uses the chord (straight-line screen distance), not the arc, and the CircosJS
    layout: records share 2*pi minus one GAP_RAD per record, in proportion to length.
    """
    angle = 2 * math.asin(min(1.0, sep_px / (2 * radius)))
    return angle * total_len / (2 * math.pi - n_records * GAP_RAD)


# --------------------------------------------------------------------------- data

def read_calls(path: str | Path) -> list[dict]:
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


def build_map_data(rows: list[dict], record_lengths: dict[str, int], query_name: str,
                   reference_name: str) -> dict:
    """Everything the page needs, as one JSON-serialisable dict."""
    missing = sorted({r["query_record"] for r in rows} - set(record_lengths))
    if missing:
        raise InputError(f"call table names query record(s) absent from the query FASTA: {missing}")
    unknown = sorted({r["call"] for r in rows} - set(RINGS))
    if unknown:
        raise InputError(f"call table has unknown call value(s): {unknown}")

    ids = list(record_lengths)                      # FASTA order (§3)
    block = {rid: f"r{i}" for i, rid in enumerate(ids)}   # CSS-safe ids for CircosJS
    total = sum(record_lengths.values())
    records = [dict(id=block[rid], name=rid, len=record_lengths[rid]) for rid in ids]

    elements = []
    for idx, r in enumerate(rows):
        s, e = int(r["element_start"]), int(r["element_end"])
        tip = [[label, r[col]] for col, label in TOOLTIP_FIELDS if r.get(col, "") != ""]
        if is_partial(r.get("element_type", "")):    # complete is the default; only flag partials
            tip.insert(3, ["Type", "partial"])
        elements.append(dict(i=idx, ring=r["call"], block=block[r["query_record"]],
                             true_pos=(s + e) / 2, pos=(s + e) / 2, tip=tip))

    # §6: nudge per (record, ring); different rings never overlap radially.
    for ring, conf in RINGS.items():
        sep = min_sep_bp(total, len(ids), conf["radius"], 2 * conf["hover_px"])
        for rid in ids:
            group = [el for el in elements if el["ring"] == ring and el["block"] == block[rid]]
            new = nudge([el["true_pos"] for el in group], sep, 0, record_lengths[rid])
            for el, p in zip(group, new):
                el["pos"] = p

    counts = {k: sum(el["ring"] == k for el in elements) for k in RINGS}
    return dict(
        meta=dict(query=query_name, reference=reference_name, version=__version__,
                  n_elements=len(elements), counts=counts),
        geometry=dict(size=SIZE, layout_inner=LAYOUT_INNER, layout_outer=LAYOUT_OUTER,
                      gap=GAP_RAD, ticks=tick_conf(total)),
        rings={k: dict(v, label=RING_LABELS[k]) for k, v in RINGS.items()},
        records=records,
        elements=elements,
    )


def tick_conf(total: int) -> dict:
    """Tick spacing ~1/100 of the genome, a 1-2-5 number; labels in Mb (or kb if small)."""
    raw = total / 100
    mag = 10 ** math.floor(math.log10(raw)) if raw > 0 else 1
    spacing = next(m * mag for m in (1, 2, 5, 10) if m * mag >= raw)
    if total >= 2_000_000:
        return dict(spacing=int(spacing), denominator=1_000_000, suffix="Mb")
    return dict(spacing=int(spacing), denominator=1_000, suffix="kb")


# --------------------------------------------------------------------------- HTML

def _resource(name: str) -> str:
    return resources.files(__package__).joinpath(name).read_text(encoding="utf-8")


def _script_safe(text: str) -> str:
    """Make text safe to inline in a <script> element."""
    return text.replace("</", "<\\/")


def render_html(data: dict) -> str:
    page = _resource("template.html")
    subs = {
        "__TITLE__": html.escape(f"MOBster map: {data['meta']['query']}"),
        "__QUERY__": html.escape(data["meta"]["query"]),
        "__REFERENCE__": html.escape(data["meta"]["reference"]),
        "__VERSION__": html.escape(data["meta"]["version"]),
        "/*__D3_JS__*/": _script_safe(_resource("vendor/d3-7.9.0.min.js")),
        "/*__CIRCOS_JS__*/": _script_safe(_resource("vendor/circos-2.1.0.min.js")),
        "/*__DATA_JSON__*/": _script_safe(json.dumps(data, separators=(",", ":"))),
    }
    for token, value in subs.items():
        if page.count(token) != 1:
            raise RuntimeError(f"template token {token} must occur exactly once")
        page = page.replace(token, value)
    return page


def write_map(rows: list[dict], query_fasta: str | Path, reference_fasta: str | Path,
              out: str | Path, record_lengths: dict[str, int] | None = None) -> Path:
    """Write the map. `record_lengths` may be passed to avoid re-reading the FASTA."""
    if record_lengths is None:
        record_lengths = {k: len(v) for k, v in load_fasta(query_fasta, "query").items()}
    data = build_map_data(rows, record_lengths, Path(query_fasta).name, Path(reference_fasta).name)
    out = Path(out)
    out.write_text(render_html(data), encoding="utf-8")
    logger.info("Wrote map %s (%d elements: %s)", out, data["meta"]["n_elements"],
                ", ".join(f"{k}={v}" for k, v in data["meta"]["counts"].items()))
    return out


def main(argv=None) -> int:
    """Regenerate a map from an existing call table: python -m MOBster.viz ..."""
    p = argparse.ArgumentParser(prog="python -m MOBster.viz",
                                description="Write the interactive circular map for a MOBster call table.")
    p.add_argument("--calls", required=True, type=Path, help="mobster_calls.csv")
    p.add_argument("--query", required=True, type=Path, help="query FASTA used for the run")
    p.add_argument("--reference", required=True, type=Path,
                   help="reference FASTA used for the run (its filename goes in the page header)")
    p.add_argument("--out", type=Path, help="output HTML (default: next to --calls)")
    args = p.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s", stream=sys.stderr)
    try:
        write_map(read_calls(args.calls), args.query, args.reference,
                  args.out or args.calls.with_name("mobster_map.html"))
    except InputError as exc:
        logger.error("%s", exc)
        return 2
    return 0
