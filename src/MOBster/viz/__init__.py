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
import re
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
# Paired map only: reference IS elements, unclassified. Drawn on the reference arcs at the
# outermost radius (the query rings never share a record with it).
REF_RING = dict(radius=312, color="#3b6ea5", shape="circle", size=36, hover_px=3.8)
REF_RING_LABEL = "Reference IS element (unclassified)"
# Record arc colours: one neutral for the single map, one per genome for the paired map.
ARC_COLOR = "#d0d7de"
GENOMES = {"reference": dict(label="Reference record", color="#9cc3e6"),
           "query": dict(label="Query record", color="#f0c08a")}
CHORD_RADIUS = 236               # ribbons end just inside the innermost (AMB) ring
CHORD_MIN_RAD = 0.07             # minimum drawn width of a ribbon end (~16 px at CHORD_RADIUS)

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


def _check_calls(rows: list[dict], record_lengths: dict[str, int]) -> None:
    missing = sorted({r["query_record"] for r in rows} - set(record_lengths))
    if missing:
        raise InputError(f"call table names query record(s) absent from the query FASTA: {missing}")
    unknown = sorted({r["call"] for r in rows} - set(RINGS))
    if unknown:
        raise InputError(f"call table has unknown call value(s): {unknown}")


def _call_elements(rows: list[dict], block: dict[str, str]) -> list[dict]:
    """One glyph per call-table row, on the ring of its call."""
    elements = []
    for idx, r in enumerate(rows):
        s, e = int(r["element_start"]), int(r["element_end"])
        tip = [[label, r[col]] for col, label in TOOLTIP_FIELDS if r.get(col, "") != ""]
        if is_partial(r.get("element_type", "")):    # complete is the default; only flag partials
            tip.insert(3, ["Type", "partial"])
        elements.append(dict(i=idx, ring=r["call"], block=block[r["query_record"]],
                             true_pos=(s + e) / 2, pos=(s + e) / 2, tip=tip))
    return elements


def _nudge_rings(elements: list[dict], rings: dict, blocks: dict[str, int]) -> None:
    """§6: nudge per (record, ring); different rings never overlap radially.
    blocks: {block id: record length}, in layout order."""
    total = sum(blocks.values())
    for ring, conf in rings.items():
        sep = min_sep_bp(total, len(blocks), conf["radius"], 2 * conf["hover_px"])
        for bid, length in blocks.items():
            group = [el for el in elements if el["ring"] == ring and el["block"] == bid]
            new = nudge([el["true_pos"] for el in group], sep, 0, length)
            for el, p in zip(group, new):
                el["pos"] = p


def build_map_data(rows: list[dict], record_lengths: dict[str, int], query_name: str,
                   reference_name: str) -> dict:
    """Everything the page needs, as one JSON-serialisable dict."""
    _check_calls(rows, record_lengths)
    ids = list(record_lengths)                      # FASTA order (§3)
    block = {rid: f"r{i}" for i, rid in enumerate(ids)}   # CSS-safe ids for CircosJS
    total = sum(record_lengths.values())
    records = [dict(id=block[rid], name=rid, len=record_lengths[rid], color=ARC_COLOR)
               for rid in ids]
    elements = _call_elements(rows, block)
    _nudge_rings(elements, RINGS, {block[rid]: record_lengths[rid] for rid in ids})

    counts = {k: sum(el["ring"] == k for el in elements) for k in RINGS}
    return dict(
        meta=dict(query=query_name, reference=reference_name, version=__version__,
                  n_elements=len(elements), counts=counts, paired=False),
        geometry=dict(size=SIZE, layout_inner=LAYOUT_INNER, layout_outer=LAYOUT_OUTER,
                      gap=GAP_RAD, ticks=tick_conf(total), chord_radius=CHORD_RADIUS, chord_min_rad=CHORD_MIN_RAD),
        rings={k: dict(v, label=RING_LABELS[k]) for k, v in RINGS.items()},
        order=list(RINGS),
        genomes=[],
        records=records,
        elements=elements,
        chords=[],
    )


# --------------------------------------------------------------------------- paired map
#
# Reference and query records share one circle (reference first), each genome's arcs in
# its own colour. Query elements keep their call rings; reference IS elements (ISEScan on
# the reference) are unclassified, on one ring. Each COMPOSITE unit with a known cargo
# source gets a ribbon from that reference span to the cargo between its two copies.

_LOCUS = re.compile(r"^(?P<ref>.+):(?P<a>\d+)-(?P<b>\d+)$")


def composite_chords(rows: list[dict]) -> list[dict]:
    """One chord per COMPOSITE unit whose `source_cargo` is filled.

    Returns dicts with the reference record and 1-based span (lo <= hi), the query record
    and the cargo's 1-based span (between the two copies), and `inverted`: the cargo is
    reverse-complemented relative to the reference (source_cargo start > end).
    """
    by_pos = {(r["query_record"], int(r["element_start"])): r for r in rows}
    chords = []
    for r in rows:
        if r["call"] != "COMPOSITE" or not r.get("source_cargo"):
            continue
        partner = r.get("composite_partner", "")
        m = re.match(r"^.*:(\d+)-(\d+)$", partner)
        if m is None:
            continue
        ps, pe = int(m.group(1)), int(m.group(2))
        s, e = int(r["element_start"]), int(r["element_end"])
        if s > ps or (r["query_record"], ps) not in by_pos:   # emit once, from the first copy
            continue
        src = _LOCUS.match(r["source_cargo"])
        if src is None:
            raise InputError(f"bad source_cargo {r['source_cargo']!r} (expected ref:start-end)")
        a, b = int(src.group("a")), int(src.group("b"))
        chords.append(dict(ref_record=src.group("ref"), ref_lo=min(a, b), ref_hi=max(a, b),
                           query_record=r["query_record"], query_lo=e + 1, query_hi=ps - 1,
                           inverted=a > b, source=r["source_cargo"],
                           copies=f"{s}-{e}, {ps}-{pe}"))
    return chords


def build_paired_data(rows: list[dict], query_lengths: dict[str, int],
                      reference_lengths: dict[str, int], reference_elements: list,
                      query_name: str, reference_name: str) -> dict:
    """The paired map: reference + query records, reference IS elements, cargo chords.

    reference_elements: inputs.Element list from ISEScan on the reference.
    """
    _check_calls(rows, query_lengths)
    missing = sorted({el.record for el in reference_elements} - set(reference_lengths))
    if missing:
        raise InputError(f"reference IS table names record(s) absent from the reference FASTA: {missing}")
    rblock = {rid: f"f{i}" for i, rid in enumerate(reference_lengths)}
    qblock = {rid: f"q{i}" for i, rid in enumerate(query_lengths)}
    records = ([dict(id=rblock[k], name=k, len=v, color=GENOMES["reference"]["color"], genome="reference")
                for k, v in reference_lengths.items()]
               + [dict(id=qblock[k], name=k, len=v, color=GENOMES["query"]["color"], genome="query")
                  for k, v in query_lengths.items()])
    blocks = {r["id"]: r["len"] for r in records}
    total = sum(blocks.values())

    elements = _call_elements(rows, qblock)
    for el in reference_elements:
        tip = [["Element", el.name], ["Family", el.family], ["Cluster", el.cluster],
               ["Genome", "reference"], ["Record", el.record], ["Start", el.s], ["End", el.e]]
        if not el.complete:
            tip.insert(3, ["Type", "partial"])
        tip = [[k, str(v)] for k, v in tip if v != ""]
        elements.append(dict(i=len(elements), ring="REF", block=rblock[el.record],
                             true_pos=(el.s + el.e) / 2, pos=(el.s + el.e) / 2, tip=tip))
    rings = {"REF": REF_RING, **RINGS}
    _nudge_rings(elements, rings, blocks)

    chords = []
    for c in composite_chords(rows):
        if c["ref_record"] not in rblock:
            raise InputError(f"source_cargo names a record absent from the reference FASTA: {c['ref_record']}")
        tip = [["Composite unit", f"{c['query_record']}:{c['copies']}"],
               ["Cargo (query)", f"{c['query_record']}:{c['query_lo']}-{c['query_hi']}"],
               ["Cargo source", c["source"]],
               ["Orientation", "inverted" if c["inverted"] else "same as reference"]]
        chords.append(dict(c, src_block=rblock[c["ref_record"]], dst_block=qblock[c["query_record"]],
                           tip=tip))

    labels = dict(RING_LABELS, REF=REF_RING_LABEL)
    counts = {k: sum(el["ring"] == k for el in elements) for k in rings}
    return dict(
        meta=dict(query=query_name, reference=reference_name, version=__version__,
                  n_elements=len(elements), counts=counts, paired=True, n_chords=len(chords)),
        geometry=dict(size=SIZE, layout_inner=LAYOUT_INNER, layout_outer=LAYOUT_OUTER,
                      gap=GAP_RAD, ticks=tick_conf(total), chord_radius=CHORD_RADIUS, chord_min_rad=CHORD_MIN_RAD),
        rings={k: dict(v, label=labels[k]) for k, v in rings.items()},
        order=list(rings),
        genomes=[dict(key=k, **v) for k, v in GENOMES.items()],
        records=records,
        elements=elements,
        chords=chords,
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


SINGLE_NOTE = """<p>Each element is drawn once, on the ring for its call, at its position in the
        <b>query</b>. Arcs are query records, in FASTA order, sized by length.</p>"""
PAIRED_NOTE = """<p>Reference records (first) and query records share the circle, each genome in
        its own colour, sized by length. Query elements sit on the ring for their call;
        reference IS elements are unclassified, on one ring.</p>
      <p>A ribbon joins each composite unit's cargo in the query to its source in the
        reference. A ribbon that crosses over itself marks cargo inserted
        <b>inverted</b> relative to the reference.</p>"""


def render_html(data: dict) -> str:
    page = _resource("template.html")
    paired = data["meta"].get("paired", False)
    subs = {
        "__TITLE__": html.escape(f"MOBster {'paired map' if paired else 'map'}: {data['meta']['query']}"),
        "__HEADING__": "MOBster: reference and query" if paired else "MOBster: IS element calls",
        "<!--__NOTE__-->": PAIRED_NOTE if paired else SINGLE_NOTE,
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


def write_paired_map(rows: list[dict], query_fasta: str | Path, reference_fasta: str | Path,
                     reference_elements: list, out: str | Path,
                     query_lengths: dict[str, int] | None = None,
                     reference_lengths: dict[str, int] | None = None) -> Path:
    """Write the paired reference + query map. Lengths may be passed to skip re-reading FASTAs."""
    if query_lengths is None:
        query_lengths = {k: len(v) for k, v in load_fasta(query_fasta, "query").items()}
    if reference_lengths is None:
        reference_lengths = {k: len(v) for k, v in load_fasta(reference_fasta, "reference").items()}
    data = build_paired_data(rows, query_lengths, reference_lengths, reference_elements,
                             Path(query_fasta).name, Path(reference_fasta).name)
    out = Path(out)
    out.write_text(render_html(data), encoding="utf-8")
    logger.info("Wrote paired map %s (%d reference IS, %d query elements, %d cargo chord(s))", out,
                data["meta"]["counts"]["REF"], len(rows), data["meta"]["n_chords"])
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
    p.add_argument("--is-table-reference", type=Path,
                   help="ISEScan table for the reference: also write the paired map")
    p.add_argument("--paired-out", type=Path,
                   help="paired map HTML (default: mobster_map_paired.html next to --calls)")
    args = p.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s", stream=sys.stderr)
    try:
        rows = read_calls(args.calls)
        write_map(rows, args.query, args.reference,
                  args.out or args.calls.with_name("mobster_map.html"))
        if args.is_table_reference is not None:
            from ..inputs import read_is_table
            write_paired_map(rows, args.query, args.reference, read_is_table(args.is_table_reference),
                             args.paired_out or args.calls.with_name("mobster_map_paired.html"))
    except InputError as exc:
        logger.error("%s", exc)
        return 2
    return 0
