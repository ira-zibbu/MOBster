"""Interactive map (docs/viz_instructions.md).

Python-level tests of the data/nudging, plus the §9 acceptance checks run on the
generated HTML in jsdom (no network). The jsdom checks need node and
`npm install --prefix test/js`; they are skipped otherwise.
"""

import json
import math
import shutil
import subprocess
from pathlib import Path

import pytest

from MOBster.inputs import Element, InputError
from MOBster.pipeline import COLUMNS, write_table
from MOBster.viz import (GAP_RAD, GENOMES, REF_RING, RINGS, TOOLTIP_FIELDS, build_map_data,
                         build_paired_data, composite_chords, nudge, read_calls, render_html,
                         write_map)

HERE = Path(__file__).parent
JSDOM_OK = shutil.which("node") and (HERE / "js" / "node_modules" / "jsdom").exists()

RECORDS = {"chrA": 600_000, "plasmid|B": 300_000, "contig_3": 100_000}   # odd id on purpose


def row(rec, s, e, call, status="OK", fam="IS1", cluster="IS1_1", **kw):
    r = {c: "" for c in COLUMNS}
    r.update(query_record=rec, element_name=cluster, element_family=fam, element_cluster=cluster,
             element_start=s, element_end=e, element_length=e - s + 1, element_type="c",
             delta=500, call=call, status=status)
    r.update(kw)
    return {k: str(v) for k, v in r.items()}


def sample_rows():
    return [
        row("chrA", 100_000, 100_768, "NEW", midpoint_gap=-9, anc_residual=-777, tsd_estimate=9,
            reference_position=99_991),
        row("chrA", 102_000, 102_768, "NEW", midpoint_gap=-9, anc_residual=-777, tsd_estimate=9,
            fam="IS150", cluster="IS150_7"),                               # 2 kb from the first
        row("chrA", 300_000, 301_342, "ANC", fam="IS186", cluster="IS186_2", midpoint_gap=1343,
            anc_residual=0),
        row("chrA", 301_500, 302_800, "ANC", fam="IS3", cluster="IS3_9", midpoint_gap=1300,
            anc_residual=-1),
        row("plasmid|B", 150_000, 150_768, "AMB", "AMB_NO_MIDPOINT_SPAN", element_type="p"),
        row("plasmid|B", 150_900, 151_700, "AMB", "AMB_MULTI_LOCUS", midpoint_gap=-9,
            anc_residual=-810),
        row("contig_3", 200, 968, "AMB", "AMB_FLANK_TRUNCATED"),
        row("contig_3", 50_000, 51_342, "ANC", fam="IS4", cluster="IS4_1", midpoint_gap=1343,
            anc_residual=0),
        row("contig_3", 80_000, 81_342, "COMPOSITE", fam="IS4", cluster="IS4_1", midpoint_gap=2,
            anc_residual=-1341, tsd_estimate=-2, composite_partner="IS4_1:95001-96342",
            source_cargo="chrR:604466-590503", n_cargo_sources=1),
        row("contig_3", 95_001, 96_342, "COMPOSITE", fam="IS4", cluster="IS4_1", midpoint_gap=2,
            anc_residual=-1340, tsd_estimate=-2, composite_partner="IS4_1:80000-81342",
            source_cargo="chrR:604466-590503", n_cargo_sources=1),
    ]


# ------------------------------------------------------------------ python-level

def test_nudge_separates_close_points_and_leaves_far_ones():
    out = nudge([100.0, 101.0, 5000.0], 10.0, 0, 10_000)
    assert out[1] - out[0] == pytest.approx(10.0)
    assert (out[0] + out[1]) / 2 == pytest.approx(100.5)        # centred on the true positions
    assert out[2] == 5000.0


def test_nudge_respects_record_bounds_and_order():
    out = nudge([0.0, 0.0, 0.0], 10.0, 0, 1000)
    assert out == pytest.approx([0.0, 10.0, 20.0])
    out = nudge([1000.0, 995.0], 10.0, 0, 1000)
    assert max(out) <= 1000 and abs(out[0] - out[1]) == pytest.approx(10.0)
    assert out[0] > out[1]                                       # input order preserved


def test_every_row_on_exactly_one_ring():
    rows = sample_rows()
    data = build_map_data(rows, RECORDS, "q.fasta", "r.fasta")
    assert sorted(e["i"] for e in data["elements"]) == list(range(len(rows)))
    assert data["meta"]["counts"] == {"NEW": 2, "COMPOSITE": 2, "ANC": 3, "AMB": 3}
    assert [r["name"] for r in data["records"]] == list(RECORDS)          # FASTA order


def test_tooltip_omits_empty_fields():
    data = build_map_data(sample_rows(), RECORDS, "q.fasta", "r.fasta")
    partial = data["elements"][4]["tip"]
    labels = [k for k, _ in partial]
    assert "midpoint_gap" not in labels and "tsd_estimate" not in labels
    assert {"Element", "Family", "Cluster", "Call", "Status"} <= set(labels)


def test_tooltip_flags_partial_elements_only():
    data = build_map_data(sample_rows(), RECORDS, "q.fasta", "r.fasta")
    tips = [dict(e["tip"]) for e in data["elements"]]
    assert tips[4]["Type"] == "partial"
    assert all("Type" not in t for i, t in enumerate(tips) if i != 4)


def test_missing_record_fails_loudly():
    with pytest.raises(InputError, match="absent from the query FASTA"):
        build_map_data(sample_rows(), {"chrA": 600_000}, "q.fasta", "r.fasta")


def test_html_is_self_contained():
    page = render_html(build_map_data(sample_rows(), RECORDS, "q.fasta", "r.fasta"))
    assert "<script src" not in page and "<link" not in page
    assert "Circos" in page and "d3js.org v7.9.0" in page
    assert page.count("</script>") == 4                          # d3, circos, data, app


def test_write_map_from_csv_and_fasta(tmp_path):
    fa = tmp_path / "q.fasta"
    fa.write_text("".join(f">{k}\n{'A' * v}\n" for k, v in RECORDS.items()))
    write_table(sample_rows(), tmp_path / "calls.csv")
    out = write_map(read_calls(tmp_path / "calls.csv"), fa, tmp_path / "ref.fasta", tmp_path / "m.html")
    page = out.read_text()
    assert "q.fasta" in page and "ref.fasta" in page and "v0.1.0" in page


# ------------------------------------------------------------------ §9 acceptance (jsdom)

@pytest.fixture(scope="module")
def report(tmp_path_factory):
    if not JSDOM_OK:
        pytest.skip("node + jsdom not installed (npm install --prefix test/js)")
    d = tmp_path_factory.mktemp("map")
    data = build_map_data(sample_rows(), RECORDS, "q.fasta", "r.fasta")
    (d / "m.html").write_text(render_html(data))
    res = subprocess.run(["node", str(HERE / "js" / "check_map.mjs"), str(d / "m.html")],
                         capture_output=True, text=True, check=True)
    return json.loads(res.stdout), data


def test_1_renders_offline(report):
    rep, _ = report
    assert rep["rendered"] and rep["errors"] == [] and rep["external_refs"] == []


def test_2_glyph_counts_match_csv(report):
    rep, data = report
    assert rep["ring_counts"] == {k: v for k, v in data["meta"]["counts"].items() if v}
    assert sorted(p["i"] for p in rep["points"]) == list(range(len(sample_rows())))


def test_3_close_elements_independently_hoverable(report):
    rep, _ = report
    by_ring = {}
    for p in rep["points"]:
        by_ring.setdefault(p["ring"], []).append(p)
    for ring, pts in by_ring.items():
        need = 2 * RINGS[ring]["hover_px"] - 1e-6
        for a in range(len(pts)):
            for b in range(a + 1, len(pts)):
                dist = math.dist((pts[a]["x"], pts[a]["y"]), (pts[b]["x"], pts[b]["y"]))
                assert dist >= need, (ring, pts[a]["i"], pts[b]["i"], dist)
    # the two NEW elements 2 kb apart each produce their own tooltip
    new = [p for p in rep["points"] if p["ring"] == "NEW"]
    assert {p["tooltip"]["Start"] for p in new} == {"100000", "102000"}


def test_4_tooltip_shows_every_nonempty_field(report):
    rep, _ = report
    rows = sample_rows()
    for p in rep["points"]:
        r = rows[p["i"]]
        expected = {label: r[col] for col, label in TOOLTIP_FIELDS if r[col] != ""}
        if r["element_type"] == "p":
            expected["Type"] = "partial"
        assert p["tooltip"] == expected


def test_5_family_in_tooltip_colour_by_ring_only(report):
    rep, _ = report
    for p in rep["points"]:
        assert p["tooltip"]["Family"] and p["tooltip"]["Cluster"]
        assert p["fill"] == RINGS[p["ring"]]["color"]           # never varies with family
    fams = {p["tooltip"]["Family"] for p in rep["points"] if p["ring"] == "ANC"}
    assert len(fams) > 1                                         # several families, one colour


def test_6_multi_record_arcs_proportional(report):
    rep, data = report
    assert rep["blocks"] == ["r0", "r1", "r2"]
    # Each glyph sits at the angle implied by proportional arc sizing (CircosJS convention:
    # theta from 12 o'clock, clockwise; arcs share 2*pi minus one gap per record).
    total, n = sum(RECORDS.values()), len(RECORDS)
    span = 2 * math.pi - n * GAP_RAD
    offsets, acc = {}, 0
    for k, rec in enumerate(data["records"]):
        offsets[rec["id"]] = (acc, k)
        acc += rec["len"]
    for p in rep["points"]:
        el = data["elements"][p["i"]]
        off, k = offsets[el["block"]]
        theta = (off + el["pos"]) / total * span + k * GAP_RAD
        got = math.atan2(p["y"], p["x"]) + math.pi / 2
        assert math.isclose(math.cos(got), math.cos(theta), abs_tol=1e-6)
        assert math.isclose(math.sin(got), math.sin(theta), abs_tol=1e-6)


def test_amb_drawn_as_tick(report):
    rep, _ = report
    assert {p["d"] for p in rep["points"] if p["ring"] == "AMB"} == {"M-1,-5h2v10h-2z"}


# ------------------------------------------------------------------ paired map

REF_RECORDS = {"chrR": 700_000, "pR": 50_000}


def paired_rows():
    """sample_rows() (one inverted unit on contig_3) plus a forward unit on chrA."""
    return sample_rows() + [
        row("chrA", 400_000, 401_342, "COMPOSITE", fam="IS4", cluster="IS4_1", midpoint_gap=2,
            composite_partner="IS4_1:420001-421342", source_cargo="chrR:100001-118658", n_cargo_sources=1),
        row("chrA", 420_001, 421_342, "COMPOSITE", fam="IS4", cluster="IS4_1", midpoint_gap=2,
            composite_partner="IS4_1:400000-401342", source_cargo="chrR:100001-118658", n_cargo_sources=1),
        # a 2 kb cargo: narrower than CHORD_MIN_RAD, so drawn widened (the larger ones may be too)
        row("chrA", 500_000, 500_768, "COMPOSITE", composite_partner="IS1_1:502769-503536",
            source_cargo="pR:22000-20001", n_cargo_sources=1),
        row("chrA", 502_769, 503_536, "COMPOSITE", composite_partner="IS1_1:500000-500768",
            source_cargo="pR:22000-20001", n_cargo_sources=1),
    ]


def ref_elements():
    return [Element(0, "chrR", "IS4_1", 50_000, 51_342, "c", "IS4", "IS4_1"),
            Element(1, "pR", "IS1_1", 1_000, 1_768, "p", "IS1", "IS1_1")]


def test_composite_chords_one_per_unit_with_orientation():
    chords = composite_chords(paired_rows())
    assert len(chords) == 3
    inv, fwd = [c for c in chords if c["ref_record"] == "chrR"]
    if not inv["inverted"]:
        inv, fwd = fwd, inv
    assert (inv["ref_record"], inv["ref_lo"], inv["ref_hi"], inv["inverted"]) == ("chrR", 590503, 604466, True)
    assert (inv["query_record"], inv["query_lo"], inv["query_hi"]) == ("contig_3", 81_343, 95_000)
    assert (fwd["ref_lo"], fwd["ref_hi"], fwd["inverted"]) == (100_001, 118_658, False)
    assert (fwd["query_lo"], fwd["query_hi"]) == (401_343, 420_000)


def test_composite_without_source_gets_no_chord():
    rows = paired_rows()
    for r in rows:
        if r["call"] == "COMPOSITE" and r["element_start"] in ("400000", "420001"):
            r["source_cargo"], r["n_cargo_sources"] = "", "2"
    assert len(composite_chords(rows)) == 2


def test_paired_layout_reference_first_and_unclassified_ring():
    data = build_paired_data(paired_rows(), RECORDS, REF_RECORDS, ref_elements(), "q.fasta", "r.fasta")
    assert [r["name"] for r in data["records"]] == list(REF_RECORDS) + list(RECORDS)
    assert {r["genome"]: r["color"] for r in data["records"]} == {
        "reference": GENOMES["reference"]["color"], "query": GENOMES["query"]["color"]}
    assert data["order"][0] == "REF" and data["meta"]["counts"]["REF"] == 2
    ref_ids = {r["id"] for r in data["records"] if r["genome"] == "reference"}
    for el in data["elements"]:
        assert (el["ring"] == "REF") == (el["block"] in ref_ids)
    tips = [dict(el["tip"]) for el in data["elements"] if el["ring"] == "REF"]
    assert tips[0] == {"Element": "IS4_1", "Family": "IS4", "Cluster": "IS4_1", "Genome": "reference",
                       "Record": "chrR", "Start": "50000", "End": "51342"}
    assert tips[1]["Type"] == "partial" and "Call" not in tips[1]
    assert len(data["chords"]) == 3


def test_paired_unknown_reference_record_fails():
    bad = ref_elements() + [Element(2, "nope", "X", 1, 10, "c")]
    with pytest.raises(InputError, match="absent from the reference FASTA"):
        build_paired_data(paired_rows(), RECORDS, REF_RECORDS, bad, "q.fasta", "r.fasta")
    rows = paired_rows()
    rows[-1]["source_cargo"] = rows[-2]["source_cargo"] = "chrX:1-20000"
    with pytest.raises(InputError, match="source_cargo names a record absent"):
        build_paired_data(rows, RECORDS, REF_RECORDS, ref_elements(), "q.fasta", "r.fasta")


@pytest.fixture(scope="module")
def paired_report(tmp_path_factory):
    if not JSDOM_OK:
        pytest.skip("node + jsdom not installed (npm install --prefix test/js)")
    d = tmp_path_factory.mktemp("paired")
    data = build_paired_data(paired_rows(), RECORDS, REF_RECORDS, ref_elements(), "q.fasta", "r.fasta")
    (d / "m.html").write_text(render_html(data))
    res = subprocess.run(["node", str(HERE / "js" / "check_map.mjs"), str(d / "m.html")],
                         capture_output=True, text=True, check=True)
    return json.loads(res.stdout), data


def test_paired_renders_with_every_glyph(paired_report):
    rep, data = paired_report
    assert rep["rendered"] and rep["errors"] == [] and rep["external_refs"] == []
    assert rep["ring_counts"] == {k: v for k, v in data["meta"]["counts"].items() if v}
    assert rep["blocks"] == ["f0", "f1", "q0", "q1", "q2"]
    for p in rep["points"]:
        if p["ring"] == "REF":
            assert p["fill"] == REF_RING["color"] and p["tooltip"]["Genome"] == "reference"


def _angle(data, block, pos):
    total = sum(r["len"] for r in data["records"])
    span = 2 * math.pi - len(data["records"]) * GAP_RAD
    off = 0
    for k, r in enumerate(data["records"]):
        if r["id"] == block:
            return (off + pos) / total * span + k * GAP_RAD
        off += r["len"]
    raise KeyError(block)


def _widen(a0, a1, m):
    if a1 - a0 >= m:
        return a0, a1
    mid = (a0 + a1) / 2
    return mid - m / 2, mid + m / 2


def _q_endpoint(d):
    """End point of the first quadratic curve in an SVG path (the source -> target leg)."""
    import re
    m = re.search(r"Q\s*[-\d.e]+\s*,\s*[-\d.e]+\s*,\s*([-\d.e]+)\s*,\s*([-\d.e]+)", d)
    return float(m.group(1)), float(m.group(2))


def test_chords_join_source_to_cargo_and_twist_when_inverted(paired_report):
    rep, data = paired_report
    assert len(rep["chords"]) == 3
    R, m = data["geometry"]["chord_radius"], data["geometry"]["chord_min_rad"]
    for ch in rep["chords"]:
        c = data["chords"][ch["i"]]
        assert ch["inverted"] == c["inverted"]
        raw_s = _angle(data, c["src_block"], c["ref_lo"] - 1), _angle(data, c["src_block"], c["ref_hi"])
        raw_t = _angle(data, c["dst_block"], c["query_lo"] - 1), _angle(data, c["dst_block"], c["query_hi"])
        s0, s1 = _widen(*raw_s, m)
        t0, t1 = _widen(*raw_t, m)
        small = raw_s[1] - raw_s[0] < m or raw_t[1] - raw_t[0] < m
        assert ch["widened"] == small and bool(ch["note"]) == small
        assert s1 - s0 >= m - 1e-9 and t1 - t0 >= m - 1e-9
        assert ch["s"] == pytest.approx([s0, s1])
        # Plain ribbon: source end -> cargo start. Twisted (inverted): source end -> cargo end.
        first = t1 if c["inverted"] else t0
        assert ch["t"] == pytest.approx([first, t0 if c["inverted"] else t1])
        x, y = _q_endpoint(ch["d"])
        assert (x, y) == pytest.approx((R * math.sin(first), -R * math.cos(first)), abs=1e-3)
        assert ch["fill"] == RINGS["COMPOSITE"]["color"]
        assert ch["shown"] == "block"
        assert ch["tooltip"]["Orientation"] == ("inverted" if c["inverted"] else "same as reference")
        assert ch["tooltip"]["Cargo source"] == c["source"]


def test_single_map_has_no_chords_or_reference(report):
    rep, data = report
    assert rep["chords"] == [] and "REF" not in data["meta"]["counts"]
