"""Integration tests with the real aligner on a simulated genome.

Includes the minus-strand verification required by spec §4.6: pysam's raw
construct coordinates must ascend with reference coordinates for reverse-strand
records, so no explicit flip is applied in classify.py.
"""

import random
import shutil
from pathlib import Path

import pysam
import pytest

from MOBster import constants as C
from MOBster.classify import Rec
from MOBster.cli import main
from MOBster.inputs import read_is_table
from MOBster.pipeline import classify_elements

from sim import insert, rand_seq, revcomp

pytestmark = pytest.mark.skipif(shutil.which("minimap2") is None, reason="minimap2 not on PATH")

IS_LEN = 1347


@pytest.fixture(scope="module")
def genomes():
    rng = random.Random(42)
    base = rand_seq(60_000, rng)
    is_a = rand_seq(IS_LEN, rng)
    is_s = rand_seq(768, rng)                              # IS1-sized: < asm5's 1 kb bandwidth
    # Reference carries ancestral copies at 10_000 (TSD 9) and 20_000 (short, TSD 9).
    ref, s_s, e_s = insert(base, 20_000, is_s, 9)
    ref, _, _ = insert(ref, 10_000, is_a, 9)
    s_s, e_s = s_s + IS_LEN + 9, e_s + IS_LEN + 9
    # Query = reference + new copies.
    q, s_c, e_c = insert(ref, 45_000, revcomp(is_a), 0)    # NEW, minus orientation, no TSD
    q, s_b, e_b = insert(q, 30_000, is_a, 5)               # NEW, TSD 5
    s_c, e_c = s_c + IS_LEN + 5, e_c + IS_LEN + 5          # B lies left of C, shifting it
    s_a, e_a = 10_009, 10_009 + IS_LEN - 1                 # ancestral copy (unchanged coords)
    assert q[s_a - 1:e_a] == is_a and q[s_b - 1:e_b] == is_a and q[s_c - 1:e_c] == revcomp(is_a)
    assert q[s_b - 6:s_b - 1] == q[e_b:e_b + 5]            # TSD flanks the element
    near_end = len(q) - 300                                # flank truncated
    elements = [("IS_A", s_a, e_a, "c"), ("IS_S", s_s, e_s, "c"), ("IS_B", s_b, e_b, "c"), ("IS_C", s_c, e_c, "c"),
                ("IS_P", s_b, e_b, "p"), ("IS_Q", s_s, e_s, "p"), ("IS_T", near_end - 100, near_end, "c")]
    expect = {"IS_A": ("ANC", C.OK), "IS_S": ("ANC", C.OK), "IS_B": ("NEW", C.OK), "IS_C": ("NEW", C.OK),
              "IS_P": ("NEW", C.OK), "IS_Q": ("ANC", C.OK), "IS_T": ("AMB", C.AMB_FLANK_TRUNCATED)}
    # Ancestral copies' TSD-excluded spans in the REFERENCE (the query is identical there).
    refspan = {"IS_A": (s_a, e_a), "IS_S": (s_s, e_s), "IS_Q": (s_s, e_s)}
    assert ref[s_a - 1:e_a] == is_a and ref[s_s - 1:e_s] == is_s
    return ref, q, elements, expect, {"IS_B": -5, "IS_C": 0, "IS_P": -5}, refspan


def write_inputs(tmp: Path, ref: str, q: str, elements, qname="query"):
    (tmp / "ref.fasta").write_text(f">ref\n{ref}\n")
    (tmp / "q.fasta").write_text(f">{qname}\n{q}\n")
    lines = ["seqID,family,cluster,isBegin,isEnd,isLen,type"]
    lines += [f"{qname},X,{n},{s},{e},{e - s + 1},{t}" for n, s, e, t in elements]
    (tmp / "q.csv").write_text("\n".join(lines) + "\n")


def run_case(tmp: Path, ref, q, elements):
    write_inputs(tmp, ref, q, elements)
    els = read_is_table(tmp / "q.csv")
    rows = classify_elements(els, {"query": q}, tmp / "ref.fasta", tmp / "work", 1)
    return {r["element_name"]: r for r in rows}


def check(rows, expect, tsd, refspan):
    for name, (call, status) in expect.items():
        assert (rows[name]["call"], rows[name]["status"]) == (call, status), rows[name]
    assert rows["IS_A"]["midpoint_gap"] == IS_LEN and rows["IS_A"]["anc_residual"] == 0
    assert rows["IS_S"]["midpoint_gap"] == 768 and rows["IS_S"]["anc_residual"] == 0
    for name, gap in tsd.items():
        assert rows[name]["midpoint_gap"] == gap
        assert rows[name]["tsd_estimate"] == -gap
    # Partial elements are classified like complete ones (README deviation 9).
    assert rows["IS_Q"]["midpoint_gap"] == 768 and rows["IS_Q"]["anc_residual"] == 0
    assert rows["IS_T"]["midpoint_gap"] == ""
    # reference_position/reference_end: the ancestral copy's span in the reference ...
    for name, (s, e) in refspan.items():
        assert (rows[name]["reference_position"], rows[name]["reference_end"]) == (s, e), rows[name]
    # ... and for NEW, the insertion point, no end.
    for name in tsd:
        assert rows[name]["reference_position"] != ""
        assert rows[name]["reference_end"] == ""


def test_plus_strand(tmp_path, genomes):
    ref, q, elements, expect, tsd, refspan = genomes
    rows = run_case(tmp_path, ref, q, elements)
    check(rows, expect, tsd, refspan)
    assert rows["IS_B"]["best_locus"].endswith(":+")
    # insertion site = reference base aligned to construct position M (1-based). With
    # minimap2's left-aligned indels this is the first duplicated base: the GD MOB position.
    assert rows["IS_B"]["reference_position"] == 30_000
    assert rows["IS_C"]["reference_position"] == 45_001   # dup 0: element follows base 45000


def test_minus_strand(tmp_path, genomes):
    """Reverse-complemented query: every construct aligns on '-'; calls must not change."""
    ref, q, elements, expect, tsd, refspan = genomes
    n = len(q)
    rc_elements = [(nm, n - e + 1, n - s + 1, t) for nm, s, e, t in elements]
    rows = run_case(tmp_path, ref, revcomp(q), rc_elements)
    check(rows, expect, tsd, refspan)
    assert rows["IS_A"]["best_locus"].endswith(":-")


def test_minus_strand_ordering_predicate_on_raw_pysam(tmp_path, genomes):
    """Spec §4.6: assert the chaining order holds on pysam's raw coordinates for '-' hits.

    Uses split alignments of the ancestral construct (left/right flank as separate
    queries is not enough, so we align the two flanks of the reverse-complemented
    construct as one read and inspect every '-' record pair that minimap2 returns).
    """
    ref, q, elements, *_ = genomes
    s, e = elements[0][1], elements[0][2]
    construct = q[s - 1 - C.DELTA:s - 1] + q[e:e + C.DELTA]
    (tmp_path / "ref.fasta").write_text(f">ref\n{ref}\n")
    (tmp_path / "c.fasta").write_text(f">0|x|y\n{revcomp(construct)}\n")
    from MOBster.align import run_minimap2
    sam = run_minimap2(tmp_path / "ref.fasta", tmp_path / "c.fasta", tmp_path / "c.sam", 1)
    with pysam.AlignmentFile(str(sam)) as fh:
        recs = [Rec.from_aln(a) for a in fh if not a.is_unmapped and not a.is_secondary]
    assert recs and all(r.strand == "-" for r in recs)
    recs.sort(key=lambda r: r.qs)
    if len(recs) == 1:
        # one record with the element-length deletion: construct and ref ascend together
        cols = [(qq, rr) for qq, rr in recs[0].columns() if qq is not None and rr is not None]
        assert all(b[0] > a[0] and b[1] > a[1] for a, b in zip(cols, cols[1:]))
    else:
        for a, b in zip(recs, recs[1:]):
            assert b.qs >= a.qe - C.OVL_MAX and b.rs >= a.re - C.OVL_MAX


def test_cli_end_to_end(tmp_path, genomes):
    ref, q, elements, expect, *_ = genomes
    write_inputs(tmp_path, ref, q, elements)
    rc = main(["--reference", str(tmp_path / "ref.fasta"), "--query", str(tmp_path / "q.fasta"),
               "--outdir", str(tmp_path / "out"), "--is-table-reference", str(tmp_path / "q.csv"),
               "--is-table-query", str(tmp_path / "q.csv")])
    assert rc == 0
    import csv
    with open(tmp_path / "out" / "mobster_calls.csv") as fh:
        rows = {r["element_name"]: r for r in csv.DictReader(fh)}
    assert {k: (v["call"], v["status"]) for k, v in rows.items()} == expect
    page = (tmp_path / "out" / "mobster_map.html").read_text()     # the map comes with every run
    assert "q.fasta" in page and "ref.fasta" in page


def test_swap_reports_loss(tmp_path, genomes):
    """Swapping reference and query: the new copies become absent -> the reference's
    copies are classified against a genome lacking them (NEW == present only in 'query')."""
    ref, q, elements, *_ = genomes
    write_inputs(tmp_path, ref, q, elements)
    (tmp_path / "r_as_q.fasta").write_text(f">query\n{ref}\n")
    (tmp_path / "r.csv").write_text("seqID,cluster,isBegin,isEnd,type\nquery,IS_A,10009,11355,c\n")
    els = read_is_table(tmp_path / "r.csv")
    rows = classify_elements(els, {"query": ref}, tmp_path / "q.fasta", tmp_path / "w", 1)
    assert rows[0]["call"] == "ANC"


def test_isescan_runs_on_query_only(tmp_path, genomes, monkeypatch):
    """The reference is never annotated; ISEScan runs once, on the query."""
    import MOBster.pipeline as pipeline
    ref, q, elements, *_ = genomes
    write_inputs(tmp_path, ref, q, elements)
    calls = []

    def fake_isescan(fasta, outdir, threads):
        calls.append(Path(fasta).name)
        return tmp_path / "q.csv"

    monkeypatch.setattr(pipeline, "run_isescan", fake_isescan)
    pipeline.run(tmp_path / "ref.fasta", tmp_path / "q.fasta", tmp_path / "out", 1, None, None)
    assert calls == ["q.fasta"]
    assert not (tmp_path / "out" / "isescan" / "reference").exists()
