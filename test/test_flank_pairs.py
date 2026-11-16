"""Flank-pair check (classify.flank_pairs): rescue and veto, with the real aligner.

Three repeat situations found in real data (K. sucrofermentans DSM15973), each
simulated with query == reference, so every correct call is ANC:

* (a) an element whose flank lies inside an adjacent multi-copy IS: the construct's
  repeat flank lands at another copy -> AMB_NO_MIDPOINT_SPAN. Rescued to ANC.
* (b1) the site lies in a repeat whose other (inverted) copy lacks the element: the
  construct aligns contiguously to the empty copy -> false NEW. Vetoed.
* (b2) a head-to-tail tandem pair of one IS: excising either copy leaves flanks that
  read as an uninterrupted copy -> false NEW. Vetoed.

Each test also checks the construct-only answer, so a change in minimap2 behaviour
that stops exercising the check shows up as a failure rather than a silent pass.
"""

import random
import shutil

import pytest

from MOBster import constants as C
from MOBster.classify import _classify_records
from MOBster.inputs import read_is_table
from MOBster.pipeline import classify_elements, read_alignments

from sim import insert, rand_seq, revcomp

pytestmark = pytest.mark.skipif(shutil.which("minimap2") is None, reason="minimap2 not on PATH")


def run(tmp, genome, elements, strand="+"):
    """Classify `elements` (name, s, e) with query == reference; returns
    ({name: row}, {name: construct-only result})."""
    n, query = len(genome), genome
    if strand == "-":                                 # flip the query only
        query = revcomp(genome)
        elements = [(nm, n - e + 1, n - s + 1) for nm, s, e in elements]
    (tmp / "ref.fasta").write_text(f">ref\n{genome}\n")
    lines = ["seqID,family,cluster,isBegin,isEnd,isLen,type"]
    lines += [f"query,X,{nm},{s},{e},{e - s + 1},c" for nm, s, e in elements]
    (tmp / "q.csv").write_text("\n".join(lines) + "\n")
    els = read_is_table(tmp / "q.csv")
    rows = classify_elements(els, {"query": query}, tmp / "ref.fasta", tmp / "work", 1)
    alns = read_alignments(tmp / "work" / "constructs.sam")
    construct_only = {el.name: _classify_records(alns.get(el.index, []), el.length, 2 * C.DELTA, C.DELTA)[0]
                      for el in els}
    return {r["element_name"]: r for r in rows}, construct_only


@pytest.mark.parametrize("strand", ["+", "-"])
def test_flank_inside_neighbouring_repeat_is_rescued(tmp_path, strand):
    rng = random.Random(7)
    g = rand_seq(120_000, rng)
    rep = rand_seq(1200, rng)                         # multi-copy IS, 8 identical copies
    for pos in (110_000, 95_000, 80_000, 65_000, 50_000, 35_000, 20_000):
        g, _, _ = insert(g, pos, rep, 4)
    g, bs, be = insert(g, 10_000, rep, 4)             # this copy gets a neighbour
    a_el, b_el = rand_seq(900, rng), rand_seq(700, rng)
    # A directly after the repeat copy: A's LEFT flank is entirely repeat sequence.
    g = g[:be] + a_el + g[be:]
    a_s, a_e = be + 1, be + len(a_el)
    # B directly before the same copy: B's RIGHT flank is entirely repeat sequence.
    g = g[:bs - 1] + b_el + g[bs - 1:]
    b_s, b_e = bs, bs + len(b_el) - 1
    a_s, a_e = a_s + len(b_el), a_e + len(b_el)
    rows, construct_only = run(tmp_path, g, [("A", a_s, a_e), ("B", b_s, b_e)], strand)
    assert {k: v["call"] for k, v in construct_only.items()} == {"A": "AMB", "B": "AMB"}
    for name, el_len in (("A", 900), ("B", 700)):
        r = rows[name]
        assert (r["call"], r["status"]) == ("ANC", C.OK), r
        assert r["n_flank_pairs"] == 1 and r["midpoint_gap"] == el_len and r["anc_residual"] == 0
        assert r["best_locus"].endswith(f":{strand}")


def test_empty_inverted_paralog_vetoes_new(tmp_path):
    rng = random.Random(11)
    g = rand_seq(80_000, rng)
    site = rand_seq(2400, rng)
    el = rand_seq(850, rng)
    g = g[:40_000] + revcomp(site) + g[40_000:]       # the copy WITHOUT the element
    g = g[:10_000] + site[:1200] + el + site[1200:] + g[10_000:]
    s, e = 10_000 + 1200 + 1, 10_000 + 1200 + len(el)
    rows, construct_only = run(tmp_path, g, [("E", s, e)])
    assert construct_only["E"]["call"] == "NEW"       # the false NEW this check exists for
    r = rows["E"]
    assert (r["call"], r["status"]) == ("AMB", C.AMB_MULTI_PAIRING), r
    assert r["n_flank_pairs"] == 2


def test_tandem_pair_vetoes_new(tmp_path):
    rng = random.Random(13)
    g = rand_seq(60_000, rng)
    t = rand_seq(1500, rng)
    g = g[:30_000] + t + t + g[30_000:]               # head-to-tail, no spacer
    first = ("T1", 30_001, 31_500)
    second = ("T2", 31_501, 33_000)
    rows, construct_only = run(tmp_path, g, [first, second])
    assert {k: v["call"] for k, v in construct_only.items()} == {"T1": "NEW", "T2": "NEW"}
    for name in ("T1", "T2"):
        assert (rows[name]["call"], rows[name]["status"]) == ("AMB", C.AMB_MULTI_PAIRING), rows[name]
