"""Composite transposition (classify.find_composites, CLAUDE.md deviation 11), real aligner.

A unit IS - cargo - revcomp(IS) is inserted as one into a copy of the reference. Each
copy alone is AMB: one flank sits at the insertion site, the other inside the cargo,
whose source is elsewhere. The unit is called from the two OUTER flanks, which must read
as a clean NEW insertion of the whole span; the INNER flanks place the cargo's source.

Each test also checks the per-copy answer without the composite pass, so a test that
stops exercising the pass fails instead of passing silently.
"""

import random
import shutil

import pytest

from MOBster import classify
from MOBster import constants as C
from MOBster.classify import find_composites
from MOBster.inputs import Element, read_is_table
from MOBster.pipeline import classify_elements

from sim import insert, rand_seq, revcomp

pytestmark = pytest.mark.skipif(shutil.which("minimap2") is None, reason="minimap2 not on PATH")

IS_LEN = (1500,)
CARGO = (6000, 66_000)        # (length, end): cargo source is reference [60000, 66000), 0-based
SITE = 120_000


def build(rng, ref, cargo, tsd=9, site_deletion=0):  # noqa: ARG001 (rng kept for symmetry)
    """Query = reference + unit at SITE. Returns (query, [(name, s, e)])."""
    is_seq = build.is_seq
    unit = is_seq + cargo + revcomp(is_seq)
    q, s, e = insert(ref, SITE, unit, tsd)
    if site_deletion:                                  # remove reference bases right after the unit
        q = q[:e + tsd] + q[e + tsd + site_deletion:]
    n = len(is_seq)
    return q, [("IS1", s, s + n - 1), ("IS2", e - n + 1, e)]


def run(tmp, ref, query, elements, strand="+"):
    tmp.mkdir(exist_ok=True)
    if strand == "-":                                  # flip the query only
        n = len(query)
        query = revcomp(query)
        elements = sorted(((nm, n - e + 1, n - s + 1) for nm, s, e in elements), key=lambda x: x[1])
    (tmp / "ref.fasta").write_text(f">ref\n{ref}\n")
    lines = ["seqID,family,cluster,isBegin,isEnd,isLen,type"]
    lines += [f"query,X,ISX,{s},{e},{e - s + 1},c" for _, s, e in elements]
    (tmp / "q.csv").write_text("\n".join(lines) + "\n")
    els = read_is_table(tmp / "q.csv")
    rows = classify_elements(els, {"query": query}, tmp / "ref.fasta", tmp / "work", 1)
    return rows


@pytest.fixture
def genome():
    rng = random.Random(23)
    build.is_seq = rand_seq(IS_LEN[0], rng)
    return rng, rand_seq(150_000, rng)


@pytest.mark.parametrize("strand", ["+", "-"])
@pytest.mark.parametrize("orient", ["fwd", "rev"])
def test_duplicative_composite(tmp_path, genome, monkeypatch, strand, orient):
    rng, ref = genome
    cargo = ref[CARGO[1] - CARGO[0]:CARGO[1]]
    if orient == "rev":
        cargo = revcomp(cargo)
    q, els = build(rng, ref, cargo)

    # Without the pass, both copies are AMB of a rescuable status.
    monkeypatch.setattr("MOBster.pipeline.find_composites", lambda *a: {})
    before = run(tmp_path / "off", ref, q, els, strand)
    assert all(r["call"] == "AMB" and r["status"] in classify.RESCUABLE for r in before), before
    monkeypatch.undo()

    rows = run(tmp_path, ref, q, els, strand)
    a, b = rows
    for r, partner in ((a, b), (b, a)):
        assert (r["call"], r["status"]) == (C.COMPOSITE, C.OK), r
        assert r["composite_partner"] == f"ISX:{partner['element_start']}-{partner['element_end']}"
        assert r["midpoint_gap"] == -9 and r["tsd_estimate"] == 9
        assert r["reference_position"] == SITE        # first TSD base, as for NEW
        assert r["n_cargo_sources"] == 1
    # Cargo read along the query from the first copy to the second; reversed when the
    # unit's cargo is reverse-complemented relative to the reference (flipping the query
    # flips it too).
    lo, hi = CARGO[1] - CARGO[0] + 1, CARGO[1]
    reverse = (orient == "rev") != (strand == "-")
    assert a["source_cargo"] == b["source_cargo"] == (f"ref:{hi}-{lo}" if reverse else f"ref:{lo}-{hi}")


def test_novel_cargo_has_no_source(tmp_path, genome):
    rng, ref = genome
    q, els = build(rng, ref, rand_seq(CARGO[0], rng))
    for r in run(tmp_path, ref, q, els):
        assert r["call"] == C.COMPOSITE, r
        assert r["n_cargo_sources"] == 0 and r["source_cargo"] == ""


def test_multicopy_cargo_source_left_empty(tmp_path, genome):
    rng, ref = genome
    cargo = ref[CARGO[1] - CARGO[0]:CARGO[1]]
    ref2 = ref[:100_000] + cargo + ref[100_000:]       # a second reference copy of the cargo
    q, els = build(rng, ref2, cargo)
    for r in run(tmp_path, ref2, q, els):
        assert r["call"] == C.COMPOSITE, r
        assert r["n_cargo_sources"] == 2 and r["source_cargo"] == ""


def test_outer_flanks_not_adjacent_stay_amb(tmp_path, genome):
    """Unit inserted with a 3 kb deletion at the site: the outer flanks are 3 kb apart in
    the reference, which is not a clean insertion, so neither copy is COMPOSITE."""
    rng, ref = genome
    q, els = build(rng, ref, ref[CARGO[1] - CARGO[0]:CARGO[1]], site_deletion=3000)
    for r in run(tmp_path, ref, q, els):
        assert r["call"] == "AMB" and r["status"] in classify.RESCUABLE, r
        assert r["composite_partner"] == ""


def test_copy_in_two_units_is_left_alone(monkeypatch):
    """One-partner rule: if every pair of three copies passed the outer test, no copy
    can be assigned a partner."""
    monkeypatch.setattr(classify, "outer_pair", lambda *a: object())
    els = [Element(i, "q", "ISX", 1000 * (2 * i + 1), 1000 * (2 * i + 1) + 500, "c", "X", "ISX")
           for i in range(3)]
    assert find_composites(els, {}) == {}
    # ... and with only two copies, they pair (call_from_group needs a real group, so
    # check the pairing step through the same monkeypatch on two copies).
    calls = []
    monkeypatch.setattr(classify, "outer_pair", lambda l, r, n: calls.append(n))
    find_composites(els[:2], {})
    assert calls == [els[1].e - els[0].s + 1]


def test_cargo_shorter_than_flank_has_no_source():
    assert classify.cargo_sources([], [], C.DELTA - 1) == []
