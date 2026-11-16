"""Unit tests for midpoint-gap arithmetic on hand-built pysam records (no aligner)."""

import re

import pysam
import pytest

from MOBster import constants as C
from MOBster.classify import Rec, classify_construct, evaluate_group, walk_record

D = C.DELTA
LC, M = 2 * D, D
HEADER = pysam.AlignmentHeader.from_dict({"SQ": [{"SN": "chr", "LN": 1_000_000},
                                                 {"SN": "chr2", "LN": 1_000_000}]})


def seg(cigar, rs, flag=0, ref=0, nm=0, as_=1000):
    a = pysam.AlignedSegment(HEADER)
    a.query_name = "0|q|IS1"
    a.flag = flag
    a.reference_id = ref
    a.reference_start = rs
    a.cigarstring = cigar
    stored = sum(int(n) for n, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar) if op in "MIS=X")
    a.query_sequence = "A" * stored
    a.set_tag("NM", nm)
    a.set_tag("AS", as_)
    return a


def test_one_record_ancestral_deletion_is_positive():
    r = Rec.from_aln(seg(f"{D}M1347D{D}M", 1000))
    wk = walk_record(r, M, C.MIDPOINT_WINDOW)
    assert wk.window_sum == +1347
    assert wk.max_off_indel == 0


def test_one_record_novel_tsd_is_negative():
    # 9 bp TSD duplicated in the construct: 9 insertion columns at the midpoint
    r = Rec.from_aln(seg(f"{D - 9}M9I{D}M", 1000))
    assert walk_record(r, M, C.MIDPOINT_WINDOW).window_sum == -9


def test_deletion_window_membership_uses_cursor():
    # junction 15 bp off nominal: still in window
    r = Rec.from_aln(seg(f"{D + 15}M1347D{D - 15}M", 1000))
    assert walk_record(r, M, C.MIDPOINT_WINDOW).window_sum == 1347
    # junction 40 bp off nominal: outside window -> sum 0, guard fires
    r = Rec.from_aln(seg(f"{D + 40}M1347D{D - 40}M", 1000))
    wk = walk_record(r, M, C.MIDPOINT_WINDOW)
    assert wk.window_sum == 0
    assert wk.max_off_indel == 1347


def test_off_midpoint_guard_makes_amb_not_new():
    a = seg(f"{D + 40}M1347D{D - 40}M", 1000)
    out = classify_construct([a], 1347, LC, M)
    assert out["call"] == "AMB" and out["status"] == C.AMB_INDEL_OFF_MIDPOINT
    assert out["midpoint_gap"] == 0 and out["anc_residual"] == -1347


def test_small_offwindow_indel_tolerated():
    a = seg(f"200M10D{D - 200}M5I{D - 5}M", 1000)
    out = classify_construct([a], 1300, LC, M)
    assert out["call"] == "NEW" and out["midpoint_gap"] == -5 and out["tsd_estimate"] == 5


def test_soft_clips_are_not_insertions():
    r = Rec.from_aln(seg(f"30S{D - 30}M{D - 30}M30S", 1000))
    assert (r.qs, r.qe) == (30, LC - 30)
    assert walk_record(r, M, C.MIDPOINT_WINDOW).window_sum == 0


def test_hard_clip_offsets_construct_coordinates():
    r = Rec.from_aln(seg(f"{D}H{D}M", 5000, flag=2048))
    assert (r.qs, r.qe) == (D, LC)


@pytest.mark.parametrize("ovl", [0, 7])
def test_one_and_two_record_shapes_agree(ovl):
    L = 1347
    one = classify_construct([seg(f"{D}M{L}D{D}M", 1000)], L, LC, M)
    # Same alignment as primary + supplementary with hard clips, with `ovl` bp of
    # micro-homology aligned by both records.
    a = seg(f"{D + ovl}M{D - ovl}H", 1000)
    b = seg(f"{D}H{D}M", 1000 + D + L, flag=2048)
    two = classify_construct([a, b], L, LC, M)
    assert one["call"] == two["call"] == "ANC"
    assert one["midpoint_gap"] == L
    assert two["midpoint_gap"] == L - ovl + ovl  # dr - dq == L + ovl - ovl
    assert two["anc_residual"] == 0


def test_two_records_not_bracketing_is_complex():
    # junction at construct 600, well right of the window
    a = seg(f"{D + 100}M{D - 100}H", 1000)
    b = seg(f"{D + 200}H{D - 200}M", 9000, flag=2048)
    out = classify_construct([a, b], 1347, LC, M)
    assert out["status"] == C.AMB_COMPLEX


def test_reversed_order_on_reference_is_complex():
    a = seg(f"{D}M{D}H", 9000)
    b = seg(f"{D}H{D}M", 1000, flag=2048)
    assert classify_construct([a, b], 1347, LC, M)["status"] == C.AMB_COMPLEX


def test_multi_locus_from_two_groups():
    a = seg(f"{LC}M", 1000)
    b = seg(f"{LC}M", 1000, ref=1, flag=256)
    out = classify_construct([a, b], 1347, LC, M)
    assert out["status"] == C.AMB_MULTI_LOCUS and out["n_loci"] == 2


def test_partial_secondary_does_not_count_as_locus():
    a = seg(f"{LC}M", 1000)
    b = seg(f"{D - 100}M", 50_000, flag=256)
    out = classify_construct([a, b], 1347, LC, M)
    assert out["call"] == "NEW" and out["n_loci"] == 1 and out["n_records"] == 2


def test_no_records_and_no_span():
    assert classify_construct([], 1347, LC, M)["status"] == C.AMB_NO_LOCUS
    left_only = seg(f"{D - 50}M", 1000)
    assert classify_construct([left_only], 1347, LC, M)["status"] == C.AMB_NO_MIDPOINT_SPAN


def test_low_coverage_locus():
    r = seg(f"{D - 100}S200M{D - 100}S", 1000)
    assert classify_construct([r], 1347, LC, M)["status"] == C.AMB_LOW_COVERAGE


def test_gap_unexplained():
    out = classify_construct([seg(f"{D}M300D{D}M", 1000)], 1347, LC, M)
    assert out["status"] == C.AMB_GAP_UNEXPLAINED and out["midpoint_gap"] == 300


def test_anc_tolerance_relative():
    L = 2000  # tolerance max(50, 100) = 100
    assert classify_construct([seg(f"{D}M{L + 90}D{D}M", 1000)], L, LC, M)["call"] == "ANC"
    assert classify_construct([seg(f"{D}M{L + 110}D{D}M", 1000)], L, LC, M)["call"] == "AMB"


def test_identity_ignores_indels():
    r = seg(f"{D}M1347D{D}M", 1000, nm=1347 + 10)
    g = evaluate_group("chr", "+", [Rec.from_aln(r)], LC, M)
    assert g.identity == pytest.approx((LC - 10) / LC)


def test_reference_span_same_for_both_shapes():
    """ANC reference_position/reference_end = the deleted reference stretch, 1-based."""
    L = 1347
    one = classify_construct([seg(f"{D}M{L}D{D}M", 1000)], L, LC, M)
    two = classify_construct([seg(f"{D}M{D}H", 1000), seg(f"{D}H{D}M", 1000 + D + L, flag=2048)],
                             L, LC, M)
    # left flank occupies 0-based [1000, 1500), so the copy is 1-based 1501..1500+L
    assert (one["reference_position"], one["reference_end"]) == (1501, 1500 + L)
    assert (two["reference_position"], two["reference_end"]) == (1501, 1500 + L)


def test_reference_span_with_microhomology_shifts_start():
    # 7 bp aligned by both records: the boundary is arbitrary within the overlap
    L = 1347
    two = classify_construct([seg(f"{D + 7}M{D - 7}H", 1000), seg(f"{D}H{D}M", 1000 + D + L, flag=2048)],
                             L, LC, M)
    assert two["call"] == "ANC"
    assert (two["reference_position"], two["reference_end"]) == (1508, 1500 + L)


def test_reference_position_new_and_amb():
    new = classify_construct([seg(f"{D - 9}M9I{D}M", 1000)], 1347, LC, M)
    assert new["reference_position"] == 1000 + D - 9 + 1
    assert new["reference_end"] == ""
    amb = classify_construct([seg(f"{D}M300D{D}M", 1000)], 1347, LC, M)
    assert amb["reference_position"] == amb["reference_end"] == ""


# ------------------------------------------------------------- flank-pair check (README 10)
# Construct NEW with a 9 bp TSD at chr:1491 (0-based ref_at_mid), element 1347 bp.
# Flank records are (side, record); the R flank's record is in its own 500 bp frame.

SEC, SUP = 256, 2048
EL = 1347


def new_construct():
    return [seg(f"{D - 9}M9I{D}M", 1000)]


def flanks(*hits):
    """hits: (side, rs, ref, flag)"""
    return [(side, seg(f"{D}M", rs, flag=flag, ref=ref)) for side, rs, ref, flag in hits]


def test_new_confirmed_by_one_flank_pair_at_same_site():
    out = classify_construct(new_construct(), EL, LC, M, flanks(("L", 1000, 0, 0), ("R", 1491, 0, 0)))
    assert (out["call"], out["status"], out["n_flank_pairs"]) == ("NEW", C.OK, 1)
    assert out["reference_position"] == 1492


def test_new_with_no_flank_pair_is_conflict():
    out = classify_construct(new_construct(), EL, LC, M, flanks(("L", 1000, 0, 0), ("R", 1491, 1, 0)))
    assert (out["call"], out["status"], out["n_flank_pairs"]) == ("AMB", C.AMB_FLANK_CONFLICT, 0)
    assert out["reference_position"] == "" and out["tsd_estimate"] == ""


def test_new_with_flank_pair_elsewhere_is_conflict():
    # Flanks pair as NEW, but 49 kb from where the construct put the insertion.
    out = classify_construct(new_construct(), EL, LC, M, flanks(("L", 50_000, 0, 0), ("R", 50_491, 0, 0)))
    assert (out["call"], out["status"]) == ("AMB", C.AMB_FLANK_CONFLICT)


def test_new_with_site_within_window_is_confirmed():
    # Pair's insertion point 20 bp from the construct's: inside the midpoint window.
    out = classify_construct(new_construct(), EL, LC, M, flanks(("L", 1020, 0, 0), ("R", 1511, 0, 0)))
    assert (out["call"], out["status"]) == ("NEW", C.OK)


def test_new_with_only_an_anc_pair_is_conflict():
    out = classify_construct(new_construct(), EL, LC, M,
                             flanks(("L", 1000, 0, 0), ("R", 1000 + D + EL, 0, 0)))
    assert (out["call"], out["status"], out["n_flank_pairs"]) == ("AMB", C.AMB_FLANK_CONFLICT, 1)


def test_new_with_new_and_anc_pairs_is_multi_pairing():
    out = classify_construct(new_construct(), EL, LC, M,
                             flanks(("L", 1000, 0, 0), ("R", 1491, 0, 0), ("R", 1000 + D + EL, 0, SEC)))
    assert (out["call"], out["status"], out["n_flank_pairs"]) == ("AMB", C.AMB_MULTI_PAIRING, 2)


def misplaced_construct():
    # L flank at its true place; R flank's primary piece at another copy on chr2.
    return [seg(f"{D}M{D}S", 1000), seg(f"{D}H{D}M", 7000, flag=SUP, ref=1)]


def test_rescue_when_repeat_flank_explains_the_miss():
    out = classify_construct(misplaced_construct(), EL, LC, M,
                             flanks(("L", 1000, 0, 0), ("R", 7000, 1, 0), ("R", 1000 + D + EL, 0, SEC)))
    assert (out["call"], out["status"], out["n_flank_pairs"]) == ("ANC", C.OK, 1)
    assert out["midpoint_gap"] == EL and out["anc_residual"] == 0
    assert (out["reference_position"], out["reference_end"]) == (1501, 1500 + EL)
    assert out["best_locus"] == f"chr:1001-{1000 + 2 * D + EL}:+"


def test_no_rescue_when_both_flanks_are_single_copy():
    # Construct AMB for a reason no repeat explains: leave it AMB.
    out = classify_construct([seg(f"{D}M{D}S", 1000)], EL, LC, M,
                             flanks(("L", 1000, 0, 0), ("R", 1000 + D + EL, 0, 0)))
    assert (out["call"], out["status"], out["n_flank_pairs"]) == ("AMB", C.AMB_NO_MIDPOINT_SPAN, 1)


def test_supplementary_piece_is_not_a_second_copy():
    # A flank split across a linear reference's ends (primary + supplementary) is one placement.
    from MOBster.classify import flank_copies
    hits = [("L", seg(f"{D - 200}M200S", 900_000)), ("L", seg(f"{D - 200}H200M", 0, flag=SUP)),
            ("R", seg(f"{D}M", 1000))]
    assert flank_copies(hits) == {"L": 1, "R": 1}


def test_confirmed_anc_reports_the_flank_pair_span():
    # Construct: 7 bp of micro-homology let the left record run into the copy, so the
    # construct's span starts 7 bp late. The separately aligned flanks end where they end.
    construct = [seg(f"{D + 7}M{D - 7}H", 1000), seg(f"{D}H{D}M", 1000 + D + EL, flag=SUP)]
    assert classify_construct(construct, EL, LC, M)["reference_position"] == 1508
    out = classify_construct(construct, EL, LC, M, flanks(("L", 1000, 0, 0), ("R", 1000 + D + EL, 0, 0)))
    assert (out["call"], out["status"]) == ("ANC", C.OK)
    assert (out["reference_position"], out["reference_end"]) == (1501, 1500 + EL)
