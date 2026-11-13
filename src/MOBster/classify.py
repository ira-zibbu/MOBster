"""Parse alignments, measure the signed midpoint gap, classify (spec §4.6).

Conventions
-----------
* All geometry here is 0-based half-open, as pysam exposes it. Nothing here
  converts coordinates except when filling 1-based output fields.
* "Construct coordinates" are positions on the full construct (length Lc),
  including any hard-clipped prefix. For reverse-strand records pysam reports
  the query as stored, i.e. reverse-complemented, so construct coordinates are
  already in the flipped frame and ascend with reference coordinates. This is
  asserted by test/test_strand.py; no explicit flip is applied here.
* midpoint_gap > 0: reference consumed more (deletion vs construct) -> ANC.
  midpoint_gap < 0: construct consumed more (insertion, duplicated TSD) -> NEW.
  Never take its absolute value.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass, field
from typing import Iterable

import pysam

from . import constants as C

_HARD_CLIP = 5


@dataclass
class Rec:
    """One alignment record in construct coordinates."""
    ref: str
    strand: str
    qs: int            # construct start (0-based)
    qe: int            # construct end (exclusive)
    rs: int            # reference start (0-based)
    re: int            # reference end (exclusive)
    secondary: bool
    aln: pysam.AlignedSegment = field(repr=False, compare=False)
    hclip: int = 0     # leading hard clip; offset from stored-seq to construct coords

    @classmethod
    def from_aln(cls, aln: pysam.AlignedSegment) -> "Rec":
        cig = aln.cigartuples or []
        hclip = cig[0][1] if cig and cig[0][0] == _HARD_CLIP else 0
        return cls(
            ref=aln.reference_name,
            strand="-" if aln.is_reverse else "+",
            qs=hclip + aln.query_alignment_start,
            qe=hclip + aln.query_alignment_end,
            rs=aln.reference_start,
            re=aln.reference_end,
            secondary=aln.is_secondary,
            aln=aln,
            hclip=hclip,
        )

    @property
    def aligned_len(self) -> int:
        return self.qe - self.qs

    def columns(self) -> Iterable[tuple[int | None, int | None]]:
        """Aligned columns in construct coordinates, soft clips excluded."""
        qas, qae, off = self.aln.query_alignment_start, self.aln.query_alignment_end, self.hclip
        for q, r in self.aln.get_aligned_pairs():
            if q is not None:
                if q < qas or q >= qae:   # soft-clipped base, not an alignment column
                    continue
                q += off
            yield q, r

    def match_stats(self) -> tuple[int, int]:
        """(matching bases, aligned M/=/X columns), mismatches from NM minus indel bases."""
        m = ins = dels = 0
        for op, ln in self.aln.cigartuples or []:
            if op in (0, 7, 8):
                m += ln
            elif op == 1:
                ins += ln
            elif op == 2:
                dels += ln
        nm = self.aln.get_tag("NM") if self.aln.has_tag("NM") else 0
        return max(m - (nm - ins - dels), 0), m


@dataclass
class Walk:
    window_sum: int            # net indel within [M-w, M+w], signed (+del, -ins)
    max_off_indel: int         # largest single indel event wholly outside the window
    ref_at_mid: int | None     # 0-based ref pos of first aligned construct base >= M
    del_span: tuple[int, int] | None = None  # ref span [start, end) of in-window deletion columns


def walk_record(rec: Rec, M: int, w: int) -> Walk:
    """Walk get_aligned_pairs once, accumulating the window sum and the off-midpoint guard."""
    lo, hi = M - w, M + w
    total = 0
    max_off = 0
    ref_at_mid = None
    del_lo = del_hi = None              # reference extent of deletion columns in the window
    cursor = rec.qs                     # last non-None construct pos + 1
    run_kind, run_len, run_in = None, 0, False

    def close_run():
        nonlocal max_off
        if run_kind is not None and not run_in:
            max_off = max(max_off, run_len)

    for q, r in rec.columns():
        if q is not None and r is not None:
            kind, pos = None, q
            if ref_at_mid is None and q >= M:
                ref_at_mid = r
        elif q is None:
            kind, pos = "D", cursor     # deletion: placed at the construct cursor
        else:
            kind, pos = "I", q
        in_win = lo <= pos <= hi
        if kind == "D" and in_win:
            total += 1
            del_lo = r if del_lo is None else min(del_lo, r)
            del_hi = r + 1 if del_hi is None else max(del_hi, r + 1)
        elif kind == "I" and in_win:
            total -= 1
        if kind != run_kind:
            close_run()
            run_kind, run_len, run_in = kind, 0, False
        if kind is not None:
            run_len += 1
            run_in = run_in or in_win
        if q is not None:
            cursor = q + 1
    close_run()
    span = (del_lo, del_hi) if del_lo is not None else None
    return Walk(total, max_off, ref_at_mid, span)


@dataclass
class GroupResult:
    ref: str
    strand: str
    recs: list[Rec]
    shape: str                  # "one", "two", "complex", "nospan"
    gap: int | None = None
    max_off_indel: int = 0
    coverage: float = 0.0
    ref_at_mid: int | None = None
    # Reference stretch absent from the construct at the junction, 0-based half-open:
    # the ancestral copy (TSD-excluded) when the call is ANC.
    ref_span: tuple[int, int] | None = None

    @property
    def spans(self) -> bool:
        return self.shape in ("one", "two")

    @property
    def score(self) -> int:
        return sum(r.aln.get_tag("AS") if r.aln.has_tag("AS") else 0 for r in self.recs)

    @property
    def locus_str(self) -> str:
        return (f"{self.ref}:{min(r.rs for r in self.recs) + 1}-"
                f"{max(r.re for r in self.recs)}:{self.strand}")

    @property
    def identity(self) -> float:
        num = den = 0
        for r in self.recs:
            a, b = r.match_stats()
            num += a
            den += b
        return num / den if den else 0.0


def union_coverage(intervals: Iterable[tuple[int, int]], Lc: int) -> float:
    covered, end = 0, -1
    for s, e in sorted(intervals):
        if e <= end:
            continue
        covered += e - max(s, end)
        end = e
    return covered / Lc


def evaluate_group(ref: str, strand: str, recs: list[Rec], Lc: int, M: int,
                   w: int = C.MIDPOINT_WINDOW, ovl_max: int = C.OVL_MAX) -> GroupResult:
    recs = sorted(recs, key=lambda r: (r.qs, r.rs))
    g = GroupResult(ref, strand, recs, shape="nospan")
    g.coverage = union_coverage(((r.qs, r.qe) for r in recs), Lc)

    if len(recs) == 1:
        r = recs[0]
        if r.qs <= M - w and r.qe >= M + w:
            wk = walk_record(r, M, w)
            g.shape, g.gap, g.max_off_indel, g.ref_at_mid = "one", wk.window_sum, wk.max_off_indel, wk.ref_at_mid
            g.ref_span = wk.del_span
        return g

    if len(recs) == 2:
        a, b = recs
        dq, dr = b.qs - a.qe, b.rs - a.re
        # Junction must lie in the midpoint window (same tolerance as the one-record
        # shape, so both shapes accept the same boundary error) and the records must
        # be ordered consistently on both axes, allowing ovl_max of micro-homology.
        brackets = a.qe <= M + w and b.qs >= M - w and a.qs < b.qs and a.qe < b.qe
        if brackets and dq >= -ovl_max and dr >= -ovl_max:
            wa, wb = walk_record(a, M, w), walk_record(b, M, w)
            g.shape, g.gap = "two", dr - dq
            g.max_off_indel = max(wa.max_off_indel, wb.max_off_indel)
            g.ref_at_mid = b.rs
            # Records are sorted along the construct, and on either strand construct
            # coordinates ascend with reference ones, so the gap [a.re, b.rs) is the
            # reference copy. With micro-homology overlap the boundary is arbitrary
            # by up to the overlap length.
            g.ref_span = (a.re, b.rs)
            return g
        # Two records that don't bracket: complex if either touches the midpoint,
        # otherwise neither reaches the junction at all.
        if any(r.qs <= M <= r.qe for r in recs) or (a.qe <= M <= b.qs):
            g.shape = "complex"
        return g

    g.shape = "complex"
    return g


def _secondary_loci(secondaries: list[Rec], Lc: int, M: int, w: int) -> int:
    """Distinct alternative placements that would themselves qualify as a locus."""
    seen = set()
    for r in secondaries:
        if r.qs <= M - w and r.qe >= M + w and r.aligned_len / Lc >= C.MIN_COVERAGE:
            seen.add((r.ref, r.strand, r.rs))
    return len(seen)


def classify_construct(alns: list[pysam.AlignedSegment], element_length: int,
                       Lc: int, M: int) -> dict:
    """Classify one construct from all of its alignment records."""
    w = C.MIDPOINT_WINDOW
    mapped = [a for a in alns if not a.is_unmapped]
    out = dict(n_records=len(mapped), n_loci=0, best_locus="", construct_coverage="",
               identity="", midpoint_gap="", anc_residual="", tsd_estimate="",
               reference_position="", reference_end="")
    if not mapped:
        return dict(out, call="AMB", status=C.AMB_NO_LOCUS)

    recs = [Rec.from_aln(a) for a in mapped]
    groups: dict[tuple[str, str], list[Rec]] = defaultdict(list)
    for r in recs:
        if not r.secondary:
            groups[(r.ref, r.strand)].append(r)
    results = [evaluate_group(ref, st, rs, Lc, M) for (ref, st), rs in groups.items()]

    spanning = [g for g in results if g.spans]
    loci = [g for g in spanning if g.coverage >= C.MIN_COVERAGE]
    n_loci = len(loci) + _secondary_loci([r for r in recs if r.secondary], Lc, M, w)
    out["n_loci"] = n_loci

    # Best group to report: a locus if any, else the best spanning group.
    pool = loci or spanning
    best = max(pool, key=lambda g: (g.coverage, g.score)) if pool else None
    if best is not None:
        out.update(best_locus=best.locus_str,
                   construct_coverage=round(best.coverage, 4),
                   identity=round(best.identity, 4),
                   midpoint_gap=best.gap,
                   anc_residual=best.gap - element_length)

    if not spanning:
        status = C.AMB_COMPLEX if any(g.shape == "complex" for g in results) else C.AMB_NO_MIDPOINT_SPAN
        return dict(out, call="AMB", status=status)
    if not loci:
        return dict(out, call="AMB", status=C.AMB_LOW_COVERAGE)
    if n_loci > 1:
        return dict(out, call="AMB", status=C.AMB_MULTI_LOCUS)

    assert best is not None and best.gap is not None
    if best.max_off_indel > C.MAX_OFFSET_INDEL:
        return dict(out, call="AMB", status=C.AMB_INDEL_OFF_MIDPOINT)
    if abs(best.gap) <= C.TSD_TOLERANCE:
        site = "" if best.ref_at_mid is None else best.ref_at_mid + 1   # 1-based for output
        return dict(out, call="NEW", status=C.OK, tsd_estimate=-best.gap, reference_position=site)
    if abs(best.gap - element_length) <= C.anc_tolerance(element_length):
        pos = end = ""
        if best.ref_span is not None:
            pos, end = best.ref_span[0] + 1, best.ref_span[1]   # 1-based inclusive for output
        return dict(out, call="ANC", status=C.OK, reference_position=pos, reference_end=end)
    return dict(out, call="AMB", status=C.AMB_GAP_UNEXPLAINED)
