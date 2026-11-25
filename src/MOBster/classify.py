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
    hclip: int = 0     # offset from stored-seq to construct coords (leading hard clip + offset)

    @classmethod
    def from_aln(cls, aln: pysam.AlignedSegment, offset: int = 0) -> "Rec":
        """offset: where the aligned sequence starts in construct coordinates (non-zero
        for a separately aligned flank, see flank_offset)."""
        cig = aln.cigartuples or []
        hclip = (cig[0][1] if cig and cig[0][0] == _HARD_CLIP else 0) + offset
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


def flank_offset(side: str, reverse: bool, delta: int) -> int:
    """Construct-coordinate start of a separately aligned flank.

    A '+' construct is L then R. A '-' record stores the reverse complement,
    rc(R) then rc(L), so the offsets swap; construct coordinates then ascend with
    reference coordinates on both strands, as for whole-construct records.
    """
    return delta if (side == "R") != reverse else 0


def flank_pairs(flank_alns: list[tuple[str, pysam.AlignedSegment]], element_length: int,
                Lc: int, M: int) -> list[tuple[str, GroupResult]]:
    """Every (left-flank hit, right-flank hit) pair that reads as a clean NEW or ANC.

    Hits include secondaries. A pair must be on one reference record and strand, bracket
    the midpoint like a two-record group, cover MIN_COVERAGE of the construct, pass the
    off-midpoint indel guard, and have a gap within the NEW or ANC tolerance.
    Returns [(call, group)], one per distinct pair of placements.
    """
    w, delta = C.MIDPOINT_WINDOW, M
    hits: dict[str, list[Rec]] = {"L": [], "R": []}
    for side, a in flank_alns:
        if not a.is_unmapped:
            hits[side].append(Rec.from_aln(a, flank_offset(side, a.is_reverse, delta)))
    max_dr = element_length + C.anc_tolerance(element_length) + 2 * w
    out, seen = [], set()
    for lh in hits["L"]:
        for rh in hits["R"]:
            if (lh.ref, lh.strand) != (rh.ref, rh.strand):
                continue
            a, b = sorted((lh, rh), key=lambda r: r.qs)
            if not -C.OVL_MAX - w <= b.rs - a.re <= max_dr:
                continue
            key = (lh.ref, lh.strand, lh.rs, rh.rs)
            if key in seen:
                continue
            g = evaluate_group(lh.ref, lh.strand, [lh, rh], Lc, M)
            if g.shape != "two" or g.coverage < C.MIN_COVERAGE or g.max_off_indel > C.MAX_OFFSET_INDEL:
                continue
            if abs(g.gap) <= C.TSD_TOLERANCE:
                call = "NEW"
            elif abs(g.gap - element_length) <= C.anc_tolerance(element_length):
                call = "ANC"
            else:
                continue
            seen.add(key)
            out.append((call, g))
    return out


# Construct-level AMB statuses a flank pair may resolve. All are what a repeat flank's
# primary hit at the wrong copy produces: on another record or strand
# (NO_MIDPOINT_SPAN), out of order (COMPLEX), or downstream on the same strand, where
# the two records still bracket and the "gap" is the distance to that copy
# (GAP_UNEXPLAINED). A real disruption is not rescued: its flanks pair with the same
# unexplained gap, which is not a consistent pair.
RESCUABLE = {C.AMB_NO_MIDPOINT_SPAN, C.AMB_COMPLEX, C.AMB_GAP_UNEXPLAINED}


def flank_copies(flank_alns: list[tuple[str, pysam.AlignedSegment]]) -> dict[str, int]:
    """Placements per flank: primary + secondary records. Supplementary records are
    pieces of one placement (e.g. a flank split across a linear reference's end)."""
    n = {"L": 0, "R": 0}
    for side, a in flank_alns:
        if not a.is_unmapped and not a.is_supplementary:
            n[side] += 1
    return n


def call_from_group(g: GroupResult, element_length: int) -> dict:
    """ANC/NEW fields for a clean flank-pair group (gap already within a tolerance)."""
    out = dict(best_locus=g.locus_str, construct_coverage=round(g.coverage, 4),
               identity=round(g.identity, 4), midpoint_gap=g.gap,
               anc_residual=g.gap - element_length)
    if abs(g.gap) <= C.TSD_TOLERANCE:
        site = "" if g.ref_at_mid is None else g.ref_at_mid + 1
        return dict(out, call="NEW", status=C.OK, tsd_estimate=-g.gap, reference_position=site)
    pos = end = ""
    if g.ref_span is not None:
        pos, end = g.ref_span[0] + 1, g.ref_span[1]
    return dict(out, call="ANC", status=C.OK, reference_position=pos, reference_end=end)


def same_site(a: GroupResult, b: GroupResult, w: int = C.MIDPOINT_WINDOW) -> bool:
    """Same reference record and strand, and the reference base at the construct
    midpoint within w: the insertion point for NEW, the copy's far end for ANC."""
    return ((a.ref, a.strand) == (b.ref, b.strand) and a.ref_at_mid is not None
            and b.ref_at_mid is not None and abs(a.ref_at_mid - b.ref_at_mid) <= w)


def classify_construct(alns: list[pysam.AlignedSegment], element_length: int,
                       Lc: int, M: int,
                       flank_alns: list[tuple[str, pysam.AlignedSegment]] | None = None) -> dict:
    """Classify one construct, then check it against its separately aligned flanks.

    flank_alns: [(side, record)], side "L" or "R". When given:
      * Confirmation: a construct ANC/NEW stands only if there is exactly one consistent
        flank pair, giving the same call at the same site (same_site). More than one
        pair -> AMB_MULTI_PAIRING; none, or one that disagrees -> AMB_FLANK_CONFLICT.
        A confirmed ANC reports the pair's reference span.
      * Rescue: a construct AMB in RESCUABLE is called from the flank pair when there is
        exactly one consistent pair and at least one flank has more than one placement
        (the repeat flank that explains why the construct's records missed).
    """
    out, best = _classify_records(alns, element_length, Lc, M)
    if flank_alns is None:
        return out
    pairs = flank_pairs(flank_alns, element_length, Lc, M)
    out["n_flank_pairs"] = len(pairs)
    if out["call"] in ("ANC", "NEW"):
        if len(pairs) == 1 and pairs[0][0] == out["call"] and best is not None \
                and same_site(best, pairs[0][1]):
            g = pairs[0][1]
            if out["call"] == "ANC" and g.ref_span is not None:
                # The pair's span, not the construct's: a construct record can extend
                # a flank through micro-homology into the copy and trim the span by a
                # few bp; a separately aligned flank cannot run past its own end.
                out["reference_position"], out["reference_end"] = g.ref_span[0] + 1, g.ref_span[1]
            return out
        status = C.AMB_MULTI_PAIRING if len(pairs) > 1 else C.AMB_FLANK_CONFLICT
        return dict(out, call="AMB", status=status, tsd_estimate="",
                    reference_position="", reference_end="")
    if (out["status"] in RESCUABLE and len(pairs) == 1
            and max(flank_copies(flank_alns).values()) > 1):
        return dict(out, **call_from_group(pairs[0][1], element_length))
    return out


# --------------------------------------------------------------------------- composites
#
# A newly transposed composite (IS - cargo - IS, inserted as one unit) presents each IS
# copy with one flank at the insertion site and one inside the cargo, whose source lies
# elsewhere. Each copy alone is therefore AMB (NO_MIDPOINT_SPAN, COMPLEX or
# GAP_UNEXPLAINED, depending on the cargo's orientation and where its source lies). The
# unit is found by pairing the OUTER flanks of two same-cluster copies: the first copy's
# left flank and the second copy's right flank must read as a clean NEW insertion of the
# whole span. The INNER flanks (first copy's right, second copy's left) are the cargo's
# two ends and give its source in the reference.


def outer_pair(left_hits: list[pysam.AlignedSegment], right_hits: list[pysam.AlignedSegment],
               unit_length: int) -> GroupResult | None:
    """The single NEW-consistent pairing of an IS-cargo-IS unit's outer flanks, else None.

    left_hits: alignments of the first copy's left flank; right_hits: of the second
    copy's right flank. Uses the flank-pair test with the unit's length as the element
    length, and, as the flank-pair check does, requires exactly one consistent pair.
    """
    Lc, M = 2 * C.DELTA, C.DELTA
    pairs = flank_pairs([("L", a) for a in left_hits] + [("R", a) for a in right_hits],
                        unit_length, Lc, M)
    new = [g for call, g in pairs if call == "NEW"]
    return new[0] if len(new) == 1 else None


def _flank_ends(a: pysam.AlignedSegment, flank_len: int) -> tuple[int, int]:
    """0-based reference positions of the stored sequence's first and last base,
    extrapolated over clips. Stored order ascends with the reference on both strands."""
    cig = a.cigartuples or []
    lead = cig[0][1] if cig and cig[0][0] == _HARD_CLIP else 0
    qs, qe = lead + a.query_alignment_start, lead + a.query_alignment_end
    return a.reference_start - qs, a.reference_end - 1 + (flank_len - qe)


def cargo_sources(start_hits: list[pysam.AlignedSegment], end_hits: list[pysam.AlignedSegment],
                  cargo_len: int, flank_len: int = C.DELTA) -> list[tuple[str, int, int]]:
    """Reference placements of a composite's cargo: [(ref, first, last)], 1-based.

    start_hits: alignments of the first 'flank_len' cargo bases (the first copy's right
    flank); end_hits: of the last ones (the second copy's left flank). A placement puts
    both on one record and strand, in cargo order, with an implied span within
    anc_tolerance of cargo_len. `first` is the reference base at the cargo's first base
    in query order, so first > last means the cargo is reverse-complemented.
    Cargo shorter than a flank has flanks reaching into the IS copies: no placement.
    """
    if cargo_len < flank_len:
        return []
    tol = C.anc_tolerance(cargo_len)
    out: set[tuple[str, int, int]] = set()
    for s in start_hits:
        for e in end_hits:
            if s.is_unmapped or e.is_unmapped or s.reference_name != e.reference_name \
                    or s.is_reverse != e.is_reverse:
                continue
            s_lo, s_hi = _flank_ends(s, flank_len)
            e_lo, e_hi = _flank_ends(e, flank_len)
            # '+': cargo runs s_lo .. e_hi. '-': stored is reverse-complemented, so the
            # cargo's first base is s_hi and its last is e_lo.
            first, last = (s_hi, e_lo) if s.is_reverse else (s_lo, e_hi)
            if abs(abs(last - first) + 1 - cargo_len) > tol:
                continue
            if (last < first) != s.is_reverse:
                continue
            out.add((s.reference_name, first + 1, last + 1))
    return sorted(out)


def find_composites(candidates: list, flank_alns: dict[int, list[tuple[str, pysam.AlignedSegment]]]
                    ) -> dict[int, dict]:
    """Pair AMB copies into composite units. Returns {element index: row update}.

    candidates: Elements (inputs.Element) eligible for pairing, i.e. AMB with a status in
    RESCUABLE. Copies pair within one query record and cluster. A unit (i, j), i before j,
    stands if its outer flanks pair as exactly one clean NEW (outer_pair). A copy found
    in more than one unit is left as it is.
    """
    groups: dict[tuple[str, str], list] = defaultdict(list)
    for el in candidates:
        groups[(el.record, el.cluster or el.name)].append(el)

    def hits(el, side):
        return [a for s, a in flank_alns.get(el.index, []) if s == side]

    units = []
    for els in groups.values():
        els = sorted(els, key=lambda el: el.s)
        for x, i in enumerate(els):
            for j in els[x + 1:]:
                if j.s <= i.e:
                    continue
                g = outer_pair(hits(i, "L"), hits(j, "R"), j.e - i.s + 1)
                if g is not None:
                    units.append((i, j, g))

    n_units = defaultdict(int)
    for i, j, _ in units:
        n_units[i.index] += 1
        n_units[j.index] += 1
    out: dict[int, dict] = {}
    for i, j, g in units:
        if n_units[i.index] != 1 or n_units[j.index] != 1:
            continue
        sources = cargo_sources(hits(i, "R"), hits(j, "L"), j.s - i.e - 1)
        cargo = f"{sources[0][0]}:{sources[0][1]}-{sources[0][2]}" if len(sources) == 1 else ""
        for el, partner in ((i, j), (j, i)):
            row = call_from_group(g, el.length)
            assert row["call"] == "NEW"
            out[el.index] = dict(row, call=C.COMPOSITE,
                                 composite_partner=f"{partner.name}:{partner.s}-{partner.e}",
                                 source_cargo=cargo, n_cargo_sources=len(sources))
    return out


def _classify_records(alns: list[pysam.AlignedSegment], element_length: int,
                      Lc: int, M: int) -> tuple[dict, GroupResult | None]:
    """Classify one construct from all of its alignment records (spec §4.6).
    Returns the output fields and the reported group (None if no group spans)."""
    w = C.MIDPOINT_WINDOW
    mapped = [a for a in alns if not a.is_unmapped]
    out = dict(n_records=len(mapped), n_loci=0, best_locus="", construct_coverage="",
               identity="", midpoint_gap="", anc_residual="", tsd_estimate="",
               reference_position="", reference_end="")
    if not mapped:
        return dict(out, call="AMB", status=C.AMB_NO_LOCUS), None

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
        return dict(out, call="AMB", status=status), best
    if not loci:
        return dict(out, call="AMB", status=C.AMB_LOW_COVERAGE), best
    if n_loci > 1:
        return dict(out, call="AMB", status=C.AMB_MULTI_LOCUS), best

    assert best is not None and best.gap is not None
    if best.max_off_indel > C.MAX_OFFSET_INDEL:
        return dict(out, call="AMB", status=C.AMB_INDEL_OFF_MIDPOINT), best
    if abs(best.gap) <= C.TSD_TOLERANCE:
        site = "" if best.ref_at_mid is None else best.ref_at_mid + 1   # 1-based for output
        return dict(out, call="NEW", status=C.OK, tsd_estimate=-best.gap, reference_position=site), best
    if abs(best.gap - element_length) <= C.anc_tolerance(element_length):
        pos = end = ""
        if best.ref_span is not None:
            pos, end = best.ref_span[0] + 1, best.ref_span[1]   # 1-based inclusive for output
        return dict(out, call="ANC", status=C.OK, reference_position=pos, reference_end=end), best
    return dict(out, call="AMB", status=C.AMB_GAP_UNEXPLAINED), best
