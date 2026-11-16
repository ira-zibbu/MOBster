"""Callable-element selection and construct building (spec §4.3, §4.4).

A construct is left flank [s-δ, s-1] + right flank [e+1, e+δ] (1-based inclusive),
i.e. the query locus with the element excised. Lc = 2δ, midpoint M = δ.
"""

from __future__ import annotations

import logging
from pathlib import Path

from . import constants as C
from .inputs import Element, InputError

logger = logging.getLogger("MOBster")


def construct_name(el: Element) -> str:
    # Leading index is the key; record/name are for humans (record ids may contain '|').
    return f"{el.index}|{el.record}|{el.name}"


def index_from_name(qname: str) -> int:
    return int(qname.split("|", 1)[0])


def precheck(el: Element, query: dict[str, str], delta: int = C.DELTA) -> str | None:
    """Return an AMB status if the element is not callable, else None.

    Partial elements are callable (deviation from spec §4.3): the construct excises
    exactly the annotated span, so ANC/NEW is decided against that span's length. The
    `element_type` column keeps them identifiable.
    """
    seq = query[el.record]
    if el.s - delta < 1 or el.e + delta > len(seq):
        return C.AMB_FLANK_TRUNCATED
    return None


def build_construct(el: Element, seq: str, delta: int = C.DELTA) -> str:
    # 1-based inclusive [s-δ, s-1] -> 0-based half-open [s-δ-1, s-1)
    left = seq[el.s - delta - 1: el.s - 1]
    # 1-based inclusive [e+1, e+δ] -> 0-based half-open [e, e+δ)
    right = seq[el.e: el.e + delta]
    assert len(left) == delta and len(right) == delta
    return left + right


def flank_name(el: Element, side: str) -> str:
    # Same leading index as the construct, so index_from_name works on both.
    return f"{el.index}|{side}|{el.record}|{el.name}"


def side_from_name(qname: str) -> str:
    return qname.split("|", 2)[1]


def write_constructs(elements: list[Element], query: dict[str, str], path: Path,
                     delta: int = C.DELTA, flanks_path: Path | None = None) -> dict[int, str]:
    """Write callable constructs to one FASTA. Returns {index: precheck status} for non-callable.

    With flanks_path, also write each construct's two flanks as separate sequences
    (`<index>|L|...` and `<index>|R|...`) for the flank-pair check in classify.py.
    """
    unknown = sorted({el.record for el in elements} - set(query))
    if unknown:
        raise InputError(f"IS table refers to record(s) not in the query FASTA: {unknown}")
    skipped: dict[int, str] = {}
    n = 0
    flanks = open(flanks_path, "w") if flanks_path is not None else None
    try:
        with open(path, "w") as fh:
            for el in elements:
                status = precheck(el, query, delta)
                if status:
                    skipped[el.index] = status
                    continue
                construct = build_construct(el, query[el.record], delta)
                fh.write(f">{construct_name(el)}\n{construct}\n")
                if flanks is not None:
                    flanks.write(f">{flank_name(el, 'L')}\n{construct[:delta]}\n"
                                 f">{flank_name(el, 'R')}\n{construct[delta:]}\n")
                n += 1
    finally:
        if flanks is not None:
            flanks.close()
    logger.info("Wrote %d construct(s) to %s; %d element(s) not callable", n, path, len(skipped))
    return skipped
