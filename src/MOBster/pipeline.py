"""End-to-end pipeline (spec §4): validate -> annotate -> select -> construct -> align -> classify -> write."""

from __future__ import annotations

import csv
import logging
from collections import Counter, defaultdict
from pathlib import Path

import pysam

from . import constants as C
from .align import run_minimap2
from .annotate import run_isescan
from .classify import RESCUABLE, classify_construct, find_composites
from .constructs import index_from_name, side_from_name, write_constructs
from .inputs import Element, InputError, load_fasta, read_is_table
from .viz import read_calls, write_map, write_paired_map

logger = logging.getLogger("MOBster")

COLUMNS = [
    "query_record", "element_name", "element_family", "element_cluster", "element_start", "element_end", "element_length",
    "element_type", "delta", "call", "status", "n_records", "n_loci", "n_flank_pairs", "best_locus",
    "construct_coverage", "identity", "midpoint_gap", "anc_residual", "tsd_estimate",
    "reference_position", "reference_end", "composite_partner", "source_cargo", "n_cargo_sources",
]


def read_alignments(sam: Path) -> dict[int, list[pysam.AlignedSegment]]:
    by_construct: dict[int, list[pysam.AlignedSegment]] = defaultdict(list)
    with pysam.AlignmentFile(str(sam), "r") as fh:
        for aln in fh:
            by_construct[index_from_name(aln.query_name)].append(aln)
    return by_construct


def read_flank_alignments(sam: Path) -> dict[int, list[tuple[str, pysam.AlignedSegment]]]:
    by_construct: dict[int, list[tuple[str, pysam.AlignedSegment]]] = defaultdict(list)
    with pysam.AlignmentFile(str(sam), "r") as fh:
        for aln in fh:
            by_construct[index_from_name(aln.query_name)].append((side_from_name(aln.query_name), aln))
    return by_construct


def classify_elements(elements: list[Element], query: dict[str, str], reference: Path,
                      workdir: Path, threads: int) -> list[dict]:
    workdir.mkdir(parents=True, exist_ok=True)
    constructs_fa, flanks_fa = workdir / "constructs.fasta", workdir / "flanks.fasta"
    skipped = write_constructs(elements, query, constructs_fa, flanks_path=flanks_fa)
    alignments: dict[int, list] = {}
    flank_alns: dict[int, list] = {}
    if len(skipped) < len(elements):
        sam = run_minimap2(reference, constructs_fa, workdir / "constructs.sam", threads)
        alignments = read_alignments(sam)
        fsam = run_minimap2(reference, flanks_fa, workdir / "flanks.sam", threads, flanks=True)
        flank_alns = read_flank_alignments(fsam)

    Lc, M = 2 * C.DELTA, C.DELTA
    rows = []
    for el in elements:
        row = dict(query_record=el.record, element_name=el.name, element_family=el.family,
                   element_cluster=el.cluster, element_start=el.s,
                   element_end=el.e, element_length=el.length, element_type=el.type,
                   delta=C.DELTA)
        if el.index in skipped:
            row.update(call="AMB", status=skipped[el.index])
        else:
            row.update(classify_construct(alignments.get(el.index, []), el.length, Lc, M,
                                          flank_alns.get(el.index, [])))
        rows.append({k: row.get(k, "") for k in COLUMNS})

    # Composite units: pair AMB copies whose outer flanks read as one NEW insertion.
    by_index = {el.index: row for el, row in zip(elements, rows)}
    candidates = [el for el in elements if el.index not in skipped
                  and by_index[el.index]["call"] == "AMB" and by_index[el.index]["status"] in RESCUABLE]
    for idx, update in find_composites(candidates, flank_alns).items():
        by_index[idx].update({k: v for k, v in update.items() if k in COLUMNS})
    return rows


def write_table(rows: list[dict], path: Path) -> None:
    with open(path, "w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=COLUMNS)
        wr.writeheader()
        wr.writerows(rows)


def summarise(rows: list[dict]) -> None:
    calls = Counter(r["call"] for r in rows)
    logger.info("Calls: %s", ", ".join(f"{k}={calls.get(k, 0)}" for k in ("ANC", "NEW", C.COMPOSITE, "AMB")))
    amb = Counter(r["status"] for r in rows if r["call"] == "AMB")
    for status, n in amb.most_common():
        logger.info("  %-24s %d", status, n)


def run(reference: Path, query: Path, outdir: Path, threads: int,
        is_table_reference: Path | None, is_table_query: Path | None) -> Path:
    outdir.mkdir(parents=True, exist_ok=True)
    reference_seqs = load_fasta(reference, "reference")
    query_seqs = load_fasta(query, "query")

    # Both genomes are annotated (spec §4.2), but classification never uses the reference
    # table (§4.3): it only marks the reference's IS elements on the paired map.
    if is_table_query is None:
        is_table_query = run_isescan(query, outdir / "isescan" / "query", threads)
    if is_table_reference is None:
        is_table_reference = run_isescan(reference, outdir / "isescan" / "reference", threads)
    reference_elements = read_is_table(is_table_reference) if is_table_reference else []
    unknown = sorted({el.record for el in reference_elements} - set(reference_seqs))
    if unknown:
        raise InputError(f"reference IS table refers to record(s) not in the reference FASTA: {unknown}")

    elements = read_is_table(is_table_query) if is_table_query else []
    rows = classify_elements(elements, query_seqs, reference, outdir / "work", threads)
    table = outdir / "mobster_calls.csv"
    write_table(rows, table)
    summarise(rows)
    logger.info("Wrote %s", table)
    # Interactive circular map (docs/viz_instructions.md). Uses the written rows so the
    # map shows exactly what is in the table.
    calls = read_calls(table)
    query_lengths = {k: len(v) for k, v in query_seqs.items()}
    write_map(calls, query, reference, outdir / "mobster_map.html", record_lengths=query_lengths)
    write_paired_map(calls, query, reference, reference_elements, outdir / "mobster_map_paired.html",
                     query_lengths=query_lengths,
                     reference_lengths={k: len(v) for k, v in reference_seqs.items()})
    return table
