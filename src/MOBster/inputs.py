"""FASTA validation and IS-table reading (spec §4.1, §4.2)."""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path

import pandas as pd
from Bio import SeqIO

logger = logging.getLogger("MOBster")


class InputError(Exception):
    """User-facing input problem; reported without a traceback."""


@dataclass(frozen=True)
class Element:
    """One IS element from an IS table. s/e are 1-based inclusive, TSD-excluded."""
    index: int
    record: str
    name: str
    s: int
    e: int
    type: str
    family: str = ""    # ISEScan `family` column, if present
    cluster: str = ""   # ISEScan `cluster` column, if present

    @property
    def length(self) -> int:
        return self.e - self.s + 1

    @property
    def complete(self) -> bool:
        return self.type.strip().lower() in ("c", "complete")


def load_fasta(path: str | Path, label: str) -> dict[str, str]:
    """Parse a (multi-)FASTA into {record_id: uppercase sequence}."""
    path = Path(path)
    if not path.exists():
        raise InputError(f"{label} FASTA does not exist: {path}")
    if not path.is_file():
        raise InputError(f"{label} FASTA is not a file: {path}")
    try:
        with open(path) as fh:
            head = fh.read(1)
            fh.seek(0)
            if not head:
                raise InputError(f"{label} FASTA is empty: {path}")
            if head != ">":
                raise InputError(f"{label} file is not FASTA (does not start with '>'): {path}")
            seqs: dict[str, str] = {}
            for rec in SeqIO.parse(fh, "fasta"):
                if rec.id in seqs:
                    raise InputError(f"{label} FASTA has duplicate record id '{rec.id}': {path}")
                seqs[rec.id] = str(rec.seq).upper()
    except (OSError, UnicodeDecodeError, ValueError) as exc:
        raise InputError(f"{label} FASTA could not be read ({exc}): {path}") from exc
    if not seqs:
        raise InputError(f"{label} FASTA contains no records: {path}")
    empty = [k for k, v in seqs.items() if not v]
    if empty:
        raise InputError(f"{label} FASTA has empty record(s) {empty}: {path}")
    logger.info("%s: %d record(s), %d bp total", label, len(seqs), sum(map(len, seqs.values())))
    return seqs


REQUIRED_COLUMNS = ("seqID", "isBegin", "isEnd", "type")


def read_is_table(path: str | Path) -> list[Element]:
    """Read an ISEScan 1.7.x table (CSV or TSV) into Elements, in file order.

    Coordinates stay 1-based inclusive here; conversion happens at the
    sequence-slicing boundary in constructs.py.
    """
    path = Path(path)
    if not path.exists():
        raise InputError(f"IS table does not exist: {path}")
    sep = "\t" if path.suffix.lower() in (".tsv", ".tab") else ","
    try:
        df = pd.read_csv(path, sep=sep, dtype=str, keep_default_na=False)
    except pd.errors.EmptyDataError:
        logger.warning("IS table is empty: %s", path)
        return []
    df.columns = [c.strip() for c in df.columns]
    missing = [c for c in REQUIRED_COLUMNS if c not in df.columns]
    if missing:
        raise InputError(f"IS table {path} lacks column(s) {missing}; found {list(df.columns)}")
    name_col = "cluster" if "cluster" in df.columns else ("family" if "family" in df.columns else None)

    elements = []
    for i, row in enumerate(df.itertuples(index=False)):
        r = row._asdict()
        try:
            s, e = int(r["isBegin"]), int(r["isEnd"])
        except ValueError as exc:
            raise InputError(f"IS table {path} row {i + 2}: non-integer coordinate") from exc
        if s < 1 or e < s:
            raise InputError(f"IS table {path} row {i + 2}: bad coordinates {s}..{e}")
        name = r[name_col] if name_col else f"IS_{i}"
        elements.append(Element(i, r["seqID"].strip(), name.strip(), s, e, r["type"].strip(),
                                r.get("family", "").strip(), r.get("cluster", "").strip()))
    logger.info("Read %d element(s) from %s", len(elements), path)
    return elements
