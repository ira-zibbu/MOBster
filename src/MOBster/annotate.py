"""Run ISEScan and locate its output table (spec §4.2).

ISEScan 1.7.3 (pinned in environment.yml) writes `<output>/<seqfile path>.csv`,
i.e. it nests under --output mirroring the *path given as --seqfile*. We pass a
bare filename by running from the FASTA's directory, so the table lands at
`<output>/<fasta basename>.csv`, and we still glob as a fallback.
"""

from __future__ import annotations

import logging
import shutil
import subprocess
from pathlib import Path

from .inputs import InputError

logger = logging.getLogger("MOBster")


def run_isescan(fasta: Path, outdir: Path, threads: int) -> Path | None:
    """Annotate `fasta`; return the path to the CSV table, or None if ISEScan found nothing."""
    exe = shutil.which("isescan.py")
    if exe is None:
        raise InputError("isescan.py not found on PATH (activate the 'mobster' conda env, "
                         "or pass --is-table-query)")
    fasta = fasta.resolve()
    outdir = outdir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    cmd = [exe, "--seqfile", fasta.name, "--output", str(outdir), "--nthread", str(threads)]
    logger.info("Running ISEScan on %s (slow step): %s", fasta.name, " ".join(cmd))
    proc = subprocess.Popen(cmd, cwd=fasta.parent, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, text=True)
    assert proc.stdout is not None
    for line in proc.stdout:
        logger.debug("[isescan] %s", line.rstrip())
    if proc.wait() != 0:
        raise RuntimeError(f"ISEScan failed on {fasta} (exit {proc.returncode})")

    expected = outdir / f"{fasta.name}.csv"
    if expected.exists():
        return expected
    hits = sorted(outdir.rglob(f"{fasta.name}.csv"))
    if len(hits) == 1:
        logger.warning("ISEScan table not at expected path %s; using %s", expected, hits[0])
        return hits[0]
    if len(hits) > 1:
        raise RuntimeError(f"Multiple ISEScan tables for {fasta.name} under {outdir}: {hits}")
    logger.warning("ISEScan produced no table for %s (no IS elements found?)", fasta.name)
    return None
