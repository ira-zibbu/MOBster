"""Align all constructs to the reference in one minimap2 call (spec §4.5)."""

from __future__ import annotations

import logging
import shutil
import subprocess
from pathlib import Path

from . import constants as C
from .inputs import InputError

logger = logging.getLogger("MOBster")


def minimap2_cmd(reference: Path, constructs: Path, threads: int) -> list[str]:
    # Preset FIRST (later flags override it). -a is mandatory: SAM + base-level CIGAR.
    #
    # DEVIATION from spec §4.5 ("v0.1 overrides nothing"): -r CHAIN_BANDWIDTH (both values).
    # asm5 sets -r 1k,100k. Either bandwidth lets chaining join an ancestral construct's
    # two flanks across the element-length reference gap (the 1k bandwidth for elements
    # < 1 kb such as IS1; the RMQ long-join for longer ones). Base-level alignment then
    # cannot pay for the deletion (asm5 gap cost ~81+L vs +1 per matched base, so no
    # δ=500 flank can bridge an element > ~420 bp), truncates, and emits ONE flank only
    # (e.g. 500S500M) with NO supplementary record for the other. Every ancestral
    # element became AMB_NO_MIDPOINT_SPAN. With a bandwidth below the shortest element,
    # each flank is its own chain -> primary + supplementary: the two-record shape of
    # §4.6. Regression tests: test_pipeline.py (1347 bp and 768 bp ancestral elements).
    bw = str(C.MM2_CHAIN_BANDWIDTH)
    return ["minimap2", "-x", "asm5", "-r", f"{bw},{bw}", "-a", "-t", str(threads),
            str(reference), str(constructs)]


def flank_minimap2_cmd(reference: Path, flanks: Path, threads: int) -> list[str]:
    # Each flank aligned on its own, keeping up to FLANK_MAX_SECONDARY secondary hits
    # (asm5 default -N 5 is fewer than the copy number of common IS families, so the
    # copy adjacent to the element could be dropped). -p stays at the preset default:
    # repeat copies score alike, unrelated hits do not.
    return ["minimap2", "-x", "asm5", "-N", str(C.FLANK_MAX_SECONDARY), "-a", "-t", str(threads),
            str(reference), str(flanks)]


def run_minimap2(reference: Path, constructs: Path, sam: Path, threads: int,
                 flanks: bool = False) -> Path:
    if shutil.which("minimap2") is None:
        raise InputError("minimap2 not found on PATH (activate the 'mobster' conda env)")
    cmd = (flank_minimap2_cmd if flanks else minimap2_cmd)(reference, constructs, threads)
    logger.info("Aligning: %s > %s", " ".join(cmd), sam)
    with open(sam, "w") as out:
        proc = subprocess.run(cmd, stdout=out, stderr=subprocess.PIPE, text=True)
    for line in proc.stderr.splitlines():
        logger.debug("[minimap2] %s", line)
    if proc.returncode != 0:
        raise RuntimeError(f"minimap2 failed (exit {proc.returncode}):\n{proc.stderr}")
    return sam
