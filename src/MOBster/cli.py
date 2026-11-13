"""Command-line interface (spec §7). The §5 constants are deliberately not exposed."""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

from . import __version__
from .inputs import InputError
from .pipeline import run

logger = logging.getLogger("MOBster")


FMT = "%(asctime)s %(levelname)s %(message)s"


def _setup_logging(outdir: Path) -> None:
    """INFO+ to stderr (stdout stays clean); full DEBUG log to <outdir>/mobster.log."""
    logger.setLevel(logging.DEBUG)
    try:
        import coloredlogs
        coloredlogs.install(level=logging.INFO, logger=logger, fmt=FMT, stream=sys.stderr)
    except ImportError:
        h = logging.StreamHandler(sys.stderr)
        h.setLevel(logging.INFO)
        h.setFormatter(logging.Formatter(FMT))
        logger.addHandler(h)
    outdir.mkdir(parents=True, exist_ok=True)
    fh = logging.FileHandler(outdir / "mobster.log", mode="w")
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(logging.Formatter(FMT))
    logger.addHandler(fh)


def parse_args(argv=None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        prog="mobster",
        description="Classify IS elements in a query genome as ancestral (ANC), "
                    "novel transpositions (NEW) or ambiguous (AMB) relative to a reference. "
                    "Swap --reference and --query to report losses.")
    p.add_argument("--reference", required=True, type=Path, help="reference genome FASTA")
    p.add_argument("--query", required=True, type=Path, help="query genome FASTA")
    p.add_argument("--outdir", type=Path, default=Path("mobster_out"),
                   help="output directory (default: mobster_out)")
    p.add_argument("--threads", type=int, default=1, help="threads for ISEScan/minimap2")
    p.add_argument("--is-table-reference", type=Path,
                   help="precomputed ISEScan table for the reference; accepted but not "
                        "used in v0.1 (the reference is not annotated)")
    p.add_argument("--is-table-query", type=Path,
                   help="precomputed ISEScan table for the query (skips ISEScan)")
    p.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    return p.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    _setup_logging(args.outdir)
    try:
        run(args.reference, args.query, args.outdir, max(1, args.threads),
            args.is_table_reference, args.is_table_query)
    except InputError as exc:
        logger.error("%s", exc)
        return 2
    except RuntimeError as exc:
        logger.error("%s", exc)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
