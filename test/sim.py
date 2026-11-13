"""Tiny random-genome simulator shared by the integration tests."""

import random

COMP = str.maketrans("ACGT", "TGCA")


def revcomp(s: str) -> str:
    return s.translate(COMP)[::-1]


def rand_seq(n: int, rng: random.Random) -> str:
    return "".join(rng.choice("ACGT") for _ in range(n))


def insert(genome: str, pos: int, element: str, tsd: int) -> tuple[str, int, int]:
    """Insert `element` with a `tsd`-bp target-site duplication.

    Follows gdtools APPLY's MOB convention: `pos` is the 1-based first duplicated
    base, and with tsd == 0 the element goes AFTER base `pos`. Returns the new
    genome and the element's 1-based inclusive TSD-excluded (s, e).
    """
    after = pos - 1 + max(tsd, 1)          # 1-based base the element follows
    new = genome[:after] + element + genome[pos - 1 if tsd else after:]
    s = after + 1
    return new, s, s + len(element) - 1
