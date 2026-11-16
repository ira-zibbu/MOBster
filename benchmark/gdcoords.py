"""GenomeDiff mutations -> element coordinates in the mutated genome (spec §8).

gdtools APPLY semantics, verified empirically against breseq 0.39.0
(test/test_gdcoords.py re-checks this against a real gdtools run):

    MOB <pos> <repeat> <strand> <dup>, dup > 0:
        new = ref[:pos-1+dup] + element(strand) + ref[pos-1:]
      `pos` is the first duplicated base; the TSD ref[pos..pos+dup-1] appears once
      before and once after the element: s = pos + dup, e = s + L - 1.
    MOB with dup == 0:
        new = ref[:pos] + element(strand) + ref[pos:]
      i.e. the element goes AFTER base `pos` (not before it): s = pos + 1.
    Both cases: the element follows base  pos - 1 + max(dup, 1).

All positions in a GD refer to the ORIGINAL reference, so every coordinate must
also be shifted by the net length change of mutations strictly upstream of it.
"""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass(frozen=True)
class Mut:
    kind: str           # SNP, DEL, INS, MOB, INV
    pos: int            # 1-based, original reference coordinates
    size: int = 0       # DEL/INV size; MOB element length L; INS length
    dup: int = 0        # MOB duplication_size
    strand: int = 1     # MOB strand
    name: str = ""      # MOB repeat_name
    seq: str = ""       # SNP/INS new sequence

    @property
    def net(self) -> int:
        """Net length change this mutation introduces."""
        return {"MOB": self.size + self.dup, "DEL": -self.size, "INS": len(self.seq)}.get(self.kind, 0)

    def gd_fields(self) -> list[str]:
        if self.kind == "MOB":
            return [self.name, str(self.strand), str(self.dup)]
        if self.kind in ("DEL", "INV"):
            return [str(self.size)]
        if self.kind in ("SNP", "INS"):
            return [self.seq]
        raise ValueError(self.kind)


def _mob_after(m: Mut) -> int:
    """Original-reference base that the MOB's element is inserted after."""
    return m.pos - 1 + max(m.dup, 1)


def map_pos(x: int, muts: list[Mut]) -> int:
    """Map an original-reference base x (1-based) to the mutated genome.

    For a base inside a MOB's TSD, returns the coordinate of the FIRST copy.
    Raises if x lies inside a deleted segment.
    """
    shift = 0
    for m in muts:
        if m.kind == "MOB":
            if x > _mob_after(m):           # past the element (and the TSD's first copy)
                shift += m.net
        elif m.kind == "DEL":
            if m.pos <= x < m.pos + m.size:
                raise ValueError(f"position {x} is deleted by {m}")
            if x >= m.pos + m.size:
                shift -= m.size
        elif m.kind == "INS":               # breseq INS inserts AFTER pos
            if x > m.pos:
                shift += m.net
    return x + shift


def mob_element(m: Mut, muts: list[Mut]) -> tuple[int, int]:
    """(s, e) of the element created by MOB `m`, 1-based inclusive, TSD-excluded."""
    assert m.kind == "MOB"
    others = [o for o in muts if o is not m]
    s = map_pos(_mob_after(m), others) + 1
    return s, s + m.size - 1


def existing_element(a: int, b: int, muts: list[Mut]) -> tuple[int, int]:
    """New span of an original element [a, b] (1-based incl.), allowing MOBs inside it."""
    return map_pos(a, muts), map_pos(b, muts)


def write_gd(muts: list[Mut], seq_id: str, path) -> None:
    lines = ["#=GENOME_DIFF\t1.0"]
    for i, m in enumerate(sorted(muts, key=lambda m: m.pos), 1):
        lines.append("\t".join([m.kind, str(i), ".", seq_id, str(m.pos), *m.gd_fields()]))
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")
