"""GD MOB -> (name, s, e) conversion: the harness's most error-prone step (spec §8).

Pure arithmetic tests, then the same mutations applied by real gdtools to a tiny
synthetic annotated GenBank, checking the element sequence and TSDs at the
predicted coordinates.
"""

import random
import shutil
import subprocess

import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

from benchmark.gdcoords import Mut, existing_element, map_pos, mob_element, write_gd
from sim import insert, rand_seq, revcomp


def test_single_mob_matches_simulator():
    rng = random.Random(1)
    g, el = rand_seq(5000, rng), rand_seq(700, rng)
    for dup in (0, 3, 9):
        _, s, e = insert(g, 2000, el, dup)
        m = Mut("MOB", 2000, size=700, dup=dup, name="ISX")
        start = 2000 + max(dup, 1)
        assert mob_element(m, [m]) == (s, e) == (start, start + 699)


def test_upstream_mutations_shift():
    m = Mut("MOB", 5000, size=1000, dup=5, name="ISX")
    d = Mut("DEL", 1000, size=10)
    i = Mut("INS", 2000, seq="ACGT")
    up = Mut("MOB", 3000, size=768, dup=9, name="IS1")
    s, e = mob_element(m, [m, d, i, up])
    assert s == 5000 + 5 - 10 + 4 + 768 + 9
    assert e - s + 1 == 1000
    # downstream mutations do not shift
    assert mob_element(m, [m, Mut("DEL", 9000, size=50)]) == (5005, 6004)


def test_mob_inside_existing_element_extends_it():
    inner = Mut("MOB", 1500, size=1443, dup=3, name="IS150")
    assert existing_element(1000, 1768, [inner]) == (1000, 1768 + 1443 + 3)


def test_deleted_position_raises():
    with pytest.raises(ValueError):
        map_pos(105, [Mut("DEL", 100, size=10)])


@pytest.mark.skipif(shutil.which("gdtools") is None, reason="gdtools not on PATH")
def test_against_gdtools(tmp_path):
    rng = random.Random(7)
    L = 777
    elem = rand_seq(L, rng)
    pre, post = rand_seq(3000, rng), rand_seq(16000, rng)
    seq = pre + elem + post                           # annotated copy at 3001..3777
    rec = SeqRecord(Seq(seq), id="SIM", name="SIM", description="sim",
                    annotations={"molecule_type": "DNA", "topology": "linear"})
    rec.features.append(SeqFeature(FeatureLocation(3000, 3000 + L, strand=1), type="repeat_region",
                                   qualifiers={"mobile_element": ["insertion sequence:ISX"]}))
    SeqIO.write(rec, tmp_path / "sim.gbk", "genbank")

    mobs = [Mut("MOB", 6000, size=L, dup=4, strand=1, name="ISX"),
            Mut("MOB", 9000, size=L, dup=0, strand=-1, name="ISX"),
            Mut("MOB", 3400, size=L, dup=6, strand=1, name="ISX"),   # inside the existing copy
            Mut("MOB", 9300, size=L, dup=9, strand=-1, name="ISX")]  # within δ of the previous
    other = [Mut("DEL", 5000, size=12), Mut("SNP", 5500, seq="A" if seq[5499] != "A" else "C"),
             Mut("INS", 8000, seq="GATTACA")]
    muts = mobs + other
    write_gd(muts, "SIM", tmp_path / "m.gd")
    subprocess.run(["gdtools", "APPLY", "-r", str(tmp_path / "sim.gbk"), "-f", "FASTA",
                    "-o", str(tmp_path / "q.fa"), str(tmp_path / "m.gd")],
                   check=True, capture_output=True)
    q = str(next(SeqIO.parse(tmp_path / "q.fa", "fasta")).seq).upper()
    assert len(q) == len(seq) + sum(m.net for m in muts)

    for m in mobs:
        s, e = mob_element(m, muts)
        got = q[s - 1:e]
        assert got == (elem if m.strand == 1 else revcomp(elem)), m
        if m.dup:
            assert q[s - 1 - m.dup:s - 1] == q[e:e + m.dup] == seq[m.pos - 1:m.pos - 1 + m.dup], m
    a, b = existing_element(3001, 3000 + L, muts)
    assert q[a - 1:a + 99] == elem[:100] and q[b - 100:b] == elem[-100:]
    assert b - a + 1 == L + L + 6
