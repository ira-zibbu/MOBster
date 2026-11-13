import pytest

from MOBster.constructs import build_construct, precheck
from MOBster import constants as C
from MOBster.inputs import Element, InputError, load_fasta, read_is_table


def test_missing_empty_and_nonfasta(tmp_path):
    with pytest.raises(InputError, match="does not exist"):
        load_fasta(tmp_path / "nope.fa", "query")
    (tmp_path / "e.fa").write_text("")
    with pytest.raises(InputError, match="empty"):
        load_fasta(tmp_path / "e.fa", "query")
    (tmp_path / "x.fa").write_text("hello\n")
    with pytest.raises(InputError, match="not FASTA"):
        load_fasta(tmp_path / "x.fa", "query")


def test_duplicate_ids(tmp_path):
    (tmp_path / "d.fa").write_text(">a\nACGT\n>a\nACGT\n")
    with pytest.raises(InputError, match="duplicate"):
        load_fasta(tmp_path / "d.fa", "query")


def test_multirecord(tmp_path):
    (tmp_path / "m.fa").write_text(">a desc\nacgt\n>b\nGG\n")
    assert load_fasta(tmp_path / "m.fa", "q") == {"a": "ACGT", "b": "GG"}


def test_read_isescan_csv(tmp_path):
    p = tmp_path / "t.csv"
    p.write_text(
        "seqID,family,cluster,isBegin,isEnd,isLen,ncopy4is,start1,end1,start2,end2,score,irId,"
        "irLen,nGaps,orfBegin,orfEnd,strand,orfLen,E-value,E-value4copy,type,ov,tir\n"
        "REL606,IS4,IS4_169,15383,16725,1343,8,15383,15408,16700,16725,52,26,26,0,0,0,,1,"
        "3.9e-260,-1,c,8,TCC:TCC\n"
        "REL606,IS110,IS110_236,22660,22815,156,1,0,0,0,0,0,0,0,0,22660,22815,-,156,7.8e-18,"
        "7.8e-18,p,1,-:-\n")
    els = read_is_table(p)
    assert els[0] == Element(0, "REL606", "IS4_169", 15383, 16725, "c", "IS4", "IS4_169")
    assert els[0].length == 1343 and els[0].complete and not els[1].complete


def test_construct_excises_element_exactly():
    # 1-based: left flank 1..5, element 6..10, right flank 11..15
    seq = "AAAAA" + "TTTTT" + "CCCCC"
    el = Element(0, "r", "IS", 6, 10, "c")
    assert build_construct(el, seq, delta=5) == "AAAAACCCCC"
    assert precheck(el, {"r": seq}, delta=5) is None
    assert precheck(el, {"r": seq}, delta=6) == C.AMB_FLANK_TRUNCATED
    assert precheck(Element(0, "r", "IS", 6, 10, "p"), {"r": seq}, delta=5) == C.AMB_PARTIAL_ELEMENT
