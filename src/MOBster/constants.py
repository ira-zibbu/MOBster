"""Fixed v0.1 constants (spec §5). Deliberately NOT exposed on the CLI."""

DELTA = 500             # flank length per side
MIDPOINT_WINDOW = 25    # half-width w of the midpoint window, construct bases
OVL_MAX = 20            # permitted construct/reference overlap when joining two records
MIN_COVERAGE = 0.80     # fraction of Lc a locus must cover
TSD_TOLERANCE = 25      # |midpoint_gap| <= this -> NEW
MAX_OFFSET_INDEL = 50   # indel > this outside the window -> AMB_INDEL_OFF_MIDPOINT


def anc_tolerance(element_length: int) -> float:
    """Both absolute and relative (spec §5)."""
    return max(50.0, 0.05 * element_length)


# Reason codes (spec §6). AMB_PARTIAL_ELEMENT is not emitted: partial elements are
# classified like complete ones (see constructs.precheck). AMB_LOW_COVERAGE is an addition: a group spans the
# midpoint but covers < MIN_COVERAGE of the construct, so no locus exists even
# though the midpoint was spanned (neither NO_LOCUS nor NO_MIDPOINT_SPAN is true).
OK = "OK"
AMB_FLANK_TRUNCATED = "AMB_FLANK_TRUNCATED"
AMB_NO_LOCUS = "AMB_NO_LOCUS"
AMB_NO_MIDPOINT_SPAN = "AMB_NO_MIDPOINT_SPAN"
AMB_LOW_COVERAGE = "AMB_LOW_COVERAGE"
AMB_MULTI_LOCUS = "AMB_MULTI_LOCUS"
AMB_GAP_UNEXPLAINED = "AMB_GAP_UNEXPLAINED"
AMB_INDEL_OFF_MIDPOINT = "AMB_INDEL_OFF_MIDPOINT"
AMB_COMPLEX = "AMB_COMPLEX"
# Additions (flank-pair check, classify.classify_construct). For a construct ANC/NEW,
# the flanks aligned separately pair consistently (a clean ANC or NEW reading):
#   at more than one place                                  -> AMB_MULTI_PAIRING
#   nowhere, or once but with another call or another site  -> AMB_FLANK_CONFLICT
AMB_MULTI_PAIRING = "AMB_MULTI_PAIRING"
AMB_FLANK_CONFLICT = "AMB_FLANK_CONFLICT"

# Call value for an IS copy that is one end of a newly transposed composite unit
# (IS - cargo - IS inserted as one), found by pairing the outer flanks of two AMB copies
# (classify.find_composites). An addition: the spec's calls are ANC, NEW and AMB, and it
# lists composites as a non-goal (§9).
COMPOSITE = "COMPOSITE"

# Secondary hits kept per separately aligned flank (minimap2 -N).
FLANK_MAX_SECONDARY = 100


# minimap2 chaining bandwidth (-r INT,INT), overriding asm5's 1k,100k. Not in the
# spec's §5 table: added because asm5's bandwidth joins ancestral flanks across the
# element and loses one flank (see align.py). Must be below the shortest element
# length and above the indels we expect to stay inside one record.
MM2_CHAIN_BANDWIDTH = 100
assert TSD_TOLERANCE < MAX_OFFSET_INDEL < MM2_CHAIN_BANDWIDTH
