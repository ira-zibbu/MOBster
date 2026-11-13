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


# Reason codes (spec §6). AMB_LOW_COVERAGE is an addition: a group spans the
# midpoint but covers < MIN_COVERAGE of the construct, so no locus exists even
# though the midpoint was spanned (neither NO_LOCUS nor NO_MIDPOINT_SPAN is true).
OK = "OK"
AMB_PARTIAL_ELEMENT = "AMB_PARTIAL_ELEMENT"
AMB_FLANK_TRUNCATED = "AMB_FLANK_TRUNCATED"
AMB_NO_LOCUS = "AMB_NO_LOCUS"
AMB_NO_MIDPOINT_SPAN = "AMB_NO_MIDPOINT_SPAN"
AMB_LOW_COVERAGE = "AMB_LOW_COVERAGE"
AMB_MULTI_LOCUS = "AMB_MULTI_LOCUS"
AMB_GAP_UNEXPLAINED = "AMB_GAP_UNEXPLAINED"
AMB_INDEL_OFF_MIDPOINT = "AMB_INDEL_OFF_MIDPOINT"
AMB_COMPLEX = "AMB_COMPLEX"


# minimap2 chaining bandwidth (-r INT,INT), overriding asm5's 1k,100k. Not in the
# spec's §5 table: added because asm5's bandwidth joins ancestral flanks across the
# element and loses one flank (see align.py). Must be below the shortest element
# length and above the indels we expect to stay inside one record.
MM2_CHAIN_BANDWIDTH = 100
assert TSD_TOLERANCE < MAX_OFFSET_INDEL < MM2_CHAIN_BANDWIDTH
