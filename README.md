# MOBster

Detect transpositions of IS elements between two bacterial genomes. Given a
**reference** and a **query** assembly, MOBster labels every IS element in the query
as ancestral (already present in the reference), new (inserted since), or ambiguous.

## Installation

```bash
git clone https://github.com/ira-zibbu/MOBster.git && cd MOBster
mamba env create -f environment.yml     # ISEScan 1.7.3, minimap2, pysam, breseq (gdtools)
mamba activate mobster
pip install -e .
```
## Running

```bash
mobster --reference R.fasta --query Q.fasta --outdir results/ --threads 8
```

| Option | Explanation |
|---|---|
| `--reference` | reference genome FASTA (required) |
| `--query` | query genome FASTA (required) |
| `--outdir` | output directory (default `mobster_out`) |
| `--threads` | threads for ISEScan and minimap2 (default 1) |
| `--is-table-query` | precomputed ISEScan table for the query; skips ISEScan on it |
| `--is-table-reference` | precomputed ISEScan table for the reference; skips ISEScan on it |

**Losses.** MOBster reports elements present in the query, and currently does not predict lost elements.

## Method

1. **Annotate.** ISEScan finds IS elements in both genomes. 
2. **Build constructs.** For each query element, take 500 bp of sequence on each side
   and join them with the element cut out (a 1,000 bp *construct*).
3. **Align.** Align all constructs to the reference with minimap2. Separately, align
   each 500 bp flank on its own. This provides a second source of information.
4. **Read the gap.** Where the construct aligns across its midpoint (either as a single alignment or split alignment), measure how much
   reference sequence lies between the two flanks (`midpoint_gap`):
   - **about 0** → the reference has an empty site; the element is **NEW**. A small
     negative gap is the target-site duplication (TSD) created by the insertion.
   - **about the element's length** → the reference has the element at that site;
     it is **ANC**.
   - anything else, or no clean alignment → **AMB**, with a code explaining why.
5. **Confirm with the flanks.** An ANC or NEW call stands only if the separately
   aligned flanks pair up at exactly one reference site, giving the same answer. When
   a flank is a repeat (for example, it runs into a neighbouring IS copy), the flank
   pairs can also rescue a call the construct missed. This also prevents false **NEW** calls if an ancestral element is present inside a repeat region, such that one copy of the repeat is element-free, and one copy has the element. 
6. **Composite transpositions.** When two copies of the same IS moved together with
   the DNA between them (IS – cargo – IS), each copy alone looks ambiguous. If the
   outer flanks of the pair read as one clean insertion, both copies are called
   **COMPOSITE**, and the cargo's source in the reference is reported.

Partial elements (ISEScan type `p`) are also included.

## Output

Files in `--outdir`:

| File | Contents |
|---|---|
| `mobster_calls.csv` | one row per query IS element (described below) |
| `mobster_map.html` | interactive circular map of the query and its calls |
| `mobster_map_paired.html` | reference and query on one circle, with ribbons linking composite cargo to its reference source |
| `mobster.log` | full debug log |
| `isescan/` | ISEScan output for each genome it was run on |
| `work/` | constructs, flanks and their alignments (SAM) |

The maps are single self-contained HTML files that open offline. Hover over a glyph for
its details.

### Calls (`call`)

| Call | Meaning |
|---|---|
| `NEW` | Inserted in the query; the reference has an empty site here. |
| `ANC` | Ancestral; the reference has a copy at the same site. |
| `COMPOSITE` | One of the two IS copies flanking a newly transposed composite unit (IS – cargo – IS). |
| `AMB` | Could not be decided. see `status` for why. |

### Status (`status`)

`OK` for every NEW, ANC and COMPOSITE call. For AMB, the `status` indicates why the call was ambiguous:

| Status | Meaning |
|---|---|
| `AMB_FLANK_TRUNCATED` | The element is within 500 bp of a contig end, so a full flank cannot be taken. |
| `AMB_NO_LOCUS` | The construct did not align to the reference at all. |
| `AMB_NO_MIDPOINT_SPAN` | The construct aligned, but no alignment crosses the element site. Typical of a site that was rearranged through an SV like an inversion or duplication, or of a flank that lies inside another repeat. |
| `AMB_LOW_COVERAGE` | An alignment crosses the site but covers less than 80% of the construct. Indicates that the site is highly diverged in sequence. |
| `AMB_MULTI_LOCUS` | The construct aligns cleanly at more than one reference site. |
| `AMB_GAP_UNEXPLAINED` | One clean site, but the gap between the flanks matches neither an empty site nor the element's length (e.g. a different element, or a deletion, at that site). |
| `AMB_INDEL_OFF_MIDPOINT` | A large indel (> 50 bp) in a flank, away from the element site. |
| `AMB_COMPLEX` | The alignment is split into pieces that don't simply bracket the site. |
| `AMB_MULTI_PAIRING` | The construct gave ANC/NEW, but the separately aligned flanks pair up at more than one reference site. |
| `AMB_FLANK_CONFLICT` | The construct gave ANC/NEW, but the separately aligned flanks don't agree (no pair, or a pair with another answer or site). |

### Columns of `mobster_calls.csv`

Coordinates are 1-based and inclusive.

| Column | Description |
|---|---|
| `query_record` | Query contig/chromosome |
| `element_name` | ISEScan element name |
| `element_family`, `element_cluster` | ISEScan family and cluster |
| `element_start`, `element_end` | Element position in the query (TSD excluded) |
| `element_length` | `element_end − element_start + 1` |
| `element_type` | ISEScan type: `c` complete, `p` partial |
| `delta` | Flank length used (default: 500) |
| `call` | `NEW`, `ANC`, `COMPOSITE` or `AMB` (see above) |
| `status` | `OK`, or the AMB reason (see above) |
| `n_records` | Number of alignment records for the construct |
| `n_loci` | Number of reference sites where the construct aligns across the midpoint with enough coverage |
| `n_flank_pairs` | Number of reference sites where the two separately aligned flanks pair consistently |
| `best_locus` | Reference site used: `record:start-end:strand` |
| `construct_coverage` | Fraction of the construct covered at `best_locus` |
| `identity` | Alignment identity at `best_locus` |
| `midpoint_gap` | Reference bases between the two flanks. ≈ 0 for NEW (negative = TSD length), ≈ `element_length` for ANC. Given for AMB calls too where an alignment exists. |
| `anc_residual` | `midpoint_gap − element_length`; ≈ 0 for a clean ANC |
| `tsd_estimate` | Estimated TSD length (`−midpoint_gap`); NEW and COMPOSITE only |
| `reference_position` | NEW / COMPOSITE: insertion point in the reference (first base of the TSD). ANC: start of the reference copy. |
| `reference_end` | ANC only: end of the reference copy |
| `composite_partner` | COMPOSITE only: the other IS copy of the unit, `record:start-end` |
| `n_cargo_sources` | COMPOSITE only: number of places in the reference the cargo maps to |
| `source_cargo` | COMPOSITE only, when the cargo maps to exactly one place: its reference span `record:start-end`. Start > end means the cargo is inverted relative to the reference. |
