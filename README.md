# MOBster

Detect transposition of mobile genetic elements (IS elements) in bacterial genomes.
Implements `docs/instructions.md` (v0.1).

For each IS element in a **query** genome, MOBster builds a construct (left flank +
right flank, with the element excised), aligns it to the **reference**, and reads the
signed indel at the construct midpoint:

| `midpoint_gap` | meaning | call |
|---|---|---|
| ≈ 0 (≤ 25 bp, negative = duplicated TSD) | empty site in the reference | `NEW` |
| ≈ element length | element present in the reference | `ANC` |
| anything else / unresolvable | — | `AMB` + reason code |

Swap `--reference` and `--query` to report losses.

## Install

```bash
mamba env create -f environment.yml     # ISEScan 1.7.3, minimap2 2.28, pysam, breseq (gdtools)
mamba activate mobster
pip install -e .
```

## Run

```bash
mobster --reference R.fasta --query Q.fasta --outdir results/ --threads 8
# skip ISEScan with a precomputed query table (ISEScan CSV/TSV format):
mobster --reference R.fasta --query Q.fasta --outdir results/ --is-table-query Q.fasta.csv

# element losses: swap the genomes (ISEScan then annotates R, as the query)
mobster --reference Q.fasta --query R.fasta --outdir results_losses/
```

Outputs in `--outdir`:

- `mobster_calls.csv`: one row per query element (spec §6 schema)
- `mobster_map.html`: interactive circular map of the calls (see below)
- `mobster.log`: full debug log (stderr shows INFO and above; stdout is unused)
- `work/constructs.fasta`, `work/constructs.sam`: the constructs and their alignments
- `isescan/query/`: ISEScan output for the query, when it was run (the reference is
  never annotated)

For NEW calls, `reference_position` is the 1-based reference base aligned to the construct
midpoint. Because minimap2 left-aligns indels, this is the first duplicated base of the
TSD, the same convention as a GenomeDiff `MOB` position.

## Interactive map

Every run also writes `mobster_map.html` (spec: `docs/viz_instructions.md`). It's a single
self-contained file, around 570 KB, with the data and both libraries inlined, so it opens
offline and can be emailed. One arc per query record, in FASTA order; NEW, ANC and AMB on
separate rings, with one fixed colour per ring. Hovering a glyph shows its family, cluster,
coordinates, status, `midpoint_gap`, `anc_residual` and `tsd_estimate`. Glyphs too close to
hover individually are spread slightly along their ring; the tooltip says when a glyph was
moved and by how much.

Regenerate from an existing table without re-running the pipeline:

```bash
python -m MOBster.viz --calls results/mobster_calls.csv --query Q.fasta --reference R.fasta
```

The page bundles D3 7.9.0 (ISC licence) and CircosJS 2.1.0, which is **GPL-2.0**. Licence
texts are in `src/MOBster/viz/vendor/`.

## Differences from the spec

1. **minimap2 `-r 100,100`** (spec §4.5 says "override nothing"). asm5's chaining
   bandwidth (`-r 1k,100k`) joins an ancestral construct's two flanks across the
   element gap. The aligner then cannot pay for the deletion, so it emits one flank
   only, with no supplementary record. Every ancestral element came out
   `AMB_NO_MIDPOINT_SPAN`. With a bandwidth below the shortest element, the flanks
   come back as primary + supplementary (the §4.6 two-record shape). The value is
   `MM2_CHAIN_BANDWIDTH` in `constants.py`.
2. **Two-record bracketing uses the midpoint window**: `qe_i ≤ M + w` and
   `qs_j ≥ M − w`, not strictly `qe_i ≤ M ≤ qs_j`. Otherwise one ancestral element with
   a few bases of junction micro-homology or boundary error would be ANC as a single
   record but AMB_COMPLEX as two, which breaks the requirement that both shapes give
   the same answer.
3. **New reason code `AMB_LOW_COVERAGE`**: a group spans the midpoint but covers less
   than `MIN_COVERAGE` of the construct. None of the spec's codes describes this case.
4. **Secondary records** count toward `n_loci` only if they would qualify as a locus on
   their own (span the midpoint window and cover ≥ `MIN_COVERAGE`). Counting every
   secondary hit, including single-flank repeat hits, would make everything
   `AMB_MULTI_LOCUS`.
5. No `-v` flag: the flag set is exactly as in §7, and debug output goes to `mobster.log`.
6. **`insertion_site` is replaced by `reference_position` and `reference_end`**, filled
   for ANC and NEW calls (empty for AMB):
   - **ANC**: the 1-based inclusive span of the ancestral copy in the reference, taken
     from the reference stretch absent from the construct at the junction (the gap
     between the two records, or the in-window deletion for one record).
   - **NEW**: `reference_position` is the insertion point (what the spec's
     `insertion_site` held); `reference_end` is empty.

   Two limits on the ANC span:
   - **Micro-homology:** where the bases at one end of the copy repeat just past the
     other end, the aligner extends a flank record through them, so the span can be
     trimmed by a few bp. In the tier 1 benchmark, all 383 ANC spans were exact or
     trimmed by 1–3 bp, never extended.
   - **Annotation error carries over:** the construct is defined by the query's
     element boundaries, so a query boundary off by N bp moves the reference span by
     the same N bp. The span is the query annotation translated into reference
     coordinates, not an independent boundary estimate. In tier 2, the span offsets
     equal ISEScan's boundary errors plus the micro-homology trim.
7. **Only the query is annotated** (spec §4.2 annotates both genomes). Classification
   never uses the reference IS table (§4.3), and losses come from the swapped run,
   which annotates the reference as its query. This skips the slowest step for one
   genome. `--is-table-reference` is still accepted, as §7 lists it, but ignored.
8. **Two extra output columns, `element_family` and `element_cluster`**, copied from the
   ISEScan table's `family` and `cluster` columns. The map's tooltip needs both
   (`viz_instructions.md` §5), and the §6 schema had only `element_name`.

## Tests and benchmark

```bash
npm install --prefix test/js                        # once: jsdom, for the map's browser checks
pytest                                              # unit + integration (needs minimap2, gdtools)
python -m benchmark.run_benchmark all --threads 10  # spec §8; tier 2 runs ISEScan (slow)
```

The benchmark simulates query genomes from the Barrick lab's annotated REL606 GenBank
(`data/REL606/REL606.barricklab.gbk`). The RefSeq record NC_012967.1 has no IS
annotation, so gdtools cannot place `MOB`s against it. Results go to
`benchmark/results/report.md`. The GD → element-coordinate conversion is
`benchmark/gdcoords.py`, tested against real gdtools in `test/test_gdcoords.py`.
Note that with `duplication_size 0`, gdtools inserts the element *after* `position`.
