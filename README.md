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
coordinates, status, `midpoint_gap`, `anc_residual` and `tsd_estimate`,
and marks partial elements. Glyphs too close to
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

   The ANC span comes from the separately aligned flanks that confirm the call
   (deviation 10). It is not taken from the construct's records: where the bases at one
   end of the copy repeat just past the other end, the aligner extends a construct
   record through them and trims the span by 1–3 bp, which a separately aligned flank
   cannot do. In tier 1, all 390 ANC spans are exact (from the construct records: 181).

   Limit: **annotation error carries over.** The flanks are defined by the query's
   element boundaries, so a query boundary off by N bp moves the reference span by the
   same N bp. The span is the query annotation translated into reference coordinates,
   not an independent boundary estimate. In tier 2, the span offsets equal ISEScan's
   boundary errors exactly (374/374).
7. **Only the query is annotated** (spec §4.2 annotates both genomes). Classification
   never uses the reference IS table (§4.3), and losses come from the swapped run,
   which annotates the reference as its query. This skips the slowest step for one
   genome. `--is-table-reference` is still accepted, as §7 lists it, but ignored.
8. **Two extra output columns, `element_family` and `element_cluster`**, copied from the
   ISEScan table's `family` and `cluster` columns. The map's tooltip needs both
   (`viz_instructions.md` §5), and the §6 schema had only `element_name`.
9. **Partial elements are classified, not set aside** (spec §4.3 emits them as
   `AMB_PARTIAL_ELEMENT`, which is no longer produced). The construct excises exactly the
   annotated span, so ANC/NEW is decided against the partial's own length: a fragment
   present in the reference gives a gap equal to its length (ANC), and an empty site
   gives ≈ 0 (NEW). They are identifiable by `element_type` (`p`), and the map's tooltip
   shows `Type: partial` for them. Where the annotated boundary falls inside a longer
   copy, a flank contains element sequence and the call ends up AMB through the normal
   predicates, not as a wrong ANC/NEW.
   - Tier 2: all 127 partials were classified, with no wrong ANC/NEW. 112 were ANC (all
     `anc_residual` 0), including 72 fragments with no truth match. 5 were new IS3
     insertions that ISEScan labels partial, now NEW. The other 10 were AMB. 8 of them
     are one IS3 copy that ISEScan truncated by 916 bp, which is `AMB_NO_MIDPOINT_SPAN`
     in every case.
   - Real data (barcode05): the 168 bp IS3 fragment is now ANC with residual 0.
10. **Flank-pair check: a second alignment pass and two new reason codes,
    `AMB_MULTI_PAIRING` and `AMB_FLANK_CONFLICT`.** The construct alignment assumes both
    flanks are unique in the reference. In IS-dense genomes they often are not, which
    breaks it in two ways:
    - **(a) A flank lies inside a neighbouring multi-copy IS.** It matches every copy
      equally, and minimap2 puts its primary hit at the wrong one. The result is AMB
      with `NO_MIDPOINT_SPAN`, `COMPLEX`, or `GAP_UNEXPLAINED` (when the wrong copy is
      downstream on the same strand).
    - **(b) The reference has a copy of the site without the element.** This can be a
      paralog of the site, a tandem copy of the element, or a nested same-family
      insertion. The construct aligns in one piece to the element-free copy, which
      beats the true split alignment. minimap2 then drops the split alignment, so the
      result is a false NEW.

    So each flank is also aligned on its own (`work/flanks.fasta`, `work/flanks.sam`;
    minimap2 asm5 with `-N FLANK_MAX_SECONDARY`, i.e. 100). Every (left hit, right hit)
    pair is listed that is on one record and strand, brackets the midpoint, covers
    `MIN_COVERAGE`, passes the off-midpoint indel guard, and has a gap within the NEW
    or ANC tolerance. These *consistent pairs* are counted in the new column
    `n_flank_pairs`.
    - **Confirmation:** a construct ANC/NEW stands only if there is exactly one
      consistent pair, with the same call, on the same record and strand, and with the
      reference base at the midpoint within `MIDPOINT_WINDOW` (25 bp) of the
      construct's. Otherwise the call becomes AMB: `AMB_MULTI_PAIRING` for more than
      one pair, and **`AMB_FLANK_CONFLICT`** for no pair or one that disagrees. A
      confirmed ANC reports the pair's reference span (see deviation 6).
    - **Rescue:** a construct AMB of one of the three statuses above is called from the
      flank pair when there is exactly one consistent pair **and** at least one flank
      has more than one placement (primary + secondary hits; supplementary pieces
      don't count). The repeated flank is what explains the construct's miss. A
      construct AMB with two single-copy flanks keeps its status, and so does a real
      disruption, which pairs with the same unexplained gap and so has no consistent
      pair.

    Results:
    - Before the change, every ANC and NEW call that was correct had exactly one
      consistent pair at the same site: 870 ANC and 58 NEW across the benchmark, and
      all 16 real NEWs in barcode05. So confirmation removed no correct call.
    - Tier 1: one IS1 copy with an IS-flanked side (`AMB_COMPLEX` in all 7 cases) is
      now ANC, and nothing else changed.
    - Tier 2: `anc` went from 45/48 to 48/48 and `anc_background` from 300/320 to
      320/320; every NEW class is unchanged. The ISAS1 false NEW (gap 25, an R flank
      smeared into a diverged repeat at the element's own start) in every case,
      including `unmodified`, is now `AMB_FLANK_CONFLICT`.
    - Real data (K. sucrofermentans DSM15973 against its ancestor): 11 AMB became ANC
      (residual 0), and both false NEWs became `AMB_MULTI_PAIRING`. All 11 rescues
      have a repeated flank.
    - Not handled: an element where both flanks are repeat sequence, e.g. a
      same-family IS nested inside another copy. That stays AMB.

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
