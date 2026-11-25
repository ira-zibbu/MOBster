# MOBster — dev notes

- Spec: `docs/instructions.md` (v0.1). Deviations are listed below under "Differences from the spec". Keep that list current. `README.md` is user-facing only (install, run, method overview, output guide); keep developer notes here.
- Env: `mamba activate mobster` (environment.yml pins ISEScan 1.7.3; outside the env, `isescan.py` picks up the system Python and fails on libssw).
- Package: `src/MOBster/` (`inputs` → `constructs` → `align` → `classify` → `pipeline`, `cli`). Constants live in `constants.py` and are deliberately not on the CLI.
- Coordinates: 0-based half-open inside `classify.py` (pysam). 1-based inclusive only in IS tables and output.
- Map: `src/MOBster/viz/` (generator, `template.html`, vendored D3 + CircosJS). Spec: `docs/viz_instructions.md`. CircosJS draws each record's glyphs inside a `<g class="block">` rotated to the record's start angle, so a glyph's own transform is record-relative.
- Tests: `pytest` (needs minimap2 + gdtools on PATH; map browser checks need `npm install --prefix test/js`, else skipped). Benchmark: `python -m benchmark.run_benchmark {prepare,tier1,tier2,report}`. Tier 1 runs in seconds and is the regression target.

## Output conventions

For NEW calls, `reference_position` is the 1-based reference base aligned to the construct
midpoint. Because minimap2 left-aligns indels, this is the first duplicated base of the
TSD, the same convention as a GenomeDiff `MOB` position.

## Interactive map (implementation notes)

Every run also writes `mobster_map.html` (spec: `docs/viz_instructions.md`). It's a single
self-contained file, around 570 KB, with the data and both libraries inlined, so it opens
offline and can be emailed. One arc per query record, in FASTA order; NEW, ANC and AMB on
separate rings (plus COMPOSITE, between NEW and ANC), with one fixed colour per ring.
Hovering a glyph shows its family, cluster, coordinates, status, `midpoint_gap`,
`anc_residual` and `tsd_estimate`, a COMPOSITE copy's partner and cargo source, and
marks partial elements. Glyphs too close to
hover individually are spread slightly along their ring; the tooltip says when a glyph was
moved and by how much.

**Paired map** (`mobster_map_paired.html`):
- Reference records come first, then query records, on one circle, with each genome's
  arcs in its own colour.
- Query elements sit on their call rings as above. The reference's own IS elements (from
  ISEScan on the reference) sit on one unclassified ring, and their tooltip gives only
  family, cluster and location.
- Each COMPOSITE unit with a known `source_cargo` gets a ribbon from that reference span
  to the cargo between its two copies in the query. A ribbon that crosses over itself
  means the cargo was inserted inverted relative to the reference (the Circos convention).
- Each ribbon end is drawn at least 0.07 rad (~16 px) wide, centred on the true span.
  A 14 kb cargo on a 7 Mb circle would otherwise be a hairline whose twist cannot be
  seen. The tooltip says when a ribbon was widened, and it always gives true coordinates.
- The ribbons are drawn with D3, not with CircosJS's chord track. That track uses
  `d3.ribbon`, which always runs both ends clockwise, so it cannot twist.

Regenerate from an existing table without re-running the pipeline:

```bash
python -m MOBster.viz --calls results/mobster_calls.csv --query Q.fasta --reference R.fasta \
       --is-table-reference results/isescan/reference/R.fasta.csv   # optional: also the paired map
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
7. **The reference IS table is used only for display.** As in spec §4.2, both genomes
   are annotated. `--is-table-reference` skips ISEScan on the reference. The reference
   table marks the reference's IS elements on the paired map and is not classified.
   Classification never uses it (§4.3), and losses still come from the swapped run.
   (Until 2026-11-24 the reference was not annotated at all.)
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

11. **Composite transpositions get their own call, `COMPOSITE`** (spec §9 lists
    composites as a non-goal). When a composite unit (IS – cargo – IS) transposes as
    one, each IS copy has one flank at the new site and one inside the cargo, whose
    source lies elsewhere. Each copy alone is therefore AMB, and its status reflects
    only the geometry. In barcode05 the same unit appears twice. With the cargo
    reversed, the two pieces land on opposite strands and both copies are
    `NO_MIDPOINT_SPAN`. With it forward, the first copy's pieces are out of order
    (`COMPLEX`), and the second copy's bracket the midpoint across 2.9 Mb
    (`GAP_UNEXPLAINED`). The flanks are unique, so the flank-pair rescue
    (deviation 10) does not apply.

    After per-element classification (`classify.find_composites`), AMB copies with one
    of those three statuses are paired within one query record and cluster. There is
    no cap on the distance between them.
    - **Unit test:** the first copy's left flank and the second copy's right flank
      (the unit's outer flanks, already in `work/flanks.sam`) must form **exactly one**
      consistent flank pair (deviation 10's test) with a NEW gap, taking the whole
      unit's length as the element length.
    - **One partner:** a copy that passes with more than one partner is left as it was.
    - **Rows:** both copies get `call = COMPOSITE`, `status = OK`, and the outer pair's
      `best_locus`, `midpoint_gap`, `tsd_estimate` and `reference_position` (the
      unit's insertion site, NEW convention). `anc_residual` is `midpoint_gap` minus
      the copy's own length, as everywhere else.
    - **New columns** (empty unless COMPOSITE):
      - `composite_partner`: the other copy, `name:start-end`.
      - `n_cargo_sources`: how many reference placements the inner flanks give the
        cargo. The first copy's right flank must be at the cargo start and the second
        copy's left flank at its end, on one record and strand, with an implied length
        within `anc_tolerance` of the cargo's.
      - `source_cargo`: `ref:start-end`, 1-based, filled only when `n_cargo_sources`
        is 1. Start is the reference base at the cargo's first base, reading along the
        query from the first copy to the second, so start > end means the cargo is
        reverse-complemented. It is empty for novel or rearranged cargo (0) and for
        multi-copy cargo (> 1).

    The call rests on the outer flanks only, so an unplaceable cargo still gives
    COMPOSITE.

    Results:
    - barcode05: the 4 IS4_84 AMBs are COMPOSITE, in two units of the same 13,962 bp
      cargo, `ADP1-WT:590505-604466`.
      - 195,877 ↔ 212,374 is inserted at 3,261,201 with the cargo reversed
        (`ADP1-WT:604466-590505`).
      - 466,375 ↔ 482,872 is inserted at 3,510,433 with it forward.
      - Both have the same 2 bp gap as every IS4_84 NEW.
      - Nothing else changed. The remaining AMBs are an IS21 partial whose flank
        crosses the reference's linear origin, and a copy at the contig end.
    - Tier 1 and tier 2 are unchanged, with no COMPOSITE (none are simulated). The four
      DSM15973 samples are unchanged.

    Not handled:
    - a unit whose copies are not both AMB with one of the three statuses;
    - copies of different clusters;
    - cargo shorter than δ (500 bp), which gets no `source_cargo`;
    - a unit already in the reference (that is an ANC gap, ignored).

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
