# MOBster — implementation specification (v0.1)

A tool to detect transposition of mobile genetic elements in bacterial genomes.

> **Status of this document.** This supersedes the earlier two-window specification.
> The detection primitive has changed (see *Background: what changed and why*). The
> figures `fig1.png` and `fig2.png` illustrate the **superseded** two-window scheme and
> should not be used to implement this version.

---

## 1. Goal

Given a reference genome `R` and a query genome `Q`, classify every insertion-sequence
(IS) element detected in `Q` as one of:

| Call | Meaning |
|------|---------|
| `ANC` | Ancestral — the same element occupies the same locus in `R`. Not a new event. |
| `NEW` | Novel transposition — the locus is empty in `R`. |
| `AMB` | Cannot be resolved; carries a reason code naming the predicate that failed. |

"Ancestral" here means **present in the reference**, which is presence/absence polarity,
not phylogenetic polarity. It is correct when `R` genuinely is the ancestor (e.g. an
evolution experiment against its founder). For two sister strains, `ANC` means "shared".

Running the tool a second time with `--reference` and `--query` swapped reports element
**losses**. This requires no extra code and is not a separate feature.

---

## 2. Background: what changed and why

The earlier design extracted two junction-spanning windows per element,
`A = [s-δ, s+δ]` and `B = [e-δ, e+δ]`, each half flank and half element. Because the
element half matches every copy of that family genome-wide, that design had to generate
a large set of alignments and then filter out the element-internal ones — which made it
dependent on the reference IS annotation being complete, and blind to an element
inserting *into* an existing element.

**v0.1 instead builds one construct per element containing no element sequence at all:**
the left flank concatenated directly to the right flank — the sequence the locus would
have if the element were excised. Aligned to the reference:

- a **novel** insertion still has its empty pre-insertion site in `R`, so the construct
  aligns contiguously, with at most a small insertion at the midpoint where the
  target-site duplication (TSD) is represented twice in `Q` and once in `R`;
- an **ancestral** element has the element sitting between those flanks in `R`, so the
  construct aligns with a **deletion of one element length** at the midpoint.

Classification therefore reduces to a single signed number: the indel at the construct
midpoint. No element-derived multi-mapping, no reference-annotation-dependent filter,
no adjacency heuristic.

---

## 3. Conventions

**Coordinates.** All internal geometry is **0-based, half-open** — what `pysam` exposes.
Convert to/from 1-based inclusive only when reading the IS table and when writing output.
Do not mix conventions inside the alignment logic. Note that SAM *text* stores `POS`
1-based while `pysam` converts on read, so `reference_start` is already 0-based; do not
convert it again.

**Element tuple.** Each element is `(name, s, e)` with `s`/`e` the first and last base of
the element *excluding* the TSD, as reported by ISEScan (`isBegin`/`isEnd`, which include
terminal inverted repeats and exclude the TSD). Element length is `e - s + 1` in 1-based
inclusive terms.

**Construct.** For flank size `δ`, the construct is the concatenation of
`[s-δ, s-1]` (left flank) and `[e+1, e+δ]` (right flank). In v0.1 both flanks are always
full length (see §5), so the construct length is `Lc = 2δ` and the **midpoint** `M = δ`.

**Midpoint gap.** The net indel at the midpoint, signed so that

- **positive** = reference consumed more than construct = a deletion relative to the
  construct = the ancestral signature;
- **negative** = construct consumed more than reference = an insertion = the duplicated
  TSD at a novel insertion site.

Never take an absolute value of this quantity. The sign carries the ANC/NEW distinction.

---

## 4. Pipeline

```
reference.fasta, query.fasta
  → [4.1] validate
  → [4.2] annotate (ISEScan, or supplied table)
  → [4.3] select callable elements
  → [4.4] build constructs
  → [4.5] align constructs to reference (minimap2)
  → [4.6] parse, locate midpoint, classify
  → [4.7] write table
```

### 4.1 Validate

Parse both FASTAs with Biopython `SeqIO`. Error clearly if a path does not exist, the
file is not readable as FASTA, or either file is empty. Multi-record FASTAs are
supported; record identifiers must be unique within each file.

### 4.2 Annotate

Run ISEScan on each genome:

```
isescan.py --seqfile <fasta> --output <outdir>
```

ISEScan's HMM search is the pipeline's slow step; the aligner cost is negligible by
comparison, so do not optimise alignment.

**Verify the output layout before writing the reader.** ISEScan nests output under
`--output` mirroring the input file's path rather than writing a flat
`<outdir>/<name>/<name>.csv`, and the output file's name and extension have changed
across versions. Pin the ISEScan version in the environment spec and confirm the actual
path.

Both genomes must be annotated with the **same ISEScan version and parameters**.

**`--is-table-reference` / `--is-table-query`** accept a precomputed IS table and skip
ISEScan for that genome. This is required for the benchmark (§8) and is the cheapest way
to re-run the classifier during development.

### 4.3 Select callable elements

Iterate over elements in the **query** table. An element is callable only if:

- ISEScan's `type` field marks it complete, not partial. Partial calls have unreliable
  inner boundaries by definition. **Emit them in the output** with status
  `AMB_PARTIAL_ELEMENT` rather than dropping them — their share of total AMB decides
  whether the deferred re-annotation work in §9 is worth doing.
- Both flanks fit within the record: `s - δ >= 1` and `e + δ <= len(record)`. Otherwise
  emit `AMB_FLANK_TRUNCATED`.

The reference IS table is **not used for classification**. Load it only if a later
feature needs it.

### 4.4 Build constructs

Write all constructs to a **single FASTA**, named so the element is recoverable
(e.g. `<query_record>|<index>|<element_name>`). Do not align per element — that
re-indexes the reference once per call and produces temp-file churn for no benefit.

### 4.5 Align

```
minimap2 -x asm5 -a -t <threads> reference.fasta constructs.fasta > constructs.sam
```

- `-a` is **mandatory**: it emits SAM and forces base-level alignment, so every record
  carries a CIGAR. Without it the midpoint indel cannot be located at all. (`-c` is the
  PAF equivalent and is not needed here.)
- Parse with `pysam`, not by hand. A hand-rolled CIGAR walker would be the most
  error-prone code in the tool, and an off-by-one there is invisible against a 1.4 kb
  ancestral gap while being fatal against a 3 bp TSD.
- Use **`get_aligned_pairs()`**, which returns one `(query_pos, reference_pos)` tuple per
  alignment column with `None` on whichever side consumes no base. Do **not** use
  `get_blocks()` for this: it reports gapless block boundaries in *reference* coordinates
  only, so it cannot locate a gap on the construct. Constructs are ~1 kb, so the per-column
  list is cheap.
- The preset must come **first**; later flags override it, so `-x asm5 -r …` applies an
  override while `-r … -x asm5` silently discards it. v0.1 overrides nothing.
- Do **not** filter alignments by MAPQ anywhere.

### 4.6 Parse and classify

For each construct, take its alignment records. Group them by `(reference_name, strand)`.

Records flagged secondary (`0x100`) are **not** chained or classified; count the distinct
reference positions among them toward `n_loci`, since they represent alternative
placements. Primary and supplementary records are used normally.

**Strand normalisation — expected to be a no-op under `pysam`, but verify.** The chaining
predicate requires construct coordinates that ascend with reference coordinates. For a
minus-strand hit that is false in the *original* construct orientation, which would
require replacing each interval `[qs, qe)` with `[Lc - qe, Lc - qs)`.

`pysam` reports query coordinates against the sequence as **stored in the record**, which
for a reverse-strand alignment is already the reverse complement — i.e. already the
flipped frame. So the explicit substitution above should be unnecessary.

**Do not take this on trust, and do not apply the flip twice.** Write a unit test over a
known minus-strand case that asserts the ordering predicate holds on `pysam`'s raw
coordinates; apply the substitution only if that test fails.

> Because v0.1 constructs are symmetric (`Lc = 2δ`, `M = δ`), the midpoint is invariant
> under this flip. **This is a v0.1-only simplification.** If truncated flanks are ever
> admitted, the midpoint becomes `Lr` in minus-strand groups and must be recomputed per
> group.

Then, within each group, recognise exactly two shapes:

**One record spanning the midpoint window.** The midpoint gap is the **net indel** within
construct interval `[M - w, M + w]`. Walk `get_aligned_pairs()` in order, maintaining a
construct cursor equal to the last non-`None` query position plus one, and for each column
whose construct position falls in the window accumulate:

- deletion column `(None, rpos)` → **+1** (reference consumed, construct not);
- insertion column `(qpos, None)` → **−1**;
- match column → 0.

Deletion columns have no query position of their own, so their window membership is
decided by the construct cursor. This reproduces the `dr - dq` convention exactly: an
ancestral element's 1 347 bp deletion appears as 1 347 deletion columns at a single
construct cursor position, summing to `+1347`; a novel insertion's 9 bp duplication
appears as nine insertion columns summing to `-9`.

"Net" is literal — if the window contains both the junction indel and an unrelated small
indel a few bases away, the sum is what you get. That is intended, and is part of why
`ANC_TOLERANCE` is not tight.

A window is used rather than the exact coordinate because the nominal midpoint is only as
accurate as the element boundary annotation; a junction 15 bp off nominal is routine.

**Off-midpoint indel guard (required).** While walking, also record the largest single
indel anywhere in the record. If an indel larger than `MAX_OFFSET_INDEL` lies **outside**
the midpoint window, emit `AMB_INDEL_OFF_MIDPOINT` regardless of what the window sum says.

Without this guard, an element whose boundary annotation is off by more than `w` has its
real junction indel fall outside the window, the sum returns near zero, and the element is
confidently classified `NEW` — an ancestral element silently becoming a false novel
transposition, which is the worst error this tool can make. The guard costs one extra
accumulator in a loop already being walked.

**Two records bracketing the midpoint**, i.e. `qe_i <= M <= qs_j` after sorting by
construct start, where with `dq = qs_j - qe_i` and `dr = rs_j - re_i`:

- `dq >= -ovl_max` and `dr >= -ovl_max` — ordered consistently on both axes, allowing a
  few bases of micro-homology overlap at the breakpoint;
- the midpoint gap is `dr - dq`.

Whether an ancestral element arrives as one record with a long `D` or as two records is
decided by minimap2's `-g`, `-r` and `-z`. Both shapes must yield the same number, so
that classification does not depend on preset or version.

**Anything else** — zero records, more than two in a group, or two that do not bracket
the midpoint — is `AMB` with the appropriate code from §6. Do not implement general
co-linear chaining in v0.1 (§9).

A **locus** is a group that spans the midpoint and whose records cover at least
`min_coverage` of the construct. Classification:

| Condition | Call |
|-----------|------|
| Exactly one locus, midpoint gap within `tsd_tolerance` of 0 | `NEW` |
| Exactly one locus, midpoint gap within `anc_tolerance` of `e - s + 1` | `ANC` |
| Exactly one locus, gap matching neither | `AMB_GAP_UNEXPLAINED` |
| Two or more loci | `AMB_MULTI_LOCUS` |
| No group spans the midpoint | `AMB_NO_MIDPOINT_SPAN` |
| No records at all | `AMB_NO_LOCUS` |

Coverage is a quality gate, not the decisive test: midpoint spanning is checked first,
then coverage is used to reject weak loci.

### 4.7 Output

One CSV row per **query element**, including non-callable ones.

---

## 5. Fixed constants

v0.1 exposes **no** tuning parameters on the CLI. Each of the following is a named
module-level constant with the stated value. Making them configurable, and sweeping them,
is deferred (§9) — the point of v0.1 is to get numbers on the board, and a tunable
parameter is an invitation to tune before there is anything to tune against.

| Name | Value | Notes |
|------|-------|-------|
| `DELTA` | 500 | Flank length per side. Kanai *et al.* (2025) resolved sites with 100 bp, so this is likely larger than necessary; it is not changed in v0.1. Must exceed `OVL_MAX` and `MIDPOINT_WINDOW` by a wide margin. |
| `MIDPOINT_WINDOW` | 25 | Half-width of the midpoint window, in construct bases. |
| `OVL_MAX` | 20 | Permitted construct/reference overlap when joining two records. |
| `MIN_COVERAGE` | 0.80 | Fraction of `Lc` a locus must cover. |
| `TSD_TOLERANCE` | 25 | How far the midpoint gap may sit from 0 and still be `NEW`. Must accommodate a negative gap: TSDs run ~3–9 bp for common families (IS1 8–9, IS186 6, IS2 5, IS3 3, IS150 3). |
| `MAX_OFFSET_INDEL` | 50 | An indel larger than this lying outside the midpoint window triggers `AMB_INDEL_OFF_MIDPOINT`. Must exceed the largest indel expected from ordinary flank variation, and be far below the smallest element length. |
| `ANC_TOLERANCE` | `max(50, 0.05 × element length)` | Both absolute and relative, because annotation boundary error shifts the gap by tens of bases and the query and reference copies can legitimately differ in length. |

None of these values is defended on a-priori grounds; they are reasonable starting points
held fixed so that benchmark results are attributable to the algorithm rather than to
parameter choice.

---

## 6. Output schema

One row per query element:

| Column | Description |
|--------|-------------|
| `query_record` | Query sequence identifier |
| `element_name` | ISEScan family/cluster name |
| `element_start`, `element_end` | 1-based inclusive, TSD-excluded |
| `element_length` | `e - s + 1` |
| `element_type` | ISEScan `type` (complete/partial) |
| `delta` | Flank length used |
| `call` | `ANC`, `NEW`, or `AMB` |
| `status` | Reason code (below); `OK` for `ANC`/`NEW` |
| `n_records` | Alignment records returned for this construct |
| `n_loci` | Groups spanning the midpoint and clearing `min_coverage` |
| `best_locus` | `<reference_record>:<start>-<end>:<strand>` or empty |
| `construct_coverage` | Fraction of `Lc` covered by the best locus |
| `identity` | Weighted mean identity of the best locus |
| `midpoint_gap` | **Signed**, raw, emitted for every call that reached alignment |
| `anc_residual` | `midpoint_gap - element_length`; the ANC test's residual |
| `tsd_estimate` | `-midpoint_gap` when `NEW`, else empty |
| `insertion_site` | Reference coordinate of the midpoint when `NEW`, else empty |

`midpoint_gap` and `anc_residual` must be written **for every element that reached
alignment, including AMB ones**. They turn ANC/NEW from a categorical decision into a
threshold on a continuous quantity, and a systematic offset in the `anc_residual`
distribution is the signature of a coordinate-convention bug — the failure mode most
likely to be hiding in this pipeline.

**Reason codes** (`status`):

`OK`, `AMB_PARTIAL_ELEMENT`, `AMB_FLANK_TRUNCATED`, `AMB_NO_LOCUS`,
`AMB_NO_MIDPOINT_SPAN`, `AMB_MULTI_LOCUS`, `AMB_GAP_UNEXPLAINED`,
`AMB_INDEL_OFF_MIDPOINT`, `AMB_COMPLEX` (more than two records in a group, or two that do
not bracket the midpoint).

AMB is a triage queue, not a dead end: every ambiguous call names the predicate that
failed, so failures can be grouped and counted rather than inspected one by one.

---

## 7. CLI

```
mobster --reference R.fasta --query Q.fasta --outdir results/ [options]
```

Required: `--reference`, `--query`. Also: `--outdir`, `--threads`,
`--is-table-reference`, `--is-table-query`. That is the complete flag set — the §5
constants are deliberately **not** exposed.
Log to stderr with levels; keep stdout clean for the table if it is ever piped.

---

## 8. Test data and benchmark

Simulate transposition events against *E. coli* B str. REL606 (RefSeq `NC_012967.1`).

```
gdtools APPLY -r REL606.gbk -f FASTA -o query.fasta mutations.gd
```

The **annotated** reference is required, not `REL606.fasta`: `MOB` lines resolve
`repeat_name` against annotated `repeat_region`/`mobile_element` features, so a bare
FASTA cannot place an element. In a `MOB` line, `position` is the **first duplicated
base** of the target site, `strand` is ±1, and `duplication_size` is the TSD length.

> Converting a GD `MOB` position plus `duplication_size` into the TSD-excluded
> `(name, s, e)` convention is the most error-prone step in the whole harness.
> **Give it its own unit test.**

**Cases to generate.** Unmodified REL606 copies (`ANC`); `MOB` at clean intergenic sites
in both orientations and with `duplication_size` 0 versus the family's native TSD
(`NEW`); `MOB` into a site already inside another element (`NEW`); two `MOB`s within δ of
each other; `MOB` within δ of a REP or *rhs* repeat (expect `AMB`); `MOB` plus a nearby
SNP or small indel in the flank (should stay `NEW` — tests the identity thresholds); a
local inversion or flanking deletion spanning an ancestral copy (expect `AMB`); and as a
hard negative, an ancestral element whose flank carries an unrelated small deletion,
which is the realistic source of ANC→AMB leakage.

**Run every benchmark twice**, on the same genomes:

1. with the exact GD-derived element coordinates supplied via `--is-table-query`;
2. with ISEScan annotation.

The difference between the two accuracies **is** the cost of annotation error, separated
from classifier error. Without this decomposition every misclassification is ambiguous
between "the logic is wrong" and "ISEScan put the boundary in the wrong place". Tier 1
also runs in milliseconds, so it is the regression target to iterate against.

**Report** the ANC/NEW/AMB partition per class and the reason-code breakdown of AMB, not
just overall accuracy — the interesting failure mode is AMB inflation. Run at the fixed
§5 constants only; parameter sweeps are deferred (§9) until this baseline exists.

---

## 9. Deferred — do not implement in v0.1

Each of these addresses a real failure mode. None should be built until a benchmark shows
it matters.

- **CLI exposure and sweeping of the §5 constants.** Once the tier-1 baseline exists,
  `DELTA` and `MIDPOINT_WINDOW` are the first worth varying — δ over roughly 100–500,
  watching the ANC/NEW/AMB partition rather than overall accuracy. Not before: a sweep
  run against an unvalidated classifier attributes algorithm errors to parameters.
- **General co-linear chaining.** Iterative maximum-weight chaining over arbitrary record
  counts, replacing the one-/two-record shapes of §4.6. Build only if `AMB_COMPLEX` is a
  large bucket. Note it must not be implemented as connected components: records covering
  the same construct interval at different loci are mutually exclusive alternatives, not
  parts of one chain, and must not be fused.
- **Per-flank diagnostic second pass.** Aligning each flank separately, for failed
  elements only, to distinguish an unplaceable flank from a multiply-placeable one.
- **Adaptive δ** — per element and per side, the smaller of global δ and the distance to
  the nearest annotated element on that side.
- **Truncated-flank support** — shrinking rather than skipping flanks near a record end.
  Requires per-group midpoint recomputation (§4.6) and coverage against actual rather
  than nominal construct length.
- **Family-consensus re-annotation.** Replacing ISEScan boundaries by clustering blastn
  hits of a family consensus along each genome (gaps ≤ 20 bp, a minimum continuous match
  expressed as a *fraction* of family length rather than Kanai *et al.*'s absolute
  300 bp). This would recover usable coordinates for internally deleted copies and give
  both genomes a shared boundary convention. Justified only if `AMB_PARTIAL_ELEMENT` is a
  large share of AMB.
- **Score margin reporting** — best locus score versus best score at a different locus,
  kept raw instead of thresholded.
- **Secondary-alignment control.** v0.1 leaves `-p`/`-N` at preset defaults, so a
  genuinely ambiguous repeat-driven placement may be pruned and reported as a single
  confident locus. If `AMB_MULTI_LOCUS` is implausibly rare in the benchmark, this is the
  first thing to check.

**Explicit non-goals.** Composite transposons (two elements flanking cargo are annotated
as two independent elements; a newly transposed composite presents each copy with one
cargo-internal flank, so both land in AMB). Read-level validation. Multi-genome time
series.

---

## 10. Hazard checklist

Consult when benchmark numbers disappoint.

- **Minus-strand handling, in either direction.** Getting the §4.6 orientation wrong
  leaves plus-strand loci correct and minus-strand loci apparently misplaced — roughly a
  50 % failure rate rather than a crash. Under `pysam` the likelier error is applying the
  flip *twice* (once by the format, once by hand), which looks identical.
- **Missing `-a`.** No CIGAR, no midpoint indel, every call falls through to AMB.
- **Double coordinate conversion.** `pysam` already returns 0-based coordinates;
  subtracting one again shifts every gap measurement by a base.
- **Absolute value on `midpoint_gap`.** Destroys the ANC/NEW distinction.
- **Boundary error exceeding `w`.** The junction indel falls outside the midpoint window,
  the window sum returns near zero, and an ancestral element is confidently called `NEW`.
  The §4.6 off-midpoint guard exists for this; if `AMB_INDEL_OFF_MIDPOINT` never fires in
  the benchmark, check that the guard is actually wired up rather than concluding boundary
  error is absent.
- **Coordinate conventions.** A systematic shift in `anc_residual` across many calls
  means 1-based/0-based or inclusive/half-open mixing, not a biological signal.
- **Record count as a discriminator.** Never branch on how many records came back:
  whether an ancestral element arrives as one record or two is a function of `-g`, `-r`
  and `-z`, so a preset change would silently reshuffle classification.
- **Top-hit selection.** Never select "the best alignment" and classify from it alone.
- **TSD as a gate.** `tsd_estimate` is a reported field and an aggregate QC signal — per
  family, `NEW` calls should cluster near the known TSD length — never a per-call filter.
  With δ = 500 and a 3 bp duplication, measurement noise is comparable to signal.
- **Ties at duplicated loci.** When an element inserts into a locus duplicated in the
  reference, flank alignments can tie with no principled tie-break. Kanai *et al.* (2025)
  flag this hazard and report one instance in their data, verified by hand, where the
  automatic assignment was nonetheless correct. `AMB_MULTI_LOCUS` exists so these surface
  without manual inspection of every call; expect them to be rare.

---

## 11. Reference

Kanai Y, Shibai A, Yokoi N, Tsuru S, Furusawa C. Laboratory evolution of the bacterial
genome structure through insertion sequence activation. *Nucleic Acids Research* 53(9),
2025. doi:10.1093/nar/gkaf331 — source for the homology-clustering annotation rule
(§9), the 100 bp flank precedent (§5), the IS1 TSD range, and the duplicated-locus tie
hazard (§10).
