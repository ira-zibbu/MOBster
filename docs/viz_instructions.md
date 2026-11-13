# MOBster — interactive circular map (implementation spec)

A companion to `instructions.md`. Produces an interactive circular genome map of
MOBster's calls.

> **Scope.** This is an **exploration and triage interface**, not a publication figure.
> Its purpose is to let the analyst scan a result set and to let a collaborator inspect
> their own data without reading a CSV. Do not optimise it for print, and do not add a
> static-image export path.

---

## 1. Inputs

| Input | Source |
|-------|--------|
| Call table | MOBster's output CSV (`instructions.md` §6) |
| Query record lengths | Read from the query FASTA at generation time (Biopython `SeqIO`); the CLI already has the path |
| Provenance strings | Query filename, reference filename, MOBster version |

No other inputs. In particular the reference IS table is **not** required — this view is
in query coordinates, and every element in the CSV has one.

## 2. Output

**A single self-contained HTML file.** Data inlined as JSON; D3 and CircosJS **vendored
into the file**, not referenced from a CDN. The file must render correctly with
networking disabled, because its delivery mechanism is email and a collaborator opening
it offline is the normal case, not an edge case. A few hundred kilobytes is acceptable.

No build step, no asset directory, no local server.

The page header must state the query filename, reference filename, and MOBster version,
so a file sitting in someone's downloads folder six months from now is self-describing.

## 3. Layout

One arc per query FASTA record, sized proportionally to record length, in input order,
with a small angular gap between arcs. This handles draft assemblies without special
casing; a closed genome is the single-arc case.

Rings, outermost inward:

| Ring | Contents | Glyph |
|------|----------|-------|
| Axis | Coordinate ticks with megabase labels | — |
| 1 | `call == NEW` | Filled circle, strong accent colour |
| 2 | `call == ANC` | Filled circle, muted neutral colour |
| 3 | `call == AMB` | Small light-grey tick |

Every row of the CSV appears on exactly one ring. Rows whose `status` is
`AMB_PARTIAL_ELEMENT` or `AMB_FLANK_TRUNCATED` never reached alignment but are still
`AMB` and still belong on ring 3.

**Rationale for the encoding, so it is not "simplified" later.** The call is carried by
radial position rather than glyph shape because small glyph shapes are not reliably
discriminable at these sizes, and because fill would otherwise have to encode both the
call and the family colour simultaneously. AMB is deliberately given less visual weight
than ANC and NEW: those are findings, AMB is "not yet resolved", and equal prominence
overstates it.

## 4. Colour

**Colour does not encode family.** One fixed colour per ring, reinforcing the call that
radial position already carries. Family and cluster appear in the tooltip only (§5).

This was considered and rejected rather than overlooked. Matching a small isolated glyph
against a legend elsewhere on the page is absolute colour identification, which fails past
roughly six to eight categories however the colours are chosen — so a family palette either
caps out (and then needs an overflow rule that encodes nothing real) or produces colours
that look readable and are not.

The cost is that you can no longer see at a glance whether the novel insertions are
concentrated in one family, which is a question IS biology genuinely asks. The answer is
that family belongs to **filtering** rather than colour — "show me only IS150" is better
served by a legend toggle than by a hue the reader has to decode — so that capability
arrives with the deferred filtering work in §8, not with a palette.

Three colours total: accent for NEW, muted neutral for ANC, light grey for AMB. Choose them
so NEW is the most visually prominent and AMB the least, matching their epistemic status.

Render a static legend naming the three rings and their colours. No interactivity on the
legend in this version (§8).

## 5. Tooltip

On hover, show:

- element name, **ISEScan family, and ISEScan cluster**
- query record and element coordinates (1-based, as in the CSV)
- `call` and `status`
- `midpoint_gap`, `anc_residual`, `tsd_estimate`

The evidence fields are not decoration. The reason `status` is a named code rather than a
flag is so ambiguous calls can be triaged, and the tooltip is where that triage happens —
without the gap and residual on hover, the user still has to open the CSV, and the
interface has bought nothing. Showing the residual on an ancestral call is also how a
systematic boundary offset gets noticed while browsing rather than only by plotting.

Fields that are empty in the CSV are omitted from the tooltip rather than shown blank.

## 6. Collisions

At 4.6 Mb rendered as a circle, two elements a few kilobases apart occupy the same pixel
and the upper one hides the lower one's hover target. Apply a small angular nudge to
glyphs closer than the hover radius so that every element remains individually hoverable.

This is a real case, not a hypothetical: adjacent insertions are among the mutation
classes in the benchmark.

## 7. Implementation notes

- CircosJS is D3-based. Check its current documentation for the correct track type for
  discrete positional glyphs and for the tooltip option name — do not assume an API
  surface from this document.
- Generate the HTML from Python (the same CLI), by templating the inlined JSON and
  vendored library sources into a single file. Keep the template as a data file in the
  package rather than a string literal in code.
- The generator must fail loudly if a query record named in the CSV is absent from the
  FASTA, rather than silently dropping those elements.

## 8. Out of scope

- **Legend-click filtering of rings and families.** Explicitly deferred, though it is the
  highest-value interaction to add next and is nearly free in D3. Note this is now the
  *only* route by which family becomes visually accessible (§4), so "are the novel
  insertions concentrated in one family?" is a question v1 cannot answer from the map —
  it has to be answered from the output table.
- **Click-through to a per-locus alignment view.** Requires per-element aligned blocks,
  which the output CSV does not carry; it would mean either persisting alignment records
  or re-running the aligner. Worth doing eventually — hover answers "what is this", click
  would answer "why" — but it is a larger change than it appears.
- **Density / hotspot track.**
- **Multi-sample comparison.** Would require reference coordinates to be comparable
  across genomes, which in turn needs a `reference_position` column defined for both ANC
  and NEW (currently `insertion_site` is populated for NEW only). Note that AMB elements
  have no reference position by construction, so a reference-frame view cannot place them
  at all; if built, they belong in a sidebar list per sample rather than at invented
  positions. **A tick whose position is a guess is worse than a tick that is not drawn.**
- Static image export, print styling, publication figure concerns of any kind.

## 9. Acceptance checks

1. The file opens and renders correctly in a browser with networking disabled.
2. Every row of the input CSV appears exactly once; per-ring glyph counts equal the
   CSV's value counts for `call`.
3. Two elements 2 kb apart on the same record are both independently hoverable.
4. The tooltip shows every field listed in §5 that is non-empty for that element.
5. Family and cluster are present in every tooltip, and glyph colour depends only on the
   ring — no glyph's colour varies with family.
6. A multi-record query FASTA renders as multiple arcs with correct proportional sizing.
