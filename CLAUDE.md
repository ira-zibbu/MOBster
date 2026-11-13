# MOBster — dev notes

- Spec: `docs/instructions.md` (v0.1). Deviations are listed in `README.md`. Keep that list current.
- Env: `mamba activate mobster` (environment.yml pins ISEScan 1.7.3; outside the env, `isescan.py` picks up the system Python and fails on libssw).
- Package: `src/MOBster/` (`inputs` → `constructs` → `align` → `classify` → `pipeline`, `cli`). Constants live in `constants.py` and are deliberately not on the CLI.
- Coordinates: 0-based half-open inside `classify.py` (pysam). 1-based inclusive only in IS tables and output.
- Map: `src/MOBster/viz/` (generator, `template.html`, vendored D3 + CircosJS). Spec: `docs/viz_instructions.md`. CircosJS draws each record's glyphs inside a `<g class="block">` rotated to the record's start angle, so a glyph's own transform is record-relative.
- Tests: `pytest` (needs minimap2 + gdtools on PATH; map browser checks need `npm install --prefix test/js`, else skipped). Benchmark: `python -m benchmark.run_benchmark {prepare,tier1,tier2,report}`. Tier 1 runs in seconds and is the regression target.
