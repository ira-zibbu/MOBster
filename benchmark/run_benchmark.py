"""MOBster v0.1 benchmark (spec §8).

    python -m benchmark.run_benchmark prepare     # simulate query genomes with gdtools
    python -m benchmark.run_benchmark tier1       # exact GD-derived element coordinates
    python -m benchmark.run_benchmark tier2       # ISEScan annotation (slow)
    python -m benchmark.run_benchmark report      # write results/report.md

Every element in a query genome gets a truth class and an expected call. The
expected call is the spec's stated expectation where it gives one (e.g. AMB near
REP/rhs repeats), otherwise the true ANC/NEW label.
"""

from __future__ import annotations

import argparse
import csv
import logging
import random
import statistics
import subprocess
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq

from .gdcoords import Mut, existing_element, mob_element, write_gd

ROOT = Path(__file__).resolve().parent.parent
GBK = ROOT / "data" / "REL606" / "REL606.barricklab.gbk"
OUT = ROOT / "benchmark" / "results"
REF_FASTA = OUT / "REL606.fasta"
SEED = 20261113
FLANK_CLEAR = 3000          # sites/copies kept this far from other IS, rhs, rRNA
NATIVE_TSD = {"IS1": 9, "IS150": 3, "IS186": 6, "IS3": 3, "IS2": 5}

log = logging.getLogger("bench")


# ----------------------------------------------------------------------------- reference model

@dataclass
class Copy:
    family: str
    a: int      # 1-based inclusive
    b: int
    strand: int

    @property
    def length(self):
        return self.b - self.a + 1


def load_reference():
    rec = next(SeqIO.parse(GBK, "genbank"))
    rec.id, rec.description = rec.name, ""   # gdtools keys on the LOCUS name ("REL606")
    copies, repeats, genes = [], [], []
    for f in rec.features:
        a, b = int(f.location.start) + 1, int(f.location.end)
        if f.type == "repeat_region" and "mobile_element" in f.qualifiers:
            fam = f.qualifiers["mobile_element"][0].split(":")[-1]
            copies.append(Copy(fam, a, b, f.location.strand or 1))
        elif f.type == "rRNA" or (f.type == "gene" and f.qualifiers.get("gene", [""])[0].startswith("rhs")):
            repeats.append((a, b, f.qualifiers.get("gene", f.qualifiers.get("product", ["?"]))[0]))
        if f.type in ("gene", "rRNA", "tRNA"):
            genes.append((a, b))
    return rec, copies, repeats, genes


def family_length(copies, fam):
    return Counter(c.length for c in copies if c.family == fam).most_common(1)[0][0]


def near(x, intervals, d):
    return any(a - d <= x <= b + d for a, b in intervals)


def isolated_copies(copies, repeats):
    blocks = [(c.a, c.b) for c in copies] + [(a, b) for a, b, _ in repeats]
    out = []
    for c in copies:
        others = [iv for iv in blocks if iv != (c.a, c.b)]
        if not any(a - FLANK_CLEAR <= c.b and c.a <= b + FLANK_CLEAR for a, b in others):
            out.append(c)
    return out


def clean_sites(rec, copies, repeats, genes, n, rng):
    """Midpoints of intergenic gaps >= 150 bp, far from IS/rhs/rRNA, >= 20 kb apart."""
    genes = sorted(genes)
    blocks = [(c.a, c.b) for c in copies] + [(a, b) for a, b, _ in repeats]
    cands, end = [], 0
    for a, b in genes:
        if a - end > 150:
            mid = (end + a) // 2
            if 10_000 < mid < len(rec) - 10_000 and not near(mid, blocks, FLANK_CLEAR):
                cands.append(mid)
        end = max(end, b)
    rng.shuffle(cands)
    chosen = []
    for c in cands:
        if all(abs(c - x) > 20_000 for x in chosen):
            chosen.append(c)
        if len(chosen) == n:
            return sorted(chosen)
    raise RuntimeError("not enough clean sites")


# ----------------------------------------------------------------------------- cases

@dataclass
class Truth:
    name: str
    s: int
    e: int
    cls: str
    expect: str
    family: str
    note: str = ""
    # Expected MOBster reference_position / reference_end, in REL606 coordinates:
    # ancestral copy span, or for a MOB the insertion point (GD position = first
    # duplicated base; with dup 0 the element follows `pos`, so pos + 1).
    ref_start: int | str = ""
    ref_end: int | str = ""


def build_cases(rec, copies, repeats, genes):
    rng = random.Random(SEED)
    L = {f: family_length(copies, f) for f in NATIVE_TSD}
    iso = [c for c in isolated_copies(copies, repeats) if c.family in NATIVE_TSD]
    rng.shuffle(iso)
    sites = iter(clean_sites(rec, copies, repeats, genes, 60, rng))
    seq = str(rec.seq)

    def mob(fam, pos, strand, dup):
        return Mut("MOB", pos, size=L[fam], dup=dup, strand=strand, name=fam)

    cases: dict[str, tuple[list[Mut], dict[int, tuple[str, str, str]], dict]] = {}
    # each: (muts, {id(mob): (cls, expect, note)}, {copy_index: (cls, expect, note)})

    cases["unmodified"] = ([], {}, {})

    muts, lab = [], {}
    for fam, tsd in NATIVE_TSD.items():
        for strand in (1, -1):
            for dup in (0, tsd):
                m = mob(fam, next(sites), strand, dup)
                muts.append(m)
                lab[id(m)] = ("new_clean", "NEW", f"strand={strand} dup={dup}")
    cases["new_clean"] = (muts, lab, {})

    muts, lab, anc = [], {}, {}
    hosts = [c for c in iso if c.family == "IS1"][:2] + [c for c in iso if c.family in ("IS3", "IS186")][:2]
    for c, fam in zip(hosts, ["IS150", "IS150", "IS1", "IS1"]):
        m = mob(fam, (c.a + c.b) // 2, 1, NATIVE_TSD[fam])
        muts.append(m)
        lab[id(m)] = ("new_in_element", "NEW", f"inside {c.family}@{c.a}")
        anc[copies.index(c)] = ("anc_disrupted", "AMB", f"host of {fam}")
    cases["new_in_element"] = (muts, lab, anc)

    muts, lab = [], {}
    for fam1, fam2 in (("IS1", "IS150"), ("IS186", "IS3"), ("IS150", "IS150"), ("IS2", "IS1")):
        p = next(sites)
        for fam, pos in ((fam1, p), (fam2, p + 250)):
            m = mob(fam, pos, rng.choice((1, -1)), NATIVE_TSD[fam])
            muts.append(m)
            lab[id(m)] = ("new_pair", "NEW", "two MOBs 250 bp apart")
    cases["new_pair"] = (muts, lab, {})

    muts, lab = [], {}
    rhs = [r for r in repeats if str(r[2]).startswith("rhs")]
    rrs = [r for r in repeats if str(r[2]).startswith("rrs")][:2]
    for a, b, nm in rhs + rrs:
        m = mob("IS1", (a + b) // 2, 1, 9)
        muts.append(m)
        lab[id(m)] = ("new_near_repeat", "AMB", f"inside {nm}")
    cases["new_near_repeat"] = (muts, lab, {})

    muts, lab = [], {}
    for fam in ("IS1", "IS150", "IS186", "IS3", "IS2", "IS1"):
        p = next(sites)
        m = mob(fam, p, 1, NATIVE_TSD[fam])
        snp_ref = seq[p + 80 - 1]
        muts += [m, Mut("SNP", p + 80, seq="A" if snp_ref != "A" else "G"), Mut("DEL", p - 150, size=4)]
        lab[id(m)] = ("new_flank_variation", "NEW", "SNP +80, 4bp DEL -150")
    cases["new_flank_variation"] = (muts, lab, {})

    muts, anc = [], {}
    for c in iso[4:7]:
        muts.append(Mut("INV", c.a - 200, size=c.length + 400))
        anc[copies.index(c)] = ("anc_rearranged", "AMB", "symmetric INV incl. 200 bp flanks")
    for c in iso[7:10]:
        muts.append(Mut("DEL", c.a - 350, size=300))
        anc[copies.index(c)] = ("anc_rearranged", "AMB", "300 bp DEL in left flank")
    cases["anc_rearranged"] = (muts, {}, anc)

    muts, anc = [], {}
    for c in iso[10:16]:
        muts.append(Mut("DEL", c.a - 150, size=6))
        anc[copies.index(c)] = ("anc_flank_smalldel", "ANC", "6 bp DEL 150 bp into left flank")
    cases["anc_flank_smalldel"] = (muts, {}, anc)
    return cases, L


def truth_for(case, muts, mob_lab, anc_lab, copies):
    out = []
    for i, c in enumerate(copies):
        s, e = existing_element(c.a, c.b, muts)
        cls, exp, note = anc_lab.get(i, ("anc" if case == "unmodified" else "anc_background", "ANC", ""))
        out.append(Truth(f"{c.family}_ref{c.a}", s, e, cls, exp, c.family, note, c.a, c.b))
    for m in muts:
        if m.kind == "MOB":
            s, e = mob_element(m, muts)
            cls, exp, note = mob_lab[id(m)]
            out.append(Truth(f"{m.name}_new{m.pos}", s, e, cls, exp, m.name, note,
                             m.pos if m.dup else m.pos + 1))
    return sorted(out, key=lambda t: t.s)


def prepare(args):
    rec, copies, repeats, genes = load_reference()
    OUT.mkdir(parents=True, exist_ok=True)
    SeqIO.write(rec, REF_FASTA, "fasta")
    fam_seqs = defaultdict(set)
    for c in copies:
        x = str(rec.seq[c.a - 1:c.b])
        fam_seqs[c.family] |= {x, str(Seq(x).reverse_complement())}
    cases, L = build_cases(rec, copies, repeats, genes)
    log.info("family lengths: %s", L)
    for case, (muts, mob_lab, anc_lab) in cases.items():
        d = OUT / "cases" / case
        d.mkdir(parents=True, exist_ok=True)
        write_gd(muts, rec.id, d / "query.gd")
        subprocess.run(["gdtools", "APPLY", "-r", str(GBK), "-f", "FASTA", "-o", str(d / "query.fasta"),
                        str(d / "query.gd")], check=True, capture_output=True)
        q = next(SeqIO.parse(d / "query.fasta", "fasta"))
        qs = str(q.seq).upper()
        assert len(qs) == len(rec) + sum(m.net for m in muts), case
        truth = truth_for(case, muts, mob_lab, anc_lab, copies)
        for t in truth:   # every new element must sit exactly where we claim
            if "_new" in t.name:
                assert qs[t.s - 1:t.e] in fam_seqs[t.family], (case, t)
        with open(d / "truth.csv", "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["seqID", "family", "cluster", "isBegin", "isEnd", "isLen", "type",
                        "truth_class", "expect", "note", "ref_start", "ref_end"])
            for t in truth:
                w.writerow([q.id, t.family, t.name, t.s, t.e, t.e - t.s + 1, "c", t.cls, t.expect, t.note,
                            t.ref_start, t.ref_end])
        log.info("%-22s %3d mutations, %3d elements", case, len(muts), len(truth))


# ----------------------------------------------------------------------------- running

def cases_list():
    return sorted(p.name for p in (OUT / "cases").iterdir() if p.is_dir())


def mobster(query, table, outdir, threads):
    # The reference table only feeds the paired map; it is computed once and cached, so
    # neither tier re-runs ISEScan on REL606.
    ref_table = isescan(REF_FASTA, OUT / "isescan" / "reference", threads)
    cmd = [sys.executable, "-m", "MOBster", "--reference", str(REF_FASTA), "--query", str(query),
           "--outdir", str(outdir), "--threads", str(threads), "--is-table-query", str(table),
           "--is-table-reference", str(ref_table)]
    subprocess.run(cmd, check=True, capture_output=True)


def tier1(args):
    for case in cases_list():
        d = OUT / "cases" / case
        mobster(d / "query.fasta", d / "truth.csv", OUT / "tier1" / case, args.threads)
        log.info("tier1 %s done", case)


def isescan(fasta: Path, outdir: Path, threads: int) -> Path:
    from MOBster.annotate import run_isescan
    table = outdir / f"{fasta.name}.csv"
    if table.exists():
        return table
    res = run_isescan(fasta, outdir, threads)
    if res is None:
        raise RuntimeError(f"no ISEScan table for {fasta}")
    return res


def tier2(args):
    logging.getLogger("MOBster").setLevel(logging.INFO)
    for case in cases_list():
        d = OUT / "cases" / case
        table = isescan(d / "query.fasta", OUT / "isescan" / case, args.threads)
        mobster(d / "query.fasta", table, OUT / "tier2" / case, args.threads)
        log.info("tier2 %s done", case)


# ----------------------------------------------------------------------------- scoring

def read_csv(p):
    with open(p) as fh:
        return list(csv.DictReader(fh))


def match_isescan(calls, truth):
    """Pair ISEScan elements with truth elements by reciprocal overlap >= 0.5."""
    out, used = [], set()
    for c in calls:
        cs, ce = int(c["element_start"]), int(c["element_end"])
        best, best_score = None, 0.5
        for i, t in enumerate(truth):     # best reciprocal overlap, not first hit
            ts, te = int(t["isBegin"]), int(t["isEnd"])
            ov = min(ce, te) - max(cs, ts) + 1
            score = min(ov / (ce - cs + 1), ov / (te - ts + 1)) if ov > 0 else 0
            if score >= best_score:
                best, best_score = i, score
        if best is None:
            out.append((c, None))
        else:
            used.add(best)
            out.append((c, truth[best]))
    missed = [t for i, t in enumerate(truth) if i not in used]
    return out, missed


def tier_rows(tier):
    rows, missed = [], []
    for case in cases_list():
        truth = read_csv(OUT / "cases" / case / "truth.csv")
        p = OUT / tier / case / "mobster_calls.csv"
        if not p.exists():
            continue
        calls = read_csv(p)
        if tier == "tier1":
            pairs = list(zip(calls, truth))
            miss = []
        else:
            pairs, miss = match_isescan(calls, truth)
        for c, t in pairs:
            rows.append(dict(c, case=case, truth_class=t["truth_class"] if t else "unmatched",
                             expect=t["expect"] if t else "", truth_family=t["family"] if t else "",
                             truth_ref_start=t.get("ref_start", "") if t else "",
                             truth_ref_end=t.get("ref_end", "") if t else "",
                             boundary_err=(int(c["element_start"]) - int(t["isBegin"]),
                                           int(c["element_end"]) - int(t["isEnd"])) if t else ""))
        missed += [dict(t, case=case) for t in miss]
    return rows, missed


CLASS_ORDER = ["anc", "anc_background", "anc_flank_smalldel", "anc_rearranged", "anc_disrupted",
               "new_clean", "new_in_element", "new_pair", "new_flank_variation", "new_near_repeat",
               "unmatched"]


def md_table(header, rows):
    out = ["| " + " | ".join(header) + " |", "|" + "---|" * len(header)]
    out += ["| " + " | ".join(str(x) for x in r) + " |" for r in rows]
    return "\n".join(out)


def tier_section(tier, title):
    rows, missed = tier_rows(tier)
    if not rows:
        return f"## {title}\n\n_not run_\n"
    by = defaultdict(list)
    for r in rows:
        by[r["truth_class"]].append(r)
    table = []
    for cls in [c for c in CLASS_ORDER if c in by]:
        rs = by[cls]
        cnt = Counter(r["call"] for r in rs)
        exp = Counter(r["expect"] for r in rs).most_common(1)[0][0] if rs[0]["expect"] else "-"
        ok = sum(r["call"] == r["expect"] for r in rs) if exp != "-" else "-"
        amb = Counter(r["status"] for r in rs if r["call"] == "AMB")
        table.append([cls, len(rs), exp, cnt["ANC"], cnt["NEW"], cnt["AMB"],
                      f"{ok}/{len(rs)}" if ok != "-" else "-",
                      ", ".join(f"{k.removeprefix('AMB_')}:{v}" for k, v in amb.most_common()) or ""])
    s = [f"## {title}\n",
         md_table(["class", "n", "expect", "ANC", "NEW", "AMB", "as expected", "AMB reasons"], table)]
    if missed:
        mc = Counter(m["truth_class"] for m in missed)
        s.append("\nTruth elements with no ISEScan match (reciprocal overlap ≥ 0.5): "
                 + ", ".join(f"{k}: {v}" for k, v in mc.most_common()))
    res = [int(r["anc_residual"]) for r in rows if r["call"] == "ANC" and r["anc_residual"] != ""]
    if res:
        s.append(f"\n`anc_residual` over ANC calls (n={len(res)}): median {statistics.median(res)}, "
                 f"mean {statistics.mean(res):.2f}, min {min(res)}, max {max(res)}; "
                 f"exactly 0: {sum(x == 0 for x in res)}")
    tsd = defaultdict(list)
    for r in rows:
        if r["call"] == "NEW" and r["truth_family"]:
            tsd[r["truth_family"]].append(int(r["tsd_estimate"]))
    if tsd:
        s.append("\n`tsd_estimate` for NEW calls by family (native TSD in brackets): "
                 + "; ".join(f"{f} [{NATIVE_TSD.get(f, '?')}]: {sorted(v)}" for f, v in sorted(tsd.items())))
    s.append(ref_coord_summary(rows))
    if tier == "tier2":
        errs = [r["boundary_err"] for r in rows if r["boundary_err"]]
        if errs:
            s.append(f"\nISEScan boundary error vs truth (start, end), bp: median start "
                     f"{statistics.median(e[0] for e in errs)}, median end "
                     f"{statistics.median(e[1] for e in errs)}; |err| > {25} on either side: "
                     f"{sum(max(abs(e[0]), abs(e[1])) > 25 for e in errs)}/{len(errs)}")
        # Partial elements are classified like complete ones (README deviation 9).
        part = [r for r in rows if r["element_type"].strip().lower() in ("p", "partial")]
        if part:
            pc = Counter((r["truth_class"], r["call"]) for r in part)
            s.append(f"\nISEScan partial elements (n={len(part)}), class → call: "
                     + ", ".join(f"{c} → {k}: {v}" for (c, k), v in
                                 sorted(pc.items(), key=lambda x: (CLASS_ORDER.index(x[0][0]), x[0][1]))))
    misc = [r for r in rows if r["expect"] and r["call"] != r["expect"]]
    if misc:
        s.append("\n<details><summary>Elements not called as expected</summary>\n")
        s.append(md_table(["case", "element", "class", "expect", "call", "status", "gap", "residual",
                           "n_loci", "cov"],
                          [[r["case"], r["element_name"], r["truth_class"], r["expect"], r["call"],
                            r["status"], r["midpoint_gap"], r["anc_residual"], r["n_loci"],
                            r["construct_coverage"]] for r in misc]))
        s.append("\n</details>")
    return "\n".join(s) + "\n"


def ref_coord_summary(rows):
    """Reported reference_position/reference_end vs the known REL606 coordinates."""
    def diffs(call, pairs):
        out = []
        for r in rows:
            if r["call"] != call or r["truth_ref_start"] in ("", None) or r["expect"] != call:
                continue
            out.append(tuple(int(r[a]) - int(r[b]) for a, b in pairs if r[b] != ""))
        return out

    def fmt(ds):
        c = Counter(ds)
        exact = c.get(tuple(0 for _ in ds[0]), 0)
        common = ", ".join(f"{k}: {v}" for k, v in c.most_common(4))
        return f"exact {exact}/{len(ds)}; most common offsets {common}"

    lines = []
    anc = diffs("ANC", [("reference_position", "truth_ref_start"), ("reference_end", "truth_ref_end")])
    new = diffs("NEW", [("reference_position", "truth_ref_start")])
    if anc:
        lines.append(f"\nANC `reference_position`, `reference_end` − annotated REL606 copy (start, end): {fmt(anc)}")
    if new:
        lines.append(f"\nNEW `reference_position` − GD insertion point: {fmt(new)}")
    return "\n".join(lines)


def report(args):
    md = ["# MOBster v0.1 benchmark\n",
          "Reference: E. coli B REL606 (Barrick lab annotated GenBank). Query genomes simulated with "
          "`gdtools APPLY`. Fixed §5 constants. Tier 1 uses exact GD-derived element coordinates; "
          "tier 2 uses ISEScan 1.7.3 annotation of each query. The difference is the cost of "
          "annotation error.\n",
          "`anc_background` = untouched ancestral copies in the modified genomes (the same copies "
          "re-counted once per case).\n",
          tier_section("tier1", "Tier 1 — exact coordinates"),
          tier_section("tier2", "Tier 2 — ISEScan annotation")]
    (OUT / "report.md").write_text("\n".join(md))
    print((OUT / "report.md").read_text())


def main():
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(name)s %(message)s", stream=sys.stderr)
    p = argparse.ArgumentParser()
    p.add_argument("step", choices=["prepare", "tier1", "tier2", "report", "all"])
    p.add_argument("--threads", type=int, default=8)
    args = p.parse_args()
    steps = ["prepare", "tier1", "tier2", "report"] if args.step == "all" else [args.step]
    for st in steps:
        globals()[st](args)


if __name__ == "__main__":
    main()
