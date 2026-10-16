#!/usr/bin/env python3
"""Score a Kraken2 run against the D6331 reference composition (species level, matched by NCBI taxid).

Rules fixed here (see PROJECT_STATE D9/D10, plan Section 6):
- A species is CALLED when its clade read count in the Kraken2 report is >= T (--min-reads). Strains roll up into their species.
- Truth = the 17 species in reference/d6331_expected_reads_by_species.csv. Expected reads = total reads x genomic-DNA fraction.
- A truth species is EVALUABLE when its expected reads >= --evaluable-min-expected (default 20). Recall is computed over evaluable
  species only; the others are listed separately. A called truth species is never a false positive, evaluable or not.
- Precision = called truth species / called species. FP audit lists every non-truth species called, with a category.
- Uncertainty: with --min-reads, a bootstrap over the pool's 5,000-read BLOCKS (not reads) gives 95% intervals (needs the per-read
  Kraken2 output, the sampling manifest and the database nodes.dmp).
- The threshold can be given as a fraction of ALL reads (--min-frac, e.g. 0.0005 = 0.05%) instead of a read count; T = ceil(frac x total reads).
  A fraction keeps its meaning when the depth changes (depth series), a fixed read count does not.
- Taxonomy aliases (--aliases, default reference/d6331_taxid_aliases.csv): extra taxids whose species-level reads count as a mock species,
  for strains that NCBI now files under a different species name (for example the D6331 Fusobacterium nucleatum strain, now F. animalis).
  Each alias needs recorded evidence (see that file). Reported in the output, never silent.
- Without --min-reads/--min-frac the script prints a sweep over thresholds (development pool only). The evaluation pool is scored only with
  --allow-eval AND an explicit --min-reads: no sweep, no tuning on it.
Outputs: data/score/<pool>/<label>/ species_table.tsv, sweep.tsv (sweep mode), fp_audit.tsv, score_summary.json, score_report.txt
"""
import argparse, bisect, collections, csv, glob, gzip, json, math, os, random, sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--run", required=True, help="classify run folder, e.g. data/classify/dev/baseline_c0")
ap.add_argument("--truth", default="reference/d6331_expected_reads_by_species.csv")
ap.add_argument("--min-reads", type=int, help="detection threshold T (reads at species rank); omit for a sweep")
ap.add_argument("--min-frac", type=float, help="detection threshold as a fraction of all reads (e.g. 0.0005); T = ceil(frac x total reads)")
ap.add_argument("--sweep", default="1,2,3,5,10,20,50,100,150,200,300,500,1000,2000", help="read-count thresholds for sweep mode")
ap.add_argument("--sweep-frac", help="comma-separated fractions of all reads for sweep mode (replaces --sweep), e.g. 0.00001,0.0001,0.0005,0.001")
ap.add_argument("--aliases", default="reference/d6331_taxid_aliases.csv", help="taxid alias file (see docstring); ignored if absent")
ap.add_argument("--evaluable-min-expected", type=float, default=20.0)
ap.add_argument("--bootstrap", type=int, default=1000, help="block-bootstrap replicates (with --min-reads)")
ap.add_argument("--seed", type=int, default=20261016)
ap.add_argument("--allow-eval", action="store_true")
ap.add_argument("--outdir", help="default data/score/<pool>/<label>")
a = ap.parse_args()

run = os.path.join(ROOT, a.run)
rec = json.load(open(os.path.join(run, "run.json")))
pool, label = rec["pool"], rec["label"]
if a.min_reads is not None and a.min_frac is not None:
    sys.exit("give --min-reads or --min-frac, not both")
if pool == "eval" and not (a.allow_eval and (a.min_reads or a.min_frac)):
    sys.exit("evaluation pool: needs --allow-eval and an explicit --min-reads or --min-frac chosen on the dev pool (no sweep on eval)")
outdir = a.outdir or os.path.join(ROOT, "data", "score", pool, label)


# ---------- report ----------
def parse_report(path):
    rows, stack = [], {}
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if len(f) < 6:
            continue
        raw = f[-1]
        depth = (len(raw) - len(raw.lstrip(" "))) // 2
        stack[depth] = (f[-3], int(f[-2]), raw.strip())
        for d in [k for k in stack if k > depth]:
            del stack[d]
        rows.append({"rank": f[-3], "taxid": int(f[-2]), "name": raw.strip(), "clade": int(f[1]), "direct": int(f[2]),
                     "distinct": int(f[4]) if len(f) >= 8 else None, "lineage": [stack[d] for d in sorted(stack)]})
    return rows


rows = parse_report(os.path.join(run, "report.txt"))
unclassified = sum(r["clade"] for r in rows if r["rank"] == "U")
classified = sum(r["clade"] for r in rows if r["taxid"] == 1)
N = unclassified + classified
if a.min_frac is not None:
    a.min_reads = max(1, math.ceil(a.min_frac * N))
outdir_tag = outdir if a.min_reads is None else outdir + (f"_f{a.min_frac:g}" if a.min_frac is not None else f"_T{a.min_reads}")
os.makedirs(outdir_tag, exist_ok=True)
species = {}
for r in rows:
    if r["rank"] == "S":
        names = [x[2] for x in r["lineage"]]
        dom = next((d for d in ("Bacteria", "Archaea", "Eukaryota", "Viruses") if d in names), "other")
        genus = next((x[2] for x in r["lineage"] if x[0] == "G"), "")
        species[r["taxid"]] = {"name": r["name"], "reads": r["clade"], "distinct": r["distinct"], "domain": dom, "genus": genus}

# ---------- truth ----------
truth = {}
for t in csv.DictReader(open(os.path.join(ROOT, a.truth))):
    tid = int(t["species_taxid"])
    exp = N * float(t["genomic_dna_pct"]) / 100.0
    truth[tid] = {"name": t["ncbi_name_in_kraken_db"], "dna_pct": float(t["genomic_dna_pct"]), "expected": exp,
                  "evaluable": exp >= a.evaluable_min_expected, "domain": t["domain"], "genus": t["ncbi_name_in_kraken_db"].split()[0]}
truth_genera = {v["genus"] for v in truth.values()}
aliases = []  # (truth_taxid, accepted_taxid, accepted_name, evidence)
alias_path = os.path.join(ROOT, a.aliases)
if os.path.exists(alias_path):
    for x in csv.DictReader(open(alias_path)):
        if int(x["species_taxid"]) not in truth:
            sys.exit(f"alias file: species_taxid {x['species_taxid']} is not a truth species")
        aliases.append((int(x["species_taxid"]), int(x["accepted_taxid"]), x["accepted_name"], x["evidence"]))


def apply_aliases(c):
    c = collections.Counter(c)
    for t, acc, _n, _e in aliases:
        c[t] += c.pop(acc, 0)
    return c


genus_reads = {r["name"]: r["clade"] for r in rows if r["rank"] == "G"}
evaluable = {t for t, v in truth.items() if v["evaluable"]}
dna_sum = sum(v["dna_pct"] for v in truth.values())


def category(tid, info):
    if info["genus"] in truth_genera:
        return "same genus as a mock species"
    return {"Bacteria": "other bacteria", "Archaea": "other archaea", "Eukaryota": "eukaryote (fungi, human, other)",
            "Viruses": "virus"}.get(info["domain"], "other")


def metrics(counts, T):
    called = {t for t, c in counts.items() if c >= T}
    correct = called & set(truth)
    fp = called - set(truth)
    tp_eval = correct & evaluable
    prec = len(correct) / len(called) if called else float("nan")
    rec_e = len(tp_eval) / len(evaluable) if evaluable else float("nan")
    rec_a = len(correct) / len(truth)
    # F1 is 0 (not undefined) when nothing is called or nothing called is right: recall is 0 in that case
    f1 = (0.0 if not called else (2 * prec * rec_e / (prec + rec_e) if (prec + rec_e) > 0 else 0.0)) if evaluable else float("nan")
    obs_tot = sum(counts[t] for t in called)
    l1 = bc = float("nan")
    if obs_tot:
        keys = called | set(truth)
        o = {t: (counts.get(t, 0) / obs_tot if t in called else 0.0) for t in keys}
        e = {t: (truth[t]["dna_pct"] / dna_sum if t in truth else 0.0) for t in keys}
        l1 = sum(abs(o[t] - e[t]) for t in keys)
        bc = 1 - sum(min(o[t], e[t]) for t in keys)
    return {"T": T, "called": len(called), "correct": len(correct), "fp": len(fp), "tp_evaluable": len(tp_eval),
            "fn_evaluable": len(evaluable - tp_eval), "precision": prec, "recall_evaluable": rec_e, "recall_all_truth": rec_a,
            "f1": f1, "l1": l1, "bray_curtis": bc, "fp_taxids": sorted(fp), "missed_evaluable": sorted(evaluable - tp_eval)}


counts = apply_aliases({t: v["reads"] for t, v in species.items()})
for t, acc, nm, _e in aliases:
    if acc in species and t not in species:
        species[t] = dict(species[acc], reads=0)  # keep metadata for the table
fmt = lambda x: "nan" if (isinstance(x, float) and math.isnan(x)) else (f"{x:.3f}" if isinstance(x, float) else str(x))
summary = {"run": a.run, "pool": pool, "total_reads": N, "classified": classified, "unclassified": unclassified,
           "evaluable_min_expected": a.evaluable_min_expected, "evaluable": [truth[t]["name"] for t in sorted(evaluable)],
           "not_evaluable": [f"{v['name']} (expected {v['expected']:.1f} reads)" for t, v in truth.items() if not v["evaluable"]]}
in_truth = sum(counts.get(t, 0) for t in truth)
summary["reads_in_mock_species"] = in_truth
summary["pct_reads_in_mock_species"] = round(100 * in_truth / N, 2)
summary["pct_classified_elsewhere"] = round(100 * (classified - in_truth) / N, 2)
alias_note = [f"Alias applied: {nm} (taxid {acc}) reads count as {truth[t]['name']}; observed {species[acc]['reads'] if acc in species else 0:,} reads moved" for t, acc, nm, _e in aliases]
out = [f"Score for {a.run} ({pool} pool): {N:,} reads, {classified:,} classified, {unclassified:,} unclassified",
       f"Reads in the {len(truth)} mock species: {in_truth:,} ({summary['pct_reads_in_mock_species']}% of all reads); classified elsewhere (other taxa or above species rank): {summary['pct_classified_elsewhere']}%",
       f"Evaluable species (expected >= {a.evaluable_min_expected:g} reads): {len(evaluable)} of {len(truth)}. Not evaluable: {'; '.join(summary['not_evaluable']) or 'none'}"] + alias_note + [""]
summary["aliases"] = [{"truth": truth[t]["name"], "accepted_taxid": acc, "accepted_name": nm, "evidence": ev} for t, acc, nm, ev in aliases]


def fp_table(T):
    r = []
    for t, c in sorted(counts.items(), key=lambda kv: -kv[1]):
        if c >= T and t not in truth:
            i = species[t]
            r.append({"taxid": t, "name": i["name"], "reads": c, "distinct_minimizers": i["distinct"], "domain": i["domain"],
                      "genus": i["genus"], "category": category(t, i)})
    return r


def write_tsv(path, rows_, cols):
    with open(path, "w") as o:
        o.write("\t".join(cols) + "\n")
        for r in rows_:
            o.write("\t".join(fmt(r[c]) if not isinstance(r[c], str) else r[c] for c in cols) + "\n")


if a.min_reads is None:
    grid = sorted({max(1, math.ceil(float(x) * N)) for x in a.sweep_frac.split(",")}) if a.sweep_frac else [int(x) for x in a.sweep.split(",")]
    sw = []
    for T in grid:
        m = metrics(counts, T)
        fps = fp_table(T)
        cat = collections.Counter(x["category"] for x in fps)
        m.update({"pct_of_reads": 100.0 * T / N, "fp_same_genus": cat.get("same genus as a mock species", 0), "fp_other": len(fps) - cat.get("same genus as a mock species", 0)})
        sw.append(m)
    cols = ["T", "pct_of_reads", "called", "correct", "fp", "fp_same_genus", "fp_other", "tp_evaluable", "fn_evaluable", "precision", "recall_evaluable",
            "recall_all_truth", "f1", "bray_curtis"]
    write_tsv(os.path.join(outdir_tag, "sweep.tsv"), sw, cols)
    out.append("Threshold sweep (T = minimum reads at species rank):")
    out.append(f"{'T':>5} {'%reads':>7} {'called':>6} {'FP':>4} {'FPsame':>6} {'FPother':>7} {'prec':>6} {'rec(eval)':>9} {'rec(all)':>8} {'F1':>6} {'BC':>6}")
    for m in sw:
        out.append(f"{m['T']:>5} {m['pct_of_reads']:>7.4f} {m['called']:>6} {m['fp']:>4} {m['fp_same_genus']:>6} {m['fp_other']:>7} {fmt(m['precision']):>6} "
                   f"{fmt(m['recall_evaluable']):>9} {fmt(m['recall_all_truth']):>8} {fmt(m['f1']):>6} {fmt(m['bray_curtis']):>6}")
    summary["sweep"] = [{k: v for k, v in m.items() if k not in ("fp_taxids", "missed_evaluable")} for m in sw]
    T_show = None
else:
    T_show = a.min_reads

# species table (always) with status at T_show if given
T0 = T_show if T_show is not None else 1
tab = []
for t, v in sorted(truth.items(), key=lambda kv: -kv[1]["dna_pct"]):
    obs = counts.get(t, 0)
    status = ("detected" if obs >= T0 else "missed") + ("" if v["evaluable"] else " (not evaluable)")
    tab.append({"taxid": t, "species": v["name"], "domain": v["domain"], "dna_pct": v["dna_pct"], "expected_reads": round(v["expected"], 1),
                "observed_reads": obs, "observed_over_expected": (round(obs / v["expected"], 2) if v["expected"] else float("nan")),
                "genus_clade_reads": genus_reads.get(v["genus"], 0), "distinct_minimizers": species[t]["distinct"] if t in species else 0, "evaluable": "yes" if v["evaluable"] else "no",
                f"status_T{T0}": status})
write_tsv(os.path.join(outdir_tag, "species_table.tsv"), tab, list(tab[0]))
out.append("")
out.append(f"Mock species (status at T={T0}):")
out.append(f"{'species':32s}{'DNA%':>8s}{'expected':>10s}{'observed':>10s}{'obs/exp':>8s}  status")
for r in tab:
    out.append(f"{r['species']:32s}{r['dna_pct']:>8g}{r['expected_reads']:>10.1f}{r['observed_reads']:>10d}{fmt(r['observed_over_expected']):>8s}  {r[f'status_T{T0}']}")

if T_show is not None:
    m = metrics(counts, T_show)
    fps = fp_table(T_show)
    write_tsv(os.path.join(outdir_tag, "fp_audit.tsv"), fps, ["taxid", "name", "reads", "distinct_minimizers", "domain", "genus", "category"])
    out.append("")
    g_ok = [t for t in evaluable if genus_reads.get(truth[t]["genus"], 0) >= T_show]
    out.append(f"Threshold T={T_show} reads = {100.0 * T_show / N:.4f}% of all reads" + (f" (--min-frac {a.min_frac:g})" if a.min_frac is not None else ""))
    out.append(f"Genus-level recall (secondary): {len(g_ok)}/{len(evaluable)} evaluable species have genus clade reads >= T")
    out.append(f"At T={T_show}: called {m['called']} species; correct {m['correct']}, false positives {m['fp']}; evaluable detected {m['tp_evaluable']}/{len(evaluable)}")
    out.append(f"precision {fmt(m['precision'])}  recall(evaluable) {fmt(m['recall_evaluable'])}  recall(all {len(truth)}) {fmt(m['recall_all_truth'])}  F1 {fmt(m['f1'])}  "
               f"L1 {fmt(m['l1'])}  Bray-Curtis {fmt(m['bray_curtis'])}")
    if m["missed_evaluable"]:
        out.append("missed evaluable: " + ", ".join(truth[t]["name"] for t in m["missed_evaluable"]))
    out.append("")
    out.append(f"False-positive audit (non-mock species with >= {T_show} reads), {len(fps)} species:")
    for x in fps[:25]:
        out.append(f"  {x['name']:38s}{x['reads']:>7d} reads  [{x['category']}]")
    if len(fps) > 25:
        out.append(f"  ... {len(fps) - 25} more in fp_audit.tsv")
    summary["threshold"] = {"min_reads": T_show, "min_frac": a.min_frac, "pct_of_reads": 100.0 * T_show / N, "genus_recall_evaluable": len(g_ok) / len(evaluable)}
    summary["at_T"] = {k: v for k, v in m.items() if k not in ("fp_taxids", "missed_evaluable")}
    summary["fp_by_category"] = dict(collections.Counter(x["category"] for x in fps))

    # ---------- block bootstrap ----------
    outfile = os.path.join(run, "kraken_output.txt.gz")
    mf = glob.glob(os.path.join(ROOT, "data", "subsets", "*.sra_blocks.manifest.json"))
    nodes_path = os.path.join(ROOT, rec["database"], "nodes.dmp")
    if a.bootstrap and os.path.exists(outfile) and len(mf) == 1 and os.path.exists(nodes_path) and rec.get("input_sha256_in_manifest"):
        man = json.load(open(mf[0]))
        blocks = sorted([b for b in man["blocks"] if b["pool"] == pool], key=lambda b: b["start"])
        starts = [b["start"] for b in blocks]
        parent, rank = {}, {}
        for line in open(nodes_path):
            p = line.split("\t|\t")
            parent[int(p[0])], rank[int(p[0])] = int(p[1]), p[2]
        memo = {}

        def sp_of(t):
            if t in memo:
                return memo[t]
            x, seen = t, 0
            while x in rank and rank[x] != "species" and x != 1 and seen < 60:
                x, seen = parent[x], seen + 1
            memo[t] = x if rank.get(x) == "species" else None
            return memo[t]

        per_block = [collections.Counter() for _ in blocks]
        with gzip.open(outfile, "rt") as f:
            for line in f:
                p = line.split("\t")
                if p[0] != "C":
                    continue
                s = sp_of(int(p[2]))
                if s is None:
                    continue
                n = int(p[1].rsplit(".", 1)[1])
                bi = bisect.bisect_right(starts, n) - 1
                per_block[bi][s] += 1
        per_block = [apply_aliases(c) for c in per_block]
        tot = collections.Counter()
        for c in per_block:
            tot.update(c)
        bad = [(truth[t]["name"], tot.get(t, 0), counts.get(t, 0)) for t in truth if tot.get(t, 0) != counts.get(t, 0)]
        out.append("")
        out.append("Cross-check, report vs per-read output (species read counts of the 17 mock species): " + ("all identical" if not bad else f"MISMATCH {bad}"))
        rng = random.Random(a.seed)
        res = {k: [] for k in ("precision", "recall_evaluable", "f1", "bray_curtis", "fp")}
        det = collections.Counter()
        for _ in range(a.bootstrap):
            c = collections.Counter()
            for _b in range(len(blocks)):
                c.update(per_block[rng.randrange(len(blocks))])
            mm = metrics(c, T_show)
            for k in res:
                res[k].append(mm[k])
            for t in truth:
                det[t] += c.get(t, 0) >= T_show
        q = lambda v: (lambda s: (s[int(0.025 * (len(s) - 1))], s[int(0.975 * (len(s) - 1))]) if s else (float("nan"), float("nan")))(sorted(x for x in v if not (isinstance(x, float) and math.isnan(x))))
        out.append(f"Block bootstrap ({a.bootstrap} resamples of {len(blocks)} blocks, seed {a.seed}); 95% intervals:")
        ci, means = {}, {}
        for k, lab in (("precision", "precision"), ("recall_evaluable", "recall (evaluable)"), ("f1", "F1"), ("bray_curtis", "Bray-Curtis"), ("fp", "false positives (n)")):
            lo, hi = q(res[k]); ci[k] = [lo, hi]
            vals = [x for x in res[k] if not (isinstance(x, float) and math.isnan(x))]
            mean = sum(vals) / len(vals) if vals else float("nan")
            means[k] = mean
            out.append(f"  {lab:22s}{fmt(float(lo)):>8s} to {fmt(float(hi)):>8s}   mean {fmt(float(mean)):>7s}")
        out.append("  detection frequency of mock species across resamples (T applied): " +
                   "; ".join(f"{truth[t]['name']} {det[t] / a.bootstrap:.2f}" for t in sorted(truth, key=lambda x: det[x]) if det[t] / a.bootstrap < 0.999))
        summary["bootstrap"] = {"n": a.bootstrap, "seed": a.seed, "ci95": ci, "mean": means, "crosscheck_mismatch": bad}
    else:
        out.append("")
        out.append("Block bootstrap skipped (needs --bootstrap > 0, kraken_output.txt.gz, the sampling manifest, nodes.dmp).")

json.dump(summary, open(os.path.join(outdir_tag, "score_summary.json"), "w"), indent=2, default=str)
open(os.path.join(outdir_tag, "score_report.txt"), "w").write("\n".join(out) + "\n")
print("\n".join(out))
print(f"\nwrote {os.path.relpath(outdir_tag, ROOT)}/")
