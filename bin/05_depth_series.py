#!/usr/bin/env python3
"""Depth series: random subsamples of one read pool at fixed depths (in Gb), classified and scored with the FROZEN rules.

What it does (all resumable; finished outputs are skipped):
  1. subsample : for each depth (default 0.1, 0.25, 0.5 Gb) and replicate (default 5) draw reads at random, without replacement,
                 until the bases reach the target (D9: depth is defined in Gb, not reads). Seeds are recorded in the manifest.
                 Replicates come from a finite pool, so they overlap (about half at 0.5 Gb) and measure within-pool variation only (D10).
  2. classify  : bin/03_classify.py on each subsample (frozen baseline: Kraken2 2.17.2 + PlusPF-8, confidence 0).
  3. score     : bin/04_score.py --min-frac 0.0008 --bootstrap 0 on each (D14b: threshold is a fraction of reads, frozen on the dev pool;
                 evaluable minimum 20 expected reads at THAT depth, D14a). No sweep, no tuning.
  4. summarise : per-depth mean and range over replicates, per-species detection frequency. The full pool run (baseline_c0) is
                 added as the single 1 Gb baseline if it exists.
Guard: --pool eval needs --allow-eval. Nothing here may change a rule; if a rule looks wrong at low depth, that is a finding to record.
Outputs: data/depth/<pool>/ (FASTQs, depth_manifest.json, depth_summary.tsv, depth_species.tsv, depth_report.txt)
         data/classify/custom/<pool>_d<depth>_r<rep>/ and data/score/custom/<same>_f0.0008/
Usage: bin/05_depth_series.py --pool eval --allow-eval [--dry-run]
"""
import argparse, glob, gzip, hashlib, json, math, os, random, subprocess, sys, time

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
FROZEN_FRAC = 0.0008   # D14b, chosen on the dev pool by bootstrap-mean F1; do not change here

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--pool", choices=["dev", "eval"], required=True)
ap.add_argument("--allow-eval", action="store_true")
ap.add_argument("--depths", default="0.1,0.25,0.5", help="depths in Gb")
ap.add_argument("--replicates", type=int, default=5)
ap.add_argument("--seed", type=int, default=20261016)
ap.add_argument("--min-frac", type=float, default=FROZEN_FRAC, help="detection threshold as a fraction of reads (frozen, D14b)")
ap.add_argument("--stages", default="subsample,classify,score,summarise")
ap.add_argument("--pool-fastq", help="testing only: read this FASTQ(.gz) instead of the pool in the manifest")
ap.add_argument("--dry-run", action="store_true")
a = ap.parse_args()
if a.pool == "eval" and not a.allow_eval:
    sys.exit("refusing to use the evaluation pool without --allow-eval (rules must already be frozen; see D14b)")
if a.min_frac != FROZEN_FRAC:
    print(f"WARNING: --min-frac {a.min_frac} differs from the frozen value {FROZEN_FRAC}; this run is not the registered analysis", file=sys.stderr)
stages = a.stages.split(",")
depths = [float(x) for x in a.depths.split(",")]
data = os.path.join(ROOT, "data")
outdir = os.path.join(data, "depth", a.pool)
os.makedirs(outdir, exist_ok=True)
tag = lambda d, r: f"{a.pool}_d{d:g}_r{r}"
fq_path = lambda d, r: os.path.join(outdir, f"{tag(d, r)}.fastq.gz")

if a.pool_fastq:
    src = os.path.join(ROOT, a.pool_fastq) if not os.path.isabs(a.pool_fastq) else a.pool_fastq
else:
    mf = glob.glob(os.path.join(data, "subsets", "*.sra_blocks.manifest.json"))
    if len(mf) != 1:
        sys.exit(f"expected exactly one manifest in data/subsets, found {len(mf)}")
    man = json.load(open(mf[0])); src = os.path.join(data, "subsets", man["pools"][a.pool]["file"])
if not os.path.exists(src):
    sys.exit(f"pool FASTQ not found: {src}")

# ---------------- 1. subsample ----------------
def read_fastq(path):
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as f:
        while True:
            h = f.readline()
            if not h:
                return
            s, p, q = f.readline(), f.readline(), f.readline()
            yield h, s, p, q

def subsample():
    mpath = os.path.join(outdir, "depth_manifest.json")
    todo = [(d, r) for d in depths for r in range(1, a.replicates + 1) if not os.path.exists(fq_path(d, r))]
    if not todo:
        print("subsample: all files exist, skipping"); return
    lens = [len(s.strip()) for _h, s, _p, _q in read_fastq(src)]
    n = len(lens); total = sum(lens)
    print(f"pool {a.pool}: {n:,} reads, {total / 1e9:.3f} Gb")
    sel, info = {}, {}
    for d, r in todo:
        target = d * 1e9
        if target > total:
            sys.exit(f"depth {d} Gb exceeds the pool ({total / 1e9:.3f} Gb)")
        rng = random.Random(f"{a.seed}:{a.pool}:{d:g}:{r}")
        order = list(range(n)); rng.shuffle(order)
        acc, chosen = 0, set()
        for i in order:
            chosen.add(i); acc += lens[i]
            if acc >= target:
                break
        sel[(d, r)] = chosen; info[(d, r)] = {"depth_gb": d, "replicate": r, "seed": f"{a.seed}:{a.pool}:{d:g}:{r}", "reads": len(chosen), "bases": acc}
    if a.dry_run:
        for k, v in info.items(): print("would write", tag(*k), v["reads"], "reads", v["bases"] / 1e9, "Gb")
        return
    outs = {k: gzip.open(fq_path(*k) + ".part", "wt", compresslevel=4) for k in sel}
    for i, rec in enumerate(read_fastq(src)):
        for k, ch in sel.items():
            if i in ch:
                outs[k].write("".join(rec))
    for k, o in outs.items():
        o.close(); os.replace(fq_path(*k) + ".part", fq_path(*k))
        h = hashlib.sha256(open(fq_path(*k), "rb").read()).hexdigest(); info[k]["sha256"] = h; info[k]["file"] = os.path.basename(fq_path(*k))
    old = json.load(open(mpath)) if os.path.exists(mpath) else {"entries": []}
    old["pool"] = a.pool; old["source"] = os.path.relpath(src, ROOT) if src.startswith(ROOT) else src
    old["pool_reads"] = n; old["pool_bases"] = total; old["seed"] = a.seed
    old["entries"] += list(info.values())
    json.dump(old, open(mpath, "w"), indent=2)
    print(f"subsample: wrote {len(info)} files")

# ---------------- 2/3. classify and score ----------------
def run(cmd):
    print(" ".join(cmd)); 
    if a.dry_run: return 0
    return subprocess.run(cmd).returncode

def classify():
    for d in depths:
        for r in range(1, a.replicates + 1):
            lab = tag(d, r)
            if os.path.exists(os.path.join(data, "classify", "custom", lab, "run.json")):
                continue
            if run([sys.executable, os.path.join(ROOT, "bin", "03_classify.py"), "--fastq", os.path.relpath(fq_path(d, r), ROOT), "--label", lab]) != 0:
                sys.exit(f"classification failed for {lab}")

def score():
    for d in depths:
        for r in range(1, a.replicates + 1):
            lab = tag(d, r)
            if os.path.exists(os.path.join(data, "score", "custom", f"{lab}_f{a.min_frac:g}", "score_summary.json")):
                continue
            if run([sys.executable, os.path.join(ROOT, "bin", "04_score.py"), "--run", f"data/classify/custom/{lab}",
                    "--min-frac", str(a.min_frac), "--bootstrap", "0"]) != 0:
                sys.exit(f"scoring failed for {lab}")

# ---------------- 4. summarise ----------------
def mean(v): return sum(v) / len(v) if v else float("nan")
def fx(x, nd=3): return "nan" if (isinstance(x, float) and math.isnan(x)) else f"{x:.{nd}f}"

def load_score(dirpath):
    sj = os.path.join(dirpath, "score_summary.json")
    if not os.path.exists(sj): return None
    s = json.load(open(sj))
    rows = {}
    with open(os.path.join(dirpath, "species_table.tsv")) as f:
        head = f.readline().rstrip("\n").split("\t")
        for line in f:
            v = dict(zip(head, line.rstrip("\n").split("\t")))
            rows[v["species"]] = v
    return s, rows

def summarise():
    cols = ["depth_gb", "n_runs", "reads_mean", "T_reads_mean", "evaluable_species", "precision_mean", "precision_min", "precision_max",
            "recall_eval_mean", "recall_eval_min", "recall_eval_max", "f1_mean", "f1_min", "f1_max", "fp_mean", "fp_min", "fp_max", "bray_curtis_mean"]
    table, sp_rows, out = [], [], []
    sets = [(d, [os.path.join(data, "score", "custom", f"{tag(d, r)}_f{a.min_frac:g}") for r in range(1, a.replicates + 1)]) for d in depths]
    base = os.path.join(data, "score", a.pool, f"baseline_c0_f{a.min_frac:g}")
    if os.path.exists(base): sets.append((1.0, [base]))
    out.append(f"Depth series, {a.pool} pool, frozen threshold {a.min_frac:g} of reads (D14b), evaluable minimum 20 expected reads at each depth (D14a)")
    out.append("Replicates overlap and come from one pool: spread = within-pool variation only. 1.0 Gb row = full pool, a single run (no replicates).")
    out.append(f"{'Gb':>5} {'runs':>4} {'reads':>8} {'T':>5} {'eval sp':>7} {'prec mean (min-max)':>22} {'recall(eval) mean (min-max)':>29} {'F1 mean (min-max)':>22} {'FP mean (min-max)':>19} {'BC':>6}")
    for d, dirs in sets:
        got = [load_score(p) for p in dirs]; got = [g for g in got if g]
        if not got: continue
        at = [g[0]["at_T"] for g in got]; th = [g[0]["threshold"] for g in got]
        g_ = lambda k: [x[k] for x in at]
        ev = len(got[0][0]["evaluable"]); N = mean([g[0]["total_reads"] for g in got])
        row = {"depth_gb": d, "n_runs": len(got), "reads_mean": N, "T_reads_mean": mean([t["min_reads"] for t in th]), "evaluable_species": ev,
               "precision_mean": mean(g_("precision")), "precision_min": min(g_("precision")), "precision_max": max(g_("precision")),
               "recall_eval_mean": mean(g_("recall_evaluable")), "recall_eval_min": min(g_("recall_evaluable")), "recall_eval_max": max(g_("recall_evaluable")),
               "f1_mean": mean(g_("f1")), "f1_min": min(g_("f1")), "f1_max": max(g_("f1")),
               "fp_mean": mean(g_("fp")), "fp_min": min(g_("fp")), "fp_max": max(g_("fp")), "bray_curtis_mean": mean(g_("bray_curtis"))}
        table.append(row)
        rng = lambda m, lo, hi: f"{fx(m)} ({fx(lo)}-{fx(hi)})"
        fp_s = f"{row['fp_mean']:.1f} ({row['fp_min']}-{row['fp_max']})"
        out.append(f"{d:>5g} {len(got):>4} {N:>8,.0f} {row['T_reads_mean']:>5.0f} {ev:>7} {rng(row['precision_mean'], row['precision_min'], row['precision_max']):>22} "
                   f"{rng(row['recall_eval_mean'], row['recall_eval_min'], row['recall_eval_max']):>29} {rng(row['f1_mean'], row['f1_min'], row['f1_max']):>22} "
                   f"{fp_s:>19} {fx(row['bray_curtis_mean']):>6}")
        # per-species: detection frequency and mean obs/exp, over runs where the species is evaluable at that depth
        names = list(got[0][1].keys())
        for sp in names:
            flags = [(g[1][sp]["evaluable"] == "yes", float(g[1][sp]["observed_reads"]) >= t["min_reads"], float(g[1][sp]["observed_over_expected"]))
                     for g, t in zip(got, th) if sp in g[1]]
            ev_runs = [f for f in flags if f[0]]
            sp_rows.append({"depth_gb": d, "species": sp, "evaluable_in_runs": f"{len(ev_runs)}/{len(flags)}",
                            "detected_in_evaluable_runs": f"{sum(1 for f in ev_runs if f[1])}/{len(ev_runs)}" if ev_runs else "n/a",
                            "obs_over_exp_mean": mean([f[2] for f in flags])})
    out.append("")
    out.append("Species that were evaluable but not detected in every run (detections/evaluable runs):")
    seen = False
    for r in sp_rows:
        if r["detected_in_evaluable_runs"] not in ("n/a",) and r["detected_in_evaluable_runs"].split("/")[0] != r["detected_in_evaluable_runs"].split("/")[1]:
            out.append(f"  {r['depth_gb']:g} Gb  {r['species']:32s} {r['detected_in_evaluable_runs']}  (obs/exp mean {r['obs_over_exp_mean']:.2f})"); seen = True
    if not seen: out.append("  none")
    with open(os.path.join(outdir, "depth_summary.tsv"), "w") as f:
        f.write("\t".join(cols) + "\n")
        for r in table: f.write("\t".join(fx(r[c], 4) if isinstance(r[c], float) else str(r[c]) for c in cols) + "\n")
    with open(os.path.join(outdir, "depth_species.tsv"), "w") as f:
        k = list(sp_rows[0].keys()) if sp_rows else []
        f.write("\t".join(k) + "\n")
        for r in sp_rows: f.write("\t".join(fx(r[c], 3) if isinstance(r[c], float) else str(r[c]) for c in k) + "\n")
    open(os.path.join(outdir, "depth_report.txt"), "w").write("\n".join(out) + "\n")
    print("\n".join(out)); print(f"\nwrote {os.path.relpath(outdir, ROOT)}/depth_summary.tsv, depth_species.tsv, depth_report.txt")

if "subsample" in stages: subsample()
if "classify" in stages: classify()
if "score" in stages: score()
if "summarise" in stages and not a.dry_run: summarise()
