#!/usr/bin/env python3
"""Database-size arms: does the 8 GB cap on the Kraken2 index cause the abundance bias? (D17; run on an AWS instance with >=128 GB RAM)

Registered design (docs/PROJECT_STATE.md, "Database arms"); nothing here may be changed after the dev stage has been looked at.
  Arms (everything else identical: image digest, --confidence 0.0, --minimum-hit-groups 2, same pools, same scoring code):
    pluspf8  = PlusPF-8  (7.5 GB index, 2026-06-26): baseline and CONTROL; threshold fixed at 0.0008 (D14b). Must reproduce the local baseline exactly.
    pluspf16 = PlusPF-16 (14.9 GB index, 2026-06-26)
    pluspf   = PlusPF full (109.9 GB index, 2026-06-26)
  Stage dev  : classify the DEV pool with every arm; for pluspf16 and pluspf score every fraction in --grid with the block bootstrap and choose the
               SMALLEST fraction whose bootstrap-mean F1 is within --tie (0.002) of that arm's maximum (D14b/D16 rule). Writes reports/dbarms/dev_tuning.tsv
               and dev_selection.json. The control is also checked: its dev species counts must equal reports/baseline_kraken2_pluspf8/dev_species_table.tsv.
               Then commit and tag `dbarms-frozen` BEFORE the eval stage.
  Stage eval : needs --allow-eval and the git tag dbarms-frozen. Classifies the EVAL pool with every arm and scores it ONCE per arm at its frozen fraction
               (bootstrap 1000). Also scores every arm at the baseline fraction 0.0008 for a like-for-like column. ALL arms are reported.
  Stage summarise : reports/dbarms/dbarms_summary.txt.
Outputs: data/classify/<pool>/<label>/, data/score/<pool>/<label>_f<frac>/, reports/dbarms/
"""
import argparse, hashlib, json, math, os, subprocess, sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
BASE_FRAC = 0.0008
ARMS = {"pluspf8": "data/db/pluspf_08_20260626", "pluspf16": "data/db/pluspf_16_20260626", "pluspf": "data/db/pluspf_20260626"}
ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--stage", choices=["dev", "eval", "summarise"], required=True)
ap.add_argument("--grid", default="0.00005,0.0001,0.0002,0.0005,0.0008,0.001,0.002,0.005")
ap.add_argument("--tie", type=float, default=0.002)
ap.add_argument("--bootstrap", type=int, default=1000)
ap.add_argument("--threads", type=int, default=16)
ap.add_argument("--memory", default="120g", help="docker memory limit for classification (the full index is 110 GB)")
ap.add_argument("--allow-eval", action="store_true")
ap.add_argument("--skip-tag-check", action="store_true", help="testing only")
ap.add_argument("--dry-run", action="store_true")
a = ap.parse_args()
grid = [float(x) for x in a.grid.split(",")]
outdir = os.path.join(ROOT, "reports", "dbarms"); os.makedirs(outdir, exist_ok=True)
sel_path = os.path.join(outdir, "dev_selection.json")
py = sys.executable
lab = lambda k: f"db_{k}"

def run(cmd):
    print(" ".join(cmd), flush=True)
    return 0 if a.dry_run else subprocess.run(cmd).returncode

def db_record(k):
    d = os.path.join(ROOT, ARMS[k]); rec = {"arm": k, "dir": ARMS[k], "files": {}}
    for fn in ("hash.k2d", "opts.k2d", "taxo.k2d"):
        p = os.path.join(d, fn)
        if os.path.exists(p):
            rec["files"][fn] = {"bytes": os.path.getsize(p)}
            if fn != "hash.k2d": rec["files"][fn]["md5"] = hashlib.md5(open(p, "rb").read()).hexdigest()
    return rec

def classify(pool, k):
    if os.path.exists(os.path.join(ROOT, "data", "classify", pool, lab(k), "run.json")): return
    cmd = [py, os.path.join(ROOT, "bin", "03_classify.py"), "--pool", pool, "--label", lab(k), "--db", ARMS[k], "--threads", str(a.threads), "--memory", a.memory] + (["--allow-eval"] if pool == "eval" else [])
    if run(cmd) != 0: sys.exit(f"classification failed: {pool} {k}")

def score_dir(pool, k, f): return os.path.join(ROOT, "data", "score", pool, f"{lab(k)}_f{f:g}")

def score(pool, k, f):
    if os.path.exists(os.path.join(score_dir(pool, k, f), "score_summary.json")): return
    cmd = [py, os.path.join(ROOT, "bin", "04_score.py"), "--run", f"data/classify/{pool}/{lab(k)}", "--min-frac", str(f), "--bootstrap", str(a.bootstrap)] + (["--allow-eval"] if pool == "eval" else [])
    if run(cmd) != 0: sys.exit(f"scoring failed: {pool} {k} {f}")

def load(pool, k, f):
    d = score_dir(pool, k, f); s = json.load(open(os.path.join(d, "score_summary.json")))
    p = os.path.join(d, "fp_audit.tsv"); fps = []
    if os.path.exists(p):
        fps = [f"{r[1]} ({r[2]})" for r in [l.rstrip("\n").split("\t") for l in open(p)][1:]]
    return s, fps

def species_obs(pool, k, f, names=("Roseburia hominis", "Bifidobacterium adolescentis", "Escherichia coli", "Limosilactobacillus fermentum")):
    out = {}
    for l in open(os.path.join(score_dir(pool, k, f), "species_table.tsv")):
        c = l.rstrip("\n").split("\t")
        if len(c) > 7 and c[1] in names: out[c[1]] = c[6]
    return out

def row(pool, k, f):
    s, fps = load(pool, k, f); at, bs = s["at_T"], s.get("bootstrap", {})
    ci = (bs.get("ci95") or {}).get("f1", [None, None])
    return {"T": s["threshold"]["min_reads"], "classified_pct": 100.0 * s["classified"] / s["total_reads"], "mock_pct": s["pct_reads_in_mock_species"],
            "precision": at["precision"], "recall_eval": at["recall_evaluable"], "f1": at["f1"], "f1_boot_mean": (bs.get("mean") or {}).get("f1"),
            "f1_lo": ci[0], "f1_hi": ci[1], "fp": at["fp"], "fp_names": fps, "bray_curtis": at["bray_curtis"], "genus_recall": s["threshold"].get("genus_recall_evaluable")}

def fx(x, nd=3): return "n/a" if x is None or (isinstance(x, float) and math.isnan(x)) else f"{x:.{nd}f}"

def wall(pool, k):
    p = os.path.join(ROOT, "data", "classify", pool, lab(k), "run.json")
    return json.load(open(p))["wall_seconds"] if os.path.exists(p) else None

if a.stage == "dev":
    for k in ARMS: classify("dev", k)
    tuning, selection = [], {}
    for k in ARMS:
        fr = [BASE_FRAC] if k == "pluspf8" else grid
        for f in fr: score("dev", k, f)
        if a.dry_run: continue
        means = {f: load("dev", k, f)[0]["bootstrap"]["mean"]["f1"] for f in fr}
        mx = max(means.values()); chosen = min(f for f in fr if means[f] >= mx - a.tie)
        for f in fr: tuning.append((k, f, means[f]))
        selection[k] = {"label": lab(k), "fraction": BASE_FRAC if k == "pluspf8" else chosen, "dev_max_boot_mean_f1": mx,
                        "rule": "baseline frozen by D14b" if k == "pluspf8" else f"smallest fraction with bootstrap-mean F1 >= max - {a.tie}",
                        "database": db_record(k), "dev_wall_seconds": wall("dev", k)}
    if not a.dry_run:
        ref = os.path.join(ROOT, "reports", "baseline_kraken2_pluspf8", "dev_species_table.tsv")
        mine = os.path.join(score_dir("dev", "pluspf8", BASE_FRAC), "species_table.tsv")
        cols = lambda path: [l.rstrip("\n").split("\t")[:8] for l in open(path)]   # taxid..genus_clade_reads; the status column depends on T and is left out
        ctl = os.path.exists(ref) and cols(ref) == cols(mine)
        print(f"\nCONTROL (PlusPF-8 on this machine reproduces the local dev baseline species counts exactly): {'YES' if ctl else 'NO  <-- stop and investigate before using the other arms'}")
        with open(os.path.join(outdir, "dev_tuning.tsv"), "w") as o:
            o.write("arm\tfraction\tboot_mean_f1\n")
            for k, f, m in tuning: o.write(f"{k}\t{f:g}\t{m:.4f}\n")
        json.dump({"arms": selection, "grid": grid, "tie": a.tie, "bootstrap": a.bootstrap, "control_reproduced": bool(ctl)}, open(sel_path, "w"), indent=2)
        print("\nDev tuning (bootstrap-mean F1 by arm and fraction):\narm        " + "".join(f"{f:>9g}" for f in grid) + "   chosen")
        for k in ARMS:
            ms = {f: m for kk, f, m in tuning if kk == k}
            print(f"{k:<10} " + "".join((f"{ms[f]:>9.3f}" if f in ms else f"{'-':>9}") for f in grid) + f"   {selection[k]['fraction']:g}")
        print(f"\nwrote {os.path.relpath(sel_path, ROOT)}. NEXT: git add reports && git commit && git tag dbarms-frozen, then --stage eval --allow-eval")

if a.stage == "eval":
    if not a.allow_eval: sys.exit("eval stage needs --allow-eval")
    if not os.path.exists(sel_path): sys.exit("run --stage dev first (reports/dbarms/dev_selection.json missing)")
    if not json.load(open(sel_path)).get("control_reproduced"): sys.exit("the control arm did not reproduce the local baseline; do not run eval")
    if not a.skip_tag_check:
        if subprocess.run(["git", "-C", ROOT, "rev-parse", "-q", "--verify", "refs/tags/dbarms-frozen"], capture_output=True).returncode != 0:
            sys.exit("git tag dbarms-frozen not found: commit reports/dbarms and tag it before touching the evaluation pool")
        if subprocess.run(["git", "-C", ROOT, "status", "--porcelain", "reports/dbarms"], capture_output=True, text=True).stdout.strip():
            sys.exit("reports/dbarms has uncommitted changes since the freeze")
    sel = json.load(open(sel_path))["arms"]
    for k in ARMS:
        classify("eval", k); score("eval", k, sel[k]["fraction"])
        if sel[k]["fraction"] != BASE_FRAC: score("eval", k, BASE_FRAC)

if a.stage in ("eval", "summarise") and not a.dry_run:
    sel = json.load(open(sel_path))["arms"]
    out = ["Database-size arms (registered design; thresholds chosen on dev, eval scored once per arm)", ""]
    for pool in ("dev", "eval"):
        for title, fsel in (("at each arm's frozen fraction", lambda k: sel[k]["fraction"]), ("at the baseline fraction 0.0008 (like for like)", lambda k: BASE_FRAC)):
            out.append(f"{pool.upper()} pool, {title}:")
            out.append(f"{'arm':<9} {'frac':>8} {'T':>5} {'classif%':>8} {'mock%':>6} {'prec':>6} {'rec(ev)':>7} {'F1':>6} {'F1 boot mean (95%)':>24} {'FP':>3} {'BC':>6} {'genus rec':>9} {'wall s':>6}  false positives")
            for k in ARMS:
                f = fsel(k)
                if not os.path.exists(os.path.join(score_dir(pool, k, f), "score_summary.json")):
                    out.append(f"{k:<9} {f:>8g}  (not run)"); continue
                r = row(pool, k, f); ci = f"{fx(r['f1_boot_mean'])} ({fx(r['f1_lo'])}-{fx(r['f1_hi'])})"
                out.append(f"{k:<9} {f:>8g} {r['T']:>5} {r['classified_pct']:>8.1f} {r['mock_pct']:>6.1f} {fx(r['precision']):>6} {fx(r['recall_eval']):>7} {fx(r['f1']):>6} {ci:>24} {r['fp']:>3} {fx(r['bray_curtis']):>6} {fx(r['genus_recall']):>9} {fx(wall(pool, k), 0):>6}  {'; '.join(r['fp_names'][:4])}")
            out.append("")
        out.append(f"{pool.upper()} observed/expected reads at 0.0008 (the abundance-bias species):")
        for k in ARMS:
            if os.path.exists(os.path.join(score_dir(pool, k, BASE_FRAC), "species_table.tsv")):
                out.append(f"  {k:<9} " + "; ".join(f"{n} {v}" for n, v in species_obs(pool, k, BASE_FRAC).items()))
        out.append("")
    out.append("dev tuning table: reports/dbarms/dev_tuning.tsv. pluspf8 is the frozen baseline and the control (fraction 0.0008, D14b).")
    open(os.path.join(outdir, "dbarms_summary.txt"), "w").write("\n".join(out) + "\n")
    print("\n".join(out)); print("wrote reports/dbarms/dbarms_summary.txt")
