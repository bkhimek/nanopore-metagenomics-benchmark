#!/usr/bin/env python3
"""Kraken2 --confidence arms: does a stricter per-read confidence remove false positives without losing rare species?

Registered design (see docs/PROJECT_STATE.md, "Confidence arms"); nothing here may be changed after the dev stage has been looked at.
  Arms: --confidence 0.0 (the frozen baseline, label baseline_c0, threshold stays 0.0008 per D14b) and each value in --arms (default 0.05,0.1,0.2,0.4).
  Stage dev  : classify the DEV pool with each arm; for each arm score every fraction in --grid with the block bootstrap; choose the arm's threshold
               = the SMALLEST fraction whose bootstrap-mean F1 is within --tie (0.002) of that arm's maximum (D14b rule, tie to the smaller).
               Writes reports/arms/dev_tuning.tsv and reports/arms/dev_selection.json. Then commit and tag `arms-frozen` BEFORE the eval stage.
  Stage eval : needs --allow-eval and the git tag arms-frozen. Classifies the EVAL pool with every arm and scores it ONCE at that arm's frozen fraction
               (bootstrap 1000). ALL arms are reported; there is no picking of a best arm on eval.
  Stage summarise : reports/arms/arms_summary.txt (dev tuning table and eval table side by side).
Outputs: data/classify/<pool>/c<conf>/, data/score/<pool>/c<conf>_f<frac>/, reports/arms/
Note: with human as an unavoidable false positive the F1 ceiling on this mock is 0.966 (14 true calls, 1 false), so arms will tie on F1;
the informative columns are the smallest usable threshold, reads classified, reads in the mock species, Bray-Curtis and the false-positive list.
"""
import argparse, json, math, os, subprocess, sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
BASE_LABEL, BASE_FRAC = "baseline_c0", 0.0008
ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--stage", choices=["dev", "eval", "summarise"], required=True)
ap.add_argument("--arms", default="0.05,0.1,0.2,0.4", help="confidence values besides the 0.0 baseline")
ap.add_argument("--grid", default="0.00005,0.0001,0.0002,0.0005,0.0008,0.001,0.002,0.005")
ap.add_argument("--tie", type=float, default=0.002)
ap.add_argument("--bootstrap", type=int, default=1000)
ap.add_argument("--allow-eval", action="store_true")
ap.add_argument("--skip-tag-check", action="store_true", help="testing only")
ap.add_argument("--dry-run", action="store_true")
a = ap.parse_args()
arms = [float(x) for x in a.arms.split(",")]
grid = [float(x) for x in a.grid.split(",")]
lab = lambda c: BASE_LABEL if c == 0 else f"c{c:g}"
outdir = os.path.join(ROOT, "reports", "arms"); os.makedirs(outdir, exist_ok=True)
sel_path = os.path.join(outdir, "dev_selection.json")
py = sys.executable

def run(cmd):
    print(" ".join(cmd), flush=True)
    if a.dry_run: return 0
    return subprocess.run(cmd).returncode

def classify(pool, c):
    if os.path.exists(os.path.join(ROOT, "data", "classify", pool, lab(c), "run.json")): return
    cmd = [py, os.path.join(ROOT, "bin", "03_classify.py"), "--pool", pool, "--label", lab(c), "--confidence", str(c)] + (["--allow-eval"] if pool == "eval" else [])
    if run(cmd) != 0: sys.exit(f"classification failed: {pool} {lab(c)}")

def score_dir(pool, c, f): return os.path.join(ROOT, "data", "score", pool, f"{lab(c)}_f{f:g}")

def score(pool, c, f):
    if os.path.exists(os.path.join(score_dir(pool, c, f), "score_summary.json")): return
    cmd = [py, os.path.join(ROOT, "bin", "04_score.py"), "--run", f"data/classify/{pool}/{lab(c)}", "--min-frac", str(f), "--bootstrap", str(a.bootstrap)] + (["--allow-eval"] if pool == "eval" else [])
    if run(cmd) != 0: sys.exit(f"scoring failed: {pool} {lab(c)} {f}")

def load(pool, c, f):
    d = score_dir(pool, c, f); s = json.load(open(os.path.join(d, "score_summary.json")))
    fps = []
    p = os.path.join(d, "fp_audit.tsv")
    if os.path.exists(p):
        rows = [l.rstrip("\n").split("\t") for l in open(p)][1:]
        fps = [f"{r[1]} ({r[2]})" for r in rows]
    return s, fps

def row(pool, c, f):
    s, fps = load(pool, c, f)
    at, bs = s["at_T"], s.get("bootstrap", {})
    return {"arm": c, "fraction": f, "T": s["threshold"]["min_reads"], "classified_pct": 100.0 * s["classified"] / s["total_reads"],
            "mock_pct": s["pct_reads_in_mock_species"], "evaluable": len(s["evaluable"]), "precision": at["precision"], "recall_eval": at["recall_evaluable"],
            "f1": at["f1"], "f1_boot_mean": (bs.get("mean") or {}).get("f1"), "f1_lo": (bs.get("ci95") or {}).get("f1", [None, None])[0],
            "f1_hi": (bs.get("ci95") or {}).get("f1", [None, None])[1], "fp": at["fp"], "fp_names": fps, "bray_curtis": at["bray_curtis"], "genus_recall": s["threshold"].get("genus_recall_evaluable")}

def fx(x, nd=3):
    return "n/a" if x is None or (isinstance(x, float) and math.isnan(x)) else f"{x:.{nd}f}"

# ---------------- dev ----------------
if a.stage == "dev":
    for c in [0.0] + arms: classify("dev", c)
    tuning, selection = [], {}
    for c in [0.0] + arms:
        for f in grid: score("dev", c, f)
        if a.dry_run: continue
        means = {f: load("dev", c, f)[0]["bootstrap"]["mean"]["f1"] for f in grid}
        mx = max(means.values())
        chosen = min(f for f in grid if means[f] >= mx - a.tie)
        for f in grid: tuning.append((c, f, means[f]))
        selection[f"{c:g}"] = {"label": lab(c), "fraction": BASE_FRAC if c == 0 else chosen, "dev_max_boot_mean_f1": mx,
                               "rule": "baseline frozen by D14b" if c == 0 else f"smallest fraction with bootstrap-mean F1 >= max - {a.tie}",
                               "dev_retuned_fraction_for_reference": chosen}
    if not a.dry_run:
        with open(os.path.join(outdir, "dev_tuning.tsv"), "w") as o:
            o.write("confidence\tfraction\tboot_mean_f1\n")
            for c, f, m in tuning: o.write(f"{c:g}\t{f:g}\t{m:.4f}\n")
        json.dump({"arms": selection, "grid": grid, "tie": a.tie, "bootstrap": a.bootstrap}, open(sel_path, "w"), indent=2)
        print("\nDev tuning (bootstrap-mean F1 by arm and fraction):")
        print("conf   " + "".join(f"{f:>9g}" for f in grid) + "   chosen")
        for c in [0.0] + arms:
            ms = {f: m for cc, f, m in tuning if cc == c}
            print(f"{c:<6g} " + "".join(f"{ms[f]:>9.3f}" for f in grid) + f"   {selection[f'{c:g}']['fraction']:g}" + ("  (baseline frozen; re-tuned would be %g)" % selection[f"{c:g}"]["dev_retuned_fraction_for_reference"] if c == 0 else ""))
        print(f"\nwrote {os.path.relpath(sel_path, ROOT)}. NEXT: git add reports && git commit && git tag arms-frozen, then --stage eval --allow-eval")

# ---------------- eval ----------------
if a.stage == "eval":
    if not a.allow_eval: sys.exit("eval stage needs --allow-eval")
    if not os.path.exists(sel_path): sys.exit("run --stage dev first (reports/arms/dev_selection.json missing)")
    if not a.skip_tag_check:
        if subprocess.run(["git", "-C", ROOT, "rev-parse", "-q", "--verify", "refs/tags/arms-frozen"], capture_output=True).returncode != 0:
            sys.exit("git tag arms-frozen not found: commit reports/arms and tag it before touching the evaluation pool")
        if subprocess.run(["git", "-C", ROOT, "status", "--porcelain", "reports/arms"], capture_output=True, text=True).stdout.strip():
            sys.exit("reports/arms has uncommitted changes since the freeze")
    sel = json.load(open(sel_path))["arms"]
    for c in [0.0] + arms:
        classify("eval", c); score("eval", c, sel[f"{c:g}"]["fraction"])

# ---------------- summarise ----------------
if a.stage in ("eval", "summarise") and not a.dry_run:
    sel = json.load(open(sel_path))["arms"]
    out = ["Kraken2 confidence arms (registered design; thresholds chosen on dev, eval scored once per arm)", ""]
    for pool in ("dev", "eval"):
        out.append(f"{pool.upper()} pool, at each arm's frozen fraction:")
        out.append(f"{'conf':>5} {'frac':>8} {'T':>5} {'classif%':>8} {'mock%':>6} {'prec':>6} {'rec(ev)':>7} {'F1':>6} {'F1 boot mean (95%)':>26} {'FP':>3} {'BC':>6} {'genus rec':>9}  false positives")
        for c in [0.0] + arms:
            f = sel[f"{c:g}"]["fraction"]
            if not os.path.exists(os.path.join(score_dir(pool, c, f), "score_summary.json")):
                out.append(f"{c:>5g} {f:>8g}  (not run)"); continue
            r = row(pool, c, f)
            ci = f"{fx(r['f1_boot_mean'])} ({fx(r['f1_lo'])}-{fx(r['f1_hi'])})"
            out.append(f"{c:>5g} {f:>8g} {r['T']:>5} {r['classified_pct']:>8.1f} {r['mock_pct']:>6.1f} {fx(r['precision']):>6} {fx(r['recall_eval']):>7} {fx(r['f1']):>6} {ci:>26} {r['fp']:>3} {fx(r['bray_curtis']):>6} {fx(r['genus_recall']):>9}  {'; '.join(r['fp_names'][:4])}")
        out.append("")
    out.append("dev tuning table: reports/arms/dev_tuning.tsv. The 0.0 arm is the frozen baseline (fraction 0.0008, D14b).")
    open(os.path.join(outdir, "arms_summary.txt"), "w").write("\n".join(out) + "\n")
    print("\n".join(out)); print("wrote reports/arms/arms_summary.txt")
