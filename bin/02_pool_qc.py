#!/usr/bin/env python3
"""Pool QC for the dev/eval read pools (no Docker, stdlib + numpy; scipy optional).

Reads the block files named in the sampling manifest, so each read keeps its block and pool labels. Writes to data/qc/:
  read_stats.tsv.gz   one row per read: block, pool, length, mean_q (error-probability mean, NanoPlot-style)
  block_stats.tsv     one row per block
  pool_summary.json   per pool summary and dev-vs-eval comparison
  pool_qc_report.txt  the same in readable form (copy this to OneDrive)
Read quality is the usual long-read definition: Q = -10 log10(mean(10^(-q/10))) over the read's bases (Phred+33).
Reads in a block are correlated, so the comparison uses blocks (37 vs 37) as the unit as well as reads; p-values
over reads would be too small and are not reported.
"""
import argparse, glob, gzip, json, os, sys
import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
ap = argparse.ArgumentParser()
ap.add_argument("--manifest", default=None, help="default: the only *.sra_blocks.manifest.json in data/subsets")
ap.add_argument("--subsets", default=os.path.join(ROOT, "data", "subsets"))
ap.add_argument("--outdir", default=os.path.join(ROOT, "data", "qc"))
a = ap.parse_args()

mp = a.manifest or (glob.glob(os.path.join(a.subsets, "*.sra_blocks.manifest.json")) or [None])[0]
if not mp:
    sys.exit("no manifest found; pass --manifest")
man = json.load(open(mp))
run = man["run"]
os.makedirs(a.outdir, exist_ok=True)
LUT = 10.0 ** (-(np.arange(256) - 33) / 10.0)      # error probability per quality character


def n50(lengths):
    s = np.sort(lengths)[::-1]
    c = np.cumsum(s)
    return int(s[np.searchsorted(c, c[-1] / 2.0)])


rows = []          # per-read (block_index, pool, length, meanq)
blocks = []
qmin, qmax = 255, 0
for bi, b in enumerate(man["blocks"]):
    path = os.path.join(a.subsets, "blocks", f"{run}.{b['start']}-{b['end']}.fastq.gz")
    L, Q = [], []
    with gzip.open(path, "rb") as f:
        while True:
            h = f.readline()
            if not h:
                break
            s = f.readline().rstrip(b"\n"); f.readline(); q = np.frombuffer(f.readline().rstrip(b"\n"), dtype=np.uint8)
            if len(q) != len(s):
                sys.exit(f"length mismatch in {path}")
            L.append(len(q))
            if len(q):
                qmin, qmax = min(qmin, int(q.min())), max(qmax, int(q.max()))
                Q.append(-10 * np.log10(LUT[q].mean()))
            else:
                Q.append(float("nan"))
    L = np.array(L); Q = np.array(Q)
    if len(L) != b["reads"] or int(L.sum()) != b["bases"]:
        sys.exit(f"block {b['start']}: {len(L)} reads / {int(L.sum())} bases differs from manifest ({b['reads']} / {b['bases']})")
    rows.append((bi, b["pool"], L, Q))
    blocks.append({"block": bi, "stratum": b["stratum"], "pool": b["pool"], "start": b["start"], "reads": len(L),
                   "bases": int(L.sum()), "mean_len": round(float(L.mean()), 1), "median_len": int(np.median(L)),
                   "n50": n50(L), "median_q": round(float(np.nanmedian(Q)), 2), "mean_q": round(float(np.nanmean(Q)), 2),
                   "pct_q10": round(100 * float(np.mean(Q >= 10)), 2), "pct_q15": round(100 * float(np.mean(Q >= 15)), 2),
                   "pct_q20": round(100 * float(np.mean(Q >= 20)), 2)})

with gzip.open(os.path.join(a.outdir, "read_stats.tsv.gz"), "wt") as o:
    o.write("block\tpool\tlength\tmean_q\n")
    for bi, pool, L, Q in rows:
        for l, q in zip(L, Q):
            o.write(f"{bi}\t{pool}\t{l}\t{q:.3f}\n")
cols = list(blocks[0])
with open(os.path.join(a.outdir, "block_stats.tsv"), "w") as o:
    o.write("\t".join(cols) + "\n")
    for r in blocks:
        o.write("\t".join(str(r[c]) for c in cols) + "\n")

P = {}
for pool in ("dev", "eval"):
    L = np.concatenate([r[2] for r in rows if r[1] == pool]); Q = np.concatenate([r[3] for r in rows if r[1] == pool])
    P[pool] = {"L": L, "Q": Q}
summ = {"manifest": os.path.basename(mp), "quality_char_range": [qmin, qmax],
        "phred_offset_33_plausible": bool(qmin >= 33), "pools": {}}
for pool, d in P.items():
    L, Q = d["L"], d["Q"]
    summ["pools"][pool] = {
        "reads": int(len(L)), "bases": int(L.sum()), "mean_len": round(float(L.mean()), 1),
        "len_q05_q25_q50_q75_q95": [int(x) for x in np.percentile(L, [5, 25, 50, 75, 95])], "n50": n50(L),
        "pct_reads_lt_500bp": round(100 * float(np.mean(L < 500)), 2),
        "mean_q": round(float(np.nanmean(Q)), 2), "q05_q25_q50_q75_q95": [round(float(x), 2) for x in np.nanpercentile(Q, [5, 25, 50, 75, 95])],
        "pct_q10": round(100 * float(np.mean(Q >= 10)), 2), "pct_q15": round(100 * float(np.mean(Q >= 15)), 2),
        "pct_q20": round(100 * float(np.mean(Q >= 20)), 2)}


def ks(x, y):
    x, y = np.sort(x[~np.isnan(x)]), np.sort(y[~np.isnan(y)])
    z = np.concatenate([x, y])
    return float(np.max(np.abs(np.searchsorted(x, z, side="right") / len(x) - np.searchsorted(y, z, side="right") / len(y))))


bd = [b for b in blocks if b["pool"] == "dev"]; be = [b for b in blocks if b["pool"] == "eval"]
cmp_ = {"ks_D_length_reads": round(ks(P["dev"]["L"], P["eval"]["L"]), 4), "ks_D_meanq_reads": round(ks(P["dev"]["Q"], P["eval"]["Q"]), 4),
        "note_ks": "effect size only; reads within a block are correlated so read-level p-values are not valid"}
try:
    from scipy.stats import mannwhitneyu
    for key in ("mean_len", "mean_q", "pct_q20"):
        r = mannwhitneyu([b[key] for b in bd], [b[key] for b in be], alternative="two-sided")
        cmp_[f"block_level_mannwhitney_p_{key}"] = round(float(r.pvalue), 3)
except ImportError:
    cmp_["block_level_tests"] = "scipy not available"
# block means by position along the run (drift check): Spearman of block mean vs read-number start
try:
    from scipy.stats import spearmanr
    for key in ("mean_len", "mean_q"):
        rho, p = spearmanr([b["start"] for b in blocks], [b[key] for b in blocks])
        cmp_[f"drift_spearman_{key}_vs_position"] = {"rho": round(float(rho), 3), "p": round(float(p), 3)}
except ImportError:
    pass
summ["comparison"] = cmp_
json.dump(summ, open(os.path.join(a.outdir, "pool_summary.json"), "w"), indent=2)

out = [f"Pool QC for {run} (manifest {os.path.basename(mp)})", ""]
out.append(f"Quality characters {qmin}..{qmax} (Phred+33 plausible: {summ['phred_offset_33_plausible']})")
out.append("")
out.append(f"{'':28s}{'dev':>14s}{'eval':>14s}")
def line(label, k, fmt="{}"):
    out.append(f"{label:28s}{fmt.format(summ['pools']['dev'][k]):>14s}{fmt.format(summ['pools']['eval'][k]):>14s}")
line("reads", "reads", "{:,}"); line("bases", "bases", "{:,}"); line("mean length (bp)", "mean_len", "{:,.1f}"); line("N50 (bp)", "n50", "{:,}")
line("% reads < 500 bp", "pct_reads_lt_500bp", "{:.2f}"); line("mean read Q", "mean_q", "{:.2f}")
line("% reads Q>=10", "pct_q10", "{:.2f}"); line("% reads Q>=15", "pct_q15", "{:.2f}"); line("% reads Q>=20", "pct_q20", "{:.2f}")
for k, lab in (("len_q05_q25_q50_q75_q95", "length 5/25/50/75/95 pct"), ("q05_q25_q50_q75_q95", "read Q 5/25/50/75/95 pct")):
    out.append(f"{lab:28s}{str(summ['pools']['dev'][k]):>30s}")
    out.append(f"{'':28s}{str(summ['pools']['eval'][k]):>30s}")
out.append("")
for k, v in cmp_.items():
    out.append(f"{k}: {v}")
open(os.path.join(a.outdir, "pool_qc_report.txt"), "w").write("\n".join(out) + "\n")
print("\n".join(out))
print(f"\nwrote {a.outdir}/ (read_stats.tsv.gz, block_stats.tsv, pool_summary.json, pool_qc_report.txt)")
