#!/usr/bin/env python3
"""22_cepa_check.py - EXPLORATORY, DEV pool only; part of no verdict; changes no registered rule.

Open question from bin/21 (dev): reads for B. fragilis cepA were 2 assigned and 27 AMBIGUOUS, the ambiguity being among the ResFinder families
cepA, cepA-29, cepA-44 and cepA-49 (reports/amr/dev_family_reads.tsv). Cause not established.

Hypothesis written BEFORE the first run: the four 'families' are the same gene under different ResFinder names (near-identical sequences), so reads
cannot separate them and the ambiguity is a naming artefact of the A4 family rule (allele number stripped, a '-NN' suffix kept).
  (P1) nucleotide identity between family representatives (first allele of each family in all.fsa; global alignment) is >= 90% for all six pairs.
  (P2) for the reads concerned, best and runner-up family differ by less than 10% in matching bases (median over reads).
  (P3) the group 'cepA' + 'cepA-NN' contains >= 3 reads once ambiguity among them is resolved, i.e. the number of distinct reads
       passing the A4 alignment filter on any cepA-group allele is >= 3.
Refuted if any pair is below 80% identity (genuinely different genes) or the median margin is >= 10%.
Reported only; if supported, the family rule could be changed ONLY by a dated addendum (A5) written before the evaluation pool is touched and
disclosed as made after seeing dev read results. Standard library only (a plain global alignment, no Biopython).
"""
import argparse, itertools, os, re, statistics
from collections import defaultdict

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
P = lambda *x: os.path.join(ROOT, *x)
HDR = re.compile(r"^(.+?)_(\d+)_(.+)$")
GROUP = re.compile(r"^cepA(-\d+)?$", re.I)


def fam(h):
    m = HDR.match(h)
    return m.group(1) if m else h


def read_fsa(path, pred):
    seqs, name = {}, None
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith(">"):
            name = line[1:].split()[0]
            if not pred(fam(name)):
                name = None
            else:
                seqs[name] = []
        elif name:
            seqs[name].append(line.strip().upper())
    return {k: "".join(v) for k, v in seqs.items()}


def global_identity(a, b):
    """Needleman-Wunsch, match +1, mismatch -1, gap -2 (linear); identity = matches / alignment length."""
    n, m = len(a), len(b)
    prev = [j * -2 for j in range(m + 1)]
    tb = [bytearray(m + 1) for _ in range(n + 1)]
    for j in range(1, m + 1):
        tb[0][j] = 2
    for i in range(1, n + 1):
        cur = [i * -2] + [0] * m; tb[i][0] = 1
        for j in range(1, m + 1):
            d = prev[j - 1] + (1 if a[i - 1] == b[j - 1] else -1); u = prev[j] - 2; l = cur[j - 1] - 2
            if d >= u and d >= l:
                cur[j] = d; tb[i][j] = 0
            elif u >= l:
                cur[j] = u; tb[i][j] = 1
            else:
                cur[j] = l; tb[i][j] = 2
        prev = cur
    i, j, mt, ln = n, m, 0, 0
    while i > 0 or j > 0:
        t = tb[i][j]
        if t == 0:
            mt += a[i - 1] == b[j - 1]; i -= 1; j -= 1
        elif t == 1:
            i -= 1
        else:
            j -= 1
        ln += 1
    return mt / ln


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fsa", default="data/db/resfinder_db/all.fsa")
    ap.add_argument("--paf", default="data/amr/dev/dev_vs_resfinder.paf")
    ap.add_argument("--out", default="reports/amr")
    a = ap.parse_args()
    seqs = read_fsa(P(a.fsa), lambda f: bool(GROUP.match(f)))
    byfam = defaultdict(list)
    for n in seqs:
        byfam[fam(n)].append(n)
    L = ["cepA group in ResFinder (EXPLORATORY, dev pool; hypothesis and predictions are in the script header)", ""]
    for f, ns in sorted(byfam.items()):
        L.append(f"family {f}: {len(ns)} allele(s), lengths {sorted(len(seqs[n]) for n in ns)}")
    reps = {f: sorted(ns)[0] for f, ns in byfam.items()}
    L += ["", "Representatives (first allele by name): " + ", ".join(f"{f}={n}" for f, n in sorted(reps.items())), "", "Pairwise global nucleotide identity between representatives:"]
    ids = {}
    for f1, f2 in itertools.combinations(sorted(reps), 2):
        ids[(f1, f2)] = global_identity(seqs[reps[f1]], seqs[reps[f2]])
        L.append(f"  {f1:10s} vs {f2:10s} {ids[(f1, f2)]*100:6.2f}%")
    # reads: alignments to any cepA-group allele passing the A4 filter (identity >= 0.80, aligned span >= 0.80 of the allele)
    per = defaultdict(dict)
    for line in open(P(a.paf)):
        f = line.split("\t", 12)
        t = f[5]
        if not GROUP.match(fam(t)):
            continue
        mt, bl, tl, ts, te = int(f[9]), int(f[10]), int(f[6]), int(f[7]), int(f[8])
        if bl <= 0 or mt / bl < 0.80 or (te - ts) / tl < 0.80:
            continue
        d = per[f[0]]; ff = fam(t)
        if ff not in d or mt > d[ff][0]:
            d[ff] = (mt, mt / bl)
    margins, best_id = [], []
    for q, d in per.items():
        r = sorted(d.values(), key=lambda x: -x[0])
        best_id.append(r[0][1])
        if len(r) > 1:
            margins.append(1 - r[1][0] / r[0][0])
    L += ["", f"Reads passing the A4 alignment filter on any cepA-group allele: {len(per)}",
          f"  reads hitting >= 2 families: {len(margins)}; median margin best vs runner-up (1 - second/best matching bases): "
          f"{(statistics.median(margins)*100 if margins else float('nan')):.1f}%",
          f"  median identity of the best alignment: {(statistics.median(best_id)*100 if best_id else float('nan')):.1f}%"]
    p1 = bool(ids) and min(ids.values()) >= 0.90
    p2 = bool(margins) and statistics.median(margins) < 0.10
    p3 = len(per) >= 3
    ref = (bool(ids) and min(ids.values()) < 0.80) or (bool(margins) and statistics.median(margins) >= 0.10)
    L += ["", f"P1 (all pairs >= 90%): {'met' if p1 else 'not met'}   P2 (median margin < 10%): {'met' if p2 else 'not met'}   P3 (>= 3 distinct reads in the group): {'met' if p3 else 'not met'}",
          f"Refutation criterion (a pair < 80% or median margin >= 10%): {'REFUTED' if ref else 'not triggered'}",
          "Exploratory; no registered rule is changed here. Any change of the family rule needs a dated addendum (A5) before the evaluation pool is touched."]
    txt = "\n".join(L)
    os.makedirs(P(a.out), exist_ok=True); open(P(a.out, "dev_cepa_check.txt"), "w").write(txt + "\n"); print(txt)


if __name__ == "__main__":
    main()
