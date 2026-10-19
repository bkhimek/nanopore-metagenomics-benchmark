#!/usr/bin/env python3
"""27_attribution_genefrac.py - EXPLORATORY, DEV POOL ONLY, follow-up to Addendum A8. Part of no verdict, changes no registered rule.
The script has no evaluation-pool option; it reads dev files only.

Question: why do tet(Q) reads often get only a family/order from Kraken2 (A8: B. fragilis tet(Q) 12/30 to species on dev, cepA 34/34)?
Mobile-element idea: the gene sits in sequence shared by several Bacteroidales, so a read that is mostly gene has few species-specific k-mers.
Measurable from existing files: read length and the gene fraction = aligned query span of the read on its best ResFinder family alignment / read length
(A7 PAF, same filter as bin/21: identity >= 0.80 base-level, span >= 0.80 of the gene, ambiguity rule as A7). Destination = bin/26 bucket
(species = correct; up_genus / up_family / up_higher = right lineage, not to species; wrong = sister / same_family / elsewhere). Reads with no truth label are left out.
Hypotheses and predictions, written BEFORE the first run (tet(Q) reads of both carriers pooled unless stated; 'above' = up_genus, up_family or up_higher):
  H8  the median gene fraction of 'above' reads is higher than that of 'species' reads by >= 0.05 AND AUC(above vs species) >= 0.60.
  H9  the median read length of 'above' reads is lower than that of 'species' reads by >= 500 bases.
  H10 gene fraction does NOT fully explain the tet(Q) loss: among reads with gene fraction < 0.30, the share reaching species is for tet(Q) at least 0.10 below that of cepA.
      (basis: cepA reads reached species in 34 of 34 at a similar span; the gene is shared on a mobile element, not only a larger share of the read)
  Descriptive only, no prediction: gene fraction and length of the wrong-organism reads; counts per bin. Small numbers: no tests beyond the AUC; the result is a pointer, not a proof.
A hypothesis that fails stays failed and is reported. Standard library only.
"""
import argparse, csv, importlib.util, os, statistics, sys
from collections import defaultdict

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
P = lambda *x: os.path.join(ROOT, *x)
ABOVE = ("up_genus", "up_family", "up_higher")
WRONG = ("sister", "same_family", "elsewhere")
LOW = 0.30


def auc(a, b):
    """P(a > b) + 0.5 P(a == b), a = above, b = species"""
    if not a or not b:
        return float("nan")
    s = sum((x > y) + 0.5 * (x == y) for x in a for y in b)
    return s / (len(a) * len(b))


def med(v):
    return statistics.median(v) if v else float("nan")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--reads-tsv", default="reports/amr/dev_attribution_reads_A8.tsv")
    ap.add_argument("--paf", default="data/amr/dev/dev_vs_resfinder_c.paf")
    ap.add_argument("--merge", default="reference/amr_family_merge.tsv")
    ap.add_argument("--out", default="reports/amr")
    ap.add_argument("--tag", default="dev_gene_fraction")
    a = ap.parse_args()
    spec = importlib.util.spec_from_file_location("x21", P("bin", "21_amr_reads.py"))
    x21 = importlib.util.module_from_spec(spec); spec.loader.exec_module(x21)
    for r in csv.DictReader(open(P(a.merge)), delimiter="\t"):
        x21.MERGE[r["alias"]] = r["family"]
    rows = list(csv.DictReader(open(P(a.reads_tsv)), delimiter="\t"))
    want = {r["read"]: r for r in rows}
    best = {}  # read -> family -> (matches, qlen, qspan)
    with open(P(a.paf)) as fh:
        for line in fh:
            f = line.split("\t", 12)
            q = f[0]
            if q not in want:
                continue
            qlen, qs, qe, tname, tlen, ts, te, mt, bl = int(f[1]), int(f[2]), int(f[3]), f[5], int(f[6]), int(f[7]), int(f[8]), int(f[9]), int(f[10])
            if bl <= 0 or mt / bl < x21.MIN_ID or (te - ts) / tlen < x21.MIN_SPAN:
                continue
            fam = x21.family(tname)
            d = best.setdefault(q, {})
            if fam not in d or mt > d[fam][0]:
                d[fam] = (mt, qlen, qe - qs)
    data = []
    for q, r in want.items():
        if r["destination"] == "no_truth" or q not in best or r["family"] not in best[q]:
            continue
        mt, qlen, span = best[q][r["family"]]
        data.append(dict(read=q, family=r["family"], true=r["true_species"], dest=r["destination"], qlen=qlen, span=span, frac=span / qlen))
    miss = len([1 for r in rows if r["destination"] != "no_truth" and r["family"] in ("cepA", "tet(Q)", "tet(W)")]) - len([d for d in data if d["family"] in ("cepA", "tet(Q)", "tet(W)")])
    L = ["EXPLORATORY follow-up to A8 (DEV only): gene fraction and read length against the Kraken2 destination", "",
         f"reads used {len(data)}; scored-family reads in the A8 table that could not be matched to an A7 alignment: {miss}", "",
         f"{'family':8s}{'class':10s}{'n':>4s}{'median read length':>20s}{'median gene fraction':>22s}"]
    grp = defaultdict(list)
    for d in data:
        c = "species" if d["dest"] == "species" else "above" if d["dest"] in ABOVE else "wrong" if d["dest"] in WRONG else "uncl"
        d["cls"] = c; grp[(d["family"], c)].append(d)
    for (f, c), v in sorted(grp.items()):
        L.append(f"{f:8s}{c:10s}{len(v):>4d}{med([x['qlen'] for x in v]):>20.0f}{med([x['frac'] for x in v]):>22.3f}")
    tq = [d for d in data if d["family"] == "tet(Q)"]
    sp = [d for d in tq if d["cls"] == "species"]; ab = [d for d in tq if d["cls"] == "above"]
    dfrac = med([d["frac"] for d in ab]) - med([d["frac"] for d in sp]); A = auc([d["frac"] for d in ab], [d["frac"] for d in sp])
    dlen = med([d["qlen"] for d in sp]) - med([d["qlen"] for d in ab])
    lowshare = lambda fam: (lambda v: (sum(1 for x in v if x["cls"] == "species") / len(v), len(v)) if v else (float("nan"), 0))([d for d in data if d["family"] == fam and d["frac"] < LOW])
    (ts, tn), (cs, cn) = lowshare("tet(Q)"), lowshare("cepA")
    L += ["", f"tet(Q) reads: species {len(sp)}, above species {len(ab)}",
          f"H8  median gene fraction above - species = {dfrac:+.3f}; AUC(above > species) = {A:.2f}: {'MET' if dfrac >= 0.05 and A >= 0.60 else 'NOT MET'}",
          f"H9  median read length species - above = {dlen:+.0f} bases: {'MET' if dlen >= 500 else 'NOT MET'}",
          f"H10 reads with gene fraction < {LOW}: tet(Q) species share {ts:.3f} (n {tn}), cepA {cs:.3f} (n {cn}); gap {cs - ts:+.3f}: {'MET' if cs - ts >= 0.10 else 'NOT MET'}",
          "", "Exploratory, dev pool only, small numbers; a pointer for a later registered test, not a result."]
    wr = [d for d in data if d["cls"] == "wrong"]
    if wr:
        L.append(f"descriptive: {len(wr)} wrong-organism reads, median read length {med([d['qlen'] for d in wr]):.0f}, median gene fraction {med([d['frac'] for d in wr]):.3f}")
    txt = "\n".join(L)
    os.makedirs(P(a.out), exist_ok=True)
    open(P(a.out, a.tag + ".txt"), "w").write(txt + "\n")
    with open(P(a.out, a.tag + ".tsv"), "w") as o:
        w = csv.writer(o, delimiter="\t"); w.writerow(["read", "family", "true_species", "destination", "read_len", "gene_span", "gene_fraction"])
        for d in data:
            w.writerow([d["read"], d["family"], d["true"], d["dest"], d["qlen"], d["span"], f"{d['frac']:.3f}"])
    print(txt)


if __name__ == "__main__":
    main()
