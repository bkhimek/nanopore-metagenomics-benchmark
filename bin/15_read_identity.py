#!/usr/bin/env python3
"""15_read_identity.py - EXPLORATORY, DEVELOPMENT pool only; not a registered arm, part of no verdict.

Question left open by bin/14: for Limosilactobacillus fermentum, Sylph's effective coverage (14.2) is 0.26 of the coverage implied by the minimap2 trace
(54.6); Faecalibacterium prausnitzii 0.42, Roseburia 0.61, Fusobacterium 0.64; most other species 0.86-1.4.

Hypothesis written BEFORE the first run: Sylph counts k-mers (k = 31) that appear intact in the reads. A read k-mer survives only if it has no error, so
Eff_cov is proportional to coverage x (probability that a 31-mer is error-free), roughly identity^31. If some species' reads carry more errors than others
(for example because of base modifications that the basecaller handles badly), their Eff_cov falls below the alignment-based coverage while their mapped
reads stay normal. Predictions:
  (P1) mean read identity (matches / alignment length against the species' own reference) is lower for the species with low Eff_cov / mapped coverage
       (L. fermentum, F. prausnitzii, Roseburia, Fusobacterium) than for those near 1 (B. fragilis, Prevotella, Akkermansia, C. difficile, E. coli).
  (P2) the species ranking by mean 31-mer survival (identity^31, averaged over reads) agrees with the ranking by Eff_cov / mapped coverage
       (Spearman rho reported; the hypothesis is supported by rho >= 0.6 and a lower identity for L. fermentum, refuted by rho < 0.3 or by identity
       differing between the two groups by under 1 percentage point).
  (P3) the size of the effect is checked, not only the sign: the survival ratio L. fermentum / median over species is compared with the observed
       Eff_cov / mapped ratio of L. fermentum / median (0.26 / about 1.0); a difference of identity of 1.5 points changes survival by about 30%, so a
       fourfold deficit needs an identity gap of about 4 points. If the gap is much smaller, errors alone cannot explain it.
Reads are assigned to species exactly as in bin/08 (largest summed matching bases; >= 0.5 x read length; ambiguous if runner-up >= 0.9 x best),
by reusing bin/08's functions. Identity of a read is (sum of PAF matches) / (sum of PAF alignment-block lengths) over its alignments to the
assigned species. Standard library only. Needs the PAF written by bin/08 (data/trace/dev/dev_vs_refs.paf).
"""
import argparse, csv, importlib.util, os, statistics, sys
from collections import defaultdict

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
P = lambda *x: os.path.join(ROOT, *x)


def load08():
    spec = importlib.util.spec_from_file_location("x08", P("bin", "08_trace_reads.py"))
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
    return m


def ranks(v):
    order = sorted(range(len(v)), key=lambda i: v[i]); r = [0.0] * len(v); i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
            j += 1
        for k in range(i, j + 1):
            r[order[k]] = (i + j) / 2 + 1
        i = j + 1
    return r


def spearman(x, y):
    rx, ry = ranks(x), ranks(y); mx, my = statistics.mean(rx), statistics.mean(ry)
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    den = (sum((a - mx) ** 2 for a in rx) * sum((b - my) ** 2 for b in ry)) ** 0.5
    return num / den if den else float("nan")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--paf", default="data/trace/dev/dev_vs_refs.paf")
    ap.add_argument("--diagnostic", default="reports/sylph/dev_diagnostic/diagnostic_by_species.tsv")
    ap.add_argument("--k", type=int, default=31)
    ap.add_argument("--min-cov", type=float, default=0.5)
    ap.add_argument("--amb", type=float, default=0.9)
    ap.add_argument("--out", default="reports/sylph/dev_diagnostic")
    a = ap.parse_args()
    x08 = load08()

    # per read: per species (query length, matches, block)
    per = {}
    with open(P(a.paf)) as fh:
        for line in fh:
            f = line.split("\t", 12)
            sp = x08.stem_species(f[5].split("__")[0])
            d = per.get(f[0])
            if d is None:
                d = per[f[0]] = [int(f[1]), {}]
            m = d[1].setdefault(sp, [0, 0])
            m[0] += int(f[9]); m[1] += int(f[10])
    # same assignment as bin/08: it needs {species: matching bases}
    ident = defaultdict(list); rlen = defaultdict(list); outcome = defaultdict(int)
    for rid, (qlen, sp_d) in per.items():
        agg = [qlen, {s: v[0] for s, v in sp_d.items()}]
        lab = x08.assign(agg, a.min_cov, a.amb)
        if lab in ("weak", "ambiguous"):
            outcome[lab] += 1; continue
        mt, bl = sp_d[lab]
        if bl <= 0:
            continue
        ident[lab].append(mt / bl); rlen[lab].append(qlen); outcome["assigned"] += 1

    diag = {r["mock_species"]: r for r in csv.DictReader(open(P(a.diagnostic)), delimiter="\t")}
    rows = []
    for sp, v in ident.items():
        if sp not in diag:
            continue
        v_sorted = sorted(v)
        surv = statistics.mean(i ** a.k for i in v)
        rows.append({"species": sp, "reads": len(v), "mean_read_len": statistics.mean(rlen[sp]), "ident_mean": statistics.mean(v) * 100,
                     "ident_median": statistics.median(v) * 100, "ident_p10": v_sorted[len(v) // 10] * 100,
                     "kmer_survival": surv, "eff_over_mapped": float(diag[sp]["eff_over_mapped"])})
    if not rows:
        sys.exit("no overlap between PAF species and the diagnostic table")
    ms = statistics.median(r["kmer_survival"] for r in rows)
    mo = statistics.median(r["eff_over_mapped"] for r in rows)
    for r in rows:
        r["survival_rel"] = r["kmer_survival"] / ms; r["observed_rel"] = r["eff_over_mapped"] / mo
    rows.sort(key=lambda r: r["eff_over_mapped"])
    rho = spearman([r["kmer_survival"] for r in rows], [r["eff_over_mapped"] for r in rows]) if len(rows) > 2 else float("nan")
    rho_i = spearman([r["ident_mean"] for r in rows], [r["eff_over_mapped"] for r in rows]) if len(rows) > 2 else float("nan")

    cols = ["species", "reads", "mean_read_len", "ident_mean", "ident_median", "ident_p10", "kmer_survival", "eff_over_mapped", "survival_rel", "observed_rel"]
    os.makedirs(P(a.out), exist_ok=True)
    with open(P(a.out, "read_identity_by_species.tsv"), "w") as o:
        o.write("\t".join(cols) + "\n")
        for r in rows:
            o.write("\t".join(f"{r[c]:.4f}" if isinstance(r[c], float) else str(r[c]) for c in cols) + "\n")
    L = [f"Read identity vs Sylph coverage deficit (EXPLORATORY, dev pool; k = {a.k}; reads assigned {outcome['assigned']}, weak {outcome['weak']}, ambiguous {outcome['ambiguous']})", "",
         f"{'species':32s}{'reads':>8s}{'len':>7s}{'id mean':>9s}{'id med':>8s}{'id p10':>8s}{'surv':>7s}{'surv/med':>9s}{'Eff/map':>8s}{'obs/med':>8s}"]
    for r in rows:
        L.append(f"{r['species'][:31]:32s}{r['reads']:>8d}{r['mean_read_len']:>7.0f}{r['ident_mean']:>9.2f}{r['ident_median']:>8.2f}{r['ident_p10']:>8.2f}"
                 f"{r['kmer_survival']:>7.3f}{r['survival_rel']:>9.2f}{r['eff_over_mapped']:>8.2f}{r['observed_rel']:>8.2f}")
    L += ["", f"Spearman rho, k-mer survival vs Eff_cov/mapped coverage: {rho:.2f};  mean identity vs the same: {rho_i:.2f}  (n = {len(rows)}; small n, descriptive only)",
          "Identity = matches / alignment length of the minimap2 trace (a lower bound on read accuracy for reads that align poorly; reads are assigned by their best species)."]
    txt = "\n".join(L)
    open(P(a.out, "read_identity_report.txt"), "w").write(txt + "\n")
    print(txt)


if __name__ == "__main__":
    main()
