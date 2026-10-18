#!/usr/bin/env python3
"""16_base_basis.py - EXPLORATORY, DEVELOPMENT pool only; post hoc (written after bin/15 results were seen); part of no verdict.

Why: bin/15 refuted its registered hypothesis (read errors do not explain Sylph's low coverage for L. fermentum), but its table showed that mean read length
differs a lot by species (L. fermentum 1.9 kb, F. prausnitzii 3.1 kb, Roseburia 4.3 kb, most others 6-8 kb) and that Sylph's deficit tracks it. bin/14 had assumed the
pool mean read length (5,491 bp) for every species. A manufacturer's composition is a share of DNA MASS (bases), so the like-for-like quantity is bases per species, not reads.

This script recomputes, per species, with the species' OWN mean read length (from bin/15's table, reads assigned as in bin/08):
  bases = reads x mean length;  coverage = bases / genome size;  expected coverage from DNA% (as bin/14);
  base-basis obs/exp (share of assigned bases among the listed species / expected DNA share among the same species);
  Sylph Eff_cov / coverage and True_cov / coverage (from bin/14's table); and Bray-Curtis of the base-basis shares versus the expected shares.
Read-basis obs/exp is shown beside it. Sylph's coverage is NOT altered. Standard library only.
"""
import argparse, csv, os, statistics
from collections import defaultdict
ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
P = lambda *x: os.path.join(ROOT, *x)

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--identity", default="reports/sylph/dev_diagnostic/read_identity_by_species.tsv")
    ap.add_argument("--diagnostic", default="reports/sylph/dev_diagnostic/diagnostic_by_species.tsv")
    ap.add_argument("--composition", default="reference/d6331_composition.csv")
    ap.add_argument("--reads", type=int, default=185000)
    ap.add_argument("--pool-mean-len", type=float, default=5491.0)
    ap.add_argument("--out", default="reports/sylph/dev_diagnostic")
    a = ap.parse_args()
    dna, sz = defaultdict(float), defaultdict(list)
    for r in csv.DictReader(open(P(a.composition))):
        n = r["species"].replace("Lactobacillus fermentum", "Limosilactobacillus fermentum")
        dna[n] += float(r["genomic_dna_pct"]); sz[n].append(float(r["genome_size_mb_table2"]))
    diag = {r["mock_species"]: r for r in csv.DictReader(open(P(a.diagnostic)), delimiter="\t")}
    rows = []
    for r in csv.DictReader(open(P(a.identity)), delimiter="\t"):
        n = r["species"]; g = sum(sz[n]) / len(sz[n]) * 1e6
        reads, ln = int(r["reads"]), float(r["mean_read_len"]); b = reads * ln
        d = diag[n]
        rows.append({"species": n, "reads": reads, "len": ln, "bases": b, "cov": b / g, "exp_cov": dna[n] / 100 * a.reads * a.pool_mean_len / g,
                     "dna": dna[n], "eff": float(d["eff_cov_max"]), "true": float(d["true_cov_max"])})
    tb = sum(r["bases"] for r in rows); tr = sum(r["reads"] for r in rows); td = sum(r["dna"] for r in rows)
    for r in rows:
        r["base_share"] = r["bases"] / tb; r["read_share"] = r["reads"] / tr; r["exp_share"] = r["dna"] / td
        r["obs_exp_bases"] = r["base_share"] / r["exp_share"]; r["obs_exp_reads"] = r["read_share"] / r["exp_share"]
        r["eff_over_cov"] = r["eff"] / r["cov"]; r["true_over_cov"] = r["true"] / r["cov"]
    bc_b = 1 - sum(min(r["base_share"], r["exp_share"]) for r in rows); bc_r = 1 - sum(min(r["read_share"], r["exp_share"]) for r in rows)
    rows.sort(key=lambda r: r["len"])
    cols = ["species", "reads", "len", "cov", "exp_cov", "obs_exp_reads", "obs_exp_bases", "eff_over_cov", "true_over_cov"]
    os.makedirs(P(a.out), exist_ok=True)
    with open(P(a.out, "base_basis_by_species.tsv"), "w") as o:
        o.write("\t".join(cols) + "\n")
        for r in rows: o.write("\t".join(f"{r[c]:.4f}" if isinstance(r[c], float) else str(r[c]) for c in cols) + "\n")
    L = [f"Base-basis recomputation (EXPLORATORY, dev pool, post hoc; {len(rows)} species)", "",
         f"{'species':32s}{'reads':>7s}{'len':>6s}{'cov':>6s}{'exp.cov':>8s}{'obs/exp reads':>14s}{'obs/exp bases':>14s}{'Eff/cov':>9s}{'True/cov':>9s}"]
    for r in rows:
        L.append(f"{r['species'][:31]:32s}{r['reads']:>7d}{r['len']:>6.0f}{r['cov']:>6.1f}{r['exp_cov']:>8.1f}{r['obs_exp_reads']:>14.2f}{r['obs_exp_bases']:>14.2f}{r['eff_over_cov']:>9.2f}{r['true_over_cov']:>9.2f}")
    t = [r["true_over_cov"] for r in rows]; e = [r["eff_over_cov"] for r in rows]
    L += ["", f"Bray-Curtis of shares versus expected DNA shares among these species: reads {bc_r:.3f}, bases {bc_b:.3f}",
          f"True_cov / coverage: median {statistics.median(t):.2f}, range {min(t):.2f}-{max(t):.2f};  Eff_cov / coverage: median {statistics.median(e):.2f}, range {min(e):.2f}-{max(e):.2f}",
          "Coverage = assigned reads x species mean read length / Table-2 genome size (E. coli: mean of five strains). Read length is the full read, not the aligned part."]
    s = "\n".join(L); open(P(a.out, "base_basis_report.txt"), "w").write(s + "\n"); print(s)

if __name__ == "__main__":
    main()
