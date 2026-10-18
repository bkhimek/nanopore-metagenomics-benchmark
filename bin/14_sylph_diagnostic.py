#!/usr/bin/env python3
"""14_sylph_diagnostic.py - EXPLORATORY diagnostic of Sylph abundance errors on the DEVELOPMENT pool (not registered, not part of any verdict).

Questions: why is Escherichia coli over-counted (share-basis obs/exp 2.34) and Limosilactobacillus fermentum under-counted (0.50) in arm 3a?

For every mock species it joins
  - the Sylph genome rows (data/sylph/dev/arm_3a.tsv; optionally arm_3c.tsv for True_cov), mapped to mock species through the crosswalk v2
    and the candidate list (accession -> database species);
  - the coverage expected from the manufacturer's DNA share (sum over strains of DNA% x reads x mean length / genome size);
  - the coverage implied by the read-to-genome alignment of the trace (reports/trace/trace_by_species.tsv: mapped reads x mean length / mean genome size);
and reports, per species: number of genome rows, summed and maximal effective coverage, True_cov (if 3c), median versus mean coverage, k-mers
reassigned, adjusted ANI, and the ratios of Sylph coverage to both references.
It also computes a post-hoc WHAT-IF: share-basis obs/exp and Bray-Curtis (among called mock species, the '09 convention') if each mock species were
represented by its best-covered genome row only (instead of the sum of all its rows, the registered rule). The what-if does NOT replace the
registered verdict; the evaluation pool is not used.

Sequence_abundance is proportional to Eff_cov x genome length and Taxonomic_abundance to Eff_cov (checked in the output by the ratio columns).
Standard library only.
"""
import argparse, csv, json, os, re, sys
from collections import defaultdict

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
P = lambda *x: os.path.join(ROOT, *x)
ACC = re.compile(r"GC[AF]_\d{9}\.\d+")


def fnum(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return float("nan")


def read_tsv(path):
    with open(P(path)) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def bray_curtis(obs, exp):
    return 1 - sum(min(obs.get(k, 0.0), exp.get(k, 0.0)) for k in set(obs) | set(exp))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--profile", default="data/sylph/dev/arm_3a.tsv")
    ap.add_argument("--profile-u", default="data/sylph/dev/arm_3c.tsv", help="arm 3c profile (True_cov); skipped if missing")
    ap.add_argument("--crosswalk", default="reference/d6331_gtdb_crosswalk_v2.csv")
    ap.add_argument("--candidates", default="reports/sylph/crosswalk_v2_candidates.tsv")
    ap.add_argument("--composition", default="reference/d6331_composition.csv")
    ap.add_argument("--trace", default="reports/trace/trace_by_species.tsv")
    ap.add_argument("--reads", type=int, default=185000)
    ap.add_argument("--mean-len", type=float, default=5491.0, help="mean read length of the pool (development: 5,491 bp)")
    ap.add_argument("--evaluable-min", type=float, default=20.0)
    ap.add_argument("--out", default="reports/sylph/dev_diagnostic")
    a = ap.parse_args()

    # accession -> database species; database species -> mock taxid
    acc2sp = {r["accession"]: r["database_species"] for r in read_tsv(a.candidates)}
    sp2mock, mock_name = {}, {}
    for r in csv.DictReader(open(P(a.crosswalk))):
        mock_name[int(r["mock_taxid"])] = r["mock_species"]
        if r["status"] == "ok":
            sp2mock[r["database_species"]] = int(r["mock_taxid"])

    # composition: expected DNA share and coverage per mock species
    dna, expcov, sizes = defaultdict(float), defaultdict(float), defaultdict(list)
    for r in csv.DictReader(open(P(a.composition))):
        t = int(r["species_taxid"]); d = fnum(r["genomic_dna_pct"]) or 0.0; g = fnum(r["genome_size_mb_table2"])
        dna[t] += d
        if g == g and g > 0:
            sizes[t].append(g)
            expcov[t] += d / 100 * a.reads * a.mean_len / (g * 1e6)
    mapped = {int(r["taxid"]): int(r["mapped_reads"]) for r in read_tsv(a.trace)}
    mapcov = {t: mapped.get(t, 0) * a.mean_len / (sum(s) / len(s) * 1e6) for t, s in sizes.items()}

    # profile rows
    rows = read_tsv(a.profile)
    ucol = {}
    if a.profile_u and os.path.exists(P(a.profile_u)):
        for r in read_tsv(a.profile_u):
            m = ACC.search(r["Genome_file"])
            ucol[m.group(0)] = fnum(r.get("True_cov"))
    gen = []
    for r in rows:
        m = ACC.search(r["Genome_file"])
        if not m:
            sys.exit(f"no accession in {r['Genome_file']}")
        acc = m.group(0)
        sp = acc2sp.get(acc)
        if sp is None:
            sys.exit(f"{acc} is not in the candidate list {a.candidates}; run the gap rule (bin/13) first")
        tid = sp2mock.get(sp)
        gen.append({"acc": acc, "db_species": sp, "mock": tid, "seq": fnum(r["Sequence_abundance"]), "tax": fnum(r["Taxonomic_abundance"]),
                    "ani": fnum(r["Adjusted_ANI"]), "eff": fnum(r.get("Eff_cov")), "true": ucol.get(acc, float("nan")),
                    "med": fnum(r["Median_cov"]), "mean": fnum(r["Mean_cov_geq1"]), "reass": fnum(r.get("kmers_reassigned"))})
    unmapped = [g for g in gen if g["mock"] is None]
    if unmapped:
        print("NOTE: genome rows outside the crosswalk (non-mock):", [(g["acc"], g["db_species"]) for g in unmapped])

    by = defaultdict(list)
    for g in gen:
        if g["mock"] is not None:
            by[g["mock"]].append(g)
    called = sorted(by)
    exp_share = {t: dna[t] for t in called}
    es = sum(exp_share.values()); exp_share = {t: v / es for t, v in exp_share.items()}

    def shares(rule):
        o = {}
        for t in called:
            gs = by[t]
            o[t] = sum(g["seq"] for g in gs) if rule == "sum" else max(gs, key=lambda g: g["eff"])["seq"]
        s = sum(o.values())
        return {t: v / s for t, v in o.items()}

    sh_sum, sh_best = shares("sum"), shares("best")
    bc_sum, bc_best = bray_curtis(sh_sum, exp_share), bray_curtis(sh_best, exp_share)

    os.makedirs(P(a.out), exist_ok=True)
    cols = ["mock_species", "taxid", "dna_pct", "genome_rows", "eff_cov_max", "eff_cov_sum", "true_cov_max", "median_cov_best", "mean_cov_best",
            "kmers_reassigned_sum", "adjusted_ani_best", "expected_cov", "mapped_cov", "eff_over_expected", "eff_over_mapped", "true_over_mapped",
            "obs_exp_sum_rule", "obs_exp_best_genome_rule"]
    out = []
    for t in sorted(called, key=lambda t: -dna[t]):
        gs = by[t]; b = max(gs, key=lambda g: g["eff"])
        tc = max((g["true"] for g in gs if g["true"] == g["true"]), default=float("nan"))
        out.append({"mock_species": mock_name[t], "taxid": t, "dna_pct": f"{dna[t]:g}", "genome_rows": len(gs), "eff_cov_max": f"{b['eff']:.2f}",
                    "eff_cov_sum": f"{sum(g['eff'] for g in gs):.2f}", "true_cov_max": f"{tc:.2f}", "median_cov_best": f"{b['med']:g}",
                    "mean_cov_best": f"{b['mean']:.2f}", "kmers_reassigned_sum": f"{sum(g['reass'] for g in gs):g}", "adjusted_ani_best": f"{b['ani']:.2f}",
                    "expected_cov": f"{expcov[t]:.2f}", "mapped_cov": f"{mapcov.get(t, float('nan')):.2f}",
                    "eff_over_expected": f"{b['eff'] / expcov[t]:.2f}" if expcov[t] else "nan",
                    "eff_over_mapped": f"{b['eff'] / mapcov[t]:.2f}" if mapcov.get(t) else "nan",
                    "true_over_mapped": f"{tc / mapcov[t]:.2f}" if mapcov.get(t) and tc == tc else "nan",
                    "obs_exp_sum_rule": f"{sh_sum[t] / exp_share[t]:.3f}", "obs_exp_best_genome_rule": f"{sh_best[t] / exp_share[t]:.3f}"})
    with open(P(a.out, "diagnostic_by_species.tsv"), "w") as o:
        w = csv.DictWriter(o, fieldnames=cols, delimiter="\t"); w.writeheader(); w.writerows(out)
    with open(P(a.out, "diagnostic_genomes.tsv"), "w") as o:
        o.write("accession\tdatabase_species\tmock_species\tsequence_abundance\ttaxonomic_abundance\teff_cov\ttrue_cov\tmedian_cov\tmean_cov_geq1\tadjusted_ani\tkmers_reassigned\n")
        for g in sorted(gen, key=lambda g: -g["seq"]):
            o.write(f"{g['acc']}\t{g['db_species']}\t{mock_name.get(g['mock'], 'NON-MOCK')}\t{g['seq']}\t{g['tax']}\t{g['eff']}\t{g['true']}\t{g['med']}\t{g['mean']}\t{g['ani']}\t{g['reass']:g}\n")

    # checks on how the abundances are built
    chk = []
    for g in gen:
        t = g["mock"]
        if t is not None and sizes.get(t) and len(by[t]) == 1 and g["eff"] > 0:
            chk.append((g["seq"] / (g["eff"] * sum(sizes[t]) / len(sizes[t])), g["tax"] / g["eff"]))
    L = [f"Sylph diagnostic (EXPLORATORY, development pool; profile {a.profile}; {a.reads} reads, mean length {a.mean_len:g})", "",
         f"Registered rule (sum of genome rows per mock species): Bray-Curtis among called mock species {bc_sum:.3f}",
         f"What-if (best-covered genome row only):                Bray-Curtis among called mock species {bc_best:.3f}", ""]
    if chk:
        r1 = [c[0] for c in chk]; r2 = [c[1] for c in chk]
        L.append(f"Single-row species: Sequence_abundance / (Eff_cov x genome size) ranges {min(r1):.3f} to {max(r1):.3f}; Taxonomic_abundance / Eff_cov ranges {min(r2):.4f} to {max(r2):.4f}")
        L.append("  (near-constant ratios mean: Sequence_abundance ~ Eff_cov x genome length, Taxonomic_abundance ~ Eff_cov)")
        L.append("")
    hdr = f"{'species':32s}{'rows':>5s}{'effmax':>8s}{'effsum':>8s}{'trueC':>7s}{'exp.cov':>8s}{'map.cov':>8s}{'eff/map':>8s}{'true/map':>9s}{'obs/exp sum':>12s}{'best-genome':>12s}"
    L.append(hdr)
    for r in out:
        L.append(f"{r['mock_species'][:31]:32s}{r['genome_rows']:>5d}{r['eff_cov_max']:>8s}{r['eff_cov_sum']:>8s}{r['true_cov_max']:>7s}{r['expected_cov']:>8s}{r['mapped_cov']:>8s}"
                 f"{r['eff_over_mapped']:>8s}{r['true_over_mapped']:>9s}{r['obs_exp_sum_rule']:>12s}{r['obs_exp_best_genome_rule']:>12s}")
    L += ["", "Notes: expected and mapped coverage assume the pool mean read length for every species and Table-2 genome sizes (mapped coverage uses the mean size",
          "for E. coli's five strains); both are approximate. The best-genome what-if is post hoc and exploratory; it does not change the registered verdict."]
    txt = "\n".join(L)
    open(P(a.out, "diagnostic_report.txt"), "w").write(txt + "\n")
    print(txt)


if __name__ == "__main__":
    main()
