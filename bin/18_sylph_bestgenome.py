#!/usr/bin/env python3
"""18_sylph_bestgenome.py - Sylph abundance with the 'best genome per mock species' rule (Addendum A2, part B).

bin/14 found on the DEVELOPMENT pool that several database genomes can share the reads of one mock species (E. coli: three genomes), so the registered rule (sum of
Sequence_abundance over all genome rows of a mock species) counts the same reads more than once. A2 fixes the alternative rule BEFORE it is applied to the evaluation pool:
  best-genome rule = for each mock species keep only the genome row with the largest Eff_cov (ties: the largest Sequence_abundance), drop its other rows,
  then renormalise the shares among the called mock species.
This script applies BOTH rules to one profile and reports Bray-Curtis in two conventions (bin/04: called species against all 17 truth species; bin/09/14: among called mock
species) and share-basis observed/expected per species. It never changes the registered scores. For the evaluation pool it needs --allow-eval and the git tag
bases-registered, and its outputs are never overwritten. Standard library only.
"""
import argparse, csv, json, os, re, subprocess, sys
from collections import defaultdict

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
P = lambda *x: os.path.join(ROOT, *x)
ACC = re.compile(r"GC[AF]_\d{9}\.\d+")


def fnum(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return float("nan")


def tsv(path):
    return list(csv.DictReader(open(P(path)), delimiter="\t"))


def bc(o, e):
    return 1 - sum(min(o.get(k, 0.0), e.get(k, 0.0)) for k in set(o) | set(e))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pool", choices=["dev", "eval"], required=True)
    ap.add_argument("--allow-eval", action="store_true")
    ap.add_argument("--profile", help="default data/sylph/<pool>/arm_3a.tsv")
    ap.add_argument("--crosswalk", help="default reference/d6331_gtdb_crosswalk_v2.csv (dev) or ..._v2_eval.csv (eval)")
    ap.add_argument("--candidates", default="reports/sylph/crosswalk_v2_candidates.tsv")
    ap.add_argument("--composition", default="reference/d6331_composition.csv")
    ap.add_argument("--out", help="default reports/sylph/<pool>_bestgenome")
    a = ap.parse_args()
    if a.pool == "eval":
        if not a.allow_eval:
            sys.exit("evaluation pool: needs --allow-eval")
        if subprocess.run(["git", "-C", ROOT, "rev-parse", "-q", "--verify", "refs/tags/bases-registered"], capture_output=True).returncode != 0:
            sys.exit("evaluation pool refused: git tag bases-registered not found (commit and tag Addendum A2 first)")
    profile = a.profile or f"data/sylph/{a.pool}/arm_3a.tsv"
    cwalk = a.crosswalk or ("reference/d6331_gtdb_crosswalk_v2_eval.csv" if a.pool == "eval" else "reference/d6331_gtdb_crosswalk_v2.csv")
    out = a.out or f"reports/sylph/{a.pool}_bestgenome"
    if a.pool == "eval" and os.path.exists(P(out, "bestgenome_report.txt")):
        sys.exit("evaluation outputs are never overwritten")

    acc2sp = {r["accession"]: r["database_species"] for r in tsv(a.candidates)}
    sp2mock, name = {}, {}
    for r in csv.DictReader(open(P(cwalk))):
        name[int(r["mock_taxid"])] = r["mock_species"]
        if r["status"] == "ok":
            sp2mock[r["database_species"]] = int(r["mock_taxid"])
    dna = defaultdict(float)
    for r in csv.DictReader(open(P(a.composition))):
        dna[int(r["species_taxid"])] += float(r["genomic_dna_pct"])
    by, nonmock = defaultdict(list), []
    for r in tsv(profile):
        m = ACC.search(r["Genome_file"]); acc = m.group(0) if m else None
        sp = acc2sp.get(acc)
        if sp is None:
            sys.exit(f"{acc} is not in {a.candidates}; run the gap rule (bin/13) first")
        g = {"acc": acc, "seq": fnum(r["Sequence_abundance"]), "eff": fnum(r.get("Eff_cov"))}
        t = sp2mock.get(sp)
        (by[t] if t is not None else nonmock).append(g)
    called = sorted(by)
    dsum_all = sum(dna.values()); me = sum(dna[t] for t in called)
    e09 = {t: dna[t] / me for t in called}; e04 = {t: dna[t] / dsum_all for t in dna}
    res = {}
    for rule in ("sum", "best"):
        v = {t: (sum(g["seq"] for g in by[t]) if rule == "sum" else max(by[t], key=lambda g: (g["eff"], g["seq"]))["seq"]) for t in called}
        tot = sum(v.values()); o = {t: v[t] / tot for t in called}
        res[rule] = {"o": o, "bc09": bc(o, e09), "bc04": bc(o, e04)}
    os.makedirs(P(out), exist_ok=True)
    L = [f"Sylph best-genome rule vs registered sum rule ({a.pool}; profile {profile}; {len(called)} mock species called, {len(nonmock)} genome rows outside the crosswalk)", "",
         f"Bray-Curtis, bin/04 convention (vs all 17 truth species):  sum {res['sum']['bc04']:.3f}   best-genome {res['best']['bc04']:.3f}",
         f"Bray-Curtis, among called mock species:                    sum {res['sum']['bc09']:.3f}   best-genome {res['best']['bc09']:.3f}", "",
         f"{'species':32s}{'rows':>5s}{'obs/exp sum':>13s}{'obs/exp best':>14s}"]
    rows = []
    for t in sorted(called, key=lambda t: -dna[t]):
        os_, ob = res["sum"]["o"][t] / e09[t], res["best"]["o"][t] / e09[t]
        L.append(f"{name[t][:31]:32s}{len(by[t]):>5d}{os_:>13.3f}{ob:>14.3f}"); rows.append((name[t], len(by[t]), os_, ob))
    s = "\n".join(L)
    open(P(out, "bestgenome_report.txt"), "w").write(s + "\n")
    with open(P(out, "bestgenome_by_species.tsv"), "w") as o:
        o.write("species\tgenome_rows\tobs_exp_sum\tobs_exp_best\n")
        for r in rows:
            o.write(f"{r[0]}\t{r[1]}\t{r[2]:.4f}\t{r[3]:.4f}\n")
    json.dump({"pool": a.pool, "bc04_sum": res["sum"]["bc04"], "bc04_best": res["best"]["bc04"], "bc09_sum": res["sum"]["bc09"], "bc09_best": res["best"]["bc09"]},
              open(P(out, "bestgenome_summary.json"), "w"), indent=1)
    print(s)


if __name__ == "__main__":
    main()
