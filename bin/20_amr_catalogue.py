#!/usr/bin/env python3
"""20_amr_catalogue.py - AMR reference catalogue step 2 (ADDENDUM A3 rules, applied exactly as registered; reference genomes + dev coverage table only).

Input : reports/amr/reference_hits_raw.tsv (bin/19), reports/sylph/dev_diagnostic/base_basis_by_species.tsv (bin/16; dev pool coverage per species).
Output: reports/amr/reference_catalogue.tsv (one row per AMR-relevant hit, category + carrier coverage), reports/amr/catalogue_report.txt.

Rules registered in A3 BEFORE the hits were seen (column names are those of AMRFinderPlus 4.0.23: 'Type' is the A3 'Element type'):
  high-confidence : Scope == core, Type == AMR, Subtype == AMR (not POINT), % Identity >= 90, % Coverage >= 90, Method not PARTIAL*,
                    and 'Element name' does not match /efflux|pump|transporter/i
  candidate       : every other hit with Type == AMR (partial, plus-scope, identity/coverage < 90, efflux-like) or Subtype == POINT
  not reported    : Type STRESS (metal, acid, biocide) and VIRULENCE hits are kept in the raw table, outside the AMR catalogue.
Opportunity (A3 step 3): carrier species coverage on the DEV pool from bin/16 (species' own mean read length). A hit is 'evaluable' only if the
carrier's coverage is >= --min-cov (default 3.0, a planning value; the registered minimum is fixed in A4 and this field is recomputed then).
Observation made after seeing the hits, not a rule: the A3 text did not say how INTERNAL_STOP is treated; the registered rule leaves it eligible
(it is not PARTIAL*), and it is reported as is. Standard library only.
"""
import argparse, csv, os, re
from collections import Counter

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
P = lambda *x: os.path.join(ROOT, *x)
EFFLUX = re.compile(r"efflux|pump|transporter", re.I)
ALIAS = {"Lactobacillus": "Limosilactobacillus fermentum", "Candida": "Candida albicans"}


def species_of(stem):
    g = stem.split("_")
    return ALIAS.get(g[0], " ".join(g[:2]))


def category(x):
    if x["Type"] != "AMR" and x["Subtype"] != "POINT":
        return None
    ident, cov = float(x["% Identity to reference"]), float(x["% Coverage of reference"])
    if (x["Scope"] == "core" and x["Type"] == "AMR" and x["Subtype"] == "AMR" and ident >= 90 and cov >= 90
            and not x["Method"].startswith("PARTIAL") and not EFFLUX.search(x["Element name"])):
        return "high-confidence"
    return "candidate"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--raw", default="reports/amr/reference_hits_raw.tsv")
    ap.add_argument("--coverage", default="reports/sylph/dev_diagnostic/base_basis_by_species.tsv")
    ap.add_argument("--min-cov", type=float, default=3.0)
    ap.add_argument("--out", default="reports/amr")
    a = ap.parse_args()
    cov = {r["species"]: float(r["cov"]) for r in csv.DictReader(open(P(a.coverage)), delimiter="\t")}
    raw = list(csv.DictReader(open(P(a.raw)), delimiter="\t"))
    rows = []
    for x in raw:
        c = category(x)
        if c is None:
            continue
        sp = species_of(x["genome"]); cv = cov.get(sp)
        rows.append({"genome": x["genome"], "species": sp, "category": c, "symbol": x["Element symbol"], "name": x["Element name"],
                     "class": x["Class"], "scope": x["Scope"], "subtype": x["Subtype"], "method": x["Method"], "identity": x["% Identity to reference"],
                     "ref_coverage": x["% Coverage of reference"], "contig": x["Contig id"], "start": x["Start"], "stop": x["Stop"],
                     "carrier_dev_cov": "" if cv is None else f"{cv:.2f}",
                     "evaluable": "no (carrier not in dev trace)" if cv is None else ("yes" if cv >= a.min_cov else f"no (carrier cov < {a.min_cov:g})")})
    os.makedirs(P(a.out), exist_ok=True)
    with open(P(a.out, "reference_catalogue.tsv"), "w") as o:
        w = csv.DictWriter(o, fieldnames=list(rows[0]), delimiter="\t"); w.writeheader(); w.writerows(rows)
    n = Counter(r["category"] for r in rows)
    L = [f"AMR reference catalogue (A3 rules; {len(raw)} raw hits, {len(rows)} AMR/POINT rows: {dict(n)})", ""]
    for cat in ("high-confidence", "candidate"):
        L.append(f"== {cat} ==")
        for r in rows:
            if r["category"] == cat:
                L.append(f"{r['genome'][:26]:27s}{r['symbol'][:16]:17s}{r['class'][:22]:23s}{r['scope']:5s}{r['subtype']:6s}{r['method']:14s}id {r['identity']:>6s} cov {r['ref_coverage']:>6s}  carrier {r['carrier_dev_cov'] or '-':>6s}  eval: {r['evaluable']}")
        L.append("")
    ev = Counter((r["category"], r["evaluable"].split(" ")[0]) for r in rows)
    L.append(f"Evaluable (carrier dev coverage >= {a.min_cov:g}x, planning value): {dict(ev)}")
    L.append("Species with a high-confidence reference ARG: " + ", ".join(sorted({r['species'] for r in rows if r['category'] == 'high-confidence'})))
    L.append("Reference annotation, not biological truth; acquired vs intrinsic is not separated automatically (A3).")
    txt = "\n".join(L); open(P(a.out, "catalogue_report.txt"), "w").write(txt + "\n"); print(txt)


if __name__ == "__main__":
    main()
