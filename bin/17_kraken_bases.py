#!/usr/bin/env python3
"""17_kraken_bases.py - EXPLORATORY, DEVELOPMENT pool only; post hoc follow-up to bin/16; not a registered arm, part of no verdict.

Why: bin/16 showed that species differ in mean read length (L. fermentum 1.9 kb ... 8 kb), so a read-count share is not a DNA-mass share. The manufacturer's
composition is a DNA-mass share. Every Kraken2 abundance reported so far (bin/04 and the arms) is a read share. Here the SAME Kraken2 classification is re-expressed
as a share of BASES: each read counts with its length (column 4 of the Kraken2 output). Nothing is reclassified and the call set is NOT changed.

Rules (fixed here, before the first run):
  - dev pool only (the script refuses an eval run). Species of a read = its Kraken2 taxon if that is a species, else its species-rank ancestor if any, else none.
  - taxid aliases (reference/d6331_taxid_aliases.csv) are applied as in bin/04. Called species = species with at least ceil(min_frac x all reads) READS
    (frozen 0.0008, as registered); the same called set is used for the read and the base shares.
  - Bray-Curtis in the bin/04 convention (called species + 17 truth species, expected = DNA% / sum of DNA%) and, to compare with bin/14-16 and Sylph, in the
    convention of bin/09 (only called mock species, both sides renormalised).
  - Consistency check: read counts per truth species must equal the observed_reads of the species table in reports/ (printed; mismatches are flagged).
Predictions written before the first run (from bin/16, where the same change of basis moved the mapping-based shares from Bray-Curtis 0.120 to 0.184):
  (P1) Kraken2 Bray-Curtis (bin/04 convention) rises when counted in bases, by at least 0.02 (expected roughly +0.03 to +0.08).
  (P2) L. fermentum obs/exp (share among called mock species) falls from about 1.7-1.9 (reads) to 0.55-0.80 (bases); F. prausnitzii falls from about 1.2-1.4 to
       about 0.7-0.9; Roseburia stays below 0.5.
  (P3) the Kraken2 minus Sylph gap in Bray-Curtis among called mock species (Sylph arm 3a dev: 0.240) narrows when Kraken2 is counted in bases. Reported, not tested.
Standard library only.
"""
import argparse, collections, csv, gzip, importlib.util, json, math, os, re, sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
P = lambda *x: os.path.join(ROOT, *x)


def load08():
    spec = importlib.util.spec_from_file_location("x08", P("bin", "08_trace_reads.py"))
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
    return m


def bc(o, e):
    return 1 - sum(min(o.get(k, 0.0), e.get(k, 0.0)) for k in set(o) | set(e))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run", nargs="+", default=["data/classify/dev/baseline_c0"], help="dev classify run folder(s)")
    ap.add_argument("--db", default="data/db/pluspf_08_20260626", help="folder holding taxonomy/nodes.dmp (all runs must share the taxonomy)")
    ap.add_argument("--truth", default="reference/d6331_expected_reads_by_species.csv")
    ap.add_argument("--aliases", default="reference/d6331_taxid_aliases.csv")
    ap.add_argument("--min-frac", type=float, default=0.0008)
    ap.add_argument("--table", default="reports/baseline_kraken2_pluspf8/dev_species_table.tsv", help="species table for the consistency check (first run only)")
    ap.add_argument("--out", default="reports/kraken_bases")
    a = ap.parse_args()
    x08 = load08()
    parent, rank = x08.load_nodes(x08.find_dmp(P(a.db), "nodes.dmp"))
    lineage = x08.make_lineage(parent)

    truth = {}
    for t in csv.DictReader(open(P(a.truth))):
        truth[int(t["species_taxid"])] = {"name": t["ncbi_name_in_kraken_db"], "dna": float(t["genomic_dna_pct"])}
    dsum = sum(v["dna"] for v in truth.values())
    alias = {}
    if os.path.exists(P(a.aliases)):
        for x in csv.DictReader(open(P(a.aliases))):
            alias[int(x["accepted_taxid"])] = int(x["species_taxid"])
    spcache = {}

    def species_of(k):
        if k not in spcache:
            s = next((x for x in lineage(k) if rank.get(x) == "species"), None) if k else None
            spcache[k] = alias.get(s, s) if s is not None else None
        return spcache[k]

    os.makedirs(P(a.out), exist_ok=True)
    lines, summary = [], {}
    for run in a.run:
        rec = json.load(open(P(run, "run.json")))
        if rec.get("pool") != "dev":
            sys.exit(f"{run}: pool is {rec.get('pool')}; this script is dev only")
        label = rec["label"]
        reads, bases = collections.Counter(), collections.Counter()
        N = Nb = 0
        with gzip.open(P(run, "kraken_output.txt.gz"), "rt") as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) < 4:
                    continue
                m = re.search(r"taxid (\d+)", f[2]); k = 0 if f[0] == "U" else (int(m.group(1)) if m else int(f[2]))
                L = int(f[3]); N += 1; Nb += L
                s = species_of(k)
                if s is not None:
                    reads[s] += 1; bases[s] += L
        T = max(1, math.ceil(a.min_frac * N))
        called = {s for s, c in reads.items() if c >= T}
        # bin/04 convention
        keys = called | set(truth)
        e04 = {t: (truth[t]["dna"] / dsum if t in truth else 0.0) for t in keys}
        res = {}
        for nm, cnt in (("reads", reads), ("bases", bases)):
            tot = sum(cnt[s] for s in called)
            o04 = {t: (cnt[t] / tot if t in called else 0.0) for t in keys}
            mock = sorted(called & set(truth)); mt = sum(cnt[s] for s in mock); me = sum(truth[s]["dna"] for s in mock)
            o09 = {s: cnt[s] / mt for s in mock}; e09 = {s: truth[s]["dna"] / me for s in mock}
            res[nm] = {"bc04": bc(o04, e04), "bc09": bc(o09, e09), "obs_exp": {s: o09[s] / e09[s] for s in mock}}
        mean_len = {s: bases[s] / reads[s] for s in reads if reads[s]}
        summary[label] = {"N": N, "T": T, "called": len(called), "bc04_reads": res["reads"]["bc04"], "bc04_bases": res["bases"]["bc04"],
                          "bc09_reads": res["reads"]["bc09"], "bc09_bases": res["bases"]["bc09"]}
        with open(P(a.out, f"dev_{label}_by_species.tsv"), "w") as o:
            o.write("taxid\tspecies\treads\tbases\tmean_read_len\tobs_exp_reads\tobs_exp_bases\n")
            for s in sorted(called & set(truth), key=lambda s: -truth[s]["dna"]):
                o.write(f"{s}\t{truth[s]['name']}\t{reads[s]}\t{bases[s]}\t{mean_len[s]:.1f}\t{res['reads']['obs_exp'][s]:.3f}\t{res['bases']['obs_exp'][s]:.3f}\n")
        lines += [f"=== run {label} (dev): {N} reads, {Nb/1e6:.1f} Mb, T = {T} reads, {len(called)} species called ({len(called & set(truth))} mock, {len(called - set(truth))} other)",
                  f"Bray-Curtis, bin/04 convention:           reads {res['reads']['bc04']:.3f}   bases {res['bases']['bc04']:.3f}   (change {res['bases']['bc04'] - res['reads']['bc04']:+.3f})",
                  f"Bray-Curtis, among called mock species:   reads {res['reads']['bc09']:.3f}   bases {res['bases']['bc09']:.3f}   (change {res['bases']['bc09'] - res['reads']['bc09']:+.3f})",
                  f"{'species':32s}{'reads':>8s}{'len':>7s}{'obs/exp reads':>15s}{'obs/exp bases':>15s}"]
        for s in sorted(called & set(truth), key=lambda s: mean_len[s]):
            lines.append(f"{truth[s]['name'][:31]:32s}{reads[s]:>8d}{mean_len[s]:>7.0f}{res['reads']['obs_exp'][s]:>15.2f}{res['bases']['obs_exp'][s]:>15.2f}")
        hum = reads.get(9606, 0)
        lines.append(f"(human: {hum} reads, {bases.get(9606, 0)/1e6:.2f} Mb; counted as a called non-mock species)")
        if run == a.run[0] and os.path.exists(P(a.table)):
            bad = []
            for r in csv.DictReader(open(P(a.table)), delimiter="\t"):
                t = int(r["taxid"])
                if t in truth and int(r["observed_reads"]) != reads.get(t, 0):
                    bad.append((r["species"], int(r["observed_reads"]), reads.get(t, 0)))
            lines.append("Consistency with the species table: " + ("all truth-species read counts equal" if not bad else f"MISMATCH {bad}"))
        lines.append("")
    lines.append("Sylph arm 3a (dev, among called mock species, bin/14): 0.240 (registered sum rule), 0.182 (best-genome what-if).")
    lines.append("Mean read length uses the full read; the call set is the registered read-count threshold. Exploratory.")
    s = "\n".join(lines); open(P(a.out, "kraken_bases_report.txt"), "w").write(s + "\n")
    json.dump(summary, open(P(a.out, "kraken_bases_summary.json"), "w"), indent=1); print(s)


if __name__ == "__main__":
    main()
