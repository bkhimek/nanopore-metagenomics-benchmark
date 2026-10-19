#!/usr/bin/env python3
"""23_locus_check.py - EXPLORATORY, DEV pool only; part of no verdict; changes no registered rule.

Open question (A6 result): observed support for tet(Q) was 15 reads on dev (22 on eval) against a planning value E of about 29 + 22 = 51 for its two
carriers (B. fragilis, P. corporis); tet(W) 14 against 22; cepA 29 against 32. E assumes the gene is covered as deeply as the species average.
Where do the reads go missing: (a) the gene locus is thinly covered in the sequenced sample, or (b) reads cover it but are lost in the ResFinder step?

Method (dev only): for every high-confidence, evaluable catalogue unit, take the reads that map (bin/08 alignments to the supplier genomes,
data/trace/dev/dev_vs_refs.paf) over >= 80% of the gene interval (reference coordinates from the catalogue; merged alignment intervals per read);
report their number (locus reads), the mean depth over the gene, and the ratio to the species' mean coverage (bin/16). Then follow those same reads in the
ResFinder alignments (data/amr/dev/dev_vs_resfinder.paf, A4/A5 rules) and say what happened to each: supported, ambiguous, no alignment to the gene
family at all, alignment below 0.80 identity, alignment below 0.80 of the gene.

Hypotheses written BEFORE the first run:
  H1 (shortfall arises after mapping to the genomes): locus reads are close to E (ratio locus reads / E >= 0.8 summed over carriers) and most missing
     reads fall in 'identity below 0.80' or 'no alignment': the read-level step, not the sample, loses them.
  H2 (locus thinly covered): locus reads are scarce (ratio <= 0.5) and close to the observed support: the sequenced sample has less of the gene than
     the species average, or the strain differs.
  Between 0.5 and 0.8 is reported as inconclusive. No threshold is tuned afterwards; the numbers are reported as they are.
Standard library only; reuses bin/21's family rule and merge list.
"""
import argparse, csv, importlib.util, os, re, sys
from collections import defaultdict

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
P = lambda *x: os.path.join(ROOT, *x)


def load21():
    spec = importlib.util.spec_from_file_location("x21", P("bin", "21_amr_reads.py"))
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
    return m


def contig_index(fasta):
    idx, n = {}, 0
    with open(fasta) as fh:
        for line in fh:
            if line.startswith(">"):
                n += 1; idx[line[1:].split()[0]] = n
    return idx


def merged_overlap(ivs, a, b):
    """total length of the union of intervals ivs inside [a, b)."""
    tot, cur_s, cur_e = 0, None, None
    for s, e in sorted(ivs):
        s, e = max(s, a), min(e, b)
        if e <= s:
            continue
        if cur_e is None or s > cur_e:
            if cur_e is not None:
                tot += cur_e - cur_s
            cur_s, cur_e = s, e
        else:
            cur_e = max(cur_e, e)
    if cur_e is not None:
        tot += cur_e - cur_s
    return tot


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--catalogue", default="reports/amr/reference_catalogue.tsv")
    ap.add_argument("--refdir", default="data/reference/D6331.refseq/genomes")
    ap.add_argument("--genome-paf", default="data/trace/dev/dev_vs_refs.paf")
    ap.add_argument("--resfinder-paf", default="data/amr/dev/dev_vs_resfinder.paf")
    ap.add_argument("--merge", default="reference/amr_family_merge.tsv")
    ap.add_argument("--map", default="reports/amr/family_map_A5.tsv")
    ap.add_argument("--out", default="reports/amr")
    a = ap.parse_args()
    x21 = load21()
    for r in csv.DictReader(open(P(a.merge)), delimiter="\t"):
        x21.MERGE[r["alias"]] = r["family"]
    fmap = {(r["species"], r["symbol"]): r for r in csv.DictReader(open(P(a.map)), delimiter="\t")}
    units = []  # evaluable high-confidence rows of the catalogue (one per genome row)
    for r in csv.DictReader(open(P(a.catalogue)), delimiter="\t"):
        f = fmap.get((r["species"], r["symbol"]))
        if r["category"] == "high-confidence" and f and f["status"] == "mapped" and f["evaluable"] == "yes":
            units.append(dict(r, family=f["resfinder_family"], E=float(f["E"]), cov=float(f["carrier_cov"])))
    # genome alignments per read and contig, for the units' contigs
    want = defaultdict(list)  # target name -> units
    for u in units:
        n = contig_index(P(a.refdir, u["genome"] + ".fasta"))[u["contig"]]
        u["target"] = f"{u['genome']}__{n}"; want[u["target"]].append(u)
    ivs = defaultdict(lambda: defaultdict(list))  # unit target -> read -> [(ts, te)]
    with open(P(a.genome_paf)) as fh:
        for line in fh:
            f = line.split("\t", 12)
            if f[5] in want:
                ivs[f[5]][f[0]].append((int(f[7]), int(f[8])))
    # resfinder alignments per read, A4 filter reasons by family
    reasons = defaultdict(lambda: defaultdict(lambda: ["none", 0]))  # read -> family -> [status, best matching bases]
    need = {}
    for u in units:
        for rd, v in ivs[u["target"]].items():
            if merged_overlap(v, int(u["start"]) - 1, int(u["stop"])) >= 0.8 * (int(u["stop"]) - int(u["start"]) + 1):
                need.setdefault(rd, set()).add(u["family"])
    best_pass = defaultdict(dict)
    seen = defaultdict(lambda: defaultdict(lambda: "none"))
    with open(P(a.resfinder_paf)) as fh:
        for line in fh:
            f = line.split("\t", 12)
            if f[0] not in need:
                continue
            fam = x21.family(f[5]); qlen, tlen, ts, te, mt, bl = int(f[1]), int(f[6]), int(f[7]), int(f[8]), int(f[9]), int(f[10])
            ok_id = bl > 0 and mt / bl >= 0.80; ok_span = (te - ts) / tlen >= 0.80
            st = "ok" if ok_id and ok_span else ("low_identity" if not ok_id else "low_span")
            cur = seen[f[0]][fam]
            rank = {"none": 0, "low_span": 1, "low_identity": 2, "ok": 3}
            if rank[st] > rank[cur] or (st == "low_identity" and cur == "low_span"):
                seen[f[0]][fam] = st
            if st == "ok":
                if fam not in best_pass[f[0]] or mt > best_pass[f[0]][fam]:
                    best_pass[f[0]][fam] = mt
    rows, L = [], ["Locus check (EXPLORATORY, dev pool; hypotheses in the script header)", "",
                   f"{'species':26s}{'gene':9s}{'E':>6s}{'locus reads':>12s}{'depth/mean':>11s}  {'supported':>9s}{'ambiguous':>10s}{'no aln':>8s}{'low id':>8s}{'low span':>9s}"]
    fam_tot = defaultdict(lambda: [0.0, 0, 0])
    for u in units:
        reads = [rd for rd, v in ivs[u["target"]].items()
                 if merged_overlap(v, int(u["start"]) - 1, int(u["stop"])) >= 0.8 * (int(u["stop"]) - int(u["start"]) + 1)]
        g = int(u["stop"]) - int(u["start"]) + 1
        depth = sum(merged_overlap(v, int(u["start"]) - 1, int(u["stop"])) for v in ivs[u["target"]].values()) / g
        c = {"supported": 0, "ambiguous": 0, "none": 0, "low_identity": 0, "low_span": 0}
        for rd in reads:
            bp = best_pass.get(rd, {})
            if u["family"] in bp:
                ranked = sorted(bp.items(), key=lambda kv: -kv[1])
                if ranked[0][0] == u["family"] and not (len(ranked) > 1 and ranked[1][1] >= 0.9 * ranked[0][1]):
                    c["supported"] += 1
                else:
                    c["ambiguous"] += 1
            else:
                c[seen[rd].get(u["family"], "none")] += 1
        rows.append(dict(species=u["species"], gene=u["symbol"], E=u["E"], locus_reads=len(reads), depth_over_mean=depth / u["cov"], **c))
        L.append(f"{u['species'][:25]:26s}{u['symbol'][:8]:9s}{u['E']:>6.1f}{len(reads):>12d}{depth / u['cov']:>11.2f}  {c['supported']:>9d}{c['ambiguous']:>10d}{c['none']:>8d}{c['low_identity']:>8d}{c['low_span']:>9d}")
        t = fam_tot[u["family"]]; t[0] += u["E"]; t[1] += len(reads); t[2] += c["supported"]
    L += ["", "Per gene family, summed over carriers (locus reads can be shared by two carriers only if the same read maps to both genomes, which the single-best mapping prevents):"]
    for f, (E, lr, sp) in sorted(fam_tot.items()):
        ratio = lr / E if E else float("nan")
        verdict = "H1 (lost after mapping to the genomes)" if ratio >= 0.8 else ("H2 (locus thinly covered)" if ratio <= 0.5 else "inconclusive")
        L.append(f"  {f:8s} E {E:6.1f}  locus reads {lr:4d}  ratio {ratio:5.2f}  supported {sp:4d}  -> {verdict}")
    L += ["", "Exploratory; the locus read count uses the supplier genome coordinates (the sequenced strains may differ), and 'supported' follows the A4/A5 rules unchanged."]
    txt = "\n".join(L)
    os.makedirs(P(a.out), exist_ok=True)
    open(P(a.out, "dev_locus_check.txt"), "w").write(txt + "\n")
    with open(P(a.out, "dev_locus_check.tsv"), "w") as o:
        w = csv.DictWriter(o, fieldnames=list(rows[0]), delimiter="\t"); w.writeheader(); w.writerows(rows)
    print(txt)


if __name__ == "__main__":
    main()
