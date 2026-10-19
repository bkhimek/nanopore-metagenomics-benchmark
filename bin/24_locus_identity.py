#!/usr/bin/env python3
"""24_locus_identity.py - EXPLORATORY, DEV pool only; part of no verdict; changes no registered rule.

Follow-up to bin/23 (reports/amr/dev_locus_check.txt). Result there: reads cover the evaluable genes at about the expected depth (locus reads / E
1.02 for tet(Q), 1.15 cepA, 0.88 tet(W)), but for B. fragilis tet(Q) 33 of 34 locus reads fell below the 0.80 identity filter against ResFinder
(1 supported), against 4 of 18 for P. corporis tet(Q) (14 supported). Cause not established. Two candidate explanations:

  (a) the reference: the supplier-genome gene differs at nucleotide level from every ResFinder allele of the family, so reads at typical error rates
      cannot reach 0.80 against any of them (the AMRFinderPlus call is protein-level, identity >= 90%, and tolerates this);
  (b) the reads: reads at this locus are poorer than reads elsewhere (locus-specific error, a mis-assembled or divergent strain locus).

Hypotheses written BEFORE the first run:
  H3 (reference): the best ResFinder allele of the family has nucleotide identity <= 0.88 to the B. fragilis tet(Q) gene and >= 0.95 to the
     P. corporis tet(Q) gene, AND the median identity of the locus reads to their own supplier genome is >= 0.90 for both carriers.
  H4 (reads): the median identity of the B. fragilis tet(Q) locus reads to the supplier genome is < 0.85 (reads, not the reference, are poor there).
  Anything else is reported as inconclusive. No threshold is tuned afterwards. All evaluable high-confidence units (cepA, tet(Q) x2, tet(W)) are
  reported with the same measurements; H3/H4 are stated for the B. fragilis tet(Q) unit only.

Method (dev only, standard library only):
  gene identity  - the gene is cut from the supplier genome (catalogue coordinates, forward or reverse-complement orientation, chosen by shared
                   15-mers); the five alleles of its ResFinder family with most shared 15-mers are aligned by global alignment (match +1, mismatch -1,
                   gap -2; identity = matches / alignment length); the best is reported with the allele name.
  read identity  - locus reads as in bin/23 (merged alignment >= 80% of the gene in data/trace/dev/dev_vs_refs.paf). For each, matches/block of its
                   records overlapping the gene (genome identity) and matches/block of its best record to any allele of the family in
                   data/amr/dev/dev_vs_resfinder.paf (ResFinder identity); medians, and the share of reads with ResFinder identity >= 0.80.
Reuses bin/23 helpers and bin/21's family rule and merge list.
"""
import argparse, csv, importlib.util, os, statistics
from collections import defaultdict

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
P = lambda *x: os.path.join(ROOT, *x)
COMP = str.maketrans("ACGTN", "TGCAN")
K = 15


def load(name, fn):
    spec = importlib.util.spec_from_file_location(name, P("bin", fn))
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
    return m


def rc(s):
    return s.translate(COMP)[::-1]


def read_fasta(path):
    seqs, name = {}, None
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith(">"):
            name = line[1:].split()[0]; seqs[name] = []
        elif name:
            seqs[name].append(line.strip().upper())
    return {k: "".join(v) for k, v in seqs.items()}


def kmers(s):
    return {s[i:i + K] for i in range(len(s) - K + 1)}


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
    i, j, match, alen = n, m, 0, 0
    while i > 0 or j > 0:
        t = tb[i][j]
        if t == 0:
            match += a[i - 1] == b[j - 1]; i -= 1; j -= 1
        elif t == 1:
            i -= 1
        else:
            j -= 1
        alen += 1
    return match / alen if alen else 0.0


def best_allele(gene, alleles):
    """orientation by shared k-mers; align the five alleles with most shared k-mers; return (identity, allele, orientation)."""
    kf, kr = kmers(gene), kmers(rc(gene))
    scored = []
    for name, s in alleles.items():
        ks = kmers(s)
        scored.append((max(len(ks & kf), len(ks & kr)), name, len(ks & kr) > len(ks & kf)))
    scored.sort(reverse=True)
    best = (0.0, "", False)
    for _, name, rev in scored[:5]:
        idv = global_identity(rc(gene) if rev else gene, alleles[name])
        if idv > best[0]:
            best = (idv, name, rev)
    return best


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--catalogue", default="reports/amr/reference_catalogue.tsv")
    ap.add_argument("--refdir", default="data/reference/D6331.refseq/genomes")
    ap.add_argument("--genome-paf", default="data/trace/dev/dev_vs_refs.paf")
    ap.add_argument("--resfinder-paf", default="data/amr/dev/dev_vs_resfinder.paf")
    ap.add_argument("--fsa", default="data/db/resfinder_db/all.fsa")
    ap.add_argument("--merge", default="reference/amr_family_merge.tsv")
    ap.add_argument("--map", default="reports/amr/family_map_A5.tsv")
    ap.add_argument("--out", default="reports/amr")
    a = ap.parse_args()
    x21, x23 = load("x21", "21_amr_reads.py"), load("x23", "23_locus_check.py")
    for r in csv.DictReader(open(P(a.merge)), delimiter="\t"):
        x21.MERGE[r["alias"]] = r["family"]
    fmap = {(r["species"], r["symbol"]): r for r in csv.DictReader(open(P(a.map)), delimiter="\t")}
    units = []
    for r in csv.DictReader(open(P(a.catalogue)), delimiter="\t"):
        f = fmap.get((r["species"], r["symbol"]))
        if r["category"] == "high-confidence" and f and f["status"] == "mapped" and f["evaluable"] == "yes":
            units.append(dict(r, family=f["resfinder_family"]))
    fsa = read_fasta(P(a.fsa))
    fam_alleles = defaultdict(dict)
    for h, s in fsa.items():
        fam_alleles[x21.family(h)][h] = s
    ivs = defaultdict(lambda: defaultdict(list))   # target -> read -> [(ts, te, mt, bl)]
    for u in units:
        n = x23.contig_index(P(a.refdir, u["genome"] + ".fasta"))[u["contig"]]
        u["target"] = f"{u['genome']}__{n}"
    want = {u["target"] for u in units}
    with open(P(a.genome_paf)) as fh:
        for line in fh:
            f = line.split("\t", 12)
            if f[5] in want:
                ivs[f[5]][f[0]].append((int(f[7]), int(f[8]), int(f[9]), int(f[10])))
    locus = {}
    for u in units:
        s0, e0 = int(u["start"]) - 1, int(u["stop"]); g = e0 - s0
        u["reads"] = [rd for rd, v in ivs[u["target"]].items()
                      if x23.merged_overlap([(x[0], x[1]) for x in v], s0, e0) >= 0.8 * g]
        locus[u["target"], u["symbol"]] = set(u["reads"])
    allreads = {rd for u in units for rd in u["reads"]}
    rf_best = defaultdict(dict)                    # read -> family -> (matches, identity)
    with open(P(a.resfinder_paf)) as fh:
        for line in fh:
            f = line.split("\t", 12)
            if f[0] not in allreads:
                continue
            mt, bl = int(f[9]), int(f[10])
            fam = x21.family(f[5])
            if bl and (fam not in rf_best[f[0]] or mt > rf_best[f[0]][fam][0]):
                rf_best[f[0]][fam] = (mt, mt / bl)
    rows = []
    for u in units:
        contig = read_fasta(P(a.refdir, u["genome"] + ".fasta"))[u["contig"]]
        gene = contig[int(u["start"]) - 1:int(u["stop"])]
        idv, allele, rev = best_allele(gene, fam_alleles[u["family"]])
        s0, e0 = int(u["start"]) - 1, int(u["stop"])
        gid, rid = [], []
        for rd in u["reads"]:
            rec = [x for x in ivs[u["target"]][rd] if x[0] < e0 and x[1] > s0]
            if rec and sum(x[3] for x in rec):
                gid.append(sum(x[2] for x in rec) / sum(x[3] for x in rec))
            if u["family"] in rf_best[rd]:
                rid.append(rf_best[rd][u["family"]][1])
        med = lambda v: statistics.median(v) if v else float("nan")
        rows.append(dict(species=u["species"], gene=u["symbol"], family=u["family"], gene_len=e0 - s0, best_allele=allele,
                         gene_identity=idv, locus_reads=len(u["reads"]), median_read_genome_id=med(gid),
                         reads_with_resfinder_aln=len(rid), median_read_resfinder_id=med(rid),
                         share_resfinder_id_ge_080=(sum(1 for x in rid if x >= 0.80) / len(rid)) if rid else float("nan")))
    L = ["Locus identity check (EXPLORATORY, dev pool; hypotheses H3/H4 in the script header)", "",
         f"{'species':26s}{'gene':9s}{'best ResFinder allele':27s}{'gene id':>8s}{'locus':>6s}{'read~genome':>12s}{'read~RF':>8s}{'n RF':>5s}{'>=0.80':>7s}"]
    for r in rows:
        L.append(f"{r['species'][:25]:26s}{r['gene'][:8]:9s}{r['best_allele'][:26]:27s}{r['gene_identity']:>8.3f}{r['locus_reads']:>6d}"
                 f"{r['median_read_genome_id']:>12.3f}{r['median_read_resfinder_id']:>8.3f}{r['reads_with_resfinder_aln']:>5d}{r['share_resfinder_id_ge_080']:>7.2f}")
    bf = [r for r in rows if r["gene"] == "tet(Q)" and r["species"].startswith("Bacteroides")]
    pc = [r for r in rows if r["gene"] == "tet(Q)" and r["species"].startswith("Prevotella")]
    L.append("")
    if bf and pc:
        b, p = bf[0], pc[0]
        h3 = b["gene_identity"] <= 0.88 and p["gene_identity"] >= 0.95 and b["median_read_genome_id"] >= 0.90 and p["median_read_genome_id"] >= 0.90
        h4 = b["median_read_genome_id"] < 0.85
        L.append(f"H3 (reference): B. fragilis tet(Q) gene id {b['gene_identity']:.3f} (<= 0.88?), P. corporis {p['gene_identity']:.3f} (>= 0.95?), "
                 f"read~genome medians {b['median_read_genome_id']:.3f} / {p['median_read_genome_id']:.3f} (both >= 0.90?) -> {'MET' if h3 else 'NOT MET'}")
        L.append(f"H4 (reads): B. fragilis tet(Q) read~genome median {b['median_read_genome_id']:.3f} (< 0.85?) -> {'MET' if h4 else 'NOT MET'}")
        L.append("Outcome: " + ("H3" if h3 and not h4 else "H4" if h4 and not h3 else "inconclusive (neither or both met)"))
    else:
        L.append("B. fragilis / P. corporis tet(Q) units not found: hypotheses cannot be scored.")
    L += ["", "Exploratory; gene identity is a plain global alignment of the supplier-genome gene against the five closest ResFinder alleles; "
          "read identities are matches/block from minimap2 (block includes gaps). 'read~RF' is over locus reads having any alignment to the family."]
    txt = "\n".join(L)
    os.makedirs(P(a.out), exist_ok=True)
    open(P(a.out, "dev_locus_identity.txt"), "w").write(txt + "\n")
    with open(P(a.out, "dev_locus_identity.tsv"), "w") as o:
        w = csv.DictWriter(o, fieldnames=list(rows[0]), delimiter="\t"); w.writeheader(); w.writerows(rows)
    print(txt)


if __name__ == "__main__":
    main()
