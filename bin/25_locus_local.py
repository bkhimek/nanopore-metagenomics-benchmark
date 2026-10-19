#!/usr/bin/env python3
"""25_locus_local.py - EXPLORATORY, DEV pool only; part of no verdict; changes no registered rule.

Follow-up to bin/23 and bin/24 (reports/amr/dev_locus_check.txt, dev_locus_identity.txt). Open question: the 34 B. fragilis tet(Q) locus reads reach a
median identity of 0.775 against their best ResFinder allele (3% >= 0.80; P. corporis tet(Q) 0.910, 78%), although the supplier-genome gene is 97.6%
identical to that allele and the reads match their genome at a median of 0.942. Neither a divergent reference (H3) nor poor reads (H4) was supported;
bin/24's read-to-genome control used the whole alignment, not the gene. Note also that every identity so far (bin/08, 15, 21, 23, 24) is minimap2's PAF
matches/block WITHOUT -c: an anchor-based approximation, not a base-level alignment.

Method: the locus reads of each evaluable high-confidence unit (as bin/23) are aligned at BASE LEVEL (minimap2 -x map-ont -c --eqx, container as bin/21)
to two gene-only references, separately: (i) the supplier-genome gene (catalogue coordinates) and (ii) the best ResFinder allele (bin/24 output). Per
read the best alignment (most '=' bases) to each reference; identity = '=' / ('=' + X + I + D) over the aligned columns; span = aligned reference
bases / reference length. Also: error composition (mismatch / insertion / deletion share of aligned columns) and a profile of the error rate along
the supplier gene in 100 bp windows (pooled over the unit's reads), and the registered-style approximate identity (PAF matches/block of the best
family record in data/amr/dev/dev_vs_resfinder.paf) for the same reads.

Hypotheses written BEFORE the first run (id_gene / id_allele = median base-level identity of the unit's reads to the supplier gene / the ResFinder allele):
  H5 (reads locally poor): B. fragilis tet(Q) id_gene <= 0.82 AND P. corporis tet(Q) id_gene >= 0.88.
  H6 (approximation artefact): B. fragilis tet(Q) id_gene >= 0.88 AND id_allele >= 0.88: at base level the reads match the allele as well as the
      other carriers do, so the 0.775 was produced by minimap2's approximate identity, not by the sequence.
  H7 (the allele really differs): B. fragilis tet(Q) id_gene >= 0.88 AND id_allele <= 0.82.
  Anything else (including 0.82-0.88 values) is reported as inconclusive. No threshold is tuned afterwards. The same measurements are reported for
  cepA, tet(W) and P. corporis tet(Q); the hypotheses are stated for B. fragilis tet(Q) only. The registered 0.80 filter is NOT changed by this script.
Standard library only; reuses bin/23 helpers and bin/21's family rule and merge list.
"""
import argparse, csv, gzip, importlib.util, os, re, statistics, subprocess
from collections import defaultdict

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
P = lambda *x: os.path.join(ROOT, *x)
CIG = re.compile(r"(\d+)([=XIDSHNM])")
WIN = 100


def load(name, fn):
    spec = importlib.util.spec_from_file_location(name, P("bin", fn))
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
    return m


def read_fasta(path):
    seqs, name = {}, None
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith(">"):
            name = line[1:].split()[0]; seqs[name] = []
        elif name:
            seqs[name].append(line.strip().upper())
    return {k: "".join(v) for k, v in seqs.items()}


def parse_cg(paf_line):
    f = paf_line.rstrip("\n").split("\t")
    cg = next((x[5:] for x in f[12:] if x.startswith("cg:Z:")), None)
    return f, cg


def analyse(cg, ts, tlen, nwin):
    """walk an eqx CIGAR from reference position ts; return counts and per-window [aligned, errors]."""
    eq = x = ins = dele = 0
    w = [[0, 0] for _ in range(nwin)]
    t = ts
    for n, op in ((int(a), b) for a, b in CIG.findall(cg)):
        if op in "=X":
            for k in range(n):
                i = min((t + k) // WIN, nwin - 1); w[i][0] += 1; w[i][1] += (op == "X")
            t += n
            if op == "=":
                eq += n
            else:
                x += n
        elif op == "D":
            for k in range(n):
                i = min((t + k) // WIN, nwin - 1); w[i][0] += 1; w[i][1] += 1
            t += n; dele += n
        elif op == "I":
            i = min(t // WIN, nwin - 1); w[i][1] += n; ins += n
        elif op == "M":
            raise SystemExit("CIGAR has M operations: run minimap2 with --eqx")
    return eq, x, ins, dele, w


def best_alignment(paf, nwin_by_ref):
    """paf path -> read -> (identity, span, eq, x, ins, dele, windows) of the alignment with most '=' bases."""
    best = {}
    for line in open(paf):
        f, cg = parse_cg(line)
        if cg is None:
            continue
        tlen, ts, te = int(f[6]), int(f[7]), int(f[8])
        eq, x, ins, dele, w = analyse(cg, ts, tlen, nwin_by_ref(f[5]))
        cols = eq + x + ins + dele
        if cols and (f[0] not in best or eq > best[f[0]][2]):
            best[f[0]] = (eq / cols, (te - ts) / tlen, eq, x, ins, dele, w)
    return best


def run_minimap2(image, ref, reads, paf, threads):
    for d in (ref, reads, paf):
        os.makedirs(os.path.dirname(d), exist_ok=True)
    if subprocess.run(["docker", "image", "inspect", image], capture_output=True).returncode != 0:
        subprocess.run(["docker", "pull", image], check=True)
    cmd = ["docker", "run", "--rm", "-u", f"{os.getuid()}:{os.getgid()}", "-m", "4g", "--cpus", str(threads),
           "-v", f"{os.path.dirname(ref)}:/ref:ro", "-v", f"{os.path.dirname(reads)}:/reads:ro", "-v", f"{os.path.dirname(paf)}:/out",
           image, "minimap2", "-x", "map-ont", "-c", "--eqx", "-N", "5", "-t", str(threads), "-o", f"/out/{os.path.basename(paf)}",
           f"/ref/{os.path.basename(ref)}", f"/reads/{os.path.basename(reads)}"]
    subprocess.run(cmd, check=True)
    return " ".join(cmd)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--catalogue", default="reports/amr/reference_catalogue.tsv")
    ap.add_argument("--refdir", default="data/reference/D6331.refseq/genomes")
    ap.add_argument("--genome-paf", default="data/trace/dev/dev_vs_refs.paf")
    ap.add_argument("--resfinder-paf", default="data/amr/dev/dev_vs_resfinder.paf")
    ap.add_argument("--fsa", default="data/db/resfinder_db/all.fsa")
    ap.add_argument("--merge", default="reference/amr_family_merge.tsv")
    ap.add_argument("--map", default="reports/amr/family_map_A5.tsv")
    ap.add_argument("--identity", default="reports/amr/dev_locus_identity.tsv", help="bin/24 output (best allele per unit)")
    ap.add_argument("--reads", default="data/subsets/SRR31222153.dev.seed20261012.blocks.fastq.gz")
    ap.add_argument("--workdir", default="data/amr/dev/locus")
    ap.add_argument("--out", default="reports/amr")
    ap.add_argument("--image", default="quay.io/biocontainers/minimap2:2.31--h118bc1c_0")
    ap.add_argument("--threads", type=int, default=2)
    ap.add_argument("--use-existing", action="store_true", help="use PAFs already in --workdir (tests)")
    a = ap.parse_args()
    x21, x23 = load("x21", "21_amr_reads.py"), load("x23", "23_locus_check.py")
    for r in csv.DictReader(open(P(a.merge)), delimiter="\t"):
        x21.MERGE[r["alias"]] = r["family"]
    fmap = {(r["species"], r["symbol"]): r for r in csv.DictReader(open(P(a.map)), delimiter="\t")}
    bestal = {(r["species"], r["gene"]): r["best_allele"] for r in csv.DictReader(open(P(a.identity)), delimiter="\t")}
    units = []
    for r in csv.DictReader(open(P(a.catalogue)), delimiter="\t"):
        f = fmap.get((r["species"], r["symbol"]))
        if r["category"] == "high-confidence" and f and f["status"] == "mapped" and f["evaluable"] == "yes":
            u = dict(r, family=f["resfinder_family"], allele=bestal[(r["species"], r["symbol"])])
            u["target"] = f"{u['genome']}__{x23.contig_index(P(a.refdir, u['genome'] + '.fasta'))[u['contig']]}"
            u["tag"] = re.sub(r"[^A-Za-z0-9]+", "_", f"{u['species']}_{u['symbol']}").strip("_")
            units.append(u)
    ivs = defaultdict(lambda: defaultdict(list))
    want = {u["target"] for u in units}
    with open(P(a.genome_paf)) as fh:
        for line in fh:
            f = line.split("\t", 12)
            if f[5] in want:
                ivs[f[5]][f[0]].append((int(f[7]), int(f[8])))
    fsa = read_fasta(P(a.fsa))
    wd = P(a.workdir); os.makedirs(wd, exist_ok=True)
    for u in units:
        s0, e0 = int(u["start"]) - 1, int(u["stop"]); g = e0 - s0
        u["reads"] = {rd for rd, v in ivs[u["target"]].items() if x23.merged_overlap(v, s0, e0) >= 0.8 * g}
        contig = read_fasta(P(a.refdir, u["genome"] + ".fasta"))[u["contig"]]
        u["gene_seq"], u["allele_seq"] = contig[s0:e0], fsa[u["allele"]]
    meta = {}
    if not a.use_existing:
        allrd = {rd: None for u in units for rd in u["reads"]}
        got = {}
        with gzip.open(P(a.reads), "rt") as fh:
            while True:
                h = fh.readline()
                if not h:
                    break
                rec = [h, fh.readline(), fh.readline(), fh.readline()]
                if h[1:].split()[0] in allrd:
                    got[h[1:].split()[0]] = "".join(rec)
        missing = [r for r in allrd if r not in got]
        if missing:
            raise SystemExit(f"{len(missing)} locus reads not found in {a.reads}")
        for u in units:
            open(os.path.join(wd, u["tag"] + ".fastq"), "w").write("".join(got[r] for r in sorted(u["reads"])))
            for kind, seq, name in (("gene", u["gene_seq"], u["tag"] + "_gene"), ("allele", u["allele_seq"], u["allele"])):
                open(os.path.join(wd, f"{u['tag']}_{kind}.fa"), "w").write(f">{name}\n{seq}\n")
        for u in units:
            for kind in ("gene", "allele"):
                meta[u["tag"] + kind] = run_minimap2(a.image, os.path.join(wd, f"{u['tag']}_{kind}.fa"), os.path.join(wd, u["tag"] + ".fastq"),
                                                     os.path.join(wd, f"{u['tag']}_{kind}.paf"), a.threads)
    # approximate (registered-style) identity of the same reads: best family record in the A4 ResFinder PAF
    need = {rd for u in units for rd in u["reads"]}
    appr = defaultdict(dict)
    with open(P(a.resfinder_paf)) as fh:
        for line in fh:
            f = line.split("\t", 12)
            if f[0] in need:
                fam = x21.family(f[5]); mt, bl = int(f[9]), int(f[10])
                if bl and (fam not in appr[f[0]] or mt > appr[f[0]][fam][0]):
                    appr[f[0]][fam] = (mt, mt / bl)
    med = lambda v: statistics.median(v) if v else float("nan")
    rows, wrows, L = [], [], ["Locus local identity (EXPLORATORY, dev pool; base-level alignment of the locus reads to gene-only references; hypotheses H5-H7 in the script header)", ""]
    L.append(f"{'species':26s}{'gene':9s}{'reads':>6s}  {'id gene':>8s}{'>=.80':>6s}{'span':>6s}   {'id allele':>9s}{'>=.80':>6s}{'span':>6s}   {'approx RF id':>12s}{'>=.80':>6s}   {'mism':>7s}{'ins':>7s}{'del':>7s}")
    for u in units:
        nwin = max(1, -(-len(u["gene_seq"]) // WIN))
        res = {}
        for kind, ref in (("gene", u["gene_seq"]), ("allele", u["allele_seq"])):
            nw = max(1, -(-len(ref) // WIN))
            res[kind] = best_alignment(os.path.join(wd, f"{u['tag']}_{kind}.paf"), lambda _r, nw=nw: nw)
        rid = [appr[rd][u["family"]][1] for rd in u["reads"] if u["family"] in appr[rd]]
        stat = {}
        for kind in ("gene", "allele"):
            v = list(res[kind].values())
            stat[kind] = dict(n=len(v), id=med([x[0] for x in v]), ge=(sum(1 for x in v if x[0] >= 0.80) / len(v)) if v else float("nan"), span=med([x[1] for x in v]))
        cols = [sum(x[i] for x in res["gene"].values()) for i in (2, 3, 4, 5)]
        tot = sum(cols) or 1
        mism, ins, dele = cols[1] / tot, cols[2] / tot, cols[3] / tot
        rows.append(dict(species=u["species"], gene=u["symbol"], reads=len(u["reads"]), allele=u["allele"],
                         gene_n=stat["gene"]["n"], id_gene=stat["gene"]["id"], ge080_gene=stat["gene"]["ge"], span_gene=stat["gene"]["span"],
                         allele_n=stat["allele"]["n"], id_allele=stat["allele"]["id"], ge080_allele=stat["allele"]["ge"], span_allele=stat["allele"]["span"],
                         approx_rf_id=med(rid), approx_ge080=(sum(1 for x in rid if x >= 0.80) / len(rid)) if rid else float("nan"),
                         mismatch=mism, insertion=ins, deletion=dele))
        r = rows[-1]
        L.append(f"{u['species'][:25]:26s}{u['symbol'][:8]:9s}{len(u['reads']):>6d}  {r['id_gene']:>8.3f}{r['ge080_gene']:>6.2f}{r['span_gene']:>6.2f}   "
                 f"{r['id_allele']:>9.3f}{r['ge080_allele']:>6.2f}{r['span_allele']:>6.2f}   {r['approx_rf_id']:>12.3f}{r['approx_ge080']:>6.2f}   {mism:>7.3f}{ins:>7.3f}{dele:>7.3f}")
        tw = [[0, 0] for _ in range(nwin)]
        for x in res["gene"].values():
            for i in range(nwin):
                tw[i][0] += x[6][i][0]; tw[i][1] += x[6][i][1]
        for i in range(nwin):
            wrows.append(dict(species=u["species"], gene=u["symbol"], window_start=i * WIN, aligned=tw[i][0], errors=tw[i][1],
                              error_rate=(tw[i][1] / (tw[i][0] + 1e-9)) if tw[i][0] else float("nan")))
        u["wrate"] = [(tw[i][1] / tw[i][0]) if tw[i][0] else float("nan") for i in range(nwin)]
    L += ["", "Error rate along the supplier gene (errors / aligned columns, % per 100 bp window, pooled over the unit's reads):"]
    for u in units:
        L.append(f"  {u['species'][:12]:13s}{u['symbol'][:8]:9s}" + " ".join(f"{100 * x:4.0f}" for x in u["wrate"]))
    bf = next((r for r in rows if r["gene"] == "tet(Q)" and r["species"].startswith("Bacteroides")), None)
    pc = next((r for r in rows if r["gene"] == "tet(Q)" and r["species"].startswith("Prevotella")), None)
    L.append("")
    if bf and pc:
        h5 = bf["id_gene"] <= 0.82 and pc["id_gene"] >= 0.88
        h6 = bf["id_gene"] >= 0.88 and bf["id_allele"] >= 0.88
        h7 = bf["id_gene"] >= 0.88 and bf["id_allele"] <= 0.82
        L.append(f"H5 (reads locally poor): B. fragilis tet(Q) id_gene {bf['id_gene']:.3f} (<= 0.82?), P. corporis {pc['id_gene']:.3f} (>= 0.88?) -> {'MET' if h5 else 'NOT MET'}")
        L.append(f"H6 (approximation artefact): id_gene {bf['id_gene']:.3f} (>= 0.88?), id_allele {bf['id_allele']:.3f} (>= 0.88?) -> {'MET' if h6 else 'NOT MET'}")
        L.append(f"H7 (allele really differs): id_gene {bf['id_gene']:.3f} (>= 0.88?), id_allele {bf['id_allele']:.3f} (<= 0.82?) -> {'MET' if h7 else 'NOT MET'}")
        L.append("Outcome: " + ("H5" if h5 else "H6" if h6 else "H7" if h7 else "inconclusive (no hypothesis met)"))
    else:
        L.append("B. fragilis / P. corporis tet(Q) units not found: hypotheses cannot be scored.")
    L += ["", "Exploratory. Identities are base-level (minimap2 -c --eqx, best alignment per read, gene-only references); 'approx RF id' is the PAF matches/block of the best family "
          "record in the registered-style run (no -c). Reads that did not align to a reference are not in its median. The registered 0.80 filter is unchanged."]
    txt = "\n".join(L)
    os.makedirs(P(a.out), exist_ok=True)
    open(P(a.out, "dev_locus_local.txt"), "w").write(txt + "\n")
    for name, rr in (("dev_locus_local.tsv", rows), ("dev_locus_local_windows.tsv", wrows)):
        with open(P(a.out, name), "w") as o:
            w = csv.DictWriter(o, fieldnames=list(rr[0]), delimiter="\t"); w.writeheader(); w.writerows(rr)
    print(txt)


if __name__ == "__main__":
    main()
