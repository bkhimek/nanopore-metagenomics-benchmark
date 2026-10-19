#!/usr/bin/env python3
"""31_veillonella_tetq.py - EXPLORATORY look at the tet(Q)-supporting read(s) whose true organism is NOT a catalogue carrier (A8, evaluation pool: SRR31222153.7328964, Veillonella rogosae).
Part of no verdict, changes no registered rule, tunes nothing. The read belongs to the evaluation pool, which was scored under A7 and A8; it is looked at AFTER scoring, one read, and the result changes no count.

Facts already known: the read supports tet(Q) (A7: identity 0.984 to a ResFinder tet(Q) allele, span >= 0.80); the supplier genome Veillonella_rogosae has NO AMRFinderPlus hit; Kraken2 gave the read to
Veillonella rogosae and the genome mapping (bin/26, bin/08 rules) agrees. Question: where does the tet(Q) sequence come from?
Evidence collected (reads and genomes of the mock only):
  1. the read's alignment to ResFinder tet(Q) (A7 PAF): allele, query interval, identity;
  2. the read's alignments to the 21 supplier genomes (bin/26 truth PAF): pieces on the Veillonella genome and on others, with query and genome coordinates;
  3. the tet(Q) allele it matched and the supplier tet(Q) genes (B. fragilis, P. corporis) aligned to the Veillonella genome (minimap2 -c, two sensitivities), to see whether the genome contains a tet(Q)-like gene at all.
Three explanations, decided by the rules below, written BEFORE the first run:
  H_a (the reference genome has the gene and the catalogue missed it): >= 80% of the read's tet(Q) interval lies inside the read's alignment to the Veillonella genome AND a gene-versus-genome alignment
      with identity >= 0.70 and gene coverage >= 0.80 overlaps (+/- 1,000) the read's genome interval on the same contig.
  H_b (insertion in the sampled strain, absent from the reference genome): <= 20% of the tet(Q) interval is covered by alignments to the Veillonella genome AND two alignment pieces (each block >= 500) flank it, one ending
      <= 200 bases after the interval starts and one starting >= 200 bases before it ends, on the same contig and strand, with a genome gap between them of -1,000 to 20,000 bases.
  H_c (anything else: chimeric read, spurious match, or a mapping that cannot be resolved): reported with the numbers.
PREDICTION: H_b. (Basis: no AMRFinderPlus hit in the reference, and the read gets the same organism from two independent methods, Kraken2 and genome mapping; a gene of 98.4% identity to tet(Q) is not a spurious match.)
A prediction that fails stays failed. One read; a worked example; the answer would not change any registered count. Standard library only; minimap2 in Docker as in bin/08 unless --gene-paf is given (tests).
"""
import argparse, csv, gzip, importlib.util, json, os, subprocess, sys, time
from collections import defaultdict
from pathlib import Path

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
P = lambda *x: os.path.join(ROOT, *x)
CARRIERS = ("Bacteroides fragilis", "Prevotella corporis")


def read_fasta(path):
    out, name, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    out[name] = "".join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name is not None:
        out[name] = "".join(buf)
    return out


def overlap(a, b):
    return max(0, min(a[1], b[1]) - max(a[0], b[0]))


def union_len(ivs, window):
    segs = []
    for s, e in sorted((max(s, window[0]), min(e, window[1])) for s, e in ivs if e > window[0] and s < window[1]):
        if segs and s <= segs[-1][1]:
            segs[-1][1] = max(segs[-1][1], e)
        else:
            segs.append([s, e])
    return sum(e - s for s, e in segs)


def run_minimap2(image, ref, query, paf, extra, threads, mem):
    ref, query, paf = Path(ref).resolve(), Path(query).resolve(), Path(paf).resolve()
    if subprocess.run(["docker", "image", "inspect", image], capture_output=True).returncode != 0:
        subprocess.run(["docker", "pull", image], check=True)
    cmd = ["docker", "run", "--rm", "-u", f"{os.getuid()}:{os.getgid()}", "-m", mem, "--cpus", str(threads), "-v", f"{ref.parent}:/ref:ro", "-v", f"{query.parent}:/q:ro",
           "-v", f"{paf.parent}:/out", image, "minimap2"] + extra + ["-c", "--secondary=yes", "-p", "0", "-N", "50", "-t", str(threads), "-o", f"/out/{paf.name}", f"/ref/{ref.name}", f"/q/{query.name}"]
    subprocess.run(cmd, check=True)
    return " ".join(cmd)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pool", default="eval")
    ap.add_argument("--read", default="", help="read id; default: tet(Q)-supporting labelled reads whose true organism is not a catalogue carrier")
    ap.add_argument("--reads-tsv", default="")
    ap.add_argument("--a7-paf", default="")
    ap.add_argument("--truth-paf", default="")
    ap.add_argument("--resfinder", default="data/db/resfinder_db/all.fsa")
    ap.add_argument("--refdir", default="data/reference/D6331.refseq/genomes")
    ap.add_argument("--raw", default="reports/amr/reference_hits_raw.tsv")
    ap.add_argument("--merge", default="reference/amr_family_merge.tsv")
    ap.add_argument("--workdir", default="data/amr/veillonella")
    ap.add_argument("--out", default="reports/amr")
    ap.add_argument("--gene-paf", nargs="*", default=[], help="existing PAFs of the gene queries against the genome (tests); skips minimap2")
    ap.add_argument("--image", default="quay.io/biocontainers/minimap2:2.31--h118bc1c_0")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--memory", default="8g")
    a = ap.parse_args()
    a.reads_tsv = a.reads_tsv or f"reports/amr/{a.pool}_attribution_reads_A8.tsv"
    a.a7_paf = a.a7_paf or f"data/amr/{a.pool}/{a.pool}_vs_resfinder_c.paf"
    a.truth_paf = a.truth_paf or f"data/amr/{a.pool}/attribution/{a.pool}_support_vs_refs.paf"
    s8 = importlib.util.spec_from_file_location("x08", P("bin", "08_trace_reads.py")); x08 = importlib.util.module_from_spec(s8); s8.loader.exec_module(x08)
    s21 = importlib.util.spec_from_file_location("x21", P("bin", "21_amr_reads.py")); x21 = importlib.util.module_from_spec(s21); s21.loader.exec_module(x21)
    for r in csv.DictReader(open(P(a.merge)), delimiter="\t"):
        x21.MERGE[r["alias"]] = r["family"]
    rows = list(csv.DictReader(open(P(a.reads_tsv)), delimiter="\t"))
    picks = [r for r in rows if (r["read"] == a.read if a.read else (r["family"] == "tet(Q)" and r["destination"] != "no_truth" and r["true_species"] not in CARRIERS))]
    if not picks:
        sys.exit("no read selected")
    work = Path(P(a.workdir)); work.mkdir(parents=True, exist_ok=True)
    meta = {"script": "bin/31_veillonella_tetq.py", "args": vars(a), "started": time.strftime("%Y-%m-%d %H:%M:%S")}
    raw = list(csv.DictReader(open(P(a.raw)), delimiter="\t"))
    rf = read_fasta(P(a.resfinder))
    L = ["EXPLORATORY: where does the tet(Q) sequence of a non-carrier read come from? (read examined after scoring; no count changed)", ""]
    outrows = []
    for pk in picks:
        rid, sp = pk["read"], pk["true_species"]
        L.append(f"read {rid}: A8 destination {pk['destination']}, true organism (genome mapping) {sp}, Kraken2 {pk['kraken_taxon']} ({pk['kraken_rank']})")
        # 1. ResFinder tet(Q) alignment (A7 filter)
        best = None; qlen = None
        with open(P(a.a7_paf)) as fh:
            for line in fh:
                f = line.split("\t", 12)
                if f[0] != rid:
                    continue
                qlen = int(f[1]); qs, qe, tname, tlen, ts, te, mt, bl = int(f[2]), int(f[3]), f[5], int(f[6]), int(f[7]), int(f[8]), int(f[9]), int(f[10])
                if bl <= 0 or mt / bl < x21.MIN_ID or (te - ts) / tlen < x21.MIN_SPAN or x21.family(tname) != "tet(Q)":
                    continue
                if best is None or mt > best[0]:
                    best = (mt, tname, qs, qe, mt / bl, (te - ts) / tlen, tlen)
        if best is None:
            L.append("  no tet(Q) alignment passing the A7 filter was found in the PAF; stop."); continue
        _, allele, tqs, tqe, idn, span, alen = best
        L.append(f"  read length {qlen}; tet(Q) alignment: {allele} ({alen} bases), read interval {tqs}-{tqe}, identity {idn:.3f}, gene span {span:.2f}")
        # 2. alignments to the supplier genomes
        pieces, others = [], []
        with open(P(a.truth_paf)) as fh:
            for line in fh:
                f = line.split("\t", 12)
                if f[0] != rid:
                    continue
                qs, qe, strand, tname, tlen, ts, te, mt, bl = int(f[2]), int(f[3]), f[4], f[5], int(f[6]), int(f[7]), int(f[8]), int(f[9]), int(f[10])
                stem, _, n = tname.rpartition("__")
                rec = dict(qs=qs, qe=qe, strand=strand, stem=stem, contig=int(n), tlen=tlen, ts=ts, te=te, mt=mt, bl=bl)
                (pieces if x08.stem_species(stem) == sp else others).append(rec)
        for p in sorted(pieces, key=lambda r: r["qs"]):
            L.append(f"  genome piece on {p['stem']} contig {p['contig']}: read {p['qs']}-{p['qe']} -> genome {p['ts']}-{p['te']} ({p['strand']}), identity {p['mt'] / p['bl']:.3f}, block {p['bl']}")
        for p in sorted(others, key=lambda r: -r["mt"])[:4]:
            L.append(f"  other genome piece on {p['stem']}: read {p['qs']}-{p['qe']} -> {p['ts']}-{p['te']}, identity {p['mt'] / p['bl']:.3f}, block {p['bl']}")
        tq = (tqs, tqe)
        cover = union_len([(p["qs"], p["qe"]) for p in pieces], tq) / (tqe - tqs)
        L.append(f"  share of the tet(Q) read interval covered by alignments to the {sp} genome: {cover:.2f}")
        # 3. gene queries against the genome of the true species
        stem = pieces[0]["stem"] if pieces else None
        genes = {allele: rf.get(allele) or next((v for k, v in rf.items() if k == allele or k.startswith(allele)), None)}
        for h in raw:
            if h["Element symbol"] == "tet(Q)":
                gseq = read_fasta(P(a.refdir, h["genome"] + ".fasta"))[h["Contig id"]][int(h["Start"]) - 1:int(h["Stop"])]
                genes[f"supplier_{h['genome']}"] = gseq
        genes = {k: v for k, v in genes.items() if v}
        gall = {}
        gfa = work / f"{rid}_genes.fa"
        with open(gfa, "w") as o:
            for k, v in genes.items():
                o.write(f">{k}\n{v}\n")
        pafs = list(a.gene_paf)
        if not pafs and stem:
            gen = P(a.refdir, stem + ".fasta")
            for tag, extra in (("asm20", ["-x", "asm20"]), ("sens", ["-k", "11", "-w", "4", "-x", "map-ont"])):
                pf = work / f"{rid}_genes_vs_genome_{tag}.paf"
                meta[f"minimap2_{tag}"] = run_minimap2(a.image, gen, gfa, pf, extra, a.threads, a.memory); pafs.append(str(pf))
        order = {n: i for i, n in enumerate(read_fasta(P(a.refdir, stem + ".fasta")), 1)} if stem and os.path.exists(P(a.refdir, (stem or "") + ".fasta")) else {}
        gl = []
        for pf in pafs:
            with open(pf) as fh:
                for line in fh:
                    f = line.split("\t", 12)
                    qlen_g, qs, qe, tname, tlen, ts, te, mt, bl = int(f[1]), int(f[2]), int(f[3]), f[5], int(f[6]), int(f[7]), int(f[8]), int(f[9]), int(f[10])
                    if bl < 200:
                        continue
                    ci = order.get(tname, None)
                    if ci is None and tname.rpartition("__")[2].isdigit():
                        ci = int(tname.rpartition("__")[2])
                    gl.append(dict(gene=f[0], idn=mt / bl, cov=(qe - qs) / qlen_g, contig=ci, ts=ts, te=te, bl=bl))
        gl.sort(key=lambda r: (-r["cov"], -r["idn"]))
        L.append(f"  tet(Q) gene sequences ({', '.join(genes)}) aligned to the {sp} genome: {len(gl)} alignments of block >= 200")
        for g in gl[:6]:
            L.append(f"    {g['gene']}: identity {g['idn']:.3f}, gene coverage {g['cov']:.2f}, contig {g['contig']} {g['ts']}-{g['te']}")
        # decision
        gen_iv = [(p["contig"], p["ts"], p["te"]) for p in pieces if overlap((p["qs"], p["qe"]), tq) > 0]
        v1 = cover >= 0.80
        v2 = any(g["idn"] >= 0.70 and g["cov"] >= 0.80 and any(g["contig"] == c and g["ts"] <= e + 1000 and g["te"] >= s - 1000 for c, s, e in gen_iv) for g in gl)
        ha = v1 and v2
        hb = False; gap_txt = "no flanking pair"
        if cover <= 0.20:
            for p1 in pieces:
                for p2 in pieces:
                    if p1 is p2 or p1["contig"] != p2["contig"] or p1["strand"] != p2["strand"] or p1["bl"] < 500 or p2["bl"] < 500:
                        continue
                    if p1["qe"] <= tqs + 200 and p2["qs"] >= tqe - 200 and p1["qs"] < p2["qs"]:
                        gap = (p2["ts"] - p1["te"]) if p1["strand"] == "+" else (p1["ts"] - p2["te"])
                        if -1000 <= gap <= 20000:
                            hb = True; gap_txt = f"flanking pieces read {p1['qs']}-{p1['qe']} and {p2['qs']}-{p2['qe']}; read gap {p2['qs'] - p1['qe']} bases against genome gap {gap} bases"
        verdict = "H_a" if ha else "H_b" if hb else "H_c"
        L += ["", f"  H_a (gene present in the reference, catalogue missed it): {'MET' if ha else 'NOT MET'}  (interval covered {cover:.2f}; gene-versus-genome alignment overlapping the read's genome interval with identity >= 0.70 and coverage >= 0.80: {'yes' if v2 else 'no'})",
              f"  H_b (insertion in the sampled strain, absent from the reference): {'MET' if hb else 'NOT MET'}  ({gap_txt})",
              f"  verdict: {verdict}; prediction was H_b: {'MET' if verdict == 'H_b' else 'NOT MET'}", ""]
        outrows.append([rid, sp, allele, f"{idn:.3f}", f"{cover:.2f}", verdict])
    L.append("Exploratory, one read, mock community. It decides no count; it says where one read's tet(Q) sequence probably comes from. Research use only, not clinically validated.")
    txt = "\n".join(L)
    os.makedirs(P(a.out), exist_ok=True)
    open(P(a.out, "veillonella_tetq_check.txt"), "w").write(txt + "\n")
    with open(P(a.out, "veillonella_tetq_check.tsv"), "w") as o:
        w = csv.writer(o, delimiter="\t"); w.writerow(["read", "true_species", "tetQ_allele", "identity", "interval_covered_by_genome", "verdict"]); w.writerows(outrows)
    json.dump(meta, open(P(a.out, "veillonella_tetq_check.json"), "w"), indent=1)
    print(txt)


if __name__ == "__main__":
    main()
