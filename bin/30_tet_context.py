#!/usr/bin/env python3
"""30_tet_context.py - EXPLORATORY, reference genomes only (no read pool is touched); part of no verdict; changes no registered rule.

Question (follow-up of A8 to A10): the species-specific k-mer signal is missing inside tet(Q) and tet(W) and, for tet(Q), up to at least 2,000 windows on both sides, whereas it is
present in and around cepA. Is that because tet(Q) sits inside a stretch of DNA that other species of the community also carry (a mobile element), while cepA sits in species-specific DNA?
Method: for every AMRFinderPlus hit of tet(Q), tet(W) and cepA (bin/19 raw table; the same genes as the catalogue) cut the supplier genome +/- FLANK bases around the gene (default 30,000) and align that
region with minimap2 (-x asm20 -c, base-level matches and block length) against all 21 supplier genomes. Alignments with identity >= 0.80 (matches / block) and block >= 500 are kept; the alignment to the region's
own place in its own genome is removed. Per partner genome the kept alignments are merged (gaps up to 1,000 allowed); the merged segment that overlaps the gene gives the extent of shared sequence
beyond the gene on the left and on the right (bases). Also reported: GC of the flanks (up to 10,000 bases each side, gene excluded) against the genome's GC, and every other AMRFinderPlus hit within the flank.
Hypotheses and predictions, written BEFORE the first run:
  H15 the two tet(Q) loci (B. fragilis, P. corporis) share sequence with each other that covers the gene and extends >= 2,000 bases beyond the gene on at least one side, in BOTH directions of the comparison
      (a shared element, consistent with the k-mer result that the species signal is lost at least 2,000 windows from the gene).
  H16 the cepA locus has no partner genome whose shared sequence extends >= 500 bases beyond the gene on either side (species-specific surroundings).
  H17 the GC of the tet(Q) flanks differs from the genome GC by >= 0.03 (absolute) in both genomes (composition different from the host, as often seen in mobile elements).
Descriptive, no prediction: tet(W) (F. prausnitzii), every partner found, co-located AMR hits.
Limits: this is sequence sharing inside a 21-genome mock; it cannot name the element, show that it moves, or say how common it is in other strains. A mobile-element hallmark gene search (integrase, transfer genes) is NOT done here.
A prediction that fails stays failed. Standard library only; minimap2 (Docker image as in bin/08) unless --paf is given (tests).
"""
import argparse, csv, importlib.util, json, os, subprocess, sys, time
from collections import defaultdict
from pathlib import Path

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
P = lambda *x: os.path.join(ROOT, *x)
MIN_ID, MIN_BLOCK, GAP, GENES = 0.80, 500, 1000, ("tet(Q)", "tet(W)", "cepA")
EXT_Q, EXT_C, GC_D, GC_WIN = 2000, 500, 0.03, 10000


def read_fasta(path):
    seqs, name, buf = [], None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    seqs.append((name, "".join(buf)))
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name is not None:
        seqs.append((name, "".join(buf)))
    return seqs


def gc(s):
    s = s.upper(); n = sum(s.count(c) for c in "ACGT")
    return (s.count("G") + s.count("C")) / n if n else float("nan")


def merged_extent(ivs, gs, ge, gap=GAP):
    """ivs: list of (qs, qe) alignment intervals in region coordinates; returns (left, right, seg_start, seg_end) of the merged segment overlapping the gene [gs, ge), or None"""
    ivs = sorted(ivs); segs = []
    for s, e in ivs:
        if segs and s <= segs[-1][1] + gap:
            segs[-1][1] = max(segs[-1][1], e)
        else:
            segs.append([s, e])
    for s, e in segs:
        if s < ge and e > gs:
            return max(0, gs - s), max(0, e - ge), s, e
    return None


def run_minimap2(image, fa, query, paf, threads, mem):
    fa, query, paf = Path(fa).resolve(), Path(query).resolve(), Path(paf).resolve()
    if subprocess.run(["docker", "image", "inspect", image], capture_output=True).returncode != 0:
        subprocess.run(["docker", "pull", image], check=True)
    ver = subprocess.run(["docker", "run", "--rm", image, "minimap2", "--version"], capture_output=True, text=True).stdout.strip()
    cmd = ["docker", "run", "--rm", "-u", f"{os.getuid()}:{os.getgid()}", "-m", mem, "--cpus", str(threads),
           "-v", f"{fa.parent}:/ref:ro", "-v", f"{query.parent}:/q:ro", "-v", f"{paf.parent}:/out", image,
           "minimap2", "-x", "asm20", "-c", "--secondary=yes", "-N", "50", "-t", str(threads), "-o", f"/out/{paf.name}", f"/ref/{fa.name}", f"/q/{query.name}"]
    t0 = time.time(); subprocess.run(cmd, check=True)
    return {"minimap2_version": ver, "command": " ".join(cmd), "wall_s": round(time.time() - t0, 1)}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--raw", default="reports/amr/reference_hits_raw.tsv")
    ap.add_argument("--refdir", default="data/reference/D6331.refseq/genomes")
    ap.add_argument("--workdir", default="data/amr/context")
    ap.add_argument("--flank", type=int, default=30000)
    ap.add_argument("--out", default="reports/amr")
    ap.add_argument("--paf", default="", help="use an existing PAF of the regions against the combined genomes (tests)")
    ap.add_argument("--image", default="quay.io/biocontainers/minimap2:2.31--h118bc1c_0")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--memory", default="8g")
    a = ap.parse_args()
    s8 = importlib.util.spec_from_file_location("x08", P("bin", "08_trace_reads.py")); x08 = importlib.util.module_from_spec(s8); s8.loader.exec_module(x08)
    work = Path(P(a.workdir)); work.mkdir(parents=True, exist_ok=True)
    hits = list(csv.DictReader(open(P(a.raw)), delimiter="\t"))
    genomes = {}
    def genome(stem):
        if stem not in genomes:
            genomes[stem] = read_fasta(P(a.refdir, stem + ".fasta"))
        return genomes[stem]
    regions = []
    for h in hits:
        if h["Element symbol"] not in GENES or h["Type"] != "AMR":
            continue
        stem = h["genome"]; contigs = genome(stem)
        idx = next((i for i, (n, _) in enumerate(contigs, 1) if n == h["Contig id"]), None)
        if idx is None:
            sys.exit(f"contig {h['Contig id']} not found in {stem}.fasta")
        seq = contigs[idx - 1][1]; gs0, ge0 = int(h["Start"]) - 1, int(h["Stop"])
        rs, re_ = max(0, gs0 - a.flank), min(len(seq), ge0 + a.flank)
        regions.append(dict(id=f"{stem}|{h['Element symbol']}|{h['Start']}", genome=stem, symbol=h["Element symbol"], contig=idx, abs_s=rs, abs_e=re_,
                            gs=gs0 - rs, ge=ge0 - rs, left_avail=gs0 - rs, right_avail=re_ - ge0, seq=seq[rs:re_], contig_seq=seq,
                            genome_gc=gc("".join(s for _, s in contigs)), hit=h))
    if not regions:
        sys.exit("no tet(Q), tet(W) or cepA hit in the raw table")
    qfa = work / "regions.fa"
    with open(qfa, "w") as o:
        for r in regions:
            o.write(f">{r['id']}\n{r['seq']}\n")
    meta = {"script": "bin/30_tet_context.py", "args": vars(a), "started": time.strftime("%Y-%m-%d %H:%M:%S")}
    if a.paf:
        paf = P(a.paf)
    else:
        allfa = work / "genomes_all.fa"
        meta["reference_md5"] = x08.build_reference(P(a.refdir), str(allfa))
        paf = str(work / "regions_vs_genomes.paf")
        meta["minimap2"] = run_minimap2(a.image, allfa, qfa, paf, a.threads, a.memory)
    byid = {r["id"]: r for r in regions}
    ivs = defaultdict(lambda: defaultdict(list))  # region id -> partner stem -> [(qs, qe, id)]
    with open(paf) as fh:
        for line in fh:
            f = line.split("\t", 12)
            r = byid.get(f[0])
            if r is None:
                continue
            qs, qe, tname, ts, te, mt, bl = int(f[2]), int(f[3]), f[5], int(f[7]), int(f[8]), int(f[9]), int(f[10])
            if bl < MIN_BLOCK or mt / bl < MIN_ID:
                continue
            stem, _, n = tname.rpartition("__")
            if stem == r["genome"] and int(n) == r["contig"] and ts < r["abs_e"] and te > r["abs_s"]:
                continue  # the region's own place
            ivs[r["id"]][stem].append((qs, qe, mt / bl))
    rows, L = [], ["EXPLORATORY: genomic context of tet(Q), tet(W) and cepA in the supplier genomes (reference genomes only)", "",
                   f"flank +/- {a.flank} bases; alignments identity >= {MIN_ID}, block >= {MIN_BLOCK}; shared extent = merged segment overlapping the gene, gaps <= {GAP}", ""]
    ext = {}
    for r in regions:
        L.append(f"{r['genome']}  {r['symbol']}  {r['hit']['Contig id']}:{r['hit']['Start']}-{r['hit']['Stop']}  (flank available: left {r['left_avail']}, right {r['right_avail']})")
        fl = r["seq"][max(0, r["gs"] - GC_WIN):r["gs"]] + r["seq"][r["ge"]:r["ge"] + GC_WIN]
        d = gc(fl) - r["genome_gc"]
        L.append(f"  GC of the flanks (<= {GC_WIN} each side, gene excluded) {gc(fl):.3f}; genome {r['genome_gc']:.3f}; difference {d:+.3f}")
        best = 0
        for stem, lst in sorted(ivs[r["id"]].items()):
            m = merged_extent([(s, e) for s, e, _ in lst], r["gs"], r["ge"])
            if m is None:
                continue
            left, right, s0, e0 = m
            idn = sum(i for _, _, i in lst) / len(lst)
            ext[(r["id"], stem)] = (left, right)
            best = max(best, left, right)
            L.append(f"  shares sequence with {stem}: gene covered, extends {left} bases to the left and {right} to the right of the gene (identity of alignments about {idn:.3f})")
            rows.append([r["genome"], r["symbol"], r["hit"]["Start"], stem, left, right, f"{idn:.3f}"])
        if not any(k[0] == r["id"] for k in ext):
            L.append("  no other genome shares the gene or its flanks (>= 500 bases at identity >= 0.80)")
        r["best"] = best; r["gcd"] = d
        near = [x for x in hits if x["genome"] == r["genome"] and x["Contig id"] == r["hit"]["Contig id"] and x is not r["hit"] and x["Element symbol"] != "NA"
                and abs((int(x["Start"]) + int(x["Stop"])) // 2 - (int(r["hit"]["Start"]) + int(r["hit"]["Stop"])) // 2) <= a.flank]
        for x in near:
            L.append(f"  other annotated element within the flank: {x['Element symbol']} ({x['Type']}, {x['Class']}) at {x['Start']}-{x['Stop']}")
        L.append("")
    tq = [r for r in regions if r["symbol"] == "tet(Q)"]
    h15 = False
    if len(tq) == 2:
        a_, b_ = tq
        e1 = ext.get((a_["id"], b_["genome"]), (0, 0)); e2 = ext.get((b_["id"], a_["genome"]), (0, 0))
        h15 = max(e1) >= EXT_Q and max(e2) >= EXT_Q
        v15 = f"{a_['genome']} -> {b_['genome']}: left {e1[0]}, right {e1[1]}; reverse: left {e2[0]}, right {e2[1]}"
    else:
        v15 = f"{len(tq)} tet(Q) loci found (2 expected)"
    ce = [r for r in regions if r["symbol"] == "cepA"]
    h16 = bool(ce) and all(r["best"] <= EXT_C for r in ce)
    h17 = len(tq) == 2 and all(abs(r["gcd"]) >= GC_D for r in tq)
    L += ["Registered predictions:",
          f"  H15 tet(Q) loci share >= {EXT_Q} bases beyond the gene with each other, both directions: {'MET' if h15 else 'NOT MET'}  ({v15})",
          f"  H16 cepA: no partner shares >= {EXT_C} bases beyond the gene: {'MET' if h16 else 'NOT MET'}  (largest extent over partners: " + ", ".join(str(r['best']) for r in ce) + ")",
          f"  H17 tet(Q) flank GC differs from the genome GC by >= {GC_D}, both genomes: {'MET' if h17 else 'NOT MET'}  (" + ", ".join(f"{r['genome']} {r['gcd']:+.3f}" for r in tq) + ")",
          "", "Exploratory, reference genomes of a mock only. Sequence sharing does not show that the DNA moves or that it is a known mobile element; no hallmark gene search was made. Research use only, not clinically validated."]
    txt = "\n".join(L)
    os.makedirs(P(a.out), exist_ok=True)
    open(P(a.out, "tet_context.txt"), "w").write(txt + "\n")
    with open(P(a.out, "tet_context.tsv"), "w") as o:
        w = csv.writer(o, delimiter="\t"); w.writerow(["genome", "gene", "gene_start", "partner_genome", "shared_left_of_gene", "shared_right_of_gene", "mean_identity"]); w.writerows(rows)
    json.dump(meta, open(P(a.out, "tet_context.json"), "w"), indent=1)
    print(txt)


if __name__ == "__main__":
    main()
