#!/usr/bin/env python3
"""32_veillonella_junction.py - EXPLORATORY follow-up to bin/31. Part of no verdict, changes no registered rule, tunes nothing, changes no count.
Both pools (dev and eval) are read; the eval pool was already scored under A7/A8, and the one read of bin/31 was examined after scoring; this look also comes after scoring.
Nothing from the evaluation pool is written to an existing evaluation output; all outputs have new names (reports/amr/veillonella_junction_check.*, data/amr/veillonella_junction/).

Background (bin/31): read SRR31222153.7328964 carries tet(Q) (98.4% identity to a ResFinder allele) but its true organism, Veillonella rogosae, has no tet(Q) in the supplier genome.
The read is two pieces joined with a 1-base gap: read 34-12,128 = Prevotella corporis (genome 59,381-71,458, identity 0.879, the tet(Q) locus region), read 12,130-25,325 = Veillonella contig 1,
genome 1,599,996-1,613,203 (+). So the junction sits at Veillonella contig 1, position J = 1,599,996, and the extra DNA is on the genome-left side of J.
Two explanations were left open: (i) the sampled strain carries an element inserted at J that the supplier genome lacks, (ii) the read is a chimera joined in library preparation.
Question: do OTHER reads of both pools that cover J show the same break?

Method. 1. Map both pools to the Veillonella rogosae genome alone (minimap2 -x map-ont -c --secondary=no, supplementary pieces kept).
2. For each read with an alignment touching J (target interval within +/-300 of J at one end, or covering J), classify by the alignments of that read on the site contig:
   span          the union of the read's alignments (gaps <= 100 merged) covers J-500..J+500 (the read continues through J as the reference does);
   split_right   an alignment starts at J +/- 300 (genome side right of J) AND the read has >= 1,000 bases beyond it on the genome-left side not aligned there;
   split_left    an alignment ends at J +/- 300 AND >= 1,000 read bases beyond it on the genome-right side;
   short_end     an alignment ends at J +/- 300 with < 1,000 read bases beyond (the read stops, or the break is too short to judge): reported, not used.
   Records used: identity >= 0.80 and block >= 300 on the site contig.
3. The read bases beyond the break (the 'tail') of every split read are mapped to all 21 supplier genomes (-x map-ont -c -p 0 -N 50) and to ResFinder.
   Tail class: carrier (best block >= 1,000 and identity >= 0.80 on B. fragilis or P. corporis), other_mock (another species), veillonella_elsewhere (the Veillonella genome >= 20 kb from J or another contig), none.
   tail_tetQ: a ResFinder tet(Q) alignment with identity >= 0.80 and gene span >= 0.80 (the A7 filter). A tail is ELEMENT-LIKE if its class is carrier or it has tet(Q).
Positive control (validity V1), decided before the run: the read SRR31222153.7328964 must come out as split_right with a carrier tail; if not, the run is VOID and reported as such with no outcome.
Registered hypotheses, written BEFORE the first run; counts exclude the control read:
  K1 (the control read is a chimera, no insertion at J in the sampled DNA): no other read is split at J AND at least 10 other reads span J.
  K2 (insertion in the sampled DNA): at least 2 other reads are split at J AND at least one of them has an element-like tail.
  OPEN (anything else, including fewer than 10 spanning reads with no split read): reported with the numbers.
PREDICTION: K1. (Basis: the supplier assembled its genome from the mock's own strain and the reference has no insertion at J; Veillonella rogosae is an abundant member of the mock, so
many reads should cover J; ONT library preparation makes a few percent chimeric reads; one chimera among a hundred or more covering reads fits that.)
A prediction that fails stays failed. Descriptive only: depth at J, the split share, tail classes. Standard library only; minimap2 in Docker as in bin/31 unless existing PAFs are given (tests).
"""
import argparse, csv, gzip, importlib.util, json, os, subprocess, sys, time
from collections import defaultdict
from pathlib import Path

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
P = lambda *x: os.path.join(ROOT, *x)
CARRIERS = ("Bacteroides fragilis", "Prevotella corporis")
CONTROL = "SRR31222153.7328964"
TOL, MINTAIL, SPANPAD, MERGE_GAP = 300, 1000, 500, 100
MIN_ID, MIN_BL = 0.80, 300
TAIL_BL, TAIL_ID, FAR = 1000, 0.80, 20000


def run_mm2(image, ref, query, paf, extra, threads, mem):
    ref, query, paf = Path(ref).resolve(), Path(query).resolve(), Path(paf).resolve()
    if subprocess.run(["docker", "image", "inspect", image], capture_output=True).returncode != 0:
        subprocess.run(["docker", "pull", image], check=True)
    cmd = ["docker", "run", "--rm", "-u", f"{os.getuid()}:{os.getgid()}", "-m", mem, "--cpus", str(threads), "-v", f"{ref.parent}:/ref:ro", "-v", f"{query.parent}:/q:ro",
           "-v", f"{paf.parent}:/out", image, "minimap2"] + extra + ["-c", "-t", str(threads), "-o", f"/out/{paf.name}", f"/ref/{ref.name}", f"/q/{query.name}"]
    subprocess.run(cmd, check=True)
    return " ".join(cmd)


def read_fasta_names(path):
    return [l[1:].split()[0] for l in open(path) if l.startswith(">")]


def union_cover(ivs, lo, hi):
    """True if the intervals (merged across gaps <= MERGE_GAP) cover lo..hi"""
    segs = []
    for s, e in sorted(ivs):
        if segs and s <= segs[-1][1] + MERGE_GAP:
            segs[-1][1] = max(segs[-1][1], e)
        else:
            segs.append([s, e])
    return any(s <= lo and e >= hi for s, e in segs)


def classify_reads(paf, contig_names, site_contig, J):
    """read -> dict(cls, rec, tail=(qstart, qend) or None, ...)"""
    recs = defaultdict(list)
    with open(paf) as fh:
        for line in fh:
            f = line.split("\t", 12)
            qlen, qs, qe, strand, tname, tlen, ts, te, mt, bl = int(f[1]), int(f[2]), int(f[3]), f[4], f[5], int(f[6]), int(f[7]), int(f[8]), int(f[9]), int(f[10])
            if bl < MIN_BL or mt / bl < MIN_ID or contig_names.get(tname) != site_contig:
                continue
            if te < J - 5000 or ts > J + 5000:
                continue
            recs[f[0]].append(dict(qlen=qlen, qs=qs, qe=qe, strand=strand, ts=ts, te=te, idn=mt / bl, bl=bl))
    out = {}
    for rid, rs in recs.items():
        qlen = rs[0]["qlen"]
        if union_cover([(r["ts"], r["te"]) for r in rs], J - SPANPAD, J + SPANPAD):
            out[rid] = dict(cls="span", qlen=qlen, ts="", te="", strand="", tail=None, tail_len=0)
            continue
        got = None
        for r in rs:
            left_un = r["qs"] if r["strand"] == "+" else qlen - r["qe"]     # unaligned read bases on the genome-left side of this alignment
            right_un = qlen - r["qe"] if r["strand"] == "+" else r["qs"]    # ... on the genome-right side
            if abs(r["ts"] - J) <= TOL and r["te"] >= J + SPANPAD:
                # alignment starts at J, runs to the right; the break faces left
                tail = (0, r["qs"]) if r["strand"] == "+" else (r["qe"], qlen)
                got = dict(cls="split_right" if left_un >= MINTAIL else "short_end", tail=tail if left_un >= MINTAIL else None, tail_len=left_un)
            elif abs(r["te"] - J) <= TOL and r["ts"] <= J - SPANPAD:
                tail = (r["qe"], qlen) if r["strand"] == "+" else (0, r["qs"])
                got = dict(cls="split_left" if right_un >= MINTAIL else "short_end", tail=tail if right_un >= MINTAIL else None, tail_len=right_un)
            if got:
                got.update(qlen=qlen, ts=r["ts"], te=r["te"], strand=r["strand"], idn=r["idn"])
                if got["cls"] != "short_end":
                    break
        if got:
            out[rid] = got
    return out


def fastq_get(path, want):
    """id -> sequence for the wanted ids (streams the file once)"""
    op = gzip.open if str(path).endswith(".gz") else open
    seqs = {}
    if not want:
        return seqs
    with op(path, "rt") as fh:
        while True:
            h = fh.readline()
            if not h:
                break
            s = fh.readline().strip(); fh.readline(); fh.readline()
            rid = h[1:].split()[0]
            if rid in want:
                seqs[rid] = s
                if len(seqs) == len(want):
                    break
    return seqs


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pools", nargs="+", default=["dev", "eval"])
    ap.add_argument("--reads", nargs="*", default=[], help="fastq per pool in the same order; default data/subsets/SRR31222153.<pool>.seed20261012.blocks.fastq.gz")
    ap.add_argument("--species", default="Veillonella_rogosae")
    ap.add_argument("--site-contig", type=int, default=1)
    ap.add_argument("--site", type=int, default=1599996)
    ap.add_argument("--refdir", default="data/reference/D6331.refseq/genomes")
    ap.add_argument("--resfinder", default="data/db/resfinder_db/all.fsa")
    ap.add_argument("--raw", default="reports/amr/reference_hits_raw.tsv")
    ap.add_argument("--merge", default="reference/amr_family_merge.tsv")
    ap.add_argument("--workdir", default="data/amr/veillonella_junction")
    ap.add_argument("--out", default="reports/amr")
    ap.add_argument("--site-paf", nargs="*", default=[], help="existing PAFs of the pools against the genome (tests), one per pool")
    ap.add_argument("--tail-paf", default="", help="existing PAF of the tails against the 21 genomes (tests)")
    ap.add_argument("--tail-rf-paf", default="", help="existing PAF of the tails against ResFinder (tests)")
    ap.add_argument("--image", default="quay.io/biocontainers/minimap2:2.31--h118bc1c_0")
    ap.add_argument("--threads", type=int, default=6)
    ap.add_argument("--memory", default="8g")
    a = ap.parse_args()
    reads = a.reads or [f"data/subsets/SRR31222153.{p}.seed20261012.blocks.fastq.gz" for p in a.pools]
    s8 = importlib.util.spec_from_file_location("x08", P("bin", "08_trace_reads.py")); x08 = importlib.util.module_from_spec(s8); s8.loader.exec_module(x08)
    s21 = importlib.util.spec_from_file_location("x21", P("bin", "21_amr_reads.py")); x21 = importlib.util.module_from_spec(s21); s21.loader.exec_module(x21)
    for r in csv.DictReader(open(P(a.merge)), delimiter="\t"):
        x21.MERGE[r["alias"]] = r["family"]
    work = Path(P(a.workdir)); work.mkdir(parents=True, exist_ok=True)
    meta = {"script": "bin/32_veillonella_junction.py", "args": vars(a), "started": time.strftime("%Y-%m-%d %H:%M:%S")}
    gen = P(a.refdir, a.species + ".fasta")
    names = read_fasta_names(gen)
    cn = {n: i for i, n in enumerate(names, 1)}
    # 1. pools against the Veillonella genome; classify
    calls = {}
    for i, pool in enumerate(a.pools):
        pf = Path(a.site_paf[i]) if a.site_paf else work / f"{pool}_vs_{a.species}.paf"
        if not a.site_paf:
            meta[f"minimap2_site_{pool}"] = run_mm2(a.image, gen, P(reads[i]), pf, ["-x", "map-ont", "--secondary=no"], a.threads, a.memory)
        for rid, d in classify_reads(str(pf), cn, a.site_contig, a.site).items():
            d["pool"] = pool; calls[rid] = d
    # 2. tails of split reads
    split = {r: d for r, d in calls.items() if d["cls"] in ("split_right", "split_left")}
    tails = work / "tails.fa"
    with open(tails, "w") as o:
        for i, pool in enumerate(a.pools):
            want = {r for r, d in split.items() if d["pool"] == pool}
            for rid, s in fastq_get(P(reads[i]), want).items():
                qs, qe = split[rid]["tail"]
                o.write(f">{rid}\n{s[qs:qe]}\n")
    tc = defaultdict(list); tq = defaultdict(bool)
    if split:
        if a.tail_paf:
            tpaf = a.tail_paf
        else:
            fa = work / "ref_combined.fa"; meta["reference_md5"] = x08.build_reference(a.refdir, fa)
            tpaf = work / "tails_vs_refs.paf"
            meta["minimap2_tails"] = run_mm2(a.image, fa, tails, tpaf, ["-x", "map-ont", "--secondary=yes", "-p", "0", "-N", "50"], a.threads, a.memory)
        with open(tpaf) as fh:
            for line in fh:
                f = line.split("\t", 12)
                stem, _, n = f[5].rpartition("__")
                mt, bl = int(f[9]), int(f[10])
                if bl >= TAIL_BL and mt / bl >= TAIL_ID:
                    tc[f[0]].append(dict(stem=stem, contig=int(n) if n.isdigit() else 0, ts=int(f[7]), te=int(f[8]), mt=mt, bl=bl, idn=mt / bl, species=x08.stem_species(stem)))
        if a.tail_rf_paf:
            rpaf = a.tail_rf_paf
        else:
            rpaf = work / "tails_vs_resfinder.paf"
            meta["minimap2_tails_resfinder"] = run_mm2(a.image, P(a.resfinder), tails, rpaf, ["-x", "map-ont", "--secondary=yes", "-p", "0", "-N", "50"], a.threads, a.memory)
        with open(rpaf) as fh:
            for line in fh:
                f = line.split("\t", 12)
                tlen, ts, te, mt, bl = int(f[6]), int(f[7]), int(f[8]), int(f[9]), int(f[10])
                if bl > 0 and mt / bl >= x21.MIN_ID and (te - ts) / tlen >= x21.MIN_SPAN and x21.family(f[5]) == "tet(Q)":
                    tq[f[0]] = True
    ref_tq = []
    if os.path.exists(P(a.raw)):
        ref_tq = [(h["genome"], int(h["Start"]), int(h["Stop"])) for h in csv.DictReader(open(P(a.raw)), delimiter="\t") if h["Element symbol"] == "tet(Q)"]
    sp_name = x08.stem_species(a.species)
    rows = []
    for rid, d in sorted(calls.items()):
        t = dict(tclass="", tspecies="", tblock="", tid="", tlocus="", tnear="")
        if rid in split:
            c = sorted(tc.get(rid, []), key=lambda r: -r["mt"])
            if not c:
                t["tclass"] = "none"
            else:
                b = c[0]
                far_v = b["species"] == sp_name and (b["contig"] != a.site_contig or abs(b["ts"] - a.site) > FAR)
                t["tclass"] = "carrier" if b["species"] in CARRIERS else "veillonella_elsewhere" if far_v else "other_mock" if b["species"] != sp_name else "none"
                t.update(tspecies=b["species"], tblock=b["bl"], tid=f"{b['idn']:.3f}", tlocus=f"{b['stem']}:{b['ts']}-{b['te']}")
                t["tnear"] = int(any(g == b["stem"] and b["te"] >= s - 30000 and b["ts"] <= e + 30000 for g, s, e in ref_tq))
        rows.append(dict(pool=d["pool"], read=rid, cls=d["cls"], qlen=d["qlen"], ts=d["ts"], te=d["te"], strand=d["strand"], tail_len=d.get("tail_len", 0), tail_tetQ=int(tq.get(rid, False)) if rid in split else "", **t))
    # decisions
    ctl = next((r for r in rows if r["read"] == CONTROL), None)
    v1 = bool(ctl and ctl["cls"] == "split_right" and ctl["tclass"] == "carrier")
    oth = [r for r in rows if r["read"] != CONTROL]
    nspan = sum(1 for r in oth if r["cls"] == "span"); nsplit = [r for r in oth if r["cls"] in ("split_right", "split_left")]
    nshort = sum(1 for r in oth if r["cls"] == "short_end")
    elem = [r for r in nsplit if r["tclass"] == "carrier" or r["tail_tetQ"] == 1]
    k1 = v1 and len(nsplit) == 0 and nspan >= 10
    k2 = v1 and len(nsplit) >= 2 and len(elem) >= 1
    verdict = "VOID" if not v1 else "K1" if k1 else "K2" if k2 else "OPEN"
    L = [f"EXPLORATORY: do other reads break at {a.species} contig {a.site_contig} position {a.site}? (both pools; after scoring; no count changed)", "",
         f"reads touching the site: {len(rows)}; spanning {nspan}; split {len(nsplit)} (right {sum(1 for r in nsplit if r['cls'] == 'split_right')}, left {sum(1 for r in nsplit if r['cls'] == 'split_left')}); short end {nshort}; control read counted separately",
         ""]
    if ctl:
        L.append(f"control read {CONTROL}: {ctl['cls']}, tail {ctl['tail_len']} bases, tail class {ctl['tclass']} ({ctl['tspecies']} {ctl['tlocus']}, identity {ctl['tid']}), tet(Q) in tail: {ctl['tail_tetQ']}")
    else:
        L.append(f"control read {CONTROL}: NOT FOUND among the reads touching the site")
    L.append(f"V1 (control read is split at the site with a carrier tail): {'MET' if v1 else 'NOT MET -> the run is VOID, no outcome is reported'}")
    if nsplit:
        L += ["", "other split reads:"]
        for r in nsplit:
            L.append(f"  {r['pool']} {r['read']} {r['cls']} length {r['qlen']} tail {r['tail_len']} -> {r['tclass']} {r['tspecies']} {r['tlocus']} id {r['tid']}; tet(Q) in tail {r['tail_tetQ']}; near reference tet(Q) {r['tnear']}")
    share = len(nsplit) / (len(nsplit) + nspan) if (len(nsplit) + nspan) else float("nan")
    L += ["", f"K1 (control read is a chimera; no other split read, >= 10 spanning): {'MET' if k1 else 'NOT MET'}  (split {len(nsplit)}, spanning {nspan})",
          f"K2 (insertion in the sampled DNA; >= 2 other split reads, >= 1 element-like tail): {'MET' if k2 else 'NOT MET'}  (split {len(nsplit)}, element-like {len(elem)})",
          f"outcome: {verdict}; prediction was K1: {'MET' if verdict == 'K1' else 'NOT MET' if verdict != 'VOID' else 'not decided (void)'}",
          f"descriptive: split share among covering reads (split / (split + spanning)) = {share:.3f}", "",
          "Exploratory, mock community, one site. It decides no count. Research use only, not clinically validated."]
    txt = "\n".join(L)
    os.makedirs(P(a.out), exist_ok=True)
    open(P(a.out, "veillonella_junction_check.txt"), "w").write(txt + "\n")
    cols = ["pool", "read", "cls", "qlen", "ts", "te", "strand", "tail_len", "tclass", "tspecies", "tblock", "tid", "tlocus", "tnear", "tail_tetQ"]
    with open(P(a.out, "veillonella_junction_check.tsv"), "w") as o:
        w = csv.writer(o, delimiter="\t"); w.writerow(cols); [w.writerow([r[c] for c in cols]) for r in rows]
    meta.update(outcome=verdict, n_span=nspan, n_split=len(nsplit), n_short=nshort)
    json.dump(meta, open(P(a.out, "veillonella_junction_check.json"), "w"), indent=1)
    print(txt)


if __name__ == "__main__":
    main()
