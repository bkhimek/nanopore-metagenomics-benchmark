#!/usr/bin/env python3
"""28_amr_kmer_position.py - AMR species attribution, ADDENDUM A9: where along a tet(Q) read does the species signal go missing?
Dev pool by default; the evaluation pool needs the guards below. Follow-up to A8 and to the exploratory bin/27 (which pointed to short or gene-heavy reads).

Idea: Kraken2 writes, for every read, the taxon of each consecutive k-mer window ('taxid:count' pieces, in read order; 0 = no hit, A = ambiguous base).
The gene sits at a known place on the read (A7 alignment to ResFinder, query start and end). So the windows can be split into: INSIDE the gene, NEAR the gene (within 2,000
windows of either gene edge), FAR (everything else). For each region S = share of windows whose taxon is the TRUE species (known in this mock; A8 labels) or a descendant.
Windows: window i starts at read position i, so a window is inside the gene if gene_start <= i <= gene_end - k (k = 35, Kraken2 default; checked per read, see V1).
Reads: the labelled supporting reads of cepA and tet(Q) in reports/amr/<pool>_attribution_reads_A8.tsv (A7 supporting reads with a true-organism label; tet(W) is reported, not tested).
A region enters a statistic only if it has >= 200 windows in that read.
VALIDITY V1: for each read the sum of window counts must equal read length - k + 1; if fewer than 95% of the reads used pass, the run is reported as invalid and no hypothesis is scored.
Hypotheses and predictions for DEV (H11 to H14) and EVALUATION (same four, scored identically), fixed in Addendum A9 BEFORE any run:
  H11 tet(Q): median over reads of (S_out - S_in) >= 0.10, S_out = windows outside the gene (near + far). The gene region supports the species less than the rest of the read.
  H12 the median (S_out - S_in) of tet(Q) reads exceeds that of cepA reads by >= 0.10. The effect is specific to tet(Q), not a feature of every gene.
  H13 the median S_out of tet(Q) reads is >= 0.10 BELOW that of cepA reads. The sequence around tet(Q) also carries less species signal (mobile element extending beyond the gene).
  H14 tet(Q) reads with a far region: median S_far >= median S_near + 0.10. The signal comes back with distance from the gene (an element boundary).
Descriptive, no prediction: what the in-gene windows map to (true species, ancestor of the species, other taxon, no hit) for each family and carrier.
Evaluation pool: needs --allow-eval, --suffix, git tag amr-a9-registered, the committed dev outputs of the same suffix, and refuses to overwrite eval_*<suffix>.* files.
Standard library only.
"""
import argparse, csv, gzip, importlib.util, json, os, statistics, subprocess, sys, time
from collections import Counter, defaultdict

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
P = lambda *x: os.path.join(ROOT, *x)
K, NEAR, MINWIN, VALID = 35, 2000, 200, 0.95
GAP = 0.10


def med(v):
    return statistics.median(v) if v else float("nan")


def windows(mapping, n):
    """'taxid:count ...' -> list of taxids (0 = no hit or ambiguous) of length sum(counts)"""
    out = []
    for tok in mapping.split():
        t, c = tok.split(":")
        if t in ("|", "|:|"):
            continue
        out.extend([0 if t in ("A", "0") else int(t)] * int(c))
    return out


def region_stats(win, gs, ge, k, cls):
    """share of windows per class for inside / near / far; cls(taxid) -> 'sp','anc','other','none'"""
    reg = {"in": [], "near": [], "far": []}
    last_in = ge - k
    for i, t in enumerate(win):
        if gs <= i <= last_in:
            r = "in"
        else:
            d = gs - i if i < gs else i - last_in
            r = "near" if d <= NEAR else "far"
        reg[r].append(cls(t))
    res = {}
    for r, v in reg.items():
        n = len(v); c = Counter(v)
        res[r] = dict(n=n, **{x: (c[x] / n if n else float("nan")) for x in ("sp", "anc", "other", "none")})
    return res


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pool", choices=["dev", "eval"], default="dev")
    ap.add_argument("--allow-eval", action="store_true")
    ap.add_argument("--suffix", default="_A9")
    ap.add_argument("--merge", default="reference/amr_family_merge.tsv")
    ap.add_argument("--reads-tsv", default="")
    ap.add_argument("--paf", default="")
    ap.add_argument("--kraken-run", default="")
    ap.add_argument("--db", default="data/db/pluspf_08_20260626")
    ap.add_argument("--out", default="reports/amr")
    a = ap.parse_args()
    a.reads_tsv = a.reads_tsv or f"reports/amr/{a.pool}_attribution_reads_A8.tsv"
    a.paf = a.paf or f"data/amr/{a.pool}/{a.pool}_vs_resfinder_c.paf"
    a.kraken_run = a.kraken_run or f"data/classify/{a.pool}/baseline_c0"
    if a.pool == "eval":
        if not a.allow_eval:
            sys.exit("evaluation pool: needs --allow-eval")
        if subprocess.run(["git", "-C", ROOT, "rev-parse", "-q", "--verify", "refs/tags/amr-a9-registered"], capture_output=True).returncode != 0:
            sys.exit("evaluation pool refused: git tag amr-a9-registered not found (commit and tag Addendum A9 first)")
        if not os.path.exists(P(a.out, f"dev_kmer_position{a.suffix}.txt")):
            sys.exit(f"evaluation pool refused: the dev run ({a.out}/dev_kmer_position{a.suffix}.txt) is not there yet")
        if any(f.startswith("eval_kmer_position") and a.suffix in f for f in os.listdir(P(a.out))):
            sys.exit("evaluation outputs are never overwritten (an eval_kmer_position*" + a.suffix + ".* file exists)")
    for need in (a.reads_tsv, a.paf, os.path.join(a.kraken_run, "kraken_output.txt.gz")):
        if not os.path.exists(P(need)):
            sys.exit(f"not found: {need}")
    sp08 = importlib.util.spec_from_file_location("x08", P("bin", "08_trace_reads.py")); x08 = importlib.util.module_from_spec(sp08); sp08.loader.exec_module(x08)
    sp21 = importlib.util.spec_from_file_location("x21", P("bin", "21_amr_reads.py")); x21 = importlib.util.module_from_spec(sp21); sp21.loader.exec_module(x21)
    for r in csv.DictReader(open(P(a.merge)), delimiter="\t"):
        x21.MERGE[r["alias"]] = r["family"]
    rows = [r for r in csv.DictReader(open(P(a.reads_tsv)), delimiter="\t") if r["destination"] != "no_truth" and r["family"] in ("cepA", "tet(Q)", "tet(W)")]
    want = {r["read"]: r for r in rows}
    # gene position on the read (A7 alignment, best alignment of the read's family)
    best = {}
    with open(P(a.paf)) as fh:
        for line in fh:
            f = line.split("\t", 12)
            q = f[0]
            if q not in want:
                continue
            qlen, qs, qe, tname, tlen, ts, te, mt, bl = int(f[1]), int(f[2]), int(f[3]), f[5], int(f[6]), int(f[7]), int(f[8]), int(f[9]), int(f[10])
            if bl <= 0 or mt / bl < x21.MIN_ID or (te - ts) / tlen < x21.MIN_SPAN:
                continue
            fam = x21.family(tname)
            d = best.setdefault(q, {})
            if fam not in d or mt > d[fam][0]:
                d[fam] = (mt, qlen, qs, qe)
    # Kraken2 window strings of these reads
    maps = {}
    with gzip.open(P(a.kraken_run, "kraken_output.txt.gz"), "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 5 and f[1] in want:
                maps[f[1]] = (int(f[3].split("|")[0]), f[4])
    parent, ranks = x08.load_nodes(x08.find_dmp(P(a.db), "nodes.dmp"))
    names, rev = x08.load_names(x08.find_dmp(P(a.db), "names.dmp"))
    lineage = x08.make_lineage(parent)
    infos = {}
    out, bad, missing = [], 0, 0
    for q, r in sorted(want.items()):
        if q not in best or r["family"] not in best[q] or q not in maps:
            missing += 1; continue
        _, qlen, qs, qe = best[q][r["family"]]
        klen, mp = maps[q]
        win = windows(mp, klen)
        if abs(len(win) - (qlen - K + 1)) > 0:
            bad += 1; continue
        sp = r["true_species"]
        if sp not in infos:
            infos[sp] = x08.SpInfo(sp, rev[sp], lineage, ranks, rev)
        i = infos[sp]

        def cls(t, i=i):
            if t == 0:
                return "none"
            lin = lineage(t)
            if i.t in lin or any(x in lin for x in i.alias):
                return "sp"
            return "anc" if t in i.idx else "other"
        st = region_stats(win, qs, qe, K, cls)
        out.append(dict(read=q, family=r["family"], true_species=sp, destination=r["destination"], qlen=qlen, gs=qs, ge=qe, st=st))
    n_ok = len(out); frac_ok = n_ok / (n_ok + bad) if (n_ok + bad) else 0.0
    valid = frac_ok >= VALID
    L = [f"AMR k-mer position check, {a.pool.upper()} pool (Addendum A9): where along the read is the species signal?", "",
         f"reads used {n_ok}; k-mer count inconsistent with read length - {K} + 1: {bad}; no alignment or no Kraken2 line: {missing}; V1 (>= {VALID:.0%} consistent): {'PASSED' if valid else 'FAILED'} ({frac_ok:.3f})", ""]
    ok = lambda v, r: v["st"][r]["n"] >= MINWIN
    tab = defaultdict(list)
    for v in out:
        tab[v["family"]].append(v)
    L.append(f"{'family':8s}{'reads':>6s}  " + "  ".join(f"S_{r} (n)" .rjust(14) for r in ("in", "near", "far")) + "    in-gene windows: species / ancestor / other / none")
    summ = {}
    for fam in ("cepA", "tet(Q)", "tet(W)"):
        vs = tab.get(fam, [])
        cells = []
        for r in ("in", "near", "far"):
            xs = [v["st"][r]["sp"] for v in vs if ok(v, r)]
            cells.append(f"{med(xs):.3f} ({len(xs)})".rjust(14))
        tot = {x: sum(v["st"]["in"][x] * v["st"]["in"]["n"] for v in vs if v["st"]["in"]["n"]) for x in ("sp", "anc", "other", "none")}
        tn = sum(tot.values()) or float("nan")
        L.append(f"{fam:8s}{len(vs):>6d}  " + "  ".join(cells) + f"    {tot['sp'] / tn:.2f} / {tot['anc'] / tn:.2f} / {tot['other'] / tn:.2f} / {tot['none'] / tn:.2f}")
        summ[fam] = vs
    def sout_safe(v):
        parts = [(v["st"][r]["sp"], v["st"][r]["n"]) for r in ("near", "far") if v["st"][r]["n"]]
        n = sum(p[1] for p in parts)
        return sum(p[0] * p[1] for p in parts) / n if n else float("nan")
    def outn(v):
        return v["st"]["near"]["n"] + v["st"]["far"]["n"]
    d_tq = [sout_safe(v) - v["st"]["in"]["sp"] for v in summ.get("tet(Q)", []) if ok(v, "in") and outn(v) >= MINWIN]
    d_ce = [sout_safe(v) - v["st"]["in"]["sp"] for v in summ.get("cepA", []) if ok(v, "in") and outn(v) >= MINWIN]
    so_tq = [sout_safe(v) for v in summ.get("tet(Q)", []) if outn(v) >= MINWIN]
    so_ce = [sout_safe(v) for v in summ.get("cepA", []) if outn(v) >= MINWIN]
    nf = [v for v in summ.get("tet(Q)", []) if ok(v, "far") and ok(v, "near")]
    h11 = med(d_tq); h12 = med(d_tq) - med(d_ce); h13 = med(so_ce) - med(so_tq)
    h14 = med([v["st"]["far"]["sp"] for v in nf]) - med([v["st"]["near"]["sp"] for v in nf])
    res = [("H11 tet(Q): median (S_out - S_in) >= 0.10", valid and h11 >= GAP, f"{h11:+.3f}, n {len(d_tq)}"),
           ("H12 tet(Q) minus cepA median (S_out - S_in) >= 0.10", valid and h12 >= GAP, f"{h12:+.3f} (tet(Q) n {len(d_tq)}, cepA n {len(d_ce)})"),
           ("H13 tet(Q) median S_out at least 0.10 below cepA", valid and h13 >= GAP, f"cepA {med(so_ce):.3f}, tet(Q) {med(so_tq):.3f}, gap {h13:+.3f}"),
           ("H14 tet(Q) median S_far >= S_near + 0.10", valid and h14 >= GAP, f"{h14:+.3f}, n {len(nf)}")]
    L += ["", "Registered predictions (A9)" + ("" if valid else " - NOT SCORED, V1 failed") + ":"]
    for nme, okk, val in res:
        L.append(f"  {nme}: {('MET' if okk else 'NOT MET') if valid else 'not scored'}  ({val})")
    L += ["", "Mock community only (true species known from the supplier genomes). S is the share of windows whose Kraken2 taxon is the true species or below; windows with no hit count against it. "
          "Small numbers, one pool per run. Research use only, not clinically validated."]
    txt = "\n".join(L)
    os.makedirs(P(a.out), exist_ok=True)
    tag = f"{a.pool}_kmer_position{a.suffix}"
    open(P(a.out, tag + ".txt"), "w").write(txt + "\n")
    with open(P(a.out, tag + ".tsv"), "w") as o:
        w = csv.writer(o, delimiter="\t")
        w.writerow(["read", "family", "true_species", "destination", "read_len", "gene_start", "gene_end"] + [f"{r}_{x}" for r in ("in", "near", "far") for x in ("n", "sp", "anc", "other", "none")])
        for v in out:
            w.writerow([v["read"], v["family"], v["true_species"], v["destination"], v["qlen"], v["gs"], v["ge"]] +
                       [(f"{v['st'][r][x]:.3f}" if x != "n" else v["st"][r][x]) for r in ("in", "near", "far") for x in ("n", "sp", "anc", "other", "none")])
    json.dump({"script": "bin/28_amr_kmer_position.py", "args": vars(a), "k": K, "reads_used": n_ok, "inconsistent": bad, "valid": valid,
               "run": time.strftime("%Y-%m-%d %H:%M:%S")}, open(P(a.out, tag + ".json"), "w"), indent=1)
    print(txt)


if __name__ == "__main__":
    main()
