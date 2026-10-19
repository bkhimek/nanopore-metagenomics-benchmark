#!/usr/bin/env python3
"""29_amr_species_windows.py - AMR species attribution, ADDENDUM A10: does a read reach its species only if it holds species-specific k-mer windows?
Dev pool by default (descriptive, informed by dev: see the addendum); the evaluation pool is the test and needs the guards below. Follows A8, bin/27 and A9 (bin/28).

Background (A9, dev, exploratory): 93 to 96% of the Kraken2 windows of a nanopore read have no hit. Among the windows that do hit, those inside the tet(Q) and tet(W)
genes all map to a shared ancestor (species share 0.0%), those inside cepA all map to the species (100%); and tet(Q) reads that stopped above species had
almost no species-specific window anywhere in the read (28 of 29 had none). Kraken2 adds the weight of a window to its whole path to the root, so a read needs only a few species-specific
windows to reach its species, and none at all means it stops at a shared ancestor.
Definitions: window = consecutive k-mer (k = 35) of the read, taxon as in the Kraken2 output (bin/28 rules, same V1 check, same gene place from the A7 alignment, same NEAR = 2,000).
A window is SPECIFIC if its taxon is the TRUE species (A8 label) or a descendant. n_sp(read) = number of specific windows in the whole read.
Pooled shares: over all reads of a family, specific windows / windows that hit (species, ancestor or other taxon), per region (inside the gene, near, far).
Reads: labelled A7 supporting reads of cepA, tet(Q), tet(W) (reports/amr/<pool>_attribution_reads_A8.tsv). destination 'species' = correct species in A8.
Predictions (A10), fixed BEFORE the evaluation pool's positional data were looked at; scored on the evaluation pool and also reported on dev, where they are informed by the data:
  M1 rule 'destination is species if and only if n_sp >= 1' is right for >= 0.90 of the reads (all three families pooled).
  M2 tet(Q): >= 0.50 of reads have n_sp = 0; cepA: <= 0.05 of reads have n_sp = 0.
  M3 inside the gene, the species share of hit windows is <= 0.02 for tet(Q) and >= 0.90 for cepA.
  M4 tet(Q), near the gene (within 2,000 windows of either edge): species share of hit windows <= 0.02.
  M5 tet(Q): species share of hit windows far from the gene >= the near share + 0.03.
Descriptive, no prediction: tet(W) shares; median n_sp per family and destination; reads that break the M1 rule.
V1 as in A9 (>= 95% of the reads with window counts equal to read length - k + 1); if it fails nothing is scored.
Evaluation pool: needs --allow-eval, --suffix, git tag amr-a10-registered, the committed dev outputs of the same suffix, and refuses to overwrite eval_species_windows*<suffix>.* files.
Standard library only.
"""
import argparse, csv, gzip, importlib.util, json, os, statistics, subprocess, sys, time
from collections import Counter, defaultdict

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
P = lambda *x: os.path.join(ROOT, *x)
K, NEAR, VALID = 35, 2000, 0.95


def med(v):
    return statistics.median(v) if v else float("nan")


def windows(mapping):
    out = []
    for tok in mapping.split():
        t, c = tok.split(":")
        if t in ("|", "|:|"):
            continue
        out.extend([0 if t in ("A", "0") else int(t)] * int(c))
    return out


def regions(win, gs, ge, k=K):
    last_in = ge - k
    reg = []
    for i in range(len(win)):
        if gs <= i <= last_in:
            reg.append("in")
        else:
            d = gs - i if i < gs else i - last_in
            reg.append("near" if d <= NEAR else "far")
    return reg


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pool", choices=["dev", "eval"], default="dev")
    ap.add_argument("--allow-eval", action="store_true")
    ap.add_argument("--suffix", default="_A10")
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
        if subprocess.run(["git", "-C", ROOT, "rev-parse", "-q", "--verify", "refs/tags/amr-a10-registered"], capture_output=True).returncode != 0:
            sys.exit("evaluation pool refused: git tag amr-a10-registered not found (commit and tag Addendum A10 first)")
        if not os.path.exists(P(a.out, f"dev_species_windows{a.suffix}.txt")):
            sys.exit(f"evaluation pool refused: the dev run ({a.out}/dev_species_windows{a.suffix}.txt) is not there yet")
        if any(f.startswith("eval_species_windows") and a.suffix in f for f in os.listdir(P(a.out))):
            sys.exit("evaluation outputs are never overwritten (an eval_species_windows*" + a.suffix + ".* file exists)")
    for need in (a.reads_tsv, a.paf, os.path.join(a.kraken_run, "kraken_output.txt.gz")):
        if not os.path.exists(P(need)):
            sys.exit(f"not found: {need}")
    s8 = importlib.util.spec_from_file_location("x08", P("bin", "08_trace_reads.py")); x08 = importlib.util.module_from_spec(s8); s8.loader.exec_module(x08)
    s21 = importlib.util.spec_from_file_location("x21", P("bin", "21_amr_reads.py")); x21 = importlib.util.module_from_spec(s21); s21.loader.exec_module(x21)
    for r in csv.DictReader(open(P(a.merge)), delimiter="\t"):
        x21.MERGE[r["alias"]] = r["family"]
    rows = [r for r in csv.DictReader(open(P(a.reads_tsv)), delimiter="\t") if r["destination"] != "no_truth" and r["family"] in ("cepA", "tet(Q)", "tet(W)")]
    want = {r["read"]: r for r in rows}
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
        win = windows(maps[q][1])
        if len(win) != qlen - K + 1:
            bad += 1; continue
        sp = r["true_species"]
        if sp not in infos:
            infos[sp] = x08.SpInfo(sp, rev[sp], lineage, ranks, rev)
        i = infos[sp]
        reg = regions(win, qs, qe)
        cnt = {rg: Counter() for rg in ("in", "near", "far")}
        for t, rg in zip(win, reg):
            if t == 0:
                continue
            lin = lineage(t)
            cnt[rg]["sp" if (i.t in lin or any(x in lin for x in i.alias)) else "anc" if t in i.idx else "other"] += 1
        nsp = sum(cnt[rg]["sp"] for rg in cnt)
        out.append(dict(read=r["read"], family=r["family"], true=sp, dest=r["destination"], qlen=qlen, nsp=nsp, cnt=cnt))
    n_ok = len(out); frac = n_ok / (n_ok + bad) if (n_ok + bad) else 0.0
    valid = frac >= VALID
    fam = defaultdict(list)
    for v in out:
        fam[v["family"]].append(v)

    def share(vs, rg):
        s = sum(v["cnt"][rg]["sp"] for v in vs); h = sum(sum(v["cnt"][rg].values()) for v in vs)
        return (s / h if h else float("nan")), h
    m1_ok = sum(1 for v in out if (v["dest"] == "species") == (v["nsp"] >= 1)); m1 = m1_ok / n_ok if n_ok else float("nan")
    z = lambda f: (sum(1 for v in fam.get(f, []) if v["nsp"] == 0) / len(fam[f])) if fam.get(f) else float("nan")
    tq_in, tq_near, tq_far = share(fam.get("tet(Q)", []), "in"), share(fam.get("tet(Q)", []), "near"), share(fam.get("tet(Q)", []), "far")
    ce_in = share(fam.get("cepA", []), "in")
    L = [f"AMR species-specific windows, {a.pool.upper()} pool (Addendum A10)", "",
         f"reads used {n_ok}; window count inconsistent with read length - {K} + 1: {bad}; no alignment or no Kraken2 line: {missing}; V1 (>= {VALID:.0%} consistent): {'PASSED' if valid else 'FAILED'} ({frac:.3f})", "",
         f"{'family':8s}{'reads':>6s}{'n_sp=0':>8s}   median n_sp (species | above | wrong)    species share of hit windows: in gene / near / far (hit windows in gene)"]
    for f in ("cepA", "tet(Q)", "tet(W)"):
        vs = fam.get(f, [])
        grp = lambda d: [v["nsp"] for v in vs if (v["dest"] == "species" if d == "species" else v["dest"].startswith("up_") if d == "above" else v["dest"] in ("sister", "same_family", "elsewhere"))]
        m = {d: med(grp(d)) for d in ("species", "above", "wrong")}
        sh = [share(vs, rg) for rg in ("in", "near", "far")]
        L.append(f"{f:8s}{len(vs):>6d}{sum(1 for v in vs if v['nsp'] == 0):>8d}   {m['species']:>8.0f} | {m['above']:>5.0f} | {m['wrong']:>5.0f}                    "
                 f"{sh[0][0]:.3f} / {sh[1][0]:.3f} / {sh[2][0]:.3f}  ({sh[0][1]})")
    brk = [v for v in out if (v["dest"] == "species") != (v["nsp"] >= 1)]
    L += ["", f"reads that break the M1 rule: {len(brk)}" + "".join(f"\n  {v['read']} {v['family']} destination {v['dest']} n_sp {v['nsp']}" for v in brk), ""]
    res = [("M1 rule 'species iff n_sp >= 1' right for >= 0.90 of reads", m1 >= 0.90, f"{m1:.3f} ({m1_ok} of {n_ok})"),
           ("M2 tet(Q) share with n_sp = 0 >= 0.50 and cepA <= 0.05", z("tet(Q)") >= 0.50 and z("cepA") <= 0.05, f"tet(Q) {z('tet(Q)'):.3f}, cepA {z('cepA'):.3f}"),
           ("M3 in gene: tet(Q) species share <= 0.02 and cepA >= 0.90", tq_in[0] <= 0.02 and ce_in[0] >= 0.90, f"tet(Q) {tq_in[0]:.3f} ({tq_in[1]} hit windows), cepA {ce_in[0]:.3f} ({ce_in[1]})"),
           ("M4 tet(Q) near the gene: species share <= 0.02", tq_near[0] <= 0.02, f"{tq_near[0]:.3f} ({tq_near[1]} hit windows)"),
           ("M5 tet(Q) far share >= near share + 0.03", tq_far[0] >= tq_near[0] + 0.03, f"far {tq_far[0]:.3f} ({tq_far[1]}), near {tq_near[0]:.3f}")]
    L.append("Registered predictions (A10)" + ("" if valid else " - NOT SCORED, V1 failed") + (" (dev: informed by these data, descriptive)" if a.pool == "dev" else "") + ":")
    for nme, okk, val in res:
        L.append(f"  {nme}: {('MET' if okk else 'NOT MET') if valid else 'not scored'}  ({val})")
    L += ["", "Mock community only (true species from the supplier genomes). A window is specific if its Kraken2 taxon is the true species or below. Small numbers, one pool per run. Research use only, not clinically validated."]
    txt = "\n".join(L)
    os.makedirs(P(a.out), exist_ok=True)
    tag = f"{a.pool}_species_windows{a.suffix}"
    open(P(a.out, tag + ".txt"), "w").write(txt + "\n")
    with open(P(a.out, tag + ".tsv"), "w") as o:
        w = csv.writer(o, delimiter="\t")
        w.writerow(["read", "family", "true_species", "destination", "read_len", "n_specific_windows"] + [f"{rg}_{x}" for rg in ("in", "near", "far") for x in ("sp", "anc", "other")])
        for v in out:
            w.writerow([v["read"], v["family"], v["true"], v["dest"], v["qlen"], v["nsp"]] + [v["cnt"][rg][x] for rg in ("in", "near", "far") for x in ("sp", "anc", "other")])
    json.dump({"script": "bin/29_amr_species_windows.py", "args": vars(a), "k": K, "reads_used": n_ok, "inconsistent": bad, "valid": valid, "run": time.strftime("%Y-%m-%d %H:%M:%S")},
              open(P(a.out, tag + ".json"), "w"), indent=1)
    print(txt)


if __name__ == "__main__":
    main()
