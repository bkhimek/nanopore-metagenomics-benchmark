#!/usr/bin/env python3
"""11_sylph_arm.py - Phase 2, arm 3 (Sylph), registered as D18 in docs/PROJECT_STATE.md (tag sylph-registered).

Subcommands (all paths relative to the project root):
  sketch     sylph sketch of a pool (default deduplication = registered arm 3a; --no-dedup only for the dedup check).
  profile    sylph profile of that sketch against the two databases JOINTLY (arm 3a; --arm 3c adds -u, exploratory).
  bootstrap  DEV ONLY. Block bootstrap: B resamples of the 37 dev blocks (seed 20261017) are built by concatenating block files,
             sketched with --no-dedup (a block drawn twice would otherwise be thrown away as duplicate reads) and profiled.
  score      Score one profile (arm 3a, 3b with --cutoff or from dev_selection.json, 3c) with the conventions of bin/04_score.py.
  tune       DEV ONLY. Arm 3b: score every effective-coverage cutoff on the bootstrap profiles, choose the SMALLEST cutoff whose
             bootstrap-mean F1 is within 0.002 of the maximum, write reports/sylph/dev_selection.json.
  diff       Species-level comparison of two profiles (dedup check).
The evaluation pool needs --allow-eval, reports/sylph/dev_selection.json and the git tag sylph-frozen; its outputs are never overwritten.

Scoring (same conventions as 04_score.py): truth = 17 species of reference/d6331_expected_reads_by_species.csv; expected reads = pool reads x
genomic-DNA fraction; EVALUABLE = expected >= 20 (14 species); a called truth species is never a false positive; precision = called mock / called;
recall over the evaluable species; F1 = 0 when nothing right is called; Bray-Curtis = 1 - sum min(o, e) with o = abundance share among ALL called
species and e = reference share over all 17 truth species (so uncalled truth species and false positives both count).
Sylph genome rows -> species through reference/d6331_gtdb_crosswalk.csv (rows summed per species); a species outside the crosswalk is a
non-mock false positive. The effective-coverage filter of arm 3b is applied to genome rows (>= cutoff) before summing and renormalising.
obs/exp: Sequence_abundance has no absolute read count without -u, so obs/exp = (share among CALLED MOCK species) / (expected DNA share among the
same species).  Kraken2 comparators are put on the same basis with --kraken label=species_table.tsv.
"""
import argparse, collections, csv, glob, gzip, importlib.util, json, math, os, random, shutil, subprocess, sys, threading, time

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
IMAGE = "quay.io/biocontainers/sylph:1.0.0--hb42e459_0"
DBS = ["data/db/sylph/gtdb-r232-c200-dbv2.syl2db", "data/db/sylph/fungi-refseq-2024-07-25-c200-v0.3.syldb"]
GTDB_TAX = ["data/db/gtdb_tax/bac120_taxonomy_r232.tsv.gz", "data/db/gtdb_tax/ar53_taxonomy_r232.tsv.gz"]
GRID = [0, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0]
SEED = 20261017
DRY = False

_spec = importlib.util.spec_from_file_location("cw10", os.path.join(ROOT, "bin", "10_sylph_crosswalk.py"))
cw = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(cw)
P = lambda *x: os.path.join(ROOT, *x)
fmt = lambda x: "nan" if (isinstance(x, float) and math.isnan(x)) else (f"{x:.3f}" if isinstance(x, float) else str(x))


# ---------------------------------------------------------------- inputs
def manifest():
    mf = glob.glob(P("data", "subsets", "*.sra_blocks.manifest.json"))
    if len(mf) != 1:
        sys.exit(f"expected exactly one manifest in data/subsets, found {len(mf)}")
    return json.load(open(mf[0]))


def pool_fastq(pool):
    return P("data", "subsets", manifest()["pools"][pool]["file"])


def rel_data(path):
    return "/data/" + os.path.relpath(path, P("data"))


def guard(a):
    if a.pool != "eval":
        return
    if not getattr(a, "allow_eval", False):
        sys.exit("evaluation pool: needs --allow-eval, reports/sylph/dev_selection.json and the git tag sylph-frozen")
    if not os.path.exists(P("reports", "sylph", "dev_selection.json")):
        sys.exit("evaluation pool refused: reports/sylph/dev_selection.json does not exist (tune on dev first)")
    tag = subprocess.run(["git", "-C", ROOT, "tag", "-l", "sylph-frozen"], capture_output=True, text=True).stdout.strip()
    if not tag:
        sys.exit("evaluation pool refused: git tag sylph-frozen not found (commit dev_selection.json and tag it first)")


def load_truth(truth_path, N, min_expected):
    truth = {}
    for t in csv.DictReader(open(P(truth_path))):
        tid = int(t["species_taxid"])
        exp = N * float(t["genomic_dna_pct"]) / 100.0
        truth[tid] = {"name": t["ncbi_name_in_kraken_db"], "dna_pct": float(t["genomic_dna_pct"]), "gc_pct": float(t["genome_copy_pct"]),
                      "expected": exp, "evaluable": exp >= min_expected, "domain": t["domain"], "genus": t["ncbi_name_in_kraken_db"].split()[0]}
    return truth


def load_crosswalk(path):
    """-> ({database species: mock taxid}, set of all mock taxids listed incl. 'not_represented' ones)"""
    m, listed = {}, set()
    for r in csv.DictReader(open(P(path))):
        tid = int(r["mock_taxid"])
        listed.add(tid)
        if r["status"] != "ok":
            continue
        if m.setdefault(r["database_species"], tid) != tid:
            sys.exit(f"crosswalk: database species {r['database_species']} maps to two mock species")
    return m, listed


def read_profile(path, gtdb):
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return []
    head, rows = cw.read_sylph(path)
    if not rows:
        return []
    gc, ac = cw.find_col(head, "Genome_file"), cw.find_col(head, "Adjusted_ANI")
    sc, tc = cw.find_col(head, "Sequence_abundance"), cw.find_col(head, "Taxonomic_abundance")
    ec = next((h for h in head if h.lower() == "eff_cov"), None) or cw.find_col(head, "Median_cov")
    cc = next((h for h in head if h.lower() == "contig_name"), None)
    out = []
    for r in rows:
        sp, src, acc = cw.species_of(r, gc, cc, gtdb)
        out.append({"species": sp or "unresolved:" + (acc or r[gc]), "seq": float(r[sc]), "tax": float(r[tc]), "eff": float(r[ec]),
                    "ani": float(r[ac]), "src": src})
    return out


def species_abund(rows, cutoff, cwmap, which="seq"):
    """genome rows -> {key: abundance}; key = mock taxid (int) or 'db:<species>' for a non-mock species. Filter: eff >= cutoff."""
    ab, info = collections.defaultdict(float), {}
    for r in rows:
        if r["eff"] < cutoff:
            continue
        key = cwmap.get(r["species"], "db:" + r["species"])
        ab[key] += r[which]
        d = info.setdefault(key, {"species": [], "genomes": 0, "eff": 0.0, "seq": 0.0, "tax": 0.0})
        if r["species"] not in d["species"]:
            d["species"].append(r["species"])
        d["genomes"] += 1; d["eff"] = max(d["eff"], r["eff"]); d["seq"] += r["seq"]; d["tax"] += r["tax"]
    return dict(ab), info


# ---------------------------------------------------------------- metrics (04_score.py conventions)
def metrics(ab, truth, ref="dna_pct"):
    evaluable = {t for t, v in truth.items() if v["evaluable"]}
    called = {k for k, v in ab.items() if v > 0}
    correct = {k for k in called if k in truth}
    fp = called - correct
    tp_eval = correct & evaluable
    prec = len(correct) / len(called) if called else float("nan")
    rec_e = len(tp_eval) / len(evaluable)
    f1 = 0.0 if not called else (2 * prec * rec_e / (prec + rec_e) if (prec + rec_e) > 0 else 0.0)
    tot = sum(ab[k] for k in called)
    rsum = sum(v[ref] for v in truth.values())
    bc = l1 = float("nan")
    if tot:
        keys = called | set(truth)
        o = {k: (ab.get(k, 0.0) / tot if k in called else 0.0) for k in keys}
        e = {k: (truth[k][ref] / rsum if k in truth else 0.0) for k in keys}
        l1 = sum(abs(o[k] - e[k]) for k in keys)
        bc = 1 - sum(min(o[k], e[k]) for k in keys)
    return {"called": len(called), "correct": len(correct), "fp": len(fp), "tp_evaluable": len(tp_eval), "precision": prec,
            "recall_evaluable": rec_e, "recall_all_truth": len(correct) / len(truth), "f1": f1, "l1": l1, "bray_curtis": bc,
            "missed_evaluable": sorted(evaluable - tp_eval), "fp_keys": sorted(fp, key=str)}


def share_stats(ab, reference):
    """Among CALLED MOCK species: observed share, expected share, obs/exp, half-L1 (BC among called mock species). reference = {taxid: value}."""
    mock = {k: v for k, v in ab.items() if isinstance(k, int) and v > 0 and k in reference}
    ot, rt = sum(mock.values()), sum(reference[k] for k in mock)
    if not mock or ot == 0 or rt == 0:
        return {}, float("nan")
    rows = {k: {"obs_share": v / ot, "exp_share": reference[k] / rt, "obs_over_exp": (v / ot) / (reference[k] / rt) if reference[k] else float("nan")}
            for k, v in mock.items()}
    return rows, 0.5 * sum(abs(r["obs_share"] - r["exp_share"]) for r in rows.values())


# ---------------------------------------------------------------- docker
def image_digest():
    r = subprocess.run(["docker", "image", "inspect", IMAGE, "--format", "{{json .RepoDigests}}"], capture_output=True, text=True)
    if r.returncode != 0:
        subprocess.run(["docker", "pull", IMAGE], check=True)
        r = subprocess.run(["docker", "image", "inspect", IMAGE, "--format", "{{json .RepoDigests}}"], capture_output=True, text=True)
    return r.stdout.strip()


def run_sylph(args_, tag, mem, threads):
    name = f"sylph_{tag}_{os.getpid()}"
    cmd = ["docker", "run", "--rm", "--name", name, "-u", f"{os.getuid()}:{os.getgid()}", "-m", mem, "--memory-swap", mem, "--cpus", str(threads),
           "-v", f"{P('data')}:/data", IMAGE, "sylph"] + args_
    if DRY:
        print(" ".join(cmd)); return {"returncode": 0, "wall_s": 0, "peak_gib": 0, "command": " ".join(cmd), "dry_run": True}
    peak, stop = [0.0], threading.Event()

    def poll():
        while not stop.is_set():
            r = subprocess.run(["docker", "stats", "--no-stream", "--format", "{{.MemUsage}}", name], capture_output=True, text=True)
            if r.returncode == 0 and r.stdout.strip():
                peak[0] = max(peak[0], cw.mem_to_gib(r.stdout.split("/")[0]))
            stop.wait(2)
    t0 = time.time(); th = threading.Thread(target=poll, daemon=True); th.start()
    rc = subprocess.run(cmd).returncode
    stop.set(); th.join(timeout=5)
    return {"returncode": rc, "wall_s": round(time.time() - t0, 1), "peak_gib": round(peak[0], 2), "command": " ".join(cmd)}


def record(pool, name, info, extra=None):
    d = {"script": "bin/11_sylph_arm.py", "time": time.strftime("%Y-%m-%d %H:%M:%S"), "image": IMAGE,
         "repo_digests": "dry-run" if DRY else image_digest(), **info, **(extra or {})}
    os.makedirs(P("reports", "sylph", "runs"), exist_ok=True)
    json.dump(d, open(P("reports", "sylph", "runs", f"{pool}_{name}.json"), "w"), indent=2)
    return d


def need_ok(info, what):
    if info["returncode"] != 0:
        sys.exit(f"{what} failed with exit code {info['returncode']} (137 = memory cap hit)")


# ---------------------------------------------------------------- sketch / profile / bootstrap
def cmd_sketch(a):
    guard(a)
    out = P("data", "sylph", a.pool, "sketch_nodedup" if a.no_dedup else "sketch")
    if glob.glob(os.path.join(out, "*.sylsp")):
        sys.exit(f"{out} already has a sketch; not overwritten")
    os.makedirs(out, exist_ok=True)
    args_ = ["sketch", "-t", str(a.threads), "-d", rel_data(out), "-r", rel_data(pool_fastq(a.pool))] + (["--no-dedup"] if a.no_dedup else [])
    info = run_sylph(args_, "sketch", a.memory, a.threads)
    need_ok(info, "sylph sketch")
    record(a.pool, "sketch_nodedup" if a.no_dedup else "sketch", info)
    print(f"sketch done: {info['wall_s']} s, polled peak {info['peak_gib']} GiB; files: {glob.glob(os.path.join(out, '*'))}")


def cmd_profile(a):
    guard(a)
    sk = glob.glob(os.path.join(P(a.sketch_dir) if a.sketch_dir else P("data", "sylph", a.pool, "sketch"), "*.sylsp"))
    if len(sk) != 1 and not DRY:
        sys.exit(f"expected exactly one .sylsp, found {sk}")
    out = P("data", "sylph", a.pool, f"arm_{a.arm}{a.suffix}.tsv")
    if os.path.exists(out):
        sys.exit(f"{out} exists; not overwritten")
    args_ = ["profile", rel_data(sk[0] if sk else "x.sylsp"), "-d"] + ["/data/" + d[len("data/"):] for d in DBS] + \
            ["-t", str(a.threads), "-o", rel_data(out)] + (["-u"] if a.arm == "3c" else [])
    info = run_sylph(args_, "profile", a.memory, a.threads)
    need_ok(info, "sylph profile")
    record(a.pool, f"profile_{a.arm}{a.suffix}", info)
    print(f"profile done: {info['wall_s']} s, polled peak {info['peak_gib']} GiB -> {os.path.relpath(out, ROOT)}")


def cmd_bootstrap(a):
    if a.pool != "dev":
        sys.exit("bootstrap is development-pool only")
    man = manifest()
    blocks = sorted([b for b in man["blocks"] if b["pool"] == "dev"], key=lambda b: b["start"])
    rng = random.Random(SEED)
    bmax = a.B if a.B else 100
    draws = [[rng.randrange(len(blocks)) for _ in blocks] for _ in range(bmax)]
    os.makedirs(P("reports", "sylph"), exist_ok=True)
    json.dump({"seed": SEED, "n_blocks": len(blocks), "draws": draws}, open(P("reports", "sylph", "bootstrap_draws_dev.json"), "w"))
    bdir = P("data", "sylph", "dev", "boot"); os.makedirs(bdir, exist_ok=True)
    tfile = os.path.join(bdir, "times.json")
    times = json.load(open(tfile)) if os.path.exists(tfile) else {}
    B, reason = bmax, "as requested" if a.B else "default 100"
    for b in range(1, bmax + 1):
        if b > B:
            break
        tsv = os.path.join(bdir, f"b{b:03d}.tsv")
        if os.path.exists(tsv):
            continue
        fq, skd = os.path.join(bdir, f"b{b:03d}.fastq.gz"), os.path.join(bdir, f"b{b:03d}_sk")
        with open(fq, "wb") as o:
            for i in draws[b - 1]:
                with open(P("data", "subsets", "blocks", blocks[i]["file"]), "rb") as f:
                    shutil.copyfileobj(f, o)
        shutil.rmtree(skd, ignore_errors=True)
        os.makedirs(skd, exist_ok=True)
        s1 = run_sylph(["sketch", "-t", str(a.threads), "-d", rel_data(skd), "-r", rel_data(fq), "--no-dedup"], "bsk", a.memory, a.threads)
        need_ok(s1, f"sketch b{b}")
        sk = glob.glob(os.path.join(skd, "*.sylsp"))
        s2 = run_sylph(["profile", rel_data(sk[0]) if sk else "x", "-d"] + ["/data/" + d[len("data/"):] for d in DBS] +
                       ["-t", str(a.threads), "-o", rel_data(tsv)], "bpr", a.memory, a.threads)
        need_ok(s2, f"profile b{b}")
        times[f"b{b:03d}"] = {"sketch_s": s1["wall_s"], "profile_s": s2["wall_s"], "peak_gib": max(s1["peak_gib"], s2["peak_gib"])}
        json.dump(times, open(tfile, "w"), indent=1)
        if not DRY:
            os.remove(fq); shutil.rmtree(skd, ignore_errors=True)
        print(f"b{b:03d}: sketch {s1['wall_s']} s, profile {s2['wall_s']} s, peak {times[f'b{b:03d}']['peak_gib']} GiB")
        if b == 1 and not a.B and s1["wall_s"] > 180:
            B, reason = 50, f"registered rule: sketching one resample took {s1['wall_s']} s (> about 180 s), B lowered to 50"
            print(reason)
    json.dump({"B": B, "reason": reason, "seed": SEED, "no_dedup_in_resamples": True}, open(P("reports", "sylph", "bootstrap_dev_run.json"), "w"), indent=1)


# ---------------------------------------------------------------- score / tune / diff
def context(a):
    man = manifest()
    N = man["pools"][a.pool]["reads"]
    truth = load_truth(a.truth, N, a.evaluable_min_expected)
    cwmap, listed = load_crosswalk(a.crosswalk)
    miss = set(truth) - listed
    if miss:
        sys.exit(f"crosswalk does not list truth taxids {sorted(miss)}")
    unrep = sorted(set(truth) - set(cwmap.values()))
    if unrep:
        print("NOTE: not represented in the database (cannot be detected): " + ", ".join(truth[t]["name"] for t in unrep), file=sys.stderr)
    gtdb = cw.load_gtdb([P(x) for x in a.gtdb_taxonomy])
    return N, truth, cwmap, gtdb


def category(key, species, truth_genera):
    return "same genus as a mock species" if species.split()[0] in truth_genera else "other"


def cmd_score(a):
    guard(a)
    N, truth, cwmap, gtdb = context(a)
    cutoff = a.cutoff
    if a.arm == "3b":
        if cutoff is None:
            sel = P("reports", "sylph", "dev_selection.json")
            if not os.path.exists(sel):
                sys.exit("arm 3b needs --cutoff (dev) or reports/sylph/dev_selection.json")
            cutoff = json.load(open(sel))["chosen_cutoff"]
    else:
        cutoff = 0
    prof = P(a.profile) if a.profile else P("data", "sylph", a.pool, f"arm_{'3c' if a.arm == '3c' else '3a'}.tsv")
    rows = read_profile(prof, gtdb)
    if not rows:
        sys.exit(f"no genome rows in {prof}")
    name = f"{a.pool}_{a.arm}" + (f"_cut{cutoff:g}" if a.arm == "3b" else "")
    out = P("reports", "sylph", name)
    if a.pool == "eval" and os.path.exists(out):
        sys.exit(f"{out} exists: the evaluation pool is scored once per arm")
    os.makedirs(out, exist_ok=True)

    ab, info = species_abund(rows, cutoff, cwmap, "seq")
    abt, _ = species_abund(rows, cutoff, cwmap, "tax")
    m = metrics(ab, truth, "dna_pct")
    shares, bc_mock = share_stats(ab, {t: v["dna_pct"] for t, v in truth.items()})
    m2 = metrics(abt, truth, "gc_pct")
    _, bc_mock_gc = share_stats(abt, {t: v["gc_pct"] for t, v in truth.items()})
    truth_genera = {v["genus"] for v in truth.values()} | {s.split()[0] for s in cwmap}
    tot = sum(v for v in ab.values())
    L = [f"Sylph arm {a.arm} on the {a.pool} pool: {len(rows)} genome rows -> {len(ab)} species after filter eff_cov >= {cutoff:g} "
         f"(profile {os.path.relpath(prof, ROOT)}, pool reads {N})",
         f"Evaluable species (expected >= {a.evaluable_min_expected:g} reads): {sum(v['evaluable'] for v in truth.values())} of {len(truth)}", "",
         f"PRIMARY (Sequence_abundance vs genomic-DNA composition): called {m['called']}; correct {m['correct']}; false positives {m['fp']}; "
         f"evaluable detected {m['tp_evaluable']}/{sum(v['evaluable'] for v in truth.values())}",
         f"  precision {fmt(m['precision'])}  recall(evaluable) {fmt(m['recall_evaluable'])}  recall(all {len(truth)}) {fmt(m['recall_all_truth'])}  F1 {fmt(m['f1'])}  "
         f"L1 {fmt(m['l1'])}  Bray-Curtis (04 convention) {fmt(m['bray_curtis'])}",
         f"  Bray-Curtis among called mock species only (09 convention) {fmt(bc_mock)}",
         "  missed evaluable: " + (", ".join(truth[t]["name"] for t in m["missed_evaluable"]) or "none"), "",
         f"SECONDARY (Taxonomic_abundance vs genome-copy composition): Bray-Curtis (04 convention) {fmt(m2['bray_curtis'])}; among called mock species {fmt(bc_mock_gc)}; "
         f"F1 {fmt(m2['f1'])}", ""]
    tab = []
    L.append(f"{'species':32s}{'DNA%':>7s}{'exp.reads':>10s}{'seq_ab%':>9s}{'tax_ab%':>9s}{'effcov':>8s}{'obs/exp*':>9s}  status")
    for t, v in sorted(truth.items(), key=lambda kv: -kv[1]["dna_pct"]):
        d = info.get(t)
        oe = shares.get(t, {}).get("obs_over_exp", float("nan"))
        row = {"taxid": t, "species": v["name"], "dna_pct": v["dna_pct"], "expected_reads": round(v["expected"], 1), "evaluable": "yes" if v["evaluable"] else "no",
               "database_species": ";".join(d["species"]) if d else "", "genomes": d["genomes"] if d else 0,
               "seq_abundance_pct": round(ab.get(t, 0.0), 4), "tax_abundance_pct": round(abt.get(t, 0.0), 4), "eff_cov_max": round(d["eff"], 3) if d else 0,
               "obs_over_exp_share": round(oe, 3) if oe == oe else "nan", "status": ("detected" if d else "missed") + ("" if v["evaluable"] else " (not evaluable)")}
        tab.append(row)
        L.append(f"{v['name'][:31]:32s}{v['dna_pct']:>7g}{row['expected_reads']:>10.1f}{row['seq_abundance_pct']:>9.3f}{row['tax_abundance_pct']:>9.3f}"
                 f"{row['eff_cov_max']:>8.3f}{fmt(oe):>9s}  {row['status']}")
    L.append("  * obs/exp = share among called mock species / expected DNA share among the same species (no absolute read counts without -u)")
    fps = sorted([(k, v) for k, v in ab.items() if not isinstance(k, int) and v > 0], key=lambda kv: -kv[1])
    L.append(""); L.append(f"False-positive audit ({len(fps)} non-mock species):")
    fpt = []
    for k, v in fps:
        sp = "; ".join(info[k]["species"]); cat = category(k, sp, truth_genera)
        fpt.append({"species": sp, "seq_abundance_pct": round(v, 4), "genomes": info[k]["genomes"], "eff_cov_max": round(info[k]["eff"], 3), "category": cat})
        L.append(f"  {sp:40s}{v:>8.3f} %  effcov {info[k]['eff']:.3f}  [{cat}]")
    n_same = sum(1 for x in fpt if x["category"].startswith("same genus"))
    L.append(f"  non-mock false positives by the REGISTERED rule: {len(fpt)} (of which {n_same} in the genus of a mock species or of a crosswalk species)")

    mp = P("reports", "trace", "trace_by_species.tsv")
    summary = {"pool": a.pool, "arm": a.arm, "cutoff": cutoff, "profile": os.path.relpath(prof, ROOT), "primary": {k: v for k, v in m.items() if k != "fp_keys"},
               "bray_curtis_among_called_mock": bc_mock, "secondary": {k: v for k, v in m2.items() if k != "fp_keys"},
               "bray_curtis_taxonomic_among_called_mock": bc_mock_gc, "false_positives": fpt,
               "b_adolescentis_obs_over_exp": shares.get(1680, {}).get("obs_over_exp"), "n_nonmock_fp_registered": len(fpt)}
    if a.pool == "dev" and os.path.exists(mp):
        mapped = {int(r["taxid"]): float(r["mapped_reads"]) for r in csv.DictReader(open(mp), delimiter="\t")}
        _, bcm = share_stats(ab, mapped)
        L.append(""); L.append(f"EXPLORATORY (dev only): Bray-Curtis among called mock species vs the MAPPING-derived composition: {fmt(bcm)}")
        summary["bray_curtis_vs_mapping_exploratory"] = bcm
    for lab, path in (a.kraken or []):
        with open(P(path)) as fh:
            head = fh.readline().rstrip("\n").split("\t")
            kr = [dict(zip(head, l.rstrip("\n").split("\t"))) for l in fh if l.strip()]
        scol = next(h for h in head if h.startswith("status_"))
        kab = {int(r["taxid"]): float(r["observed_reads"]) for r in kr if r[scol].startswith("detected") and float(r["dna_pct"]) > 0}
        ks, kbc = share_stats(kab, {t: v["dna_pct"] for t, v in truth.items()})
        L.append(f"Kraken2 comparator {lab}: B. adolescentis obs/exp on the SAME share basis {fmt(ks.get(1680, {}).get('obs_over_exp', float('nan')))}; "
                 f"Bray-Curtis among called mock species {fmt(kbc)}")
        summary.setdefault("kraken_comparators", {})[lab] = {"b_adolescentis_share_obs_over_exp": ks.get(1680, {}).get("obs_over_exp"), "bc_among_called_mock": kbc}
    with open(os.path.join(out, "species_table.tsv"), "w") as o:
        o.write("\t".join(tab[0]) + "\n")
        for r in tab:
            o.write("\t".join(str(r[c]) for c in tab[0]) + "\n")
    with open(os.path.join(out, "fp_audit.tsv"), "w") as o:
        o.write("species\tseq_abundance_pct\tgenomes\teff_cov_max\tcategory\n")
        for x in fpt:
            o.write("\t".join(str(x[c]) for c in ("species", "seq_abundance_pct", "genomes", "eff_cov_max", "category")) + "\n")
    json.dump(summary, open(os.path.join(out, "score_summary.json"), "w"), indent=2, default=str)
    open(os.path.join(out, "score_report.txt"), "w").write("\n".join(L) + "\n")
    print("\n".join(L)); print(f"\nwrote {os.path.relpath(out, ROOT)}/")


def pct(v, q):
    s = sorted(x for x in v if x == x)
    return s[int(q * (len(s) - 1))] if s else float("nan")


def cmd_tune(a):
    if a.pool != "dev":
        sys.exit("tune is development-pool only")
    N, truth, cwmap, gtdb = context(a)
    files = sorted(glob.glob(P("data", "sylph", "dev", "boot", "b*.tsv")))
    runf = P("reports", "sylph", "bootstrap_dev_run.json")
    B = json.load(open(runf))["B"] if os.path.exists(runf) else len(files)
    files = files[:B]
    if len(files) < B:
        sys.exit(f"only {len(files)} of {B} bootstrap profiles present")
    res = {c: {"f1": [], "precision": [], "recall_evaluable": [], "bray_curtis": [], "fp": [], "badol": []} for c in GRID}
    for f in files:
        rows = read_profile(f, gtdb)
        for c in GRID:
            ab, _ = species_abund(rows, c, cwmap, "seq")
            m = metrics(ab, truth)
            sh, _ = share_stats(ab, {t: v["dna_pct"] for t, v in truth.items()})
            for k, v in (("f1", m["f1"]), ("precision", m["precision"]), ("recall_evaluable", m["recall_evaluable"]), ("bray_curtis", m["bray_curtis"]),
                         ("fp", m["fp"]), ("badol", sh.get(1680, {}).get("obs_over_exp", float("nan")))):
                res[c][k].append(v)
    mean = lambda v: (lambda s: sum(s) / len(s) if s else float("nan"))([x for x in v if x == x])
    table = []
    point = None
    pp = P("data", "sylph", "dev", "arm_3a.tsv")
    if os.path.exists(pp):
        prow = read_profile(pp, gtdb)
        point = {c: metrics(species_abund(prow, c, cwmap, "seq")[0], truth)["f1"] for c in GRID}
    for c in GRID:
        r = res[c]
        table.append({"cutoff": c, "mean_f1": mean(r["f1"]), "f1_lo": pct(r["f1"], 0.025), "f1_hi": pct(r["f1"], 0.975), "mean_precision": mean(r["precision"]),
                      "mean_recall_evaluable": mean(r["recall_evaluable"]), "mean_bc": mean(r["bray_curtis"]), "mean_fp": mean(r["fp"]),
                      "mean_b_adol_obs_over_exp": mean(r["badol"]), "point_f1": point[c] if point else float("nan")})
    best = max(t["mean_f1"] for t in table)
    chosen = min(t["cutoff"] for t in table if t["mean_f1"] >= best - 0.002)
    L = [f"Arm 3b tuning on the dev pool: {len(files)} bootstrap profiles (seed {SEED}); rule: smallest cutoff with mean F1 within 0.002 of the maximum ({best:.4f})",
         f"{'cutoff':>7s}{'meanF1':>8s}{'lo':>7s}{'hi':>7s}{'prec':>7s}{'rec(ev)':>8s}{'BC':>7s}{'FP':>6s}{'Badol':>7s}{'pointF1':>8s}"]
    for t in table:
        L.append(f"{t['cutoff']:>7g}{fmt(t['mean_f1']):>8s}{fmt(t['f1_lo']):>7s}{fmt(t['f1_hi']):>7s}{fmt(t['mean_precision']):>7s}"
                 f"{fmt(t['mean_recall_evaluable']):>8s}{fmt(t['mean_bc']):>7s}{fmt(t['mean_fp']):>6s}{fmt(t['mean_b_adol_obs_over_exp']):>7s}{fmt(t['point_f1']):>8s}")
    L.append(f"CHOSEN cutoff: {chosen:g}")
    os.makedirs(P("reports", "sylph"), exist_ok=True)
    json.dump({"rule": "smallest cutoff whose bootstrap-mean F1 is within 0.002 of the maximum (D16 rule)", "B": len(files), "seed": SEED, "grid": GRID,
               "chosen_cutoff": chosen, "table": table}, open(P("reports", "sylph", "dev_selection.json"), "w"), indent=2, default=str)
    open(P("reports", "sylph", "dev_tuning.txt"), "w").write("\n".join(L) + "\n")
    print("\n".join(L)); print("\nwrote reports/sylph/dev_selection.json and dev_tuning.txt")


def cmd_diff(a):
    gtdb = cw.load_gtdb([P(x) for x in a.gtdb_taxonomy])
    def agg(p):
        d = collections.defaultdict(float)
        for r in read_profile(P(p), gtdb):
            d[r["species"]] += r["seq"]
        return d
    x, y = agg(a.a), agg(a.b)
    keys = sorted(set(x) | set(y), key=lambda k: -max(x.get(k, 0), y.get(k, 0)))
    print(f"{'species':42s}{a.a[-24:]:>26s}{a.b[-24:]:>26s}")
    for k in keys:
        print(f"{k[:41]:42s}{x.get(k, 0):>26.4f}{y.get(k, 0):>26.4f}")
    print(f"species only in first: {sorted(set(x) - set(y))}; only in second: {sorted(set(y) - set(x))}; max abs difference {max(abs(x.get(k, 0) - y.get(k, 0)) for k in keys):.4f} percentage points")


def main():
    global DRY
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dry-run", action="store_true", help="print the docker commands without running them")
    sub = ap.add_subparsers(dest="cmd", required=True)

    def common(p, docker=False):
        p.add_argument("--pool", choices=["dev", "eval"], default="dev")
        p.add_argument("--allow-eval", action="store_true")
        if docker:
            p.add_argument("--threads", type=int, default=6)
            p.add_argument("--memory", default="10g")
        else:
            p.add_argument("--truth", default="reference/d6331_expected_reads_by_species.csv")
            p.add_argument("--crosswalk", default="reference/d6331_gtdb_crosswalk.csv")
            p.add_argument("--gtdb-taxonomy", nargs="+", default=GTDB_TAX)
            p.add_argument("--evaluable-min-expected", type=float, default=20.0)
    p = sub.add_parser("sketch"); common(p, True); p.add_argument("--no-dedup", action="store_true"); p.set_defaults(f=cmd_sketch)
    p = sub.add_parser("profile"); common(p, True); p.add_argument("--arm", choices=["3a", "3c"], default="3a")
    p.add_argument("--sketch-dir"); p.add_argument("--suffix", default=""); p.set_defaults(f=cmd_profile)
    p = sub.add_parser("bootstrap"); common(p, True); p.add_argument("--B", type=int, help="number of resamples (default 100; 50 if the first sketch takes > 180 s)"); p.set_defaults(f=cmd_bootstrap)
    p = sub.add_parser("score"); common(p); p.add_argument("--arm", choices=["3a", "3b", "3c"], required=True); p.add_argument("--cutoff", type=float)
    p.add_argument("--profile"); p.add_argument("--kraken", nargs="*", type=lambda s: tuple(s.split("=", 1)), help="label=path/to/species_table.tsv"); p.set_defaults(f=cmd_score)
    p = sub.add_parser("tune"); common(p); p.set_defaults(f=cmd_tune)
    p = sub.add_parser("diff"); p.add_argument("a"); p.add_argument("b"); p.add_argument("--gtdb-taxonomy", nargs="+", default=GTDB_TAX); p.set_defaults(f=cmd_diff)
    a = ap.parse_args()
    DRY = a.dry_run
    a.f(a)


if __name__ == "__main__":
    main()
