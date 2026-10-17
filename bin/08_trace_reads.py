#!/usr/bin/env python3
"""08_trace_reads.py - where do the reads of each mock species go in Kraken2?  (EXPLORATORY, dev pool only)

Idea: label every dev read by what it MAPS to (minimap2 against the D6331 reference genomes, independent of Kraken2),
then cross that label with the taxon Kraken2 gave the same read.  For each mock species this gives
  - how many reads truly come from it (mapped share) versus its theoretical DNA share  -> tests "the tube has less of it"
  - where its reads ended up in Kraken2: correct species / genus / family / higher / sister species / same family / elsewhere / unclassified
  -> tests "Kraken2 loses them" and says where.

Rules (fixed here, before the first run; none is tuned on results):
  - dev pool only; there is no option to touch the eval pool.
  - minimap2 -x map-ont --secondary=no against all 21 genome FASTA files (contigs renamed STEM__n).
  - A read is assigned to the species with the largest summed number of matching bases over its alignments.
    It needs matching bases >= --min-cov (default 0.5) x read length, otherwise 'weak' (unmapped).
    If the runner-up species has >= --amb (default 0.9) x the best, the read is 'ambiguous'.
  - Strains of E. coli are one species.  Fusobacterium: Kraken2 calls under F. animalis count as correct (alias D13).
  - Destination of a read, given its true species S and the Kraken2 taxon k:
      species      k is S or below it (or, for Fusobacterium, F. animalis or below)
      up_genus     k is an ancestor of S at or below the genus node
      up_family    k is an ancestor above genus, up to the family node
      up_higher    k is an ancestor above the family
      sister       k is inside the same genus, another species
      same_family  k is inside the same family, another genus
      elsewhere    anything else (for example Homo sapiens)
      unclassified Kraken2 gave taxid 0
Caveats: mapping needs the sequenced strain to be close to the reference (a divergent strain shows as unmapped, not as a
Kraken2 loss); the two yeast references are draft assemblies; mapping shares are read shares (follow DNA mass), as expected_reads do.
Standard library only; minimap2 runs in a Docker image (recorded in run.json).
"""
import argparse, gzip, hashlib, json, os, re, subprocess, sys, time
from collections import Counter, defaultdict
from pathlib import Path

BUCKETS = ["species", "up_genus", "up_family", "up_higher", "sister", "same_family", "elsewhere", "unclassified"]
ALIAS_NAMES = {"Fusobacterium nucleatum": ["Fusobacterium animalis"]}  # D13, reference/d6331_taxid_aliases.csv


def stem_species(stem):
    if stem.startswith("Escherichia_coli"):
        return "Escherichia coli"
    special = {"Candida_albican": "Candida albicans",
               "Lactobacillus_fermentum": "Limosilactobacillus fermentum"}
    return special.get(stem, stem.replace("_", " "))


def opener(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "rt")


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def find_dmp(db, name):
    db = Path(db)
    for cand in (db / "taxonomy" / name, db / name):
        if cand.exists():
            return cand
    hits = sorted(db.rglob(name))
    if not hits:
        sys.exit(f"ERROR: {name} not found under {db}")
    return hits[0]


def load_nodes(path):
    parent, rank = {}, {}
    with open(path) as fh:
        for line in fh:
            f = line.split("\t|\t")
            t = int(f[0]); parent[t] = int(f[1]); rank[t] = f[2]
    return parent, rank


def load_names(path):
    names, rev = {}, {}
    with open(path) as fh:
        for line in fh:
            f = line.split("\t|\t")
            if f[3].startswith("scientific name"):
                t = int(f[0]); names[t] = f[1]; rev.setdefault(f[1], t)
    return names, rev


def make_lineage(parent):
    cache = {}
    def lineage(t):
        r = cache.get(t)
        if r is None:
            r, x, seen = [], t, set()
            while x in parent and x not in seen:
                r.append(x); seen.add(x)
                if parent[x] == x:
                    break
                x = parent[x]
            r = tuple(r)
            cache[t] = r
        return r
    return lineage


def read_kraken(path):
    out = {}
    with opener(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 3:
                continue
            if f[0] == "U":
                out[f[1]] = 0
            else:
                m = re.search(r"taxid (\d+)", f[2])
                out[f[1]] = int(m.group(1)) if m else int(f[2])
    return out


def read_paf(path):
    agg = {}
    with open(path) as fh:
        for line in fh:
            f = line.split("\t", 12)
            sp = stem_species(f[5].split("__")[0])
            d = agg.get(f[0])
            if d is None:
                d = agg[f[0]] = [int(f[1]), {}]
            d[1][sp] = d[1].get(sp, 0) + int(f[9])
    return agg


def assign(agg_entry, min_cov, amb):
    qlen, per = agg_entry
    ranked = sorted(per.items(), key=lambda kv: -kv[1])
    best_sp, best = ranked[0]
    if best / qlen < min_cov:
        return "weak"
    if len(ranked) > 1 and ranked[1][1] >= amb * best:
        return "ambiguous"
    return best_sp


class SpInfo:
    def __init__(self, name, taxid, lineage, ranks, rev):
        self.name, self.t = name, taxid
        self.anc = lineage(taxid)
        self.idx = {x: i for i, x in enumerate(self.anc)}
        self.gi = next((i for i, x in enumerate(self.anc) if ranks.get(x) == "genus"), None)
        self.fi = next((i for i, x in enumerate(self.anc) if ranks.get(x) == "family"), None)
        self.g = self.anc[self.gi] if self.gi is not None else None
        self.f = self.anc[self.fi] if self.fi is not None else None
        self.alias = [rev[n] for n in ALIAS_NAMES.get(name, []) if n in rev]


def destination(info, k, lineage, names, ranks):
    """returns (bucket, label)"""
    if k == 0:
        return "unclassified", "unclassified"
    lin = lineage(k)
    if info.t in lin or any(a in lin for a in info.alias):
        return "species", "species (correct)"
    nm = names.get(k, str(k))
    if k in info.idx:
        i = info.idx[k]
        if info.gi is not None and i <= info.gi:
            b = "up_genus"
        elif info.fi is not None and i <= info.fi:
            b = "up_family"
        else:
            b = "up_higher"
        return b, f"{b}: {nm} [{ranks.get(k, '?')}]"
    sp = next((x for x in lin if ranks.get(x) == "species"), k)
    spn = names.get(sp, str(sp))
    if info.g is not None and info.g in lin:
        return "sister", f"sister: {spn}"
    if info.f is not None and info.f in lin:
        return "same_family", f"same_family: {spn}"
    return "elsewhere", f"elsewhere: {spn}"


def build_reference(refdir, outfa):
    files = sorted(Path(refdir).glob("*.fasta"))
    if not files:
        sys.exit(f"ERROR: no *.fasta in {refdir}")
    sums = {}
    with open(outfa, "w") as out:
        for p in files:
            stem, n = p.stem, 0
            sums[p.name] = md5(p)
            with open(p) as fh:
                for line in fh:
                    if line.startswith(">"):
                        n += 1
                        out.write(f">{stem}__{n}\n")
                    else:
                        out.write(line if line.endswith("\n") else line + "\n")
    return sums


def docker_out(args):
    return subprocess.run(args, check=True, capture_output=True, text=True).stdout.strip()


def run_minimap2(image, fa, reads, paf, threads, mem):
    if subprocess.run(["docker", "image", "inspect", image], capture_output=True).returncode != 0:
        subprocess.run(["docker", "pull", image], check=True)
    digests = docker_out(["docker", "image", "inspect", image, "--format", "{{json .RepoDigests}}"])
    version = docker_out(["docker", "run", "--rm", image, "minimap2", "--version"])
    fa, reads, paf = Path(fa).resolve(), Path(reads).resolve(), Path(paf).resolve()
    cmd = ["docker", "run", "--rm", "-u", f"{os.getuid()}:{os.getgid()}", "-m", mem, "--cpus", str(threads),
           "-v", f"{fa.parent}:/ref:ro", "-v", f"{reads.parent}:/reads:ro", "-v", f"{paf.parent}:/out",
           image, "minimap2", "-x", "map-ont", "--secondary=no", "-t", str(threads),
           "-o", f"/out/{paf.name}", f"/ref/{fa.name}", f"/reads/{reads.name}"]
    t0 = time.time()
    subprocess.run(cmd, check=True)
    return {"image": image, "repo_digests": digests, "minimap2_version": version,
            "command": " ".join(cmd), "wall_s": round(time.time() - t0, 1)}


def load_species_table(path):
    exp = {}
    if path and Path(path).exists():
        with open(path) as fh:
            head = fh.readline().rstrip("\n").split("\t")
            ti, di = head.index("taxid"), head.index("dna_pct")
            for line in fh:
                f = line.rstrip("\n").split("\t")
                exp[int(f[ti])] = float(f[di])
    return exp


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pool", choices=["dev"], default="dev", help="dev only (exploratory)")
    ap.add_argument("--reads", default="data/subsets/SRR31222153.dev.seed20261012.blocks.fastq.gz")
    ap.add_argument("--refdir", default="data/reference/D6331.refseq/genomes")
    ap.add_argument("--kraken-run", default="data/classify/dev/baseline_c0")
    ap.add_argument("--db", default="data/db/pluspf_08_20260626", help="for nodes.dmp / names.dmp")
    ap.add_argument("--species-table", default="reports/baseline_kraken2_pluspf8/dev_species_table.tsv")
    ap.add_argument("--workdir", default="data/trace/dev")
    ap.add_argument("--outdir", default="reports/trace")
    ap.add_argument("--image", default="quay.io/biocontainers/minimap2:2.31--h118bc1c_0")
    ap.add_argument("--threads", type=int, default=6)
    ap.add_argument("--memory", default="8g")
    ap.add_argument("--paf", help="use an existing PAF (skips reference build and mapping)")
    ap.add_argument("--min-cov", type=float, default=0.5)
    ap.add_argument("--amb", type=float, default=0.9)
    a = ap.parse_args()

    work, out = Path(a.workdir), Path(a.outdir)
    work.mkdir(parents=True, exist_ok=True); out.mkdir(parents=True, exist_ok=True)
    run = {"script": "bin/08_trace_reads.py", "pool": a.pool, "args": vars(a),
           "started": time.strftime("%Y-%m-%d %H:%M:%S")}

    if a.paf:
        paf = Path(a.paf)
    else:
        fa = work / "ref_combined.fa"
        run["reference_md5"] = build_reference(a.refdir, fa)
        paf = work / "dev_vs_refs.paf"
        run["minimap2"] = run_minimap2(a.image, fa, a.reads, paf, a.threads, a.memory)

    kraken = read_kraken(Path(a.kraken_run) / "kraken_output.txt.gz")
    agg = read_paf(paf)
    extra = set(agg) - set(kraken)
    if extra:
        sys.exit(f"ERROR: {len(extra)} PAF read ids are not in the Kraken2 output (example {next(iter(extra))}); wrong pair of files?")
    parent, ranks = load_nodes(find_dmp(a.db, "nodes.dmp"))
    names, rev = load_names(find_dmp(a.db, "names.dmp"))
    lineage = make_lineage(parent)
    exp = load_species_table(a.species_table)

    labels = {}
    for rid in kraken:
        labels[rid] = assign(agg[rid], a.min_cov, a.amb) if rid in agg else "none"
    sp_names = sorted({l for l in labels.values() if l not in ("weak", "ambiguous", "none")})
    infos = {}
    for n in sp_names:
        if n not in rev:
            sys.exit(f"ERROR: species name '{n}' not found in names.dmp")
        infos[n] = SpInfo(n, rev[n], lineage, ranks, rev)

    N = len(kraken)
    status = Counter(l if l in ("weak", "ambiguous", "none") else "mapped" for l in labels.values())
    mat = {n: Counter() for n in sp_names}
    dests = {n: Counter() for n in sp_names}
    unm = Counter(); unm_class = Counter()
    for rid, k in kraken.items():
        l = labels[rid]
        if l in infos:
            b, lab = destination(infos[l], k, lineage, names, ranks)
            mat[l][b] += 1; dests[l][lab] += 1
        else:
            unm_class["unclassified" if k == 0 else "classified"] += 1
            if k:
                lin = lineage(k)
                sp = next((x for x in lin if ranks.get(x) == "species"), k)
                unm[names.get(sp, str(sp))] += 1
    mapped_total = sum(sum(c.values()) for c in mat.values())

    rows = []
    for n in sp_names:
        c = mat[n]; tot = sum(c.values()); info = infos[n]
        e = exp.get(info.t)
        share_all = 100.0 * tot / N
        share_map = 100.0 * tot / mapped_total if mapped_total else 0.0
        rows.append({"species": n, "taxid": info.t, "expected_dna_pct": e, "mapped_reads": tot,
                     "share_of_all_reads_pct": share_all, "share_of_mapped_pct": share_map,
                     "mapped_over_expected": (share_map / e) if e else None,
                     **{"to_" + b: c[b] for b in BUCKETS},
                     "pct_species": 100.0 * c["species"] / tot if tot else 0.0,
                     "pct_above_species": 100.0 * (c["up_genus"] + c["up_family"] + c["up_higher"]) / tot if tot else 0.0,
                     "pct_unclassified": 100.0 * c["unclassified"] / tot if tot else 0.0})
    rows.sort(key=lambda r: -(r["expected_dna_pct"] or 0))

    cols = ["species", "taxid", "expected_dna_pct", "mapped_reads", "share_of_all_reads_pct", "share_of_mapped_pct",
            "mapped_over_expected"] + ["to_" + b for b in BUCKETS] + ["pct_species", "pct_above_species", "pct_unclassified"]
    def fmt(v):
        return "" if v is None else (f"{v:.3f}" if isinstance(v, float) else str(v))
    with open(out / "trace_by_species.tsv", "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(fmt(r[c]) for c in cols) + "\n")
    with open(out / "trace_destinations.tsv", "w") as fh:
        fh.write("true_species\tdestination\treads\n")
        for n in sp_names:
            for lab, cnt in dests[n].most_common():
                if cnt >= 5:
                    fh.write(f"{n}\t{lab}\t{cnt}\n")
    with open(out / "trace_unmapped_kraken.tsv", "w") as fh:
        fh.write("kraken_species_of_unmapped_reads\treads\n")
        for lab, cnt in unm.most_common(40):
            fh.write(f"{lab}\t{cnt}\n")

    lines = [f"Trace of dev reads (minimap2 truth vs Kraken2 {a.kraken_run}); N = {N} reads",
             f"mapped to one species {status['mapped']} ({100*status['mapped']/N:.2f}%), ambiguous {status['ambiguous']}, "
             f"weak/low coverage {status['weak']}, no alignment {status['none']}  (min-cov {a.min_cov}, amb {a.amb})",
             f"unmapped/ambiguous reads: Kraken2 classified {unm_class['classified']}, unclassified {unm_class['unclassified']}", "",
             f"{'species':32s} {'exp%':>6s} {'map%':>6s} {'map/exp':>7s} {'reads':>7s} | {'species':>7s} {'genus':>6s} {'family':>6s} {'higher':>6s} {'sister':>6s} {'sameF':>6s} {'else':>6s} {'unclas':>6s}"]
    for r in rows:
        lines.append(f"{r['species'][:32]:32s} {fmt(r['expected_dna_pct']):>6s} {r['share_of_mapped_pct']:6.2f} "
                     f"{('%.2f' % r['mapped_over_expected']) if r['mapped_over_expected'] else '':>7s} {r['mapped_reads']:7d} | "
                     + " ".join(f"{r['to_' + b]:{w}d}" for b, w in zip(BUCKETS, (7, 6, 6, 6, 6, 6, 6, 6))))
    lines += ["", "Top Kraken2 species among unmapped/ambiguous reads: " + "; ".join(f"{k} {v}" for k, v in unm.most_common(8))]
    text = "\n".join(lines)
    (out / "trace_summary.txt").write_text(text + "\n")
    run["counts"] = dict(status); run["finished"] = time.strftime("%Y-%m-%d %H:%M:%S")
    (out / "run.json").write_text(json.dumps(run, indent=2, default=str) + "\n")
    print(text)


if __name__ == "__main__":
    main()
