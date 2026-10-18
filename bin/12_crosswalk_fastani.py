#!/usr/bin/env python3
"""12_crosswalk_fastani.py - crosswalk v2 for the Sylph arm (addendum to D18; v1 = registered rule, commit 58a48a0, tag sylph-crosswalk).

Why: v1 profiled each reference genome AS READS, which gives an effective coverage of only ~0.001-0.08, so several database species tie at
adjusted ANI 100 and 'highest ANI' picked one arbitrarily. v2 uses the criterion that DEFINES GTDB species: fastANI >= 95% with an aligned
fraction >= 0.5 between the mock reference genome and the database genome.

Method (reference genomes only decide; reads never do):
 1. Candidate database genomes = every accession that appears in the Genome_file column of (a) the reference-as-reads Sylph output
    (all hits, not only the top 5) and (b) the dev read profile (candidates only; it does not decide anything).
 2. Download them with NCBI `datasets` (data/crosswalk_v2/genomes/<accession>.fna).
 3. fastANI of each of the 21 reference genomes against each candidate.
 4. A database species is mapped to a mock species when ANY of its candidate genomes has ANI >= --min-ani (95) and aligned fraction >=
    --min-af (0.5) to ANY reference genome of that mock species. A database species may not map to two mock species (abort).
    A mock species with no such species is 'not_represented'.
Outputs: reference/d6331_gtdb_crosswalk_v2.csv (same first columns as v1), reports/sylph/crosswalk_v2_pairs.tsv (all fastANI pairs),
         reports/sylph/crosswalk_v2_candidates.tsv, reports/sylph/crosswalk_v2_run.json.  Standard library only (fastANI and datasets are run as tools).
"""
import argparse, csv, glob, json, os, re, shutil, subprocess, sys, time, zipfile
from collections import defaultdict
import importlib.util

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
P = lambda *x: os.path.join(ROOT, *x)
_spec = importlib.util.spec_from_file_location("cw10", P("bin", "10_sylph_crosswalk.py"))
cw = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(cw)
ACC = re.compile(r"GC[AF]_\d{9}\.\d+")


def candidates(hit_files, gtdb):
    """accession -> {'species':..., 'src':..., 'from': set(files)}"""
    out = {}
    for f in hit_files:
        head, rows = cw.read_sylph(P(f))
        gc = cw.find_col(head, "Genome_file")
        cc = next((h for h in head if h.lower() == "contig_name"), None)
        for r in rows:
            m = ACC.search(r[gc])
            if not m:
                continue
            sp, src, _ = cw.species_of(r, gc, cc, gtdb)
            d = out.setdefault(m.group(0), {"species": sp, "src": src, "from": set()})
            if sp and not d["species"]:
                d["species"], d["src"] = sp, src
            d["from"].add(os.path.basename(f))
    return out


def download(accs, gdir):
    os.makedirs(gdir, exist_ok=True)
    need = [a for a in accs if not os.path.exists(os.path.join(gdir, a + ".fna"))]
    if not need:
        return []
    lst = os.path.join(gdir, "accessions.txt"); open(lst, "w").write("\n".join(need) + "\n")
    zp = os.path.join(gdir, "download.zip")
    r = subprocess.run(["datasets", "download", "genome", "accession", "--inputfile", lst, "--include", "genome", "--filename", zp, "--no-progressbar"],
                       capture_output=True, text=True)
    if r.returncode != 0 or not os.path.exists(zp):
        sys.exit(f"datasets failed: {r.stderr[-800:]}")
    with zipfile.ZipFile(zp) as z:
        by_acc = defaultdict(list)
        for n in z.namelist():
            m = ACC.search(n)
            if m and n.endswith(".fna"):
                by_acc[m.group(0)].append(n)
        for a, names in by_acc.items():
            with open(os.path.join(gdir, a + ".fna"), "wb") as o:
                for n in sorted(names):
                    o.write(z.read(n))
    os.remove(zp)
    return [a for a in need if not os.path.exists(os.path.join(gdir, a + ".fna"))]


def run_fastani(refs, cands, outp, threads):
    ql, rl = outp + ".ql", outp + ".rl"
    open(ql, "w").write("\n".join(refs) + "\n"); open(rl, "w").write("\n".join(cands) + "\n")
    r = subprocess.run(["fastANI", "--ql", ql, "--rl", rl, "-t", str(threads), "-o", outp], capture_output=True, text=True)
    if r.returncode != 0:
        sys.exit(f"fastANI failed: {r.stderr[-800:]}")
    ver = subprocess.run(["fastANI", "--version"], capture_output=True, text=True)
    return (ver.stdout + ver.stderr).strip().splitlines()[0]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--refdir", default="data/reference/D6331.refseq/genomes")
    ap.add_argument("--hits", nargs="+", default=["data/sylph/crosswalk/refs_vs_dbs.tsv", "data/sylph/dev/arm_3a.tsv"])
    ap.add_argument("--gtdb-taxonomy", nargs="+", default=["data/db/gtdb_tax/bac120_taxonomy_r232.tsv.gz", "data/db/gtdb_tax/ar53_taxonomy_r232.tsv.gz"])
    ap.add_argument("--truth", default="reference/d6331_expected_reads_by_species.csv")
    ap.add_argument("--workdir", default="data/crosswalk_v2")
    ap.add_argument("--out", default="reference/d6331_gtdb_crosswalk_v2.csv")
    ap.add_argument("--min-ani", type=float, default=95.0)
    ap.add_argument("--min-af", type=float, default=0.5)
    ap.add_argument("--threads", type=int, default=6)
    ap.add_argument("--fastani-tsv", help="use an existing fastANI output (skips download and fastANI; for tests)")
    a = ap.parse_args()

    gtdb = cw.load_gtdb([P(x) for x in a.gtdb_taxonomy])
    cand = candidates(a.hits, gtdb)
    names, tids = {}, {}
    for t in csv.DictReader(open(P(a.truth))):
        rec = (t["ncbi_name_in_kraken_db"].strip(), int(t["species_taxid"]))
        tids[rec[1]] = rec[0]
        for k in (t["species"], t["ncbi_name_in_kraken_db"]):
            names[k.strip().lower()] = rec

    def mock_of(stem):
        for c in (cw.stem_species(stem), stem.replace("_", " ")):
            if c.lower() in names:
                return names[c.lower()]
        sys.exit(f"reference file {stem} matches no truth species")

    run = {"script": "bin/12_crosswalk_fastani.py", "args": vars(a), "started": time.strftime("%Y-%m-%d %H:%M:%S"), "n_candidates": len(cand)}
    refs = sorted(glob.glob(P(a.refdir, "*.fasta")))
    if not refs:
        sys.exit(f"no *.fasta in {a.refdir}")
    work = P(a.workdir); os.makedirs(work, exist_ok=True)
    if a.fastani_tsv:
        pairs_file = P(a.fastani_tsv)
    else:
        missing = download(sorted(cand), os.path.join(work, "genomes"))
        run["not_downloadable"] = missing
        if missing:
            print(f"WARNING: {len(missing)} candidate genomes could not be downloaded: {missing}", file=sys.stderr)
        pairs_file = os.path.join(work, "fastani_out.tsv")
        run["fastani_version"] = run_fastani(refs, [os.path.join(work, "genomes", x + ".fna") for x in sorted(cand) if os.path.exists(os.path.join(work, "genomes", x + ".fna"))],
                                             pairs_file, a.threads)
        run["datasets_version"] = (subprocess.run(["datasets", "--version"], capture_output=True, text=True).stdout or "").strip()

    best = {}   # (mock_taxid, db_species) -> record
    allpairs = []
    for line in open(pairs_file):
        f = line.rstrip("\n").split("\t")
        if len(f) < 5:
            continue
        q, rfile, ani, mm, tot = os.path.basename(f[0])[:-len(".fasta")], os.path.basename(f[1]), float(f[2]), int(f[3]), int(f[4])
        m = ACC.search(rfile)
        if not m or m.group(0) not in cand:
            sys.exit(f"fastANI reference {rfile} is not a known candidate accession")
        acc, af = m.group(0), mm / tot
        mock, tid = mock_of(q)
        sp = cand[acc]["species"] or "unresolved:" + acc
        allpairs.append((q, mock, acc, sp, ani, af))
        if ani >= a.min_ani and af >= a.min_af:
            d = best.setdefault((tid, sp), {"ani": 0.0, "af": 0.0, "acc": acc, "files": set()})
            d["files"].add(q)
            if ani > d["ani"]:
                d.update(ani=ani, af=af, acc=acc)
    owner = {}
    for (tid, sp) in best:
        if owner.setdefault(sp, tid) != tid:
            sys.exit(f"database species {sp} reaches >= {a.min_ani} ANI to two mock species ({tids[owner[sp]]}, {tids[tid]}); decide explicitly, not silently")

    os.makedirs(os.path.dirname(P(a.out)), exist_ok=True)
    with open(P(a.out), "w", newline="") as o:
        w = csv.writer(o)
        w.writerow(["mock_species", "mock_taxid", "database_species", "status", "fastani", "aligned_fraction", "genome_accession", "reference_files"])
        for tid in sorted(tids, key=lambda t: tids[t]):
            rows = sorted(((sp, d) for (t, sp), d in best.items() if t == tid), key=lambda x: -x[1]["ani"])
            if not rows:
                w.writerow([tids[tid], tid, "", "not_represented", "", "", "", ""])
            for sp, d in rows:
                w.writerow([tids[tid], tid, sp, "ok", f"{d['ani']:.2f}", f"{d['af']:.2f}", d["acc"], ";".join(sorted(d["files"]))])
    os.makedirs(P("reports", "sylph"), exist_ok=True)
    with open(P("reports", "sylph", "crosswalk_v2_pairs.tsv"), "w") as o:
        o.write("reference_file\tmock_species\tcandidate_accession\tdatabase_species\tfastani\taligned_fraction\n")
        for q, mock, acc, sp, ani, af in sorted(allpairs, key=lambda x: (x[1], -x[4])):
            o.write(f"{q}\t{mock}\t{acc}\t{sp}\t{ani:.2f}\t{af:.2f}\n")
    with open(P("reports", "sylph", "crosswalk_v2_candidates.tsv"), "w") as o:
        o.write("accession\tdatabase_species\tspecies_source\tfound_in\n")
        for acc, d in sorted(cand.items()):
            o.write(f"{acc}\t{d['species']}\t{d['src']}\t{','.join(sorted(d['from']))}\n")
    run["finished"] = time.strftime("%Y-%m-%d %H:%M:%S")
    run["not_represented"] = [tids[t] for t in tids if not any(tt == t for tt, _ in best)]
    json.dump(run, open(P("reports", "sylph", "crosswalk_v2_run.json"), "w"), indent=2, default=list)

    print(f"crosswalk v2: {len(cand)} candidate genomes, {len(allpairs)} fastANI pairs reported; thresholds ANI >= {a.min_ani}, aligned fraction >= {a.min_af}")
    for tid in sorted(tids, key=lambda t: tids[t]):
        rows = sorted(((sp, d) for (t, sp), d in best.items() if t == tid), key=lambda x: -x[1]["ani"])
        print(f"  {tids[tid]:32s} -> " + ("; ".join(f"{sp} (ANI {d['ani']:.2f}, AF {d['af']:.2f})" for sp, d in rows) if rows else "NOT REPRESENTED"))


if __name__ == "__main__":
    main()
