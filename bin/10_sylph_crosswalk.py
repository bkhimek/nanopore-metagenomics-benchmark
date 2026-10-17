#!/usr/bin/env python3
"""10_sylph_crosswalk.py - build the taxonomy crosswalk for the Sylph arm (D18) from the REFERENCE GENOMES ONLY (no reads).

Registered rule (PROJECT_STATE.md, 'Sylph arm'): each of the 21 D6331 genome FASTA files is profiled (as if it were a read set) against the same
databases that will be used for the reads; the mock species is mapped to the GTDB (or fungal) species of its highest-ANI database genome
(ties: higher effective coverage).  A reference file with no hit at adjusted ANI >= --min-ani (default 95) is 'not_represented'.  GTDB species names
come from the GTDB taxonomy files (accession -> s__ name); for hits that are not in GTDB (the fungal database) the binomial is taken from the
Sylph 'Contig_name' column.  The result is committed BEFORE any read is profiled and is not changed after seeing read results (same logic as D13).

Outputs:  reference/d6331_gtdb_crosswalk.csv   (one row per mock species and database species)
          reports/sylph/crosswalk_evidence.tsv (top hits per reference file)
          reports/sylph/crosswalk_run.json     (command, image digest, peak container memory, wall time)
Standard library only; Sylph runs in a Docker image.
"""
import argparse, csv, gzip, json, os, re, subprocess, sys, threading, time
from collections import defaultdict
from pathlib import Path

ACC = re.compile(r"(GC[AF])_(\d{9})\.(\d+)")


def stem_species(stem):
    if stem.startswith("Escherichia_coli"):
        return "Escherichia coli"
    special = {"Candida_albican": "Candida albicans", "Lactobacillus_fermentum": "Limosilactobacillus fermentum"}
    return special.get(stem, stem.replace("_", " "))


def acc_of(text):
    m = ACC.search(text or "")
    return f"{m.group(1)}_{m.group(2)}.{m.group(3)}" if m else None


def load_gtdb(paths):
    m = {}
    for p in paths:
        op = gzip.open if str(p).endswith(".gz") else open
        with op(p, "rt") as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) < 2:
                    continue
                acc = f[0][3:] if f[0][:3] in ("RS_", "GB_") else f[0]
                sp = [x for x in f[1].split(";") if x.startswith("s__")]
                if sp:
                    m[acc] = sp[0][3:]
    return m


def read_sylph(path):
    with open(path) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        rows = [dict(zip(head, l.rstrip("\n").split("\t"))) for l in fh if l.strip()]
    return head, rows


def find_col(head, *names):
    low = {h.lower(): h for h in head}
    for n in names:
        if n.lower() in low:
            return low[n.lower()]
    sys.exit(f"ERROR: none of the columns {names} in the Sylph output header: {head}")


def sample_stem(name):
    b = os.path.basename(name)
    for ext in (".gz", ".fasta", ".fa", ".fna", ".fastq", ".fq"):
        if b.endswith(ext):
            b = b[: -len(ext)]
    return b


def species_of(row, genome_col, contig_col, gtdb):
    g = row.get(genome_col, "")
    acc = acc_of(g)
    if acc and acc in gtdb:
        return gtdb[acc], "gtdb", acc
    if contig_col and row.get(contig_col):
        words = re.sub(r"^\S+\s+", "", row[contig_col]).split()
        if len(words) >= 2:
            return " ".join(words[:2]), "contig_name", acc or g
    return None, "unresolved", acc or g


def mem_to_gib(txt):
    m = re.match(r"([\d.]+)\s*([KMGT]?i?B)", txt.strip())
    if not m:
        return 0.0
    v, u = float(m.group(1)), m.group(2)
    return v * {"B": 1 / 2**30, "KiB": 1 / 2**20, "MiB": 1 / 2**10, "GiB": 1, "TiB": 1024,
                "kB": 1e3 / 2**30, "MB": 1e6 / 2**30, "GB": 1e9 / 2**30}.get(u, 0)


def run_sylph(image, dbs, refs, out_tsv, threads, mem, sample_threads):
    if subprocess.run(["docker", "image", "inspect", image], capture_output=True).returncode != 0:
        subprocess.run(["docker", "pull", image], check=True)
    digests = subprocess.run(["docker", "image", "inspect", image, "--format", "{{json .RepoDigests}}"],
                             check=True, capture_output=True, text=True).stdout.strip()
    dbdir = {Path(d).resolve().parent for d in dbs}
    if len(dbdir) != 1:
        sys.exit("ERROR: all database files must be in one directory")
    dbdir = dbdir.pop()
    refdir = Path(refs[0]).resolve().parent
    out_tsv = Path(out_tsv).resolve(); out_tsv.parent.mkdir(parents=True, exist_ok=True)
    name = f"sylph_cw_{os.getpid()}"
    cmd = ["docker", "run", "--rm", "--name", name, "-u", f"{os.getuid()}:{os.getgid()}", "-m", mem, "--cpus", str(threads),
           "-v", f"{dbdir}:/db:ro", "-v", f"{refdir}:/ref:ro", "-v", f"{out_tsv.parent}:/out", image,
           "sylph", "profile", "-t", str(threads), "-s", str(sample_threads),
           "-d"] + [f"/db/{Path(d).name}" for d in dbs] + ["-r"] + [f"/ref/{Path(r).name}" for r in refs] + ["-o", f"/out/{out_tsv.name}"]
    peak, stop = [0.0], threading.Event()

    def poll():
        while not stop.is_set():
            r = subprocess.run(["docker", "stats", "--no-stream", "--format", "{{.MemUsage}}", name], capture_output=True, text=True)
            if r.returncode == 0 and r.stdout.strip():
                peak[0] = max(peak[0], mem_to_gib(r.stdout.split("/")[0]))
            stop.wait(2)
    t0 = time.time()
    th = threading.Thread(target=poll, daemon=True); th.start()
    rc = subprocess.run(cmd).returncode
    stop.set(); th.join(timeout=5)
    info = {"image": image, "repo_digests": digests, "command": " ".join(cmd), "returncode": rc,
            "wall_s": round(time.time() - t0, 1), "peak_container_mem_gib_polled": round(peak[0], 2), "memory_cap": mem}
    if rc != 0:
        print("sylph exit code", rc, "(137 usually means the memory cap was hit)", file=sys.stderr)
    return info, rc


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--refdir", default="data/reference/D6331.refseq/genomes")
    ap.add_argument("--dbs", nargs="+", default=["data/db/sylph/gtdb-r232-c200-dbv2.syl2db", "data/db/sylph/fungi-refseq-2024-07-25-c200-v0.3.syldb"])
    ap.add_argument("--gtdb-taxonomy", nargs="+", required=True, help="GTDB r232 bac120/ar53 taxonomy TSV(.gz) files")
    ap.add_argument("--image", default="quay.io/biocontainers/sylph:1.0.0--hb42e459_0")
    ap.add_argument("--tsv", help="use an existing Sylph output instead of running Sylph")
    ap.add_argument("--workdir", default="data/sylph/crosswalk")
    ap.add_argument("--outdir", default="reports/sylph")
    ap.add_argument("--crosswalk", default="reference/d6331_gtdb_crosswalk.csv")
    ap.add_argument("--truth", default="reference/d6331_expected_reads_by_species.csv")
    ap.add_argument("--threads", type=int, default=6)
    ap.add_argument("--sample-threads", type=int, default=2)
    ap.add_argument("--memory", default="10g")
    ap.add_argument("--min-ani", type=float, default=95.0)
    a = ap.parse_args()

    refs = sorted(Path(a.refdir).glob("*.fasta"))
    if not refs:
        sys.exit(f"ERROR: no *.fasta in {a.refdir}")
    Path(a.outdir).mkdir(parents=True, exist_ok=True)
    run = {"script": "bin/10_sylph_crosswalk.py", "args": vars(a), "started": time.strftime("%Y-%m-%d %H:%M:%S"),
           "n_reference_files": len(refs)}
    if a.tsv:
        tsv = Path(a.tsv)
    else:
        tsv = Path(a.workdir) / "refs_vs_dbs.tsv"
        info, rc = run_sylph(a.image, a.dbs, refs, tsv, a.threads, a.memory, a.sample_threads)
        run["sylph"] = info
        if rc != 0:
            (Path(a.outdir) / "crosswalk_run.json").write_text(json.dumps(run, indent=2) + "\n")
            sys.exit(f"ERROR: sylph failed with exit code {rc}; see {a.outdir}/crosswalk_run.json")
    head, rows = read_sylph(tsv)
    sc = find_col(head, "Sample_file"); gc = find_col(head, "Genome_file")
    ac = find_col(head, "Adjusted_ANI"); ec = find_col(head, "Eff_cov")
    cc = next((h for h in head if h.lower() == "contig_name"), None)
    run["sylph_header"] = head
    gtdb = load_gtdb(a.gtdb_taxonomy)
    if not gtdb:
        sys.exit("ERROR: no species names read from the GTDB taxonomy files")

    # truth file = the registered list of mock species; its names and taxids are authoritative
    names, truth_ids = {}, {}
    with open(a.truth) as fh:
        rd = csv.DictReader(fh)
        tcol = next(c for c in rd.fieldnames if "taxid" in c.lower())
        ncol = next(c for c in rd.fieldnames if c.lower() in ("species", "species_name", "name"))
        kcol = next((c for c in rd.fieldnames if "ncbi_name" in c.lower()), None)
        for r in rd:
            rec = (r[ncol].strip(), r[tcol].strip())
            truth_ids[rec[1]] = rec[0]
            for k in (r[ncol], r[kcol] if kcol else ""):
                if k.strip():
                    names[k.strip().lower()] = rec

    def resolve(stem):
        for cand in (stem_species(stem), stem.replace("_", " ")):
            if cand.lower() in names:
                return names[cand.lower()]
        sys.exit(f"ERROR: reference file {stem} does not match any species in {a.truth}")

    by_sample = {}
    for r in rows:
        by_sample.setdefault(sample_stem(r[sc]), []).append(r)
    cw = defaultdict(dict)   # mock species -> db species -> record
    notrep, evidence = [], []
    for ref in refs:
        stem = ref.stem
        mock, mock_tid = resolve(stem)
        hits = sorted(by_sample.get(stem, []), key=lambda r: (-float(r[ac]), -float(r[ec])))
        for r in hits[:5]:
            sp, src, acc = species_of(r, gc, cc, gtdb)
            evidence.append([stem, mock, sp or "", src, acc or "", r[ac], r[ec], r.get(cc, "") if cc else ""])
        top = hits[0] if hits and float(hits[0][ac]) >= a.min_ani else None
        if top is None:
            notrep.append((mock, stem))
            continue
        sp, src, acc = species_of(top, gc, cc, gtdb)
        if sp is None:
            sys.exit(f"ERROR: cannot resolve a species name for the top hit of {stem}: {top}")
        d = cw[mock].setdefault(sp, {"ani": float(top[ac]), "eff_cov": float(top[ec]), "acc": acc, "src": src, "files": []})
        d["files"].append(stem)
        d["ani"] = max(d["ani"], float(top[ac]))

    taxid = {v: k for k, v in truth_ids.items()}
    covered = {taxid[m] for m in list(cw) + [m for m, _ in notrep]}
    if covered != set(truth_ids):
        sys.exit(f"ERROR: crosswalk covers {len(covered)} of {len(truth_ids)} truth species; missing: "
                 f"{sorted(truth_ids[t] for t in set(truth_ids) - covered)}")

    Path(a.crosswalk).parent.mkdir(parents=True, exist_ok=True)
    with open(a.crosswalk, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["mock_species", "mock_taxid", "database_species", "status", "adjusted_ani", "eff_cov", "genome_accession", "source", "reference_files"])
        for mock in sorted(set(list(cw) + [m for m, _ in notrep])):
            tid = taxid[mock]
            if mock in cw:
                for sp, d in sorted(cw[mock].items()):
                    w.writerow([mock, tid, sp, "ok", f"{d['ani']:.2f}", f"{d['eff_cov']:.3f}", d["acc"], d["src"], ";".join(d["files"])])
            else:
                w.writerow([mock, tid, "", "not_represented", "", "", "", "", ";".join(s for m, s in notrep if m == mock)])
    with open(Path(a.outdir) / "crosswalk_evidence.tsv", "w") as fh:
        fh.write("reference_file\tmock_species\tdatabase_species\tsource\taccession\tadjusted_ani\teff_cov\tcontig_name\n")
        for e in evidence:
            fh.write("\t".join(e) + "\n")
    run["finished"] = time.strftime("%Y-%m-%d %H:%M:%S")
    run["not_represented"] = [m for m, _ in notrep]
    (Path(a.outdir) / "crosswalk_run.json").write_text(json.dumps(run, indent=2) + "\n")

    print(f"crosswalk: {len(cw)} mock species represented, {len(set(m for m, _ in notrep) - set(cw))} not represented  (min ANI {a.min_ani})")
    for mock in sorted(cw):
        print(f"  {mock:32s} -> " + "; ".join(f"{sp} (ANI {d['ani']:.1f}, {d['src']}, {len(d['files'])} ref file(s))" for sp, d in sorted(cw[mock].items())))
    for m, s in notrep:
        print(f"  {m:32s} -> NOT REPRESENTED ({s})")


if __name__ == "__main__":
    main()
