#!/usr/bin/env python3
"""19_amr_reference.py - AMR reference catalogue step 1 (ADDENDUM A3; reference genomes only, no reads, no pool is touched).

Runs AMRFinderPlus (nucleotide mode, --plus) once on each manufacturer reference genome of the Zymo D6331 archive and writes
  - data/amr/reference/<genome>.amrfinder.tsv   (raw AMRFinderPlus output, not committed)
  - reports/amr/reference_hits_raw.tsv          (all genomes combined, one row per hit, genome and organism flag added)
  - reports/amr/reference_run.json              (tool and database versions, command lines, input SHA-256, organism flags used)
No filtering and no category is applied here; categories are assigned by bin/20 using the rules registered in A3.

Nucleotide-only limitation (documented, not hidden): without a protein FASTA and GFF, AMRFinderPlus translates the genome in six frames
(BLASTX-style), so partial and fragmented hits can be reported differently from a protein-based run.
Point-mutation search needs --organism; it is passed only for organisms listed by `amrfinder --list_organisms` (checked at run time,
not hard-coded), otherwise no --organism is given and point mutations are NOT searched for that genome.
Standard library only. Needs amrfinder on PATH (conda install -c bioconda -c conda-forge ncbi-amrfinderplus; amrfinder -u).
"""
import argparse, csv, hashlib, json, os, re, subprocess, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
# genome file stem (first word) -> AMRFinderPlus organism keyword, used only if the keyword is in `amrfinder --list_organisms`
ORG = {"Escherichia": "Escherichia", "Salmonella": "Salmonella", "Enterococcus": "Enterococcus_faecalis",
       "Clostridioides": "Clostridioides_difficile"}
SKIP = ("Candida", "Saccharomyces")  # yeasts: not a bacterial AMR target (draft assemblies; AMRFinderPlus is bacterial)


def sh(cmd):
    return subprocess.run(cmd, capture_output=True, text=True)


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--refdir", default="data/reference/D6331.refseq/genomes")
    ap.add_argument("--out-raw", default="data/amr/reference")
    ap.add_argument("--out", default="reports/amr")
    ap.add_argument("--threads", type=int, default=4)
    a = ap.parse_args()

    v = sh(["amrfinder", "--version"])
    if v.returncode:
        sys.exit("amrfinder not found: conda install -c bioconda -c conda-forge ncbi-amrfinderplus; amrfinder -u")
    ver = sh(["amrfinder", "--database_version"]); orgs = set(re.findall(r"[A-Za-z_]+", sh(["amrfinder", "--list_organisms"]).stdout))
    files = sorted((ROOT / a.refdir).glob("*.fasta"))
    if not files:
        sys.exit(f"no *.fasta in {a.refdir}")
    (ROOT / a.out_raw).mkdir(parents=True, exist_ok=True); (ROOT / a.out).mkdir(parents=True, exist_ok=True)
    meta = {"amrfinder_version": v.stdout.strip(), "database_version": (ver.stdout + ver.stderr).strip(), "genomes": {}}
    rows, header = [], None
    for f in files:
        stem = f.stem
        if stem.startswith(SKIP):
            meta["genomes"][stem] = {"skipped": "yeast", "sha256": sha(f)}; continue
        key = ORG.get(stem.split("_")[0]); org = key if key in orgs else None
        out = ROOT / a.out_raw / f"{stem}.amrfinder.tsv"
        cmd = ["amrfinder", "-n", str(f), "--plus", "--threads", str(a.threads), "-o", str(out)] + (["--organism", org] if org else [])
        r = sh(cmd)
        if r.returncode:
            sys.exit(f"amrfinder failed on {f.name}: {r.stderr[-500:]}")
        meta["genomes"][stem] = {"sha256": sha(f), "organism_flag": org, "point_mutations_searched": bool(org), "cmd": " ".join(cmd)}
        with open(out) as fh:
            rd = csv.DictReader(fh, delimiter="\t"); header = header or rd.fieldnames
            for x in rd:
                x["genome"] = stem; x["organism_flag"] = org or ""; rows.append(x)
    with open(ROOT / a.out / "reference_hits_raw.tsv", "w") as o:
        w = csv.DictWriter(o, fieldnames=["genome", "organism_flag"] + (header or []), delimiter="\t", extrasaction="ignore")
        w.writeheader(); w.writerows(rows)
    json.dump(meta, open(ROOT / a.out / "reference_run.json", "w"), indent=1)
    print(f"{len(files)} genomes, {len(rows)} raw hits; {meta['amrfinder_version']}; {meta['database_version']}")


if __name__ == "__main__":
    main()
