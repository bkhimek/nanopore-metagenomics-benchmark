#!/usr/bin/env python3
"""Run Kraken2 on one read pool (or any FASTQ) in the pinned Docker image, with a memory cap and a run record.

Design rules encoded here:
- The evaluation pool is refused unless --allow-eval is given (thresholds must be frozen first; see PROJECT_STATE D5/D9).
- The image is pinned by digest, the database directory is recorded with the md5 of its archive, and every run writes run.json.
- Memory is capped with `docker --memory` so a run that does not fit the 12 GB budget fails instead of silently swapping (D4).
Outputs go to data/classify/<pool>/<label>/: report.txt, kraken_output.txt.gz, run.json, run.log
"""
import argparse, datetime, glob, gzip, json, os, shutil, subprocess, sys, time

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
DEFAULT_IMAGE = "quay.io/biocontainers/kraken2@sha256:db5d772d9df323cf0ddc644fb9793c0fca25b61f86469be42f1a9f5b903b9fdf"
DEFAULT_DB = "data/db/pluspf_08_20260626"
DB_ARCHIVE_MD5 = "79a153b99f045bc2ae95e6d57c17a02d"   # k2_pluspf_08_GB_20260626.tar.gz, verified 2026-10-15

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--pool", choices=["dev", "eval"], help="read pool from the sampling manifest")
ap.add_argument("--fastq", help="alternative input FASTQ(.gz), relative to the project root; overrides --pool")
ap.add_argument("--label", default="full", help="name of this run (default 'full'); used for the output folder")
ap.add_argument("--confidence", type=float, default=0.0, help="Kraken2 --confidence (default 0.0)")
ap.add_argument("--min-hit-groups", type=int, default=2, help="Kraken2 --minimum-hit-groups (default 2, the Kraken2 default)")
ap.add_argument("--threads", type=int, default=6)
ap.add_argument("--memory", default="10g", help="docker memory limit (default 10g; WSL cap is 12 GB)")
ap.add_argument("--db", default=DEFAULT_DB, help="database dir relative to project root (default %(default)s)")
ap.add_argument("--image", default=DEFAULT_IMAGE)
ap.add_argument("--allow-eval", action="store_true", help="required for the evaluation pool; use only after thresholds are frozen")
ap.add_argument("--dry-run", action="store_true", help="print the command and exit")
a = ap.parse_args()

if not a.pool and not a.fastq:
    sys.exit("give --pool dev|eval or --fastq PATH")
if a.pool == "eval" and not a.allow_eval:
    sys.exit("refusing to run on the evaluation pool without --allow-eval. Freeze thresholds on the dev pool first.")

data = os.path.join(ROOT, "data")
if a.fastq:
    inp = os.path.join(ROOT, a.fastq); pool_name = "custom"
else:
    mf = glob.glob(os.path.join(data, "subsets", "*.sra_blocks.manifest.json"))
    if len(mf) != 1:
        sys.exit(f"expected exactly one manifest in data/subsets, found {len(mf)}")
    man = json.load(open(mf[0])); inp = os.path.join(data, "subsets", man["pools"][a.pool]["file"]); pool_name = a.pool
if not os.path.exists(inp):
    sys.exit(f"input not found: {inp}")
if not inp.startswith(data + os.sep):
    sys.exit("input must be inside the project data/ folder (it is mounted into the container)")
dbdir = os.path.join(ROOT, a.db)
if not os.path.exists(os.path.join(dbdir, "hash.k2d")):
    sys.exit(f"database not found: {dbdir}/hash.k2d")

outdir = os.path.join(data, "classify", pool_name, a.label)
if os.path.exists(os.path.join(outdir, "run.json")):
    sys.exit(f"{outdir} already has a completed run; choose a different --label (results are not overwritten)")
os.makedirs(outdir, exist_ok=True)
rel = lambda p: "/data/" + os.path.relpath(p, data)
cmd = ["docker", "run", "--rm", "--memory", a.memory, "--memory-swap", a.memory, "-v", f"{data}:/data", a.image,
       "kraken2", "--db", rel(dbdir), "--threads", str(a.threads), "--gzip-compressed",
       "--confidence", str(a.confidence), "--minimum-hit-groups", str(a.min_hit_groups), "--report-minimizer-data",
       "--report", rel(os.path.join(outdir, "report.txt")), "--output", rel(os.path.join(outdir, "kraken_output.txt")), rel(inp)]
print(" ".join(cmd))
if a.dry_run:
    sys.exit(0)

t0 = time.time(); start = datetime.datetime.now().isoformat(timespec="seconds")
with open(os.path.join(outdir, "run.log"), "w") as lg:
    rc = subprocess.run(cmd, stdout=lg, stderr=subprocess.STDOUT).returncode
wall = time.time() - t0
log = open(os.path.join(outdir, "run.log")).read()
if rc != 0:
    print(log[-2000:]); sys.exit(f"kraken2 failed (exit {rc}); 137 usually means the memory limit was hit")
ko = os.path.join(outdir, "kraken_output.txt")
with open(ko, "rb") as f, gzip.open(ko + ".gz", "wb", compresslevel=6) as g:
    shutil.copyfileobj(f, g)
os.remove(ko)
rec = {"started": start, "wall_seconds": round(wall, 1), "pool": pool_name, "label": a.label, "input": os.path.relpath(inp, ROOT),
       "input_sha256_in_manifest": (man["pools"][a.pool]["sha256"] if not a.fastq else None),
       "image": a.image, "database": a.db, "database_archive_md5": DB_ARCHIVE_MD5 if a.db == DEFAULT_DB else None,
       "options": {"confidence": a.confidence, "minimum_hit_groups": a.min_hit_groups, "threads": a.threads, "memory_limit": a.memory,
                   "report_minimizer_data": True}, "command": " ".join(cmd), "kraken_summary": [l.strip() for l in log.splitlines() if "sequences" in l]}
json.dump(rec, open(os.path.join(outdir, "run.json"), "w"), indent=2)
print(log.strip()); print(f"wall {wall:.0f}s; outputs in {os.path.relpath(outdir, ROOT)}/")
