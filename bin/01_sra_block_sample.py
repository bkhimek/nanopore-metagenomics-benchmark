#!/usr/bin/env python3
"""Seeded, stratified block sampling of reads from an SRA run, by read-number range.

Why: the ENA FASTQ is one gzip stream (no random access) and this connection is slow (~0.5 MB/s), so
the full 120 GB cannot be streamed. SRA lets us fetch reads by number (fastq-dump -N/-X), so we take
blocks of consecutive reads spread across the whole run and download only those.

Design (all deterministic from --seed):
  * The read-number range 1..TOTAL is cut into 2*BLOCKS_PER_POOL equal strata.
  * One block of BLOCK_READS consecutive reads is placed at a random position inside each stratum.
  * Even-numbered strata -> development pool, odd-numbered strata -> evaluation pool.
    Pools are disjoint, both span the whole run, and neighbouring reads (which may be related, e.g.
    split reads, same pore state) always stay in the same pool.
  * Each block is downloaded to its own gzip file (resumable: finished blocks are skipped on rerun),
    scanned (read count, consecutive read IDs, FASTQ structure, bases) and then concatenated into pools
    (concatenated gzip members are a valid gzip file). A manifest records everything.

Usage:
  bin/01_sra_block_sample.py --plan-only
  bin/01_sra_block_sample.py            # run it (in tmux; safe to stop and rerun)
"""
import argparse, csv, gzip, hashlib, json, os, random, shlex, subprocess, sys, time
import signal
signal.signal(signal.SIGPIPE, signal.SIG_DFL)   # allow `... --plan-only | head` without a traceback

IMAGE = "ncbi/sra-tools@sha256:4f048ba6754246e1ed7799a2c8ed469a38427f392e403ca08b93c4255104c9ce"
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

ap = argparse.ArgumentParser()
ap.add_argument("--run", default="SRR31222153")
ap.add_argument("--project", default="PRJNA1181840")
ap.add_argument("--total-reads", type=int, help="default: read_count from metadata/<project>_read_run.tsv")
ap.add_argument("--seed", type=int, default=20261012)
ap.add_argument("--block-reads", type=int, default=5000)
ap.add_argument("--blocks-per-pool", type=int, default=37)
ap.add_argument("--outdir", default=os.path.join(ROOT, "data", "subsets"))
ap.add_argument("--image", default=IMAGE)
ap.add_argument("--fastq-dump-cmd", help="override the command (for testing); default runs fastq-dump in Docker")
ap.add_argument("--retries", type=int, default=3)
ap.add_argument("--plan-only", action="store_true")
a = ap.parse_args()

total = a.total_reads
if total is None:
    meta = os.path.join(ROOT, "metadata", f"{a.project}_read_run.tsv")
    with open(meta) as fh:
        rows = [r for r in csv.DictReader(fh, delimiter="\t") if r["run_accession"] == a.run]
    if not rows:
        sys.exit(f"{a.run} not found in {meta}")
    total = int(rows[0]["read_count"])

S = 2 * a.blocks_per_pool
L = total // S
if L < a.block_reads:
    sys.exit("strata smaller than block size: reduce --block-reads or --blocks-per-pool")
rng = random.Random(a.seed)
plan = []
for i in range(S):
    lo = i * L + 1
    hi = (i + 1) * L if i < S - 1 else total
    start = rng.randint(lo, hi - a.block_reads + 1)
    plan.append({"stratum": i, "pool": "dev" if i % 2 == 0 else "eval",
                 "start": start, "end": start + a.block_reads - 1})

print(f"run {a.run}: {total:,} reads; {S} strata of ~{L:,}; block {a.block_reads:,} reads; "
      f"{a.blocks_per_pool} blocks/pool = {a.blocks_per_pool * a.block_reads:,} reads/pool; seed {a.seed}")
if a.plan_only:
    for b in plan:
        print(f"{b['stratum']:3d} {b['pool']:4s} {b['start']:>10,} - {b['end']:>10,}")
    sys.exit(0)

blockdir = os.path.join(a.outdir, "blocks")
os.makedirs(blockdir, exist_ok=True)
base_cmd = shlex.split(a.fastq_dump_cmd) if a.fastq_dump_cmd else ["docker", "run", "--rm", a.image, "fastq-dump"]
logf = open(os.path.join(a.outdir, f"{a.run}.seed{a.seed}.fetch.log"), "a")


def bpath(b):
    return os.path.join(blockdir, f"{a.run}.{b['start']}-{b['end']}.fastq.gz")


def scan(path, b):
    n = bases = 0
    with gzip.open(path, "rb") as f:
        while True:
            h = f.readline()
            if not h:
                break
            s, p, q = f.readline(), f.readline(), f.readline()
            if not q or h[:1] != b"@" or p[:1] != b"+" or len(s) != len(q):
                raise ValueError(f"malformed FASTQ record {n} in {path}")
            rid = h[1:].split(None, 1)[0].decode()
            if rid != f"{a.run}.{b['start'] + n}":
                raise ValueError(f"read id {rid} != expected {a.run}.{b['start'] + n}")
            n += 1
            bases += len(s) - 1
    if n != b["end"] - b["start"] + 1:
        raise ValueError(f"{path}: {n} reads, expected {b['end'] - b['start'] + 1}")
    return n, bases


def fetch(b):
    path = bpath(b)
    if os.path.exists(path):
        return "cached"
    tmp = path + ".tmp"
    cmd = base_cmd + ["-N", str(b["start"]), "-X", str(b["end"]), "--stdout", a.run]
    with gzip.open(tmp, "wb", compresslevel=4) as out:
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=logf)
        for chunk in iter(lambda: p.stdout.read(1 << 20), b""):
            out.write(chunk)
        rc = p.wait()
    if rc != 0:
        os.remove(tmp)
        raise RuntimeError(f"fastq-dump exit {rc}")
    scan(tmp, b)               # structure + consecutive IDs + count, before accepting the block
    os.replace(tmp, path)
    return "downloaded"


t0 = time.time()
done_dl = 0
for k, b in enumerate(plan, 1):
    for attempt in range(1, a.retries + 1):
        try:
            tb = time.time()
            status = fetch(b)
            break
        except Exception as e:
            print(f"  block {k} attempt {attempt} failed: {e}", file=sys.stderr, flush=True)
            if attempt == a.retries:
                sys.exit(f"giving up on block {k}; rerun to resume (finished blocks are kept)")
            time.sleep(10 * attempt)
    if status == "downloaded":
        done_dl += 1
    el = time.time() - t0
    left = len(plan) - k
    eta = (el / done_dl * left / 60) if done_dl else 0
    print(f"[{k:3d}/{len(plan)}] {b['pool']:4s} {b['start']:>10,}-{b['end']:<10,} {status:10s} "
          f"{time.time() - tb:5.0f}s  elapsed {el / 60:5.1f} min  eta ~{eta:5.1f} min", flush=True)

# ---- assemble pools and write manifest ----
manifest = {"run": a.run, "project": a.project, "total_reads": total, "seed": a.seed,
            "block_reads": a.block_reads, "blocks_per_pool": a.blocks_per_pool, "image": a.image,
            "rng": "python random.Random(seed).randint, one start per stratum", "pools": {}, "blocks": []}
for pool in ("dev", "eval"):
    blocks = [b for b in plan if b["pool"] == pool]
    out = os.path.join(a.outdir, f"{a.run}.{pool}.seed{a.seed}.blocks.fastq.gz")
    reads = bases = 0
    with open(out + ".tmp", "wb") as o:
        for b in blocks:
            try:
                n, nb = scan(bpath(b), b)
            except ValueError as e:
                sys.exit(f"{e}\nDelete {bpath(b)} and rerun to re-fetch just that block.")
            reads += n
            bases += nb
            with open(bpath(b), "rb") as f:
                o.write(f.read())
            manifest["blocks"].append({**b, "reads": n, "bases": nb, "file": os.path.basename(bpath(b))})
    os.replace(out + ".tmp", out)
    h = hashlib.sha256()
    with open(out, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 22), b""):
            h.update(chunk)
    manifest["pools"][pool] = {"file": os.path.basename(out), "blocks": len(blocks), "reads": reads,
                               "bases": bases, "mean_read_length": round(bases / reads, 1), "sha256": h.hexdigest()}
mpath = os.path.join(a.outdir, f"{a.run}.seed{a.seed}.sra_blocks.manifest.json")
with open(mpath, "w") as f:
    json.dump(manifest, f, indent=2)
for pool, d in manifest["pools"].items():
    print(f"{pool}: {d['reads']:,} reads, {d['bases'] / 1e9:.3f} Gb, mean length {d['mean_read_length']:,.0f}")
print(f"manifest: {mpath}")
