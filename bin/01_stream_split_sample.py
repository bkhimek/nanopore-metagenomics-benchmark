#!/usr/bin/env python3
"""Single-pass, order-independent, seeded read split of a FASTQ stream.

Reads an UNCOMPRESSED FASTQ from stdin. For every read a deterministic value
u in [0,1) is derived from hash(seed + read_id). Reads with u < dev_frac go to
the development pool; reads with dev_frac <= u < dev_frac + eval_frac go to the
evaluation pool; all other reads are discarded. The pools are therefore disjoint
by construction, and the selection does not depend on read order, so it is
reproducible from (seed, fractions, read IDs) alone.

Usage:
  pigz -dc run.fastq.gz | 01_stream_split_sample.py --seed 20261012 \
      --dev-frac 0.0156 --eval-frac 0.0156 --outdir subsets --prefix SRR31222153
"""
import argparse, gzip, hashlib, json, os, sys, time

ap = argparse.ArgumentParser()
ap.add_argument("--seed", type=int, required=True)
ap.add_argument("--dev-frac", type=float, required=True)
ap.add_argument("--eval-frac", type=float, required=True)
ap.add_argument("--outdir", required=True)
ap.add_argument("--prefix", required=True)
ap.add_argument("--progress-every", type=int, default=1_000_000)
a = ap.parse_args()
if a.dev_frac < 0 or a.eval_frac < 0 or a.dev_frac + a.eval_frac > 1:
    sys.exit("fractions must be >=0 and sum to <=1")

os.makedirs(a.outdir, exist_ok=True)
paths = {k: os.path.join(a.outdir, f"{a.prefix}.{k}.seed{a.seed}.fastq.gz") for k in ("dev", "eval")}
tmp = {k: v + ".tmp" for k, v in paths.items()}
out = {k: gzip.open(v, "wb", compresslevel=4) for k, v in tmp.items()}
seed_b = str(a.seed).encode() + b"|"
two64 = float(1 << 64)
cut_dev, cut_eval = a.dev_frac, a.dev_frac + a.eval_frac
stats = {"total_reads": 0, "total_bases": 0,
         "dev": {"reads": 0, "bases": 0}, "eval": {"reads": 0, "bases": 0}}
t0 = time.time()
rd = sys.stdin.buffer.readline
while True:
    h = rd()
    if not h:
        break
    s, p, q = rd(), rd(), rd()
    if not q or h[:1] != b"@" or p[:1] != b"+":
        sys.exit(f"malformed FASTQ record after {stats['total_reads']} reads")
    n = len(s) - 1
    stats["total_reads"] += 1
    stats["total_bases"] += n
    rid = h[1:].split(None, 1)[0]
    u = int.from_bytes(hashlib.blake2b(seed_b + rid, digest_size=8).digest(), "big") / two64
    k = "dev" if u < cut_dev else ("eval" if u < cut_eval else None)
    if k:
        out[k].write(h + s + p + q)
        stats[k]["reads"] += 1
        stats[k]["bases"] += n
    if stats["total_reads"] % a.progress_every == 0:
        print(f"[{time.time()-t0:7.0f}s] reads={stats['total_reads']:,} "
              f"dev={stats['dev']['reads']:,} eval={stats['eval']['reads']:,}",
              file=sys.stderr, flush=True)
for f in out.values():
    f.close()
for k in out:
    os.replace(tmp[k], paths[k])
stats.update(seed=a.seed, dev_frac=a.dev_frac, eval_frac=a.eval_frac,
             outputs=paths, elapsed_s=round(time.time() - t0, 1))
json.dump(stats, sys.stdout, indent=2)
print()
