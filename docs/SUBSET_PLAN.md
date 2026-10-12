# Subset plan: PRJNA1181840 / SRR31222153 (updated 2026-10-12)

## Goal
Obtain reproducible read sets from a single 128.4 Gb run (23,429,059 reads, mean ~5.5 kb) without
downloading the 120 GB FASTQ, and keep development and evaluation reads separate.

## What we measured (from WSL, 2026-10-12)
| Check | Result |
|---|---|
| ENA FASTQ structure | one gzip stream, no random access; read IDs rewritten as `SRR31222153.n` |
| Transfer speed | ENA 0.47, NCBI/AWS 0.56, Cloudflare 0.62 MB/s: the connection is the limit |
| Full ENA stream | ~120 GB at 0.47 MB/s = ~3 days (rejected) |
| SRA read-range fetch | works: `fastq-dump -N a -X b`, 20,000 reads in 4 min 14 s (~83 reads/s), ~14 s fixed overhead per call |
| SRA vs ENA | first 200 reads identical (sequence + quality md5 66b4146b...) |

## Design: stratified block sampling via SRA (`bin/01_sra_block_sample.py`)
1. Cut read numbers 1..23,429,059 into 74 equal strata (~316,600 reads each).
2. In each stratum place one block of 5,000 consecutive reads at a seeded random position (seed 20261012).
3. Even strata -> development pool (37 blocks), odd strata -> evaluation pool (37 blocks).
   Each pool is 185,000 reads (~1 Gb). Both span the whole run; neighbouring reads always stay in the same pool.
4. Each block is its own gzip file in `data/subsets/blocks/`; finished blocks are skipped on rerun (resumable),
   and every block is checked (FASTQ structure, consecutive read IDs, read count) before it is accepted.
5. Pools are concatenations of blocks (valid gzip). A manifest records parameters, per-block and per-pool counts,
   bases, mean read length, SHA-256 of pools, and the pinned Docker image digest.

Docker image: `ncbi/sra-tools@sha256:4f048ba6754246e1ed7799a2c8ed469a38427f392e403ca08b93c4255104c9ce`.

## Time
74 blocks x (~60 s transfer + ~14 s overhead) = about 1.5 hours. Safe to stop and rerun.

## Running it (WSL)
```
cd ~/projects/metagenomics_1_Nanopore
bin/01_sra_block_sample.py --plan-only | head      # inspect the plan first
tmux new -s subset
bin/01_sra_block_sample.py                          # defaults: seed 20261012, 37 blocks/pool x 5,000 reads
```

## Depth series (later, in Nextflow)
From each pool, seeded downsampling to about 0.1, 0.25, 0.5 and 1 Gb with 5 replicate seeds per depth.
Development pool: tune thresholds only. Evaluation pool: report only, thresholds frozen first.

## Limits to state wherever results are reported
- Within-run held-out data only: both pools share library prep, run and chemistry. Not an independent experiment.
- Replicates drawn from one ~1 Gb pool overlap heavily at larger depths: they measure sampling variance only.
- ENA/SRA read IDs are serial numbers, so duplex sister reads cannot be identified. Keeping blocks together
  within a pool limits leakage from adjacent related reads but does not rule out all of it.
- Reads are sampled in blocks of consecutive reads, so within-block reads are not independent draws.
- Chemistry (Q20+, probably Kit 14) and basecaller are still to be confirmed from the SRA experiment or paper.

## After it runs: representativeness checks (record in PROJECT_STATE.md)
- Pool mean read length vs run mean 5,479 bp (from ENA counts); per-block mean length across the run (any drift?).
- `seqkit stats` / NanoPlot on each pool: length and quality distributions should match between pools.
