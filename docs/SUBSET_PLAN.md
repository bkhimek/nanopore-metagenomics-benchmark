# Subset plan: PRJNA1181840 / SRR31222153 (updated 2026-10-13)

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
Estimated 74 blocks x (~60 s transfer + ~14 s overhead) = about 1.5 hours. Actual: 168.6 minutes (2026-10-12). Safe to stop and rerun.

## Running it (WSL)
```
cd ~/projects/metagenomics_1_Nanopore
bin/01_sra_block_sample.py --plan-only | head      # inspect the plan first
tmux new -s subset
bin/01_sra_block_sample.py                          # defaults: seed 20261012, 37 blocks/pool x 5,000 reads
```

## Depth series (later, in Nextflow)
Depth is defined by Gb. From each pool, seeded downsampling to about 0.1, 0.25 and 0.5 Gb with 5 replicate seeds per depth.
The full pool (~1 Gb) is a single baseline without replicates: a full-pool replicate would be identical.
Development pool: tune thresholds only. Evaluation pool: report only, thresholds frozen first.

## Limits to state wherever results are reported
- Within-run held-out data only: both pools share library prep, run and chemistry. Not an independent experiment.
- Replicates drawn from one ~1 Gb pool overlap heavily (about half of each replicate overlaps another at 0.5 Gb): they
  measure within-pool sampling variation only, not run-to-run variation, and are not independent.
- ENA/SRA read IDs are serial numbers, so duplex sister reads cannot be identified. Keeping blocks together
  within a pool limits leakage from adjacent related reads but does not rule out all of it.
- Reads are sampled in blocks of consecutive reads, so within-block reads are not independent draws. Confidence
  intervals must use block-level resampling (bootstrap over the 37 blocks per pool), not read-level.
- Chemistry and basecaller are still to be confirmed from the SRA experiment or paper. "Q20+" in the study title suggests
  Kit 14 chemistry, but that is an inference, not a verified fact.

## Result (2026-10-12; details in PROJECT_STATE.md)
- Finished in 168.6 min with no failure or retry lines. Dev pool 185,000 reads, 1.016 Gb, mean 5,491 bp; eval pool
  185,000 reads, 1.012 Gb, mean 5,473 bp (run mean 5,479 bp). Mean length drift across the run about -30 bp (-0.5%).
- Read quality and pool similarity not yet checked.

## Representativeness checks (done for length; quality still to do)
- Pool mean read length vs run mean 5,479 bp (from ENA counts); per-block mean length across the run (any drift?).
- `seqkit stats` / NanoPlot on each pool: length and quality distributions should match between pools.
