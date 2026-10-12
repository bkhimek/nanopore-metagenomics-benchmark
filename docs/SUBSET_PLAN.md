# Subset plan: PRJNA1181840 / SRR31222153 (draft, 2026-10-12)

## Goal
Obtain reproducible, manageable read sets from a single 128.4 Gb run (23.43 M reads, mean ~5.5 kb,
~120 GB compressed FASTQ) without storing the full file, and keep development and evaluation reads separate.

## Design
1. **Stream once, never store.** `curl | tee(md5) | pigz -dc | 01_stream_split_sample.py`.
   Peak RAM is small; disk use is only the two output pools.
2. **Read-level seeded split.** Each read gets u = hash(seed + read_id) in [0,1). u < 0.0156 -> development pool;
   0.0156 <= u < 0.0312 -> evaluation pool; otherwise discarded. Fractions of reads, about 365,000 reads (~2.0 Gb) per pool.
   - Pools are disjoint by construction.
   - Selection does not depend on read order, so it is reproducible from the seed alone (tested: shuffled input gives identical pools).
   - Not a prefix of the file, so no pore/time-order bias from taking "the first N reads".
3. **Integrity.** The MD5 of the downloaded stream is compared with ENA's `fastq_md5`; a mismatch writes a `.FAILED` marker.
   A manifest (counts, bases, seed, fractions) and SHA-256 of the outputs are written next to the pools.
4. **Depth series (done later, in Nextflow).** From each pool, seeded downsampling (for example `seqkit sample -s <seed> -p <fraction>`)
   to about 0.25, 0.5, 1 and 2 Gb, with 5 replicate seeds per depth. Development pool: tune thresholds only.
   Evaluation pool: report results only, with thresholds frozen first.

## What this does and does not give you
- It gives a within-run, read-level held-out set. It does **not** give an independent sequencing experiment:
  both pools share library prep, run conditions and chemistry. State this wherever results are reported.
- Replicate subsamples drawn from one 2 Gb pool overlap heavily at larger depths, so replicate spread is sampling variance only.
- Duplex or split reads (possible with Q20+ data) could place sister reads in different pools; check read-ID format
  (e.g. `id1;id2`) in the first reads before running and note any leakage.
- Fractions are by reads, so pool sizes in Gb vary slightly with the read-length distribution; the manifest records the real values.

## Running it (WSL)
```
cd ~/projects/metagenomics_1_Nanopore
bin/00_ena_metadata.sh PRJNA1181840            # writes metadata/PRJNA1181840_read_run.tsv (includes fastq_md5)
curl -s https://<fastq_ftp from the TSV> | zcat | head -8   # optional: look at read IDs first (curl may report a pipe error; harmless)
tmux new -s subset                              # long job; survives a closed terminal
bin/01_stream_subset.sh SRR31222153 20261012    # defaults: seed 20261012, 0.0156 + 0.0156
```
- Keep outputs in `~/projects/metagenomics_1_Nanopore/data/subsets/` (gitignored), NOT in OneDrive (large, would sync).
- A failed or interrupted run cannot resume mid-stream; it restarts from byte 0. The output is deterministic, so a rerun gives identical pools.
- Needs about 15 GB free (script refuses otherwise; override with MIN_FREE_GB).

## Time estimate (download-bound, ~120 GB)
| Connection | Approx. time |
|---|---|
| 50 Mbit/s | 5.3 h |
| 100 Mbit/s | 2.7 h |
| 500 Mbit/s | 32 min |
Decompression and the Python pass should keep up with 100-500 Mbit/s but are single-threaded; measure on the first
million reads (progress is printed every 1 M reads) and extrapolate.

## After it runs: representativeness checks (record in PROJECT_STATE.md)
- Pool mean read length vs the run mean (5,479 bp from ENA counts); read counts vs expected 0.0156 x 23,429,059 = ~365k.
- `seqkit stats` and NanoPlot on each pool: length and quality distributions should match between pools.
- Chemistry and basecaller still need to be recorded from the SRA experiment page or the Argo paper methods.

## Open decision (for you)
Pool size: 2 Gb each is a first guess matched to 12 GB RAM and quick iteration. If the Kraken2 database proves too constrained,
depth matters less than database choice, and smaller pools would do. Revisit after the first classifier run.
