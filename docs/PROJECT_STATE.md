# PROJECT_STATE: decision log and findings

Last updated: 2026-10-15 (evening)

## Status
Phase 0, dataset feasibility audit, in progress. Run record, chemistry and basecaller verified; both pools sampled;
theoretical composition checked. Remaining: QC of both pools (script ready: bin/02_pool_qc.py), lot-specific CoA (lot unknown),
basecalling model tier, database choice (must fit 12 GB and include fungi/archaea). Plan is at revision 4 (minor update 15 October).

## Environment
- WSL2 Ubuntu, 12 GB RAM cap, 6 processors (C:\Users\krist\.wslconfig)
- Repo: ~/projects/metagenomics_1_Nanopore/ ; data (gitignored, inside the project folder): ~/projects/metagenomics_1_Nanopore/data/
- OneDrive project folder (copy of the plan, subset plan, handover, manifest and selected outputs; no scripts):
  Windows: C:\Users\krist\OneDrive\Documents\Projects\metagenomics\metagenomics_1_Nanopore
  From WSL: /mnt/c/Users/krist/OneDrive/Documents/Projects/metagenomics/metagenomics_1_Nanopore
- Nextflow from the `nextflow` conda env; tools in pinned Docker images; no new conda envs

## Decisions
- D1 Nextflow DSL2 will be used. Reasoning: replicated downsampling (3 depths x 5 replicates = 15+
  tasks, more with extra databases/classifiers), per-process containers, -resume on long steps, and
  the deliverable is a portfolio pipeline meant to be rerun by others. Counter-argument: single
  machine, single developer, so first modules are written and tested as plain scripts, then wrapped.
- D2 Taxonomy benchmark first; AMR only if an independent truth set exists (otherwise descriptive).
- D3 One frozen database for the baseline; a second database is a Phase 2 experiment.
- D4 Memory: full standard Kraken2 DB does not fit 12 GB. Use a capped-size index, memory mapping,
  or minimap2; record the sensitivity limitation.
- D5 Evaluation: development and held-out data separated; detection-threshold rule fixed before evaluation.

- D8 (2026-10-13) SRR31222153 is the PROVISIONAL primary dataset for the taxonomy benchmark, pending verification of
  the D6331 composition and pool QC. LomanLab data remain a fallback (older chemistry, so confounded).
- D9 (2026-10-13) Metrics: primary = precision/recall/F1 at a pre-specified threshold plus a false-positive audit;
  AUPRC supplementary only. Recall is computed over "evaluable" species (expected read count at that depth above a
  pre-specified minimum); the rest are listed separately. Depth is defined by Gb (0.1, 0.25, 0.5 Gb x 5 replicates);
  the full ~1 Gb pool is a single baseline without replicates.
- D10 (2026-10-13) Uncertainty: block-level resampling (reads come in 5,000-read blocks). Replicates measure
  within-pool variation only (about half overlap at 0.5 Gb), not run-to-run variation.
- D12 (2026-10-15) Baseline classifier = Kraken2 2.17.2 (image quay.io/biocontainers/kraken2@sha256:db5d772d9df323cf0ddc644fb9793c0fca25b61f86469be42f1a9f5b903b9fdf)
  with the PlusPF-8 index dated 2026-06-26 (7.5 GB, includes fungi and archaea; FROZEN as the baseline database, D3). Run with a docker memory cap
  (10g default). AWS deferred: a large-database run on a big-RAM instance is the Phase 2 experiment to measure what the 8 GB cap costs.
- D11 (2026-10-13) Database must contain fungi and archaea (standard Kraken2 indices lack fungi); otherwise those
  mock members are declared out of evaluable scope. Check coverage before running. Refines D4.

## Findings (verified from primary records)
- PRJNA1181840 (ENA, checked 2026-10-12): one run SRR31222153; PromethION; WGS/METAGENOMIC/RANDOM;
  23,429,059 reads; 128,369,700,312 bases (mean ~5.5 kb); FASTQ 119,695,816,967 bytes;
  study "Deep Nanopore Q20+ sequencing of Zymo D6331"; first public 2025-02-23.
- LomanLab mockcommunity README: even (ZRC190633) and log (ZRC190842) Zymo standards, 10 species each,
  GridION ~14-16 Gb and PromethION ~146-148 Gb; hosted on CLIMB, no ENA accessions listed there.

## Chemistry and basecaller (verified 2026-10-15 from the ENA record of experiment SRX26602699)
- DESIGN_DESCRIPTION (submitter free text, not a structured field): "SQK-LSK114 + FLO-PRO114M + Dorado v7.2.13".
  Instrument PromethION; WGS / METAGENOMIC / RANDOM; single layout. So Kit 14 ligation chemistry, PromethION flow cell, Dorado 7.2.13.
- Not stated: basecalling model tier (fast / hac / sup). Link to the Argo paper: still unverified by me.
- Results apply to this chemistry and basecaller only (plan Section 3).

## D6331 expected composition (checked 2026-10-13 against the manufacturer's manual, Tables 1 and 2)
- Source: ZymoBIOMICS Gut Microbiome Standard D6331 instruction manual (PDF SHA-256 a6fdb1b4...cc76eb). THEORETICAL composition;
  the lot-specific certificate of analysis (CoA) is not obtained: lot number of the sequenced tube unknown (BioSample has no lot).
- 21 strains: 18 bacteria (5 E. coli), 1 archaeon (M. smithii), 2 fungi (C. albicans, S. cerevisiae); 17 species.
- Five columns: genomic DNA, 16S only, 16S+18S, genome copy, cell number. Genomic DNA runs 14% to 0.0001% (sums to 100.011).
- Verified: genome-copy % = genomic-DNA % / genome size, renormalised. Reads follow DNA mass (use the DNA column for read-fraction
  outputs); coverage per genome follows the genome-copy column (use for coverage-based abundance, e.g. Sylph).
- Expected reads at 0.1/0.25/0.5/1.0 Gb (DNA basis, 5,479 bp mean): Salmonella 1.8/4.6/9.1/18; E. faecalis 0.2/0.5/0.9/1.8;
  C. perfringens ~0/0/0.1/0.2; M. smithii 18/46/91/183. Salmonella (0.01%) is at the manual's stated impurity limit (<0.01%
  foreign DNA) and is a close relative of the E. coli (14% of DNA combined): hard and interesting, but contamination cannot be
  separated from signal at that level.
- Reference genomes (D6331.refseq.zip, SHA-256 4795b7df...97d3d): bacterial/archaeal lengths match Table 2; the two yeast files are
  draft assemblies (C. albicans 13,062 contigs, 12.58 Mb vs 14.68 in Table 2; S. cerevisiae 51 contigs, 12.84 vs 13.30).
- Files: reference/d6331_composition.csv, reference/d6331_expected_reads_by_species.csv, reference/README.md (in the repo).
- Manual states relative abundance deviation <15% on average and impurity <0.01%.

## Baseline database and classifier (2026-10-15)
- Index: Kraken2 PlusPF-8, 2026-06-26 (genome-idx.s3.amazonaws.com/kraken/k2_pluspf_08_GB_20260626.tar.gz). Archive md5 79a153b99f045bc2ae95e6d57c17a02d (verified
  against the published md5 list). Archive 5.5 GB; index hash.k2d 7.5 GB; in data/db/pluspf_08_20260626/ (gitignored). Contents: RefSeq archaea, bacteria, viral, plasmid, human,
  UniVec_Core plus protozoa and fungi (per the index listing page, read 2026-10-15).
- Coverage check (inspect.txt): all 17 D6331 species present at rank S. Zymo's "Lactobacillus fermentum" is "Limosilactobacillus fermentum" (taxid 1613) in the current taxonomy:
  score by taxid. Minimizer counts per species are uneven (E. coli 1.79 M; M. smithii 38 k; Fusobacterium 101 k; Candida albicans 252 k): species with few retained minimizers
  may be detected less sensitively in a capped index. To be measured, not assumed.
- Smoke test (one 5,000-read block, start 10,336,400 = stratum 32, dev pool): 94.38% classified, 5.62% unclassified (confidence 0); 1.45 s classification, 11 s wall incl. database load.
  Single block, not scored.
- Script: bin/03_classify.py (pinned image digest, docker memory cap, run.json record; refuses the eval pool without --allow-eval).
- Reference files now carry species_taxid (reference/d6331_*.csv).

## Unverified / open
- Argo paper link and basecalling model tier (see below)
- Lot-specific CoA for the sequenced tube: lot number unknown. Ask Zymo (tech@zymoresearch.com) or the Argo paper authors. Do not block the project on it; label results as against theoretical composition.
- Basecalling model tier (fast/hac/sup) not stated in the ENA record; look in the Argo paper methods or ask the authors.
- Argo paper link: reported as confirmed in the paper's data availability statement; not verified by me.
- Decisions 1 (dataset) and 2 (metrics, threshold rule) not yet formally settled
- LomanLab ENA accessions (if any)
- PRJNA751627, PRJEB41184, PRJNA1335824, PRJNA554461 (clinical leads, deferred)
- Disk: 664 GB free at 2026-10-12 (checked); recheck before database downloads

## Data access findings (2026-10-12)
- ENA FASTQ for SRR31222153 is one gzip stream (no complete gzip member in the first 100 MB): no random access.
  Read IDs are rewritten as SRR31222153.n. ENA MD5 b60b15eee28618650e8fc5c0b8066b06 (119,695,816,967 bytes).
- Throughput from WSL ~0.5 MB/s everywhere (ENA 0.47, NCBI/AWS 0.56, Cloudflare 0.62): the connection is the limit.
  Full ENA stream would take ~3 days: rejected.
- SRA route works: image ncbi/sra-tools@sha256:4f048ba6754246e1ed7799a2c8ed469a38427f392e403ca08b93c4255104c9ce;
  fastq-dump -N/-X by read number; 20,000 reads in 4 min 14 s (~83 reads/s; ~14 s fixed overhead per call);
  first 200 reads identical to ENA (sequence+quality md5 66b4146bba99f796f337d62a248f5b21).
  SRA run object is 109,839,515,692 bytes (not smaller than ENA, but range-readable).
- D6: stratified block sampling via SRA, 74 strata, 5,000-read blocks, 37 per pool (~185,000 reads, ~1 Gb), seed 20261012;
  even strata -> dev pool, odd -> eval pool (whole blocks stay together). Within-run held-out only. See docs/SUBSET_PLAN.md.
- D7: the ENA full-stream scripts (01_stream_*) were removed as superseded.
- Disk: 664 GB free on WSL root at 2026-10-12.

## Subset sampling result (2026-10-12)
- bin/01_sra_block_sample.py, seed 20261012: 74 strata, 5,000-read blocks, 37 per pool. Finished in 168.6 min;
  no failure or retry lines in the run output or fetch log.
- dev pool: 185,000 reads, 1,015,906,091 bases, mean length 5,491.4 bp.
  eval pool: 185,000 reads, 1,012,464,842 bases, mean length 5,472.8 bp. Whole-run mean (ENA counts): 5,479 bp.
- Mean length by run quarter: 5,489 / 5,493 / 5,481 / 5,467 bp (drift about -30 bp, -0.5%); block means 5,295-5,706 bp (sd 81).
- Pool files (data/subsets/): dev 914 MB, eval 911 MB. SHA-256 prefixes: dev cd5e6500996fb456, eval ea568fa2704b977e
  (full values in the manifest, copied to metadata/).
- Not yet checked: read quality distributions and pool similarity (seqkit/NanoPlot), chemistry and basecaller.
- Status: sampling complete; QC done (next section).

## Pool QC result (2026-10-15, bin/02_pool_qc.py; outputs in data/qc/, scipy 1.18.1 pip-installed into base)
- Read quality = per-read error-probability mean (Phred+33; quality characters 34..83). Reads: 185,000 per pool.
- dev vs eval: mean length 5,491 vs 5,473 bp; N50 9,252 vs 9,227; 8.9% vs 9.0% of reads <500 bp; mean read Q 16.95 vs 16.97; median Q 17.91 vs 17.92;
  Q>=10 89.9% vs 90.1%; Q>=15 74.0% vs 74.0%; Q>=20 27.6% vs 27.9%.
- Length percentiles 5/25/50/75/95 (dev): 382 / 1,307 / 4,061 / 8,036 / 15,554 bp. Read Q percentiles 5/25/50/75/95: 6.3 / 14.8 / 17.9 / 20.3 / 22.8.
- Pool similarity: KS D (effect size) length 0.0025, mean Q 0.0029. Block-level Mann-Whitney (37 vs 37 blocks): mean length p=0.55, mean Q p=0.22,
  % reads Q>=20 p=0.11. No detectable difference between pools.
- Drift along the run (Spearman, block mean vs read-number position): length rho -0.11 (p 0.36); mean Q rho -0.005 (p 0.96). No detectable drift.
- Interpretation: "Q20+" in the study title describes the chemistry; only ~28% of reads reach Q20 and the typical read is about Q18.
  Report as Kit 14 data, median read Q about 18. No filtering applied; any length/quality filter must be pre-specified and applied identically to both pools.
- Phase 0 audit items now closed except: lot CoA (unknown lot), basecalling model tier, Argo paper link.
