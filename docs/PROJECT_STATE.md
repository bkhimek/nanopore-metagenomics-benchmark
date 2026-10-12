# PROJECT_STATE: decision log and findings

Last updated: 2026-10-12

## Status
Phase 0, dataset feasibility audit, in progress.

## Environment
- WSL2 Ubuntu, 12 GB RAM cap, 6 processors (C:\Users\krist\.wslconfig)
- Repo: ~/projects/metagenomics_1_Nanopore/ ; data (gitignored, inside the project folder): ~/projects/metagenomics_1_Nanopore/data/
- OneDrive copy of selected outputs: Projects/metagenomics/metagenomics_1_Nanopore
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

## Findings (verified from primary records)
- PRJNA1181840 (ENA, checked 2026-10-12): one run SRR31222153; PromethION; WGS/METAGENOMIC/RANDOM;
  23,429,059 reads; 128,369,700,312 bases (mean ~5.5 kb); FASTQ 119,695,816,967 bytes;
  study "Deep Nanopore Q20+ sequencing of Zymo D6331"; first public 2025-02-23.
- LomanLab mockcommunity README: even (ZRC190633) and log (ZRC190842) Zymo standards, 10 species each,
  GridION ~14-16 Gb and PromethION ~146-148 Gb; hosted on CLIMB, no ENA accessions listed there.

## Unverified / open
- PRJNA1181840 chemistry (Kit 14?), basecaller/model, Zymo D6331 composition, link to Argo paper (Nat Commun 2025)
- LomanLab ENA accessions (if any)
- PRJNA751627, PRJEB41184, PRJNA1335824, PRJNA554461 (clinical leads, deferred)
- Free disk space for subsets and databases
