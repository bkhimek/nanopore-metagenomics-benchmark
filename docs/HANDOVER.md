# HANDOVER

## 2026-10-12 session (Cowork, planning)
- Plan reviewed and revised to rev 3 (OneDrive project folder).
- PRJNA1181840 ENA run record verified; see PROJECT_STATE.md.
- Next: run bin/00_ena_metadata.sh in WSL, fetch chemistry/basecaller from SRA/paper, audit LomanLab
  accessions, check disk space, then choose the Phase 0 dataset (Decision 1) and fix metrics (Decision 2).
- Rule: update PROJECT_STATE.md and this file before ending any session, in any tool.

## 2026-10-12 session (data access)
- Measured data access; switched to SRA block sampling (see PROJECT_STATE.md, docs/SUBSET_PLAN.md).
- Next: bin/01_sra_block_sample.py --plan-only, then run it in tmux (~1.5 h); then seqkit/NanoPlot checks on both pools;
  confirm chemistry/basecaller; pick Phase 0 dataset (Decision 1) and fix metrics (Decision 2).
- GitHub: https://github.com/bkhimek/nanopore-metagenomics-benchmark

## 2026-10-12 end of day (supersedes the earlier "Next:" lines)
State: Phase 0 in progress. Sampling is complete; QC is next.
- Done: plan rev 3 (OneDrive copy); repo scaffold in WSL; ENA run record verified for PRJNA1181840 (SRR31222153);
  data access measured (~0.5 MB/s everywhere, ENA FASTQ is one gzip stream); SRA read-range route verified
  (first 200 reads identical to ENA); bin/01_sra_block_sample.py run to completion (74 blocks, 168.6 min).
- Result: dev and eval pools, each 185,000 reads (~1 Gb), mean length 5,491 / 5,473 bp (run mean 5,479);
  no length drift across the run (-30 bp). Files in data/subsets/, manifest in metadata/. Details: docs/PROJECT_STATE.md.
- Decisions so far: D1 Nextflow (plain scripts first), D2 taxonomy first, D3 one frozen database, D4 memory-capped
  classifier setup (12 GB RAM), D5 dev/eval separation, D6 SRA block sampling, D7 stream scripts removed.
- Not yet checked: read QUALITY distributions and pool similarity; chemistry/basecaller; Zymo D6331 composition.

Next steps, in order:
1. QC both pools (seqkit stats, NanoPlot) in pinned Docker images; compare length and quality between pools.
2. Confirm chemistry/kit and basecaller (SRA experiment SRX26602699 page, or the Argo paper methods); get the
   Zymo D6331 reference composition from the manufacturer spec sheet. Record both in the dataset register.
3. Settle Decision 1 (dataset) and Decision 2 (metrics, detection threshold rule) in PROJECT_STATE.md.
4. Optional: check whether LomanLab data has ENA accessions (second dataset for a chemistry-transfer test).
5. Phase 1: choose the baseline classifier and a database that fits 12 GB RAM; freeze and record it.
6. Only then wrap steps in Nextflow (D1).

Housekeeping:
- GitHub: https://github.com/bkhimek/nanopore-metagenomics-benchmark. Check `git status`; commit and push if not done.
  Still to add: LICENSE, README line "research use only, not clinically validated".
- OneDrive folder: delete setup_project.sh and the stray cd.docx; the plan copy and SUBSET_PLAN.md are snapshots (recopy if changed).
- Rules: update this file and PROJECT_STATE.md before ending any session; verify numbers from primary records, not summaries.

## 2026-10-13 session (plan review and document update)
State: Phase 0 in progress. Pools sampled; QC still NOT done. Supersedes the 2026-10-12 "Next steps" ordering.
- Reviewed external (ChatGPT) feedback point by point; adopted the valid points. Rejected none outright; items about
  the D6331 composition are unverified until checked against the Zymo certificate.
- Plan updated to rev 4 (Nanopore_Metagenomics_Project_Plan_Rev4.docx; Rev3 kept). Updated SUBSET_PLAN.md, PROJECT_STATE.md
  (new decisions D8-D11, D6331 composition section).
- Changes in method: AUPRC supplementary; recall over evaluable species; depth by Gb with full pool as single baseline;
  block-level resampling; database must include fungi/archaea.

Next steps, in order:
1. DONE (evening): D6331 composition checked against the Zymo manual and stored as reference/d6331_composition.csv and d6331_expected_reads_by_species.csv. Still open: lot CoA (lot number unknown; email Zymo or Argo authors).
2. QC both pools (seqkit stats, NanoPlot) in pinned Docker images; compare length and quality between pools.
3. Confirm chemistry/kit and basecaller from the SRA/ENA experiment record (reported, not verified: SQK-LSK114, FLO-PRO114M, Dorado v7.2.13).
4. Settle Decisions 1 and 2; choose the baseline classifier and a database that fits 12 GB and includes fungi/archaea.
5. Then wrap in Nextflow (D1). An educational workflow guide is being written alongside (Word document in OneDrive only: Nanopore_Metagenomics_Workflow_Guide.docx; not in the repo).

## 2026-10-15 session (chemistry confirmed, QC script)
State: Phase 0 essentially complete. Pools sampled and QC done (see PROJECT_STATE.md, Pool QC result): pools indistinguishable, no drift, median read Q ~18.
- Chemistry and basecaller verified from the ENA record of SRX26602699 (WSL curl): SQK-LSK114, FLO-PRO114M, Dorado v7.2.13
  (submitter free text; model tier not stated). Recorded in PROJECT_STATE.md, SUBSET_PLAN.md and the plan (rev 4, minor update).
- Added bin/02_pool_qc.py: per-read length and mean Q (error-probability mean), per-block and per-pool summaries, dev-vs-eval
  comparison at block level, drift vs position. Output in data/qc/. Supersedes the 2026-10-13 step 2 and 3 wording
  (seqkit/NanoPlot images not needed for this step).

Evening update: PlusPF-8 (2026-06-26) downloaded and verified, all 17 species present, Kraken2 2.17.2 pinned by digest, smoke test OK (see PROJECT_STATE.md, D12),
bin/03_classify.py written. AWS idea: deferred to Phase 2 (large-database comparison).

Next steps, in order:
0. DONE (evening 2026-10-15): bin/03_classify.py --pool dev --label baseline_c0 (94.14% classified, 17 s wall, within the 10g cap). See the 2026-10-16 entry below for what comes next.
1. DONE: bin/02_pool_qc.py run on the real pools; numbers recorded in PROJECT_STATE.md.
2. Optional: lot CoA (email Zymo / Argo authors); basecalling model tier (Argo paper methods).
3. Settle Decisions 1 and 2 (dataset confirmation; metrics and threshold rule, including the evaluable-species minimum).
4. Choose the baseline classifier and a database that fits 12 GB and includes fungi/archaea; freeze and record.
5. Then wrap in Nextflow (D1).

## 2026-10-16 session (scoring script)
State: Phase 1 started. Baseline run exists on the dev pool; not yet scored. NOTE: the OneDrive copies of PROJECT_STATE.md and HANDOVER.md written on 15 Oct did not overwrite
the older versions (the repo commit 4390d57 therefore lacks the baseline-run record); this entry and the updated PROJECT_STATE.md restore it. Cause: files sent to OneDrive right after being written can arrive with stale content; send, list the folder to compare sizes, and send again if they differ.
- Wrote bin/04_score.py (species-level scoring by taxid, threshold sweep on dev, FP audit, block bootstrap, eval-pool guard) and tested it on a synthetic dataset with hand-derived answers.
- Method decisions still needed from the user (Decision 2): evaluable-species minimum (proposed 20 expected reads) and the rule for choosing T on the dev pool (proposed: maximise F1 on dev, ties to larger T).

Next steps, in order:
1. Copy bin/04_score.py into the repo, run `bin/04_score.py --run data/classify/dev/baseline_c0` (sweep; dev pool only) and paste the output.
2. Decide the evaluable minimum and the T rule (Decision 2); record as D13; then run with --min-reads T on dev for the bootstrap and FP audit.
3. Only then freeze and run bin/03_classify.py --pool eval --allow-eval --label baseline_c0, and score it once with the frozen T.
4. Phase 1 extras: depth series (0.1/0.25/0.5 Gb x 5 replicates), then Phase 2 (AWS large database, long-read-aware profiler, AMR).
5. Optional: lot CoA, basecalling model tier, Argo paper link.

## Paths (reference)
- WSL repo: ~/projects/metagenomics_1_Nanopore/ (data in its gitignored data/ subfolder)
- OneDrive project folder (Windows): C:\Users\krist\OneDrive\Documents\Projects\metagenomics\metagenomics_1_Nanopore
- OneDrive project folder (from WSL): /mnt/c/Users/krist/OneDrive/Documents/Projects/metagenomics/metagenomics_1_Nanopore
- GitHub: https://github.com/bkhimek/nanopore-metagenomics-benchmark
