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

## 2026-10-16 afternoon (first real scoring, Fusobacterium finding)
State: bin/04_score.py committed (895e437) and run on the real dev baseline (sweep). Results and numbers: PROJECT_STATE.md, Scoring. Repo docs restored (baseline-run record present).
- Finding: Fusobacterium nucleatum showed 92 of 11,100 expected reads; the reads are under F. animalis (8,597). The D6331 reference genome, classified in 5 kb pieces with the frozen baseline, is 97.1% F. animalis,
  so this is a truth-table naming mismatch, not low sensitivity. Fixed by an alias file (D13), reported in the score output. Found after seeing results; recorded as such.
- Other low recoveries are real classifier behaviour (Bifidobacterium adolescentis stops at genus level; Roseburia loses reads above genus/unclassified; E. coli shares reads with relatives).
- 04_score.py updated: --min-frac, --sweep-frac, --aliases, genus-level recall, wider default grid. Regression-tested against the previous version.
- Alias applied (commit e161f77): Fusobacterium nucleatum 8,689 reads (0.78x). Sweep rerun; F1 peaks at T=150 (0.966) one read below Methanobrevibacter smithii (151 reads): a knife edge.
- User decided: evaluable minimum stays 20 (D14a); threshold = fraction of reads chosen by bootstrap-mean F1 on dev, ties to the smaller (D14b, accepted).
- 04_score.py prints the bootstrap MEAN of each metric (commit 8f6c741). Candidate comparison done: best is --min-frac 0.0008 (T=148 on 185,000 reads), FROZEN. Remaining FP at that threshold is Homo sapiens (923 reads).
  Predictions for the eval run are written in PROJECT_STATE.md before it is run.

- EVAL RUN DONE (freeze-d14 tag precedes it; commit 8e1e2d3): 14/14 evaluable detected, 1 false positive (human, 818 reads), precision 0.933, F1 0.966, Bray-Curtis 0.168; all written predictions met. Details in PROJECT_STATE.md.
  The eval pool has now been used once: do not re-score it with other thresholds or tune on it.

Next steps, in order (supersedes the 2026-10-16 morning list):
1. DONE: candidate comparison; threshold frozen at 0.0008 (D14b).
2. DONE: eval pool classified and scored once; compared with predictions.
2a. DONE: depth series on the eval pool (tag depth-registered precedes it); all registered predictions met; results table in PROJECT_STATE.md; small outputs to be copied to reports/baseline_kraken2_pluspf8/ (depth_summary.tsv, depth_species.tsv, depth_report.txt).
2f. 04_score.py edge-case fix (no species called) and the results/ -> reports/ move (results/ is gitignored) are described in PROJECT_STATE.md; the arms-registered tag was deleted and recreated after the real commit.
2e. Chosen by the user: arm 1 = Kraken2 confidence sweep (D16, registered in PROJECT_STATE.md). Steps: copy bin/06_confidence_arms.py, commit + tag arms-registered, run `--stage dev` (about 5-10 min), paste the dev tuning table, then commit reports/arms and tag arms-frozen, then `--stage eval --allow-eval`.
2d. (original list of candidates, others still open) Next: choose Phase 2 arms. Candidates: (i) Kraken2 --confidence sweep (0.05/0.1/0.2) as a cheap check on false positives; (ii) large database on AWS to test what the 8 GB cap costs (abundance bias); (iii) a long-read-aware profiler (Sylph, minimap2-based); (iv) trace the human reads (about 0.4% in every run).
    RULE for any arm that has a tunable setting: tune it on the dev pool with the same rule (bootstrap-mean F1, ties to the smaller), register the full list of arms and predictions in advance, then score ALL arms on eval and report all of them (no picking the best arm on eval).
2c. DONE 2026-10-16: plan refreshed to revision 5 (Nanopore_Metagenomics_Project_Plan_Rev5.docx) and the workflow guide to v2 (Nanopore_Metagenomics_Workflow_Guide_v2.docx), both in OneDrive. Delete the old guide file once happy.
2b. Optional small checks: where the human reads come from (BLAST/minimap2 of a sample of them), genome accession of the Fusobacterium strain, whether L. fermentum 1.66x is relatives' reads.
3. (merged into step 2 above)
4. Phase 1 extras: depth series (0.1/0.25/0.5 Gb x 5 replicates; the fraction threshold carries over), then Phase 2 (AWS large database to test the capped index, long-read-aware profiler, AMR).
5. Optional: lot CoA, basecalling model tier, Argo paper link, genome accession for the Fusobacterium strain.

## Paths (reference)
- WSL repo: ~/projects/metagenomics_1_Nanopore/ (data in its gitignored data/ subfolder)
- OneDrive project folder (Windows): C:\Users\krist\OneDrive\Documents\Projects\metagenomics\metagenomics_1_Nanopore
- OneDrive project folder (from WSL): /mnt/c/Users/krist/OneDrive/Documents/Projects/metagenomics/metagenomics_1_Nanopore
- GitHub: https://github.com/bkhimek/nanopore-metagenomics-benchmark

## 2026-10-16 afternoon (confidence arms done)
State: Phase 2 arm 1 (Kraken2 confidence, D16) complete. Tags: arms-registered (f3c9ba7), arms-frozen (15e92f3); eval result committed a427e54 (reports/arms/arms_summary.txt, dev_tuning.tsv, dev_selection.json).
- Result: baseline c0 at 0.0008 is best (eval F1 0.966; c0.05 0.929 at T=10 but only 19% of reads in mock species at species rank and Bray-Curtis 0.34; c0.1 F1 0.60, recall 6/14; c0.2 and c0.4 call nothing). Full table and comparison with the six registered expectations (1 and 5 met, 2 partly met, 3 partly contradicted, 4 and 6 contradicted) in PROJECT_STATE.md, Confidence arms RESULT.
- Open from this arm: c0.05's optimum is at the grid edge (5e-5); the second half of expectation 2 (non-human false positives at fraction 0.0002) not checked; genus-level use of higher confidence is untested.
- Next: choose Phase 2 arm 2. Candidates: (ii) large database on AWS (tests whether the 8 GB cap causes the abundance bias: Roseburia 0.41x, B. adolescentis 0.22x, E. coli 0.57x, L. fermentum 1.66x); (iii) long-read-aware profiler (Sylph or minimap2); (iv) trace the human reads. Same rule: register arms and predictions first, tune on dev, score all arms on eval.
- Still to do: refresh the guide (depth series, confidence arms, reports/ path instead of results/) and the plan (Rev6); housekeeping (LICENSE, README line, stray cd.docx).

## 2026-10-16 late afternoon (database arms registered, AWS in progress)
State: Phase 2 arm 2 (database size, D17) registered in PROJECT_STATE.md with decision rule and six expectations; bin/07_database_arms.py written and dry-run tested (not yet run on data). Commit it and tag `dbarms-registered` BEFORE uploading any reads to AWS.
- AWS set-up done so far: region eu-north-1 (Stockholm), account budget $15 with alerts at 50% and 85%, key pair nanopore-bench (private key in ~/.ssh/nanopore-bench.pem, mode 600, not in the repo), vCPU quota 32, no leftover EC2 resources.
- Next: launch r6i.4xlarge (Ubuntu 24.04, 200 GB gp3, SSH from my IP only), run the steps in 'Run order' above, TERMINATE the instance and check Volumes = 0 afterwards.

## 2026-10-16 evening (database arms done)
State: Phase 2 arm 2 (D17) complete. Tags dbarms-registered (820490a) and dbarms-frozen (afe5577); the eval summary is reports/dbarms/dbarms_summary.txt (commit it); small outputs are in data/aws_run/dbarms_outputs.tgz (gitignored). Details and the comparison with the six registered expectations are in PROJECT_STATE.md, Database arms RESULT.
- Result: the 8 GB cap is NOT the main cause of the abundance bias (registered rule: Roseburia and B. adolescentis both <= 0.50 with PlusPF full). Larger index: classified 94.2 -> 98.0%, reads in mock species 83.9 -> 92.3%, Bray-Curtis 0.168 -> 0.137, E. coli obs/exp 0.57 -> 0.86, B. adolescentis 0.22 -> 0.36, Roseburia 0.41 -> 0.43; F1 unchanged at the ceiling 0.966; human reads at species rank rise 818 -> 1,997 on eval. Expectations: 1, 3, 4, 5, 6 met; 2 partly (Roseburia guess contradicted).
- AWS: instance i-0575eeb3dbb5c13bc (eu-north-1) was TERMINATED at about 18:40 BST on 2026-10-16 (console showed 'Terminated') after the output pack was copied home; the private key file was shredded from WSL and the host entries removed from known_hosts. Confirmed by the user at 18:50 BST: EC2 Volumes = 0 and the instance shows Terminated. The key pair nanopore-bench still exists in the console (no cost; can be deleted). Total cost roughly USD 2 (r5.4xlarge, about 1.75 h; check the AWS bill).
- Next: choose Phase 2 arm 3: (iii) long-read-aware profiler (Sylph or minimap2 against a species-level reference) is now the natural next test of whether the remaining bias is a classifier property; (iv) trace the human reads. Same rule: register arms and predictions first, tune on dev, score all arms on eval.
- DONE 2026-10-16 evening: guide v3 (Nanopore_Metagenomics_Workflow_Guide_v3.docx) and plan Rev6 (Nanopore_Metagenomics_Project_Plan_Rev6.docx) written and delivered to OneDrive; delete Guide_v2 and Plan_Rev5 once checked. Still to do: housekeeping (LICENSE, README line 'research use only, not clinically validated', stray cd.docx in OneDrive).

## 2026-10-17 session (read trace; abundance 'bias' re-interpreted)
State: Phase 2 arm 2 results unchanged and registered verdicts untouched. Two EXPLORATORY dev-only scripts added: bin/08_trace_reads.py (committed c3d3fa2; minimap2 labels vs Kraken2 destinations; outputs in reports/trace/) and bin/09_composition_vs_mapping.py (to commit with reports/trace/bc_vs_mapping.tsv). Full numbers and caveats: PROJECT_STATE.md, 'Read trace and composition check'.
- Main finding: the Roseburia (0.42), L. fermentum (1.71) and F. prausnitzii (1.26) deviations are in the mapped reads (tube composition), and Kraken2 reproduces them (obs/mapped 0.98-1.03). Classifier loss is confined to E. coli (stops at family, 0.62 -> 0.94 with the full index) and B. adolescentis (stops at genus, 0.44 -> 0.72). Bray-Curtis against the mapped composition: 0.065 / 0.048 / 0.022 for PlusPF-8 / 16 / full (against theory 0.165 / 0.150 / 0.131).
- Wording to correct: plan Rev6 and guide v3 describe the abundance deviations as a database/classifier property; that overstates it (see PROJECT_STATE addendum). Plan Rev7 and guide v4 should say: presence/absence good; abundance deviations = tube composition (mapping shows it) plus Kraken2 caution for close relatives; the cap matters for E. coli and B. adolescentis only. Not done yet.
- Data kept locally (gitignored): data/trace/dev/ (PAF, combined reference), data/aws_run/unpacked/ (AWS pack). Temporary: /tmp/pack can be deleted.

Next steps, in order:
1. Commit bin/09, reports/trace/bc_vs_mapping.tsv and the two updated docs (cp from OneDrive transfer/d3).
2. Email Zymo (tech@zymoresearch.com) for the lot CoA / or ask the Argo authors; the trace makes this the most valuable open item.
3. Correct the wording in plan (Rev7) and guide (v4).
4. Phase 2 arm 3: a long-read-aware profiler (Sylph or minimap2-based). Register arms and predictions first, tune on dev, tag, score eval once per arm. Score every arm against BOTH references (theoretical composition = registered; mapped composition = exploratory) and decide in advance which one is primary.
5. Housekeeping: LICENSE, README line 'research use only, not clinically validated', stray cd.docx, delete Guide_v2 and Plan_Rev5 once checked, optionally the unused AWS key pair; confirm the AWS bill (about USD 2 expected).
