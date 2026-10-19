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
- DONE 2026-10-17: plan Rev7 (Nanopore_Metagenomics_Project_Plan_Rev7.docx) and guide v4 (Nanopore_Metagenomics_Workflow_Guide_v4.docx) written with the corrected wording (abundance deviations = tube composition plus Kraken2 caution for close relatives; the cap matters for E. coli and B. adolescentis only; registered verdicts kept). Delete Rev6 and guide v3 (and the older Rev5 and guide v2) from OneDrive once checked.
- Data kept locally (gitignored): data/trace/dev/ (PAF, combined reference), data/aws_run/unpacked/ (AWS pack). Temporary: /tmp/pack can be deleted.

Next steps, in order:
1. DONE (commit 45c3896): bin/09, reports/trace/bc_vs_mapping.tsv and the updated docs committed and pushed.
2. Email Zymo (tech@zymoresearch.com) for the lot CoA / or ask the Argo authors; the trace makes this the most valuable open item.
3. Phase 2 arm 3: a long-read-aware profiler (Sylph or minimap2-based). Register arms and predictions first, tune on dev, tag, score eval once per arm. Score every arm against BOTH references (theoretical composition = registered; mapped composition = exploratory) and decide in advance which one is primary.
4. Housekeeping: LICENSE, README line 'research use only, not clinically validated', stray cd.docx, delete Guide_v2 and Plan_Rev5 once checked, optionally the unused AWS key pair; confirm the AWS bill (about USD 2 expected).

## 2026-10-17 afternoon (Sylph arm registered; downloads running)
State: Phase 2 arm 3 (Sylph, D18) REGISTERED in PROJECT_STATE.md (section 'Sylph arm'); nothing run on reads yet. Commit this file set and tag `sylph-registered` BEFORE any read is sketched or profiled. Zymo email sent (reply pending).
- Running (user's WSL, tmux session sylpdl): download of gtdb-r232-c200-dbv2.syl2db (20,044,737,329 bytes) and fungi-refseq-2024-07-25-c200-v0.3.syldb (710,490,375 bytes) into data/db/sylph/; done when dl_done.txt exists. Check sizes and the .headers files (server hashes).
- Sylph 1.0.0 image pulled (quay.io/biocontainers/sylph:1.0.0--hb42e459_0).
Next steps, in order:
1. Commit PROJECT_STATE.md and HANDOVER.md; tag `sylph-registered`; push (tags: git push --tags).
2. Verify the downloads (sizes, hash), then the crosswalk step (sketch the 21 D6331 genomes, profile against the two databases, map to GTDB species with the GTDB r232 taxonomy file); commit reference/d6331_gtdb_crosswalk.csv before any reads are profiled. I write bin/10 (crosswalk) and bin/11 (sketch, profile, bootstrap, score).
3. Small test of joint profiling (-d with .syl2db plus .syldb) and peak memory (/usr/bin/time -v) on one block; then dev: sketch, profile, bootstrap, tune 3b, commit, tag `sylph-frozen`; then eval once per arm.
4. Housekeeping as before; ask the study authors for lot number and basecalling model tier (short email).

## 2026-10-17 evening (Sylph: crosswalk run, first dev profile, crosswalk RULE problem)
State: databases verified; crosswalk v1 produced with the registered rule; first dev profile done (sketch 5 s, profile 4 s, 2.4 GiB: no cloud machine needed); bin/11 committed (0308610). Scoring the dev profile against crosswalk v1 gives F1 0.50 for a NAMING reason (reads hit GTDB species other than the tied one the crosswalk picked); see PROJECT_STATE 'Sylph arm: progress'. Nothing from the eval pool has been used. Running: bootstrap in tmux session sylboot (100 resamples, resumable).
Done tonight: crosswalk v1 committed unchanged (58a48a0, tag sylph-crosswalk, pushed); dev bootstrap finished and committed (3691d89: B = 100, seed 20261017, slowest sketch 6.2 s, peak 2.84 GiB). reports/sylph/dev_3a/ (scoring against crosswalk v1) is deliberately NOT committed: not interpretable, do not cite.
Next steps, in order:
1. Start of next session: `git log --oneline -3` should show 3691d89 on top. Profiles for tuning are in data/sylph/dev/boot/b001..b100.tsv (gitignored, regenerable in about 25 minutes with `bin/11_sylph_arm.py bootstrap`).
2. Fix the crosswalk rule: test `sylph query` (genome sketches `sylph sketch -g ...` of the 21 references against the databases) for ANI without a coverage model; if unusable, an evidence-based alias table (D13 style). Write the addendum (disclose: made after seeing dev read profiles; eval untouched), also deciding the obs/exp basis and the false-positive rule; commit as crosswalk v2; rescore dev (`score --arm 3a`), then `tune`, write dev_selection.json, commit, tag `sylph-frozen`, then eval once per arm (`--allow-eval`). Report v1 and v2 results side by side.
3. Update the plan (Rev8) and guide (v5): Sylph registration, the crosswalk finding, and replace 'the tube' with 'the sequenced sample' (headings such as 'The tube is not the label').
4. Email the study authors (lot number; basecalling model tier); record any Zymo reply.
5. Housekeeping: LICENSE; README 'research use only, not clinically validated'; stray cd.docx; delete Plan_Rev5/Rev6, Guide_v2/v3 and OneDrive transfer/d1..d12 once the files are in the repo; optionally the unused AWS key pair; confirm the AWS bill (about USD 2).

## 2026-10-18 morning (Sylph: crosswalk v2, dev rescored, Addendum A1 written)
State: crosswalk v2 built with fastANI (bin/12, commit 6390ab5); dev rescored: 13/14 evaluable detected, 0 false positives, F1 0.963, but abundance is worse than Kraken2 (Bray-Curtis 0.248; E. coli 2.34x over). Arm 3b cutoff = 0 (identical to 3a). Addendum A1 (PROJECT_STATE) records the crosswalk change (made after seeing dev reads), the like-for-like verdict basis and the eval gap rule. Eval pool untouched; nothing frozen yet.
Next steps, in order:
1. Dev arm 3c (exploratory): `python3 bin/11_sylph_arm.py profile --pool dev --arm 3c` then `score --arm 3c --pool dev --crosswalk reference/d6331_gtdb_crosswalk_v2.csv --kraken ...`.
2. Commit addendum, crosswalk v2 and evidence (reference/d6331_gtdb_crosswalk_v2.csv, reports/sylph/crosswalk_v2_*, docs), tag `sylph-crosswalk-v2`. Then commit reports/sylph/dev_selection.json, dev_tuning.txt, dev_3a, dev_3a_crosswalk_v1, dev_3c; tag `sylph-frozen`; push with tags.
3. Claude writes the gap-rule extension of bin/12 (A1.2) and commits it BEFORE the eval profile exists.
4. Eval ONCE per arm: `sketch --pool eval --allow-eval`, `profile --pool eval --arm 3a|3c --allow-eval`, gap rule, `score --pool eval --arm 3a|3b|3c --allow-eval --crosswalk reference/d6331_gtdb_crosswalk_v2.csv --kraken PlusPF8=... PlusPFfull=...` (eval species tables under data/aws_run/unpacked/data/score/eval/). Apply the A1.3 verdict. Report v1 beside v2.
5. Optional dev investigation (not part of the verdict): why E. coli is over- and L. fermentum under-estimated (genome sizes, shared k-mers).
6. Plan Rev8 and guide v5 (Sylph, the crosswalk finding, 'the sequenced sample' wording); study-author email; housekeeping (see earlier entries; OneDrive transfer/ folders can be deleted).

## 2026-10-18 (Sylph arm DONE: eval scored once per arm)
State: Phase 2 arm 3 complete. Gap rule found nothing new; eval arms 3a = 3b = 3c: F1 0.963, 13/14 evaluable, 0 false positives, Bray-Curtis 0.249. B. adolescentis share-basis obs/exp 0.4947 vs Kraken2 PlusPF-full 0.3926 -> A1.3 verdict 'resolves' by 0.002 (marginal; literal 0.46 also exceeded). Abundance worse than Kraken2 (E. coli 2.38x over). v1 crosswalk on eval: F1 0.50, 7 false positives. Commits: 9734019 (gap rule, tag sylph-gaprule), 9051fb0 (--label), e5f7e8e (eval reports, tag sylph-eval-scored). Tag history: sylph-frozen moved from 8832e12 to 4337e4c before any eval data existed. Details in PROJECT_STATE 'Sylph arm: EVAL RESULT'.
Next steps, in order (supersedes the 2026-10-18 morning list):
1. Commit these two docs (PROJECT_STATE.md, HANDOVER.md) and push.
2. Optional dev investigation (not part of the verdict): why E. coli is over- and L. fermentum under-estimated (genome sizes, shared k-mers).
3. Plan Rev8 and guide v5 (Sylph registration and results, the crosswalk finding, 'the sequenced sample' wording instead of 'the tube').
4. Study-author email (lot number; basecalling model tier); record any Zymo reply.
5. Housekeeping: LICENSE; README 'research use only, not clinically validated'; stray cd.docx; delete Plan_Rev5/Rev6, Guide_v2/v3 and OneDrive transfer/d1..d19 once files are in the repo; optionally the unused AWS key pair; confirm the AWS bill (about USD 2). Optional later: fresh third pool (needs its own registered addendum), AMR design, Nextflow wrapping.

## 2026-10-18 later (documents written; Argo paper checked)
State: Plan Rev8 and Guide v5 written (OneDrive project folder; Rev7/v4 can be deleted; the guide explains Bray-Curtis and all new terms for a novice). The Argo paper and its supplement were read (PROJECT_STATE 'External literature check'): they confirm PromethION, Kit 14, Dorado v7.2.13 and PRJNA1181840 but give no lot number, extraction method or model tier. Contact: corresponding author Tong Zhang (zhangt@hku.hk); sequencing was assisted by Xianghui Shi. Email drafted; sending to be confirmed by the user. ENA/BioSample lot check not yet run (user to run the command and paste the output).
Next steps, in order (replaces the previous list):
1. Confirm the email was sent; record any reply (also the Zymo reply, if any). If a lot number arrives, request its certificate of analysis from Zymo and compare it with the read-trace composition.
2. DONE 2026-10-18 (see the next entry): Sylph E. coli over-count and L. fermentum deficit diagnosed.
3. Housekeeping: LICENSE; README 'research use only, not clinically validated'; stray cd.docx; delete Plan_Rev5/Rev6/Rev7, Guide_v2/v3/v4 and OneDrive transfer/d1..d21; copy PROJECT_STATE.md and HANDOVER.md from docs/ to the OneDrive project root; optionally the unused AWS key pair; confirm the AWS bill (about USD 2).
4. Optional later: fresh third pool (needs its own registered addendum), AMR design, Nextflow wrapping.

## 2026-10-18 later still (Sylph abundance diagnostic, dev only)
State: bin/14 (E. coli double counting: three database genomes share the reads; best-genome what-if 2.34 -> 1.15, post hoc), bin/15 (registered read-error hypothesis for the L. fermentum coverage deficit: REFUTED; identity 90-94% in all species, rho 0.23) and bin/16 (post hoc: mean read length differs by species, 1.9 kb for L. fermentum to 8 kb; with each species' own length Sylph's True_cov / coverage has median 1.05, so the L. fermentum 'under-count' was a bin/14 assumption artefact). Commits e372a22, f48c60a, 67c98ae; full text in PROJECT_STATE 'Sylph abundance error: EXPLORATORY diagnostic'. Registered verdicts unchanged. In bases L. fermentum is about 0.65 of expected, not 1.7 as in reads; Bray-Curtis of shares versus DNA expectation 0.120 on reads, 0.184 on bases.
bin/17 (commit 1c6d603) re-expressed Kraken2 in bases: Bray-Curtis 0.166 -> 0.191 (reads -> bases), predictions P1/P2 met, gap to Sylph 0.075 -> 0.051; Kraken2 stays ahead on abundance by a smaller margin; Kraken2's species-level reads are longer than the mapped reads (classification length bias). Next steps: (1) decide whether to run the same base-share view on the eval pool (needs its own registered addendum; do not run it without one) and whether to correct Plan Rev8 / Guide v5 (they still say the cause was not investigated); (2) record whether the email to Dr Zhang was sent and any replies (also ask about extraction / fragmentation if useful); (3) housekeeping items as before; (4) optional: AMR design, Nextflow wrapping, fresh third pool.
Session note: GitHub was linked to the Claude session for reading only; pushes still go through the user (WSL).

## 2026-10-18 evening (Addendum A2 registered, evaluation run pending)
State: Plan Rev8 and Guide v5 amended and committed (d1e535a). Addendum A2 (D19) written: eval-pool base-share view of Kraken2 PlusPF-8 baseline_c0 (bin/17 v2) and the Sylph best-genome rule (bin/18), predictions fixed from dev, guards require git tag bases-registered. Not covered: the 16 GB and full PlusPF classifications (AWS machine deleted; re-run about USD 2) and the confidence arms.
Next steps, in order: (1) copy bin/17 (v2, 11,467 bytes) and bin/18 (5,988 bytes) into bin/, run python3 bin/18_sylph_bestgenome.py --pool dev and keep its output; (2) commit PROJECT_STATE.md, HANDOVER.md, bin/17, bin/18 and reports/sylph/dev_bestgenome; (3) git tag bases-registered and push the tag; (4) run once: python3 bin/17_kraken_bases.py --run data/classify/eval/baseline_c0 --allow-eval --min-frac 0.0008 --compare-dev reports/kraken_bases/dev_baseline_c0_by_species.tsv and python3 bin/18_sylph_bestgenome.py --pool eval --allow-eval; (5) commit the outputs, tag bases-eval-scored; (6) record results against the predictions in PROJECT_STATE, then update Plan and Guide; (7) send the email to Dr Zhang (user, 19 October) and record any reply.

## 2026-10-18 late evening (Addendum A2 scored)
State: A2 run once on the evaluation pool after tag bases-registered (de94048): all eight predictions met (Kraken2 in bases 0.168 -> 0.191 Bray-Curtis; L. fermentum 1.98 -> 0.68; read length replicates, 13 of 14 species within 10%; Sylph best-genome rule 0.242 -> 0.182, E. coli 2.38 -> 1.17; Sylph best-genome vs Kraken2 in bases 0.182 vs 0.190). Full text in PROJECT_STATE 'ADDENDUM A2 RESULT'. Registered verdicts unchanged. Pools are from one run, so this is a within-run replication.
Next steps, in order: (1) user: commit the eval outputs (reports/kraken_bases/eval_*, reports/sylph/eval_bestgenome) with the updated docs and README, tag bases-eval-scored, push tags; (2) Plan Rev8 and Guide v5: add the A2 result (short paragraph each); (3) user sends the email to Dr Zhang (19 October) and records any reply; Zymo reply pending; (4) housekeeping: GitHub About/topics, OneDrive cleanup (needs the ls output), optional AWS key pair and bill; (5) optional later: re-run the 16 GB and full PlusPF classifications in bases (about USD 2), AMR design, Nextflow wrapping, fresh third pool or another library (needs its own registered addendum).

## 2026-10-18 night (Guide v6 after an external review)
State: a ChatGPT review of the guide was checked point by point (PROJECT_STATE 'External review of the Workflow Guide'); the damaged 5.4 table and the stale Section 7 status were real and are fixed; reads-versus-DNA example, glossary, cautious composition wording, overview, worked example (Roseburia), LCA/database note and the split of section 5.5 added; guide is now v6 (replaces v5). Decided against re-running the 16 GB and full-index Kraken2 classifications in bases (cloud re-run). Not done: appendix for the timeline and six-section restructure (larger change, later).
Next steps: (1) user: copy Guide v6, README, PROJECT_STATE, HANDOVER, git rm the v5 guide, commit; (2) email to Dr Zhang (19 October), record any reply; Zymo reply pending; (3) stray file 'main' in the repo root: look at it, then delete; (4) AMR design when time allows.


## 2026-10-18 late night (AMR design registered as Addendum A3)
- PROJECT_STATE.md now holds ADDENDUM A3 (AMR design, decision D20): reference catalogue from AMRFinderPlus on the 15 bacterial Zymo reference genomes (annotation, not biological truth), three categories, categorisation rules fixed before the hits are seen, detectability from carrier coverage before any false negative is declared, one read-level pilot arm (ResFinder planned) on the dev pool, assembly and Argo deferred. A second-agent review was adopted.
- bin/19_amr_reference.py written (step 1; not yet committed or run). Needs on the user's WSL: conda install -c bioconda -c conda-forge ncbi-amrfinderplus; amrfinder -u. The shell in my session has no network access to NCBI.
- NEXT: user commits docs + bin/19, tags amr-registered, runs bin/19, sends reports/amr/reference_hits_raw.tsv and reference_run.json; then bin/20 (categories) and addendum A4 (read-level method, thresholds, evaluable minimum).
- Still open from before: email to Dr Zhang (19 Oct), Zymo reply, stray empty file `main`.

## 2026-10-18 night II (AMR reference catalogue done; Addendum A4 registered)
- bd61f30: AMRFinderPlus 4.0.23 / db 2025-07-16.1 on the reference genomes: 8 high-confidence and 30 candidate reference ARGs (see PROJECT_STATE 'AMR reference catalogue RESULT'). amr-registered tag is on 99629f3 (A3 + bin/19).
- A4 (this commit; tag amr-a4-registered): read-level pilot design with ResFinder, thresholds (identity 0.80, gene coverage 0.80, >= 3 reads, evaluable E >= 10), family mapping rule, protocol. Nothing run.
- NEXT: user obtains ResFinder DB, records commit/date in reports/amr/resfinder_version.txt; I write bin/21_amr_reads.py (family map first, then dev run); eval only after a freeze addendum.
- Still open: email to Dr Zhang (19 Oct), Zymo reply, stray `main` file.

## 2026-10-18 night III (AMR read-level pilot, dev result)
- 7c712a4 (tag amr-map-registered): bin/21 and family map committed before the read run. The dev run followed: 2 of 3 evaluable high-confidence families detected (tet(Q), tet(W)); cepA missed at 2 reads with 27 ambiguous reads (cause not established); see PROJECT_STATE 'AMR read-level pilot: DEV RESULT'. 6 of 8 high-confidence genes are in ResFinder; most candidate genes are not.
- NEXT: look at reports/amr/dev_family_reads.tsv for the families sharing the ambiguous cepA reads (user sends it); decide whether an addendum A5 is needed (family rule or ambiguity rule) BEFORE any freeze; then a freeze addendum and one eval run. The eval pool has not been used for AMR.
- Still open: email to Dr Zhang (19 Oct), Zymo reply, stray `main` file, AMR dev outputs to commit if not yet done (reports/amr/dev_*).

## 2026-10-19 (cepA ambiguity explained; Addendum A5 registered)
- Email to Dr Zhang sent 19 October (no reply yet; Zymo reply pending).
- bin/22 (89f249f, hypothesis committed before the run): the four ResFinder names cepA, cepA-29, -44, -49 are 98.8-99.7% identical; 29 reads in the group; P1-P3 met. The A4 cepA miss was a naming artefact of the family rule.
- A5 (D22) registered: explicit merge list reference/amr_family_merge.tsv; bin/21 gains --merge/--suffix; A4-rule outputs kept; predictions for the A5 dev run fixed beforehand. The A5 dev run is NOT an independent test (the change was motivated by dev).
- NEXT: user commits docs + bin/21 + reference/amr_family_merge.tsv, tags amr-a5-registered, runs the two commands in A5; then a freeze addendum (A6) and ONE eval pass (needs an eval AMR path in bin/21 with guards: --allow-eval, tag check, never overwrite; I write it); then Plan Rev9 / Guide v7 with the AMR work.
