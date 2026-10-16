# PROJECT_STATE: decision log and findings

Last updated: 2026-10-16 (late afternoon, database-size arms REGISTERED, AWS set-up in progress)

## Status
Phase 0 (dataset audit) essentially complete; Phase 1 (baseline benchmark) started. Done: run record, chemistry/basecaller, pools sampled and QC'd,
theoretical composition checked, baseline database and classifier frozen (D12), baseline run on the dev pool, first scoring done (sweep on dev, see Scoring). Truth-table naming issue found and fixed by alias (D13), sweep rerun with it. Evaluable minimum fixed at 20 (D14a). Threshold FROZEN at 0.08% of reads (D14b), chosen on the dev pool by bootstrap-mean F1.
Eval pool classified and scored ONCE with the frozen threshold (result below, within the written predictions). Depth series run on the eval pool and compared with the registered predictions (all met). Phase 2 arm 1 (Kraken2 confidence, D16) done: registered (tag arms-registered, f3c9ba7), dev-tuned and frozen (tag arms-frozen, 15e92f3), eval scored once per arm (result committed a427e54); the baseline (confidence 0) is best; result and comparison with the registered expectations below.
Still open, non-blocking: lot-specific CoA, basecalling model tier, Argo paper link. Plan is at revision 5 (16 October; Rev4 kept). Workflow guide is at version 2 (16 October; OneDrive only; new section 4 explains thresholds and scoring for a novice).

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
- D13 (2026-10-16) Truth-table alias: reads classified as Fusobacterium animalis (taxid 76859) count as the D6331 "Fusobacterium nucleatum" strain (taxid 851).
  Found after the first dev scoring (F. nucleatum 92 reads vs 11,100 expected; F. animalis 8,597; genus Fusobacterium 9,300). Evidence independent of the
  reads: the D6331 reference genome cut into 489 x 5 kb pieces and classified with the frozen baseline gave 475 (97.1%) F. animalis, 3 F. nucleatum. Stored in
  reference/d6331_taxid_aliases.csv with the evidence; 04_score.py reports every alias it applies. Caveat: found after seeing results, strain not checked against an
  assembly accession (the FASTA header is only the file name). Not a classifier error: a taxonomy-name mismatch in the truth table, same family as the Lactobacillus rename.
- D14a (2026-10-16, user) Evaluable-species minimum = 20 expected reads (set before results, kept). Salmonella (18.5 expected) is "detected, not evaluable".
- D14b (2026-10-16, FROZEN, user accepted the rule) Detection threshold = a FRACTION of all reads, chosen on the dev pool as the candidate with the highest block-bootstrap MEAN F1 (not the
  single-pool point value); ties go to the SMALLER threshold (keeps real low-abundance species). Candidates 0.0002, 0.0005, 0.0008, 0.001, 0.002, 0.005 (the list includes values picked after
  seeing the sweep, all on dev). Result: 0.0008 (mean F1 0.948; 0.005 is 0.947, a tie, and the tie rule agrees). FROZEN: --min-frac 0.0008, i.e. T = ceil(0.0008 x reads) = 148 on a 185,000-read pool.
  No further tuning on the eval pool. Strict truth: any species not in the 17 (including human) counts as a false positive; no exclusions after the fact.
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
- Baseline run on the dev pool (label baseline_c0; confidence 0.0, minimum-hit-groups 2, 6 threads, docker memory cap 10g): 185,000 reads / 1,015.91 Mbp,
  174,153 classified (94.14%), 10,847 unclassified (5.86%); Kraken2 classification 14.0 s, wall 17 s including database load; no out-of-memory kill (supports D4).
  Outputs: data/classify/dev/baseline_c0/ (report.txt, kraken_output.txt.gz, run.json, run.log). Evaluation pool untouched.

## Scoring (2026-10-16) - bin/04_score.py
- Truth matched by species taxid (reference/d6331_expected_reads_by_species.csv); strains roll up to species; expected reads = total reads x genomic-DNA fraction.
- Called species = clade reads at species rank >= T. Precision = called mock species / called species (a called mock species is never a false positive). Recall is over
  EVALUABLE species only (expected reads >= --evaluable-min-expected, default 20, PROPOSED not yet decided); recall over all 17 is also printed. F1 from precision and evaluable recall.
  Abundance: L1 and Bray-Curtis of observed read fractions among called species vs expected DNA fractions.
- False-positive audit lists every non-mock species called, categorised (same genus as a mock species / other bacteria / archaea / eukaryote / virus).
- Uncertainty: bootstrap over the pool's 5,000-read blocks (37 per pool), using the per-read Kraken2 output, the manifest and nodes.dmp; cross-checks report vs per-read species counts.
- Guard: eval pool only with --allow-eval and an explicit --min-reads (no sweep on eval). Sweep mode runs on dev only.
- Tested on a synthetic dataset with hand-derived answers (precision 4/7, 4/5, 1.0 at T=1, 3, 5; evaluable recall 1.0; recall over all 0.75 at T>=5). Real-data run: see below.
- Added 2026-10-16: --min-frac / --sweep-frac (threshold as a fraction of ALL reads, T = ceil(frac x N), so it keeps its meaning across the depth series); --aliases (D13);
  genus-level recall as a secondary figure (genus clade reads >= T); sweep grid widened to 2000. Regression-tested: identical metrics and bootstrap to the previous version when no alias file is present.

### First real scoring (dev pool, baseline_c0, sweep, BEFORE the alias; 2026-10-16)
- 185,000 reads; 174,153 classified; 146,319 reads (79.09%) in the 17 mock species; 15.05% classified elsewhere. 14 of 17 species evaluable at the proposed minimum 20
  (not evaluable: Salmonella enterica 18.5 expected, Enterococcus faecalis 1.9, Clostridium perfringens 0.2).
- Sweep (T = min reads): called / FP / precision / evaluable recall / F1 = T1 497/481/0.032/1.000/0.062; T5 71/56/0.211/1.000/0.349; T20 36/21/0.417/1.000/0.588;
  T50 27/12/0.556/1.000/0.714; T100 17/4/0.765/0.929/0.839 (loses Fusobacterium, 92 reads). Bray-Curtis about 0.22 at every T (composition error is not a thresholding effect).
  F1 was still rising at the edge of the grid (T=100), so the optimum had not been reached.
- Observed/expected (T=1): F. prausnitzii 1.23, V. rogosae 0.87, R. hominis 0.42, B. fragilis 0.95, E. coli 0.57, P. corporis 0.94, B. adolescentis 0.22, F. nucleatum 0.01 (see D13),
  L. fermentum 1.67, C. difficile 0.96, A. muciniphila 0.97, C. albicans 0.95, S. cerevisiae 0.75, M. smithii 0.82, Salmonella 3.03, E. faecalis 0 (missed), C. perfringens 5.4.
- Where reads went (genus clade vs species clade): Roseburia 10,914 vs 10,882 (shortfall is above genus or unclassified); Bifidobacterium 5,737 vs 2,421 for adolescentis
  (about 3,000 reads stay at genus level: shared between close relatives); Escherichia 15,956 vs 14,646 with close relatives (Shigella, E. albertii, E. fergusonii) taking some.
  Limosilactobacillus genus 18,560 vs fermentum 18,529 (1.67x: possible absorption from relatives or composition difference; not investigated).
- Note: fp_audit.tsv is written only in --min-reads / --min-frac mode, not in sweep mode.

### Sweep after the D13 alias (dev pool, baseline_c0; 2026-10-16)
- Reads in the 17 mock species: 154,916 (83.74%); classified elsewhere 10.4%; Fusobacterium nucleatum now 8,689 (obs/exp 0.78); Bray-Curtis about 0.17 at every T.
- T (reads, % of reads): called / FP / precision / evaluable recall / F1: 50 (0.027%) 26/11/0.577/1.000/0.732; 100 (0.054%) 17/3/0.824/1.000/0.903;
  150 (0.081%) 15/1/0.933/1.000/0.966 (maximum); 200-500 14/1/0.929/0.929/0.929; 1000 (0.54%) 13/0/1.000/0.929/0.963; 2000 (1.08%) 12/0/1.000/0.857/0.923.
- Fraction grid run (pre-listed): 0.0001 -> T 19, F1 0.577; 0.0005 -> 93, 0.875; 0.001 -> 185, 0.929; 0.005 -> 925, 0.963; 0.01 -> 1850, 0.963.
- KNIFE EDGE: the F1 maximum (T=150) lies one read below Methanobrevibacter smithii (151 observed reads, 0.1% DNA). Poisson noise alone is about +/-12 reads, so on another pool this
  species would be missed about half the time. The two near-equal options: T about 150 (all 14 evaluable found, 1 false positive, fragile) or T about 925 or more (0 false positives,
  loses M. smithii by design; F1 0.963 vs 0.966). The F1 difference is far smaller than the sampling noise. Ties-to-larger would push the choice to 1% of reads, which is silly for a mock whose
  lowest evaluable species is 0.1%, hence the proposal to break ties towards the smaller threshold and to compare candidates by bootstrap mean F1.
- The 1 remaining false positive at T=148 is Homo sapiens: 923 reads (0.5% of reads), 8,053 distinct minimizers. The mock has no human DNA, so this is host/lab contamination or a
  misassignment; not resolved (many distinct minimizers lean towards real human reads, but this is not shown). Kept as a false positive (strict truth). Note the 0.5% option (T=925) is zero-FP
  only because human has 923 reads: both candidate thresholds sit within a few reads of an edge (M. smithii 151 vs 148; human 923 vs 925).

### Eval-pool result (2026-10-16, scored once, --min-frac 0.0008 = T 148, no tuning; data/score/eval/baseline_c0_f0.0008/)
- 185,000 reads; 174,323 classified (94.23%); 155,279 (83.93%) in the 17 mock species; 10.29% classified elsewhere; Fusobacterium alias moved 8,388 reads.
- Called 15 species, 14 correct, 1 false positive (Homo sapiens, 818 reads; dev had 923). Evaluable recall 14/14, recall over all 17 = 0.824. Precision 0.933, F1 0.966 (bootstrap mean 0.962, 95% 0.929-0.966), Bray-Curtis 0.168.
  Report vs per-read cross-check identical. Methanobrevibacter smithii 163 reads (dev 151), detected in 91% of resamples. Genus-level recall 14/14.
- Obs/exp vs dev: F. prausnitzii 1.24, Roseburia 0.41, E. coli 0.57, B. adolescentis 0.22, F. nucleatum 0.76, L. fermentum 1.66, M. smithii 0.88; Salmonella 59 reads (3.2x, below T, not evaluable).
- Against the written predictions: recall 14/14 (predicted 13-14) OK; 1 FP (1-2) OK, human is the FP as expected; F1 0.966 (0.90-0.97) OK, at the top of the range; Fusobacterium 0.76x (about 0.78) OK;
  Roseburia, Bifidobacterium, E. coli the same as dev OK. The result is a within-run held-out check only (same sample, run, chemistry): it shows the threshold was not over-fitted to the dev blocks, not that it transfers
  to other runs or chemistries. Both pools put M. smithii within 3 to 15 reads of the threshold (151, 163 vs 148), so this species stays fragile at this depth. Human at about 0.45-0.5% of reads in both pools
  looks systematic (library or database), not random; origin not resolved.
- Baseline conclusions for Kraken2 PlusPF-8 (confidence 0): presence/absence is good at a sensible threshold (14/14 evaluable, 1 FP), but abundance is biased: Bray-Curtis 0.17; Roseburia 0.41x,
  Bifidobacterium adolescentis 0.22x, E. coli 0.57x low (reads stop above species or go unclassified) and Limosilactobacillus fermentum 1.66x high. Whether the 8 GB capped index causes the low values is the Phase 2 question (large database on AWS).

### Threshold choice on dev (bootstrap, 1000 resamples of 37 blocks, seed 20261016)
- frac (T) : precision / F1 (point; mean; 95% interval) / false positives (mean)
  0.0002 (37): 0.536 / 0.698; 0.690; 0.667-0.714 / 13.5
  0.0005 (93): 0.778 / 0.875; 0.874; 0.824-0.903 / 4.1
  0.0008 (148): 0.933 / 0.966; 0.948; 0.897-0.966 / 1.07  (14/14 evaluable detected; M. smithii detected in 58% of resamples)
  0.001 (185): 0.929 / 0.929; 0.929; 0.929-0.929 / 1.0  (M. smithii lost: 13/14)
  0.002 (370): same as 0.001
  0.005 (925): 1.000 / 0.963; 0.947; 0.929-0.963 / 0.46  (M. smithii lost)
- Predictions written BEFORE the eval run (so they can be checked): with 0.0008, expect evaluable recall 13/14 or 14/14 (M. smithii about a coin flip), 1 to 2 false positives (human likely among them),
  F1 roughly 0.90 to 0.97, Fusobacterium about 0.78x, Roseburia/Bifidobacterium/E. coli recoveries similar to dev. A result far outside this is a finding, not something to tune away.
- Depth note: the fraction transfers across depth, but at low depth the evaluable-species rule (expected reads >= 20 at that depth) removes the rarest species, so recall is reported over fewer species.
- Open (D14), proposals not yet accepted by the user: (a) keep the evaluable minimum at 20 expected reads (pre-specified before results; Salmonella, 18.5 expected, is then "detected, not evaluable");
  (b) choose the threshold as a FRACTION of reads (not a read count) on the dev pool, maximising F1 over a wide grid, ties to the larger value, after the D13 alias is applied; report the whole curve with
  bootstrap intervals; freeze before touching eval.

## Depth series (design registered 2026-10-16, BEFORE running; bin/05_depth_series.py)
- D15 (2026-10-16, proposed by Claude, default unless the user changes it): the depth series is drawn from the EVALUATION pool. Reason: the rules are frozen (D14), nothing will be changed in response, and drawing from dev would be
  optimistic because the threshold was tuned on dev blocks. The full-pool eval result (baseline_c0) is the 1.0 Gb baseline. A second series on dev can be added later as a labelled "tuning pool, optimistic" comparison.
- Design: depths 0.1, 0.25, 0.5 Gb (D9), 5 replicates each, random reads without replacement until bases reach the target, seed 20261016 with per-(pool, depth, replicate) seeds recorded in data/depth/<pool>/depth_manifest.json.
  Replicates overlap (about half at 0.5 Gb) and come from one pool: spread = within-pool variation only (D10). Classified with the frozen baseline (03_classify.py defaults), scored with --min-frac 0.0008 --bootstrap 0
  (block bootstrap does not apply to random-read subsamples). Evaluable minimum 20 expected reads at EACH depth (D14a). No sweep, no retuning, no rule changes.
- Expected sizes: about 18,500 / 46,000 / 92,500 reads; threshold T about 15 / 37 / 74 reads; evaluable species 13 / 14 / 14 (Methanobrevibacter smithii expects about 18.5 reads at 0.1 Gb, so it drops out there).
- Predictions written before the run: (1) evaluable recall 13/13 in all replicates at 0.1 Gb; at 0.25 and 0.5 Gb M. smithii (observed 0.088% of reads vs threshold 0.08%) detected in roughly 70-80% of replicates, so recall 13/14 or 14/14;
  (2) false positives 1 to 2.5 on average (human always present at about 0.45%), at most about 4 in any run, possibly a little higher at 0.1 Gb because counts are small; (3) F1 mean between 0.90 and 0.97 at every depth;
  (4) Bray-Curtis about 0.17 (+/- 0.01) at every depth; (5) per-species obs/exp the same as at full depth (for example F. nucleatum about 0.76, B. adolescentis about 0.22). Anything far outside is a finding, not something to tune away.
- Cost: subsampling reads the pool once (a minute or two); 15 classifications at about 17 s each.

### Depth series result (eval pool, 2026-10-16; data/depth/eval/; frozen threshold 0.0008 of reads; each replicate classified in 4-25 s)
| Gb | runs | reads | T | evaluable | precision mean (min-max) | evaluable recall mean (min-max) | F1 mean (min-max) | FP mean (min-max) | Bray-Curtis |
| 0.1 | 5 | 18,298 | 15 | 13 | 0.910 (0.824-0.933) | 1.000 (1.000-1.000) | 0.953 (0.903-0.966) | 1.4 (1-3) | 0.171 |
| 0.25 | 5 | 45,890 | 37 | 14 | 0.932 (0.929-0.933) | 0.986 (0.929-1.000) | 0.958 (0.929-0.966) | 1.0 (1-1) | 0.169 |
| 0.5 | 5 | 91,318 | 74 | 14 | 0.933 (0.933-0.933) | 1.000 (1.000-1.000) | 0.966 (0.966-0.966) | 1.0 (1-1) | 0.168 |
| 1.0 | 1 (full pool) | 185,000 | 148 | 14 | 0.933 | 1.000 | 0.966 | 1 | 0.168 |
- Predictions vs result: all five met. (1) 0.1 Gb recall 13/13 in every replicate OK; M. smithii at 0.25 Gb detected in 4/5 (predicted about 70-80%) OK, at 0.5 Gb in 5/5 (a little above the predicted 70-80%, but 5 replicates cannot tell 80% from 100%; the
  lowest count was 74 against T=74, so it was detected by one read); (2) FP mean 1.4 / 1.0 / 1.0 (predicted 1-2.5, max about 4) OK, maximum 3; (3) F1 means 0.953 / 0.958 / 0.966 (predicted 0.90-0.97) OK; (4) Bray-Curtis 0.171 / 0.169 / 0.168 OK;
  (5) obs/exp unchanged (F. nucleatum 0.73-0.78, B. adolescentis 0.20-0.23, Roseburia 0.40-0.43, E. coli 0.56-0.59, L. fermentum 1.65-1.69) OK.
- Findings: (a) with a threshold that is a fraction of reads, presence/absence performance is flat from 1 Gb down to 0.1 Gb (about 18,000 reads) for this mock; the cost of lower depth is only that the rarest species drop out of
  evaluability (M. smithii at 0.1 Gb: 6-20 reads against T=15; detected in 4 of 5 although not evaluable). (b) 14 of 15 runs have exactly one false positive, and it is human in every run inspected (about 0.35-0.5% of reads at every depth: 62-75 reads at 0.1 Gb,
  190-232 at 0.25 Gb, about 400 at 0.5 Gb): it is systematic. (c) The extra false positives (Escherichia albertii 17 reads, Shigella flexneri 17 reads) occurred only in one 0.1 Gb replicate, where T=15 reads is a low absolute bar:
  relatives of E. coli cross a small absolute threshold by chance. (d) Abundance bias (Roseburia, Bifidobacterium, E. coli low; Limosilactobacillus high) does not depend on depth, so it is a database/classifier property, not noise.
- Caveats: replicates are drawn from one 1 Gb pool and overlap (about half at 0.5 Gb), so the very small spreads at 0.5 Gb are partly an artefact of overlap; they measure within-pool variation only (D10). The mock has 14 species at >=0.1% DNA;
  a clinical sample with a dominant host or few abundant species would behave differently. One run, one chemistry.

## Scoring script fix (2026-10-16, 04_score.py)
- Edge case found by the confidence arms: when NO species reaches the threshold in any bootstrap resample, the interval code crashed (empty list). Fixed: F1 is now 0 (not undefined) when nothing is called or nothing called is right (recall is 0);
  precision stays undefined (nan) in that case; intervals of empty lists are nan. Regression-tested: identical metrics, bootstrap and sweep output to the previous version whenever at least one species is called (all earlier results unaffected).

## Confidence arms (Phase 2, arm 1; REGISTERED 2026-10-16 before running; bin/06_confidence_arms.py)
- D16 (2026-10-16, user chose this as the next experiment): test Kraken2 --confidence as a cheap lever against false positives. Arms: 0.0 (the frozen baseline, threshold stays 0.0008) plus 0.05, 0.1, 0.2, 0.4. Everything else as in D12.
- Registration note (honest record): the first registration commit FAILED because results/ is gitignored (Nextflow scaffold), yet the tag arms-registered was created on the previous commit; the tag was deleted and recreated after the real
  commit, and tracked outputs now live in reports/ (not results/). Before that commit the dev stage had already been started. It crashed at confidence 0.1, fraction 0.002 (threshold 370 reads) because NO species reached 370 reads at species
  rank, while scoring at 0.001 (threshold 185) had completed. So the registration is partly informed by one glimpse: species-level counts collapse at confidence 0.1, which probably breaks expectation (4). The expectation text below was
  written before any run and is unchanged.
- Rule (same as D14b, made explicit): for each non-baseline arm, score the dev pool with the block bootstrap at fractions 0.00005, 0.0001, 0.0002, 0.0005, 0.0008, 0.001, 0.002, 0.005 and choose the SMALLEST fraction whose bootstrap-mean F1 is within 0.002 of that arm's maximum.
  Selections are written to reports/arms/dev_selection.json, committed and tagged `arms-frozen` BEFORE the eval pool is touched (the script refuses otherwise). Then each arm is scored ONCE on eval at its own frozen fraction. ALL arms are reported;
  no arm is chosen on eval. The 0.0 arm is not re-tuned (its dev optimum on the wider grid is shown for reference only).
- Why F1 alone will not separate the arms: with human DNA as an unavoidable false positive, the F1 ceiling on this mock is 0.966 (14 true calls, 1 false). The informative columns are therefore: the smallest usable fraction (the lowest
  abundance the arm can call without false positives), reads classified, reads in the 17 mock species, Bray-Curtis, per-species obs/exp (Roseburia, Bifidobacterium adolescentis) and the false-positive list.
- Expectations written before the run (low confidence, mechanism: Kraken2 confidence is the share of a read's k-mers supporting a clade; nanopore reads at about Q18 have many k-mers with errors, so high values push reads up to genus or leave them unclassified):
  (1) reads classified (any rank) and reads in the mock species at species rank both fall monotonically with confidence, c0.4 clearly lowest (below about 80% classified);
  (2) the dev-chosen fraction for c0.05 to c0.2 is at or below 0.0008 for at least one arm (noise falls faster than signal), and false-positive species other than human at a fixed low fraction (0.0002) fall monotonically with confidence;
  (3) human stays a false positive in every arm whose chosen fraction is below about 0.4%;
  (4) evaluable recall stays 14/14 for c <= 0.1 and may drop at c >= 0.2 (rare species lose species-level reads); (5) the already-low species (Bifidobacterium adolescentis 0.22, Roseburia hominis 0.42) fall further as confidence rises;
  (6) Bray-Curtis changes by less than 0.03 for c <= 0.1. A result far outside this is a finding. The eval pool is looked at once per arm; this is a second look at it overall, so it is reported as a pre-registered comparison, not as fresh held-out data.

### Confidence arms RESULT (2026-10-16; dev tuned, frozen at tag arms-frozen = 15e92f3, eval scored once per arm; summary in reports/arms/arms_summary.txt, commit a427e54)
Dev tuning, bootstrap-mean F1 (grid fractions 5e-5 ... 0.005; chosen = smallest within 0.002 of the arm's maximum):
- c0 (baseline, frozen 0.0008): 0.479, 0.584, 0.690, 0.874, 0.948, 0.929, 0.929, 0.947 (re-tuned would also be 0.0008)
- c0.05: 0.924, 0.901, 0.923, 0.880, 0.839, 0.833, 0.783, 0.600; chosen 5e-5 (T=10), the EDGE of the registered grid, so the true optimum may lie lower (not explored; would need a labelled exploratory run on dev only)
- c0.1: 0.604, 0.496, 0.251, 0.133, 0.131, 0.040, 0, 0; chosen 5e-5. c0.2 and c0.4: 0 everywhere (no species reaches even 10 reads at species rank); chosen 5e-5 by the tie rule (meaningless).
EVAL pool at each arm's frozen fraction (all arms reported):
| conf | frac (T) | classified % | reads in mock species % | precision | recall (evaluable) | F1 | F1 boot mean (95%) | FP | Bray-Curtis | genus recall |
| 0 | 0.0008 (148) | 94.2 | 83.9 | 0.933 | 14/14 (1.000) | 0.966 | 0.962 (0.929-0.966) | 1 (human 818) | 0.168 | 1.000 |
| 0.05 | 5e-5 (10) | 57.5 | 19.4 | 0.929 | 13/14 (0.929) | 0.929 | 0.943 (0.929-0.966) | 1 (human 21) | 0.336 | 1.000 |
| 0.1 | 5e-5 (10) | 0.4 | 0.2 | 1.000 | 6/14 (0.429) | 0.600 | 0.596 (0.571-0.600) | 0 | 0.516 | 0.500 |
| 0.2 | 5e-5 (10) | 0.0 | 0.0 | n/a | 0 | 0 | 0 | 0 | n/a | 0 |
| 0.4 | 5e-5 (10) | 0.0 | 0.0 | n/a | 0 | 0 | 0 | 0 | n/a | 0 |
Dev (same columns, abbreviated): c0 F1 0.966 (boot mean 0.948), c0.05 0.929 (0.924), c0.1 0.600 (0.604), c0.2/c0.4 0. Dev and eval agree for every arm.
Reading: confidence does not help at species rank on this data. c0.05 reaches a similar F1 only because its threshold falls to 10 reads; it moves 37% of reads above species rank (classified 94 -> 57%), cuts reads at species rank in mock species from 84% to 19%
and doubles Bray-Curtis (0.17 -> 0.34), so abundance estimates get much worse for no gain in presence/absence. c0.1 and above are unusable at species rank. Genus recall is 1.000 for c0 and c0.05 and 0.500 on eval for c0.1; whether a high confidence is useful at GENUS level is open (not tested as a separate arm).
Human DNA is a false positive only where anything is classified at species rank with a threshold below its read count (c0, c0.05); it is not called in c0.1 to c0.4 because almost nothing is classified.
Comparison with the registered expectations (text above unchanged):
- (1) MET: classified 94.2, 57.5, 0.4, 0.0, 0.0 and mock species at species rank 83.9, 19.4, 0.2, 0.0, 0.0 (monotone); c0.4 far below 80%.
- (2) PARTLY MET: chosen fraction is at or below 0.0008 in every arm (all 5e-5), but at the grid edge. The second half (non-human false positives at fraction 0.0002 fall monotonically) is NOT checked from the summary; the summary shows no non-human false positive at any chosen fraction.
- (3) PARTLY CONTRADICTED: human is a false positive for c0 and c0.05, but not for c0.1, c0.2, c0.4 (chosen fraction below 0.4% in each) because almost nothing is classified at species rank.
- (4) CONTRADICTED: evaluable recall is 14/14 only for c0; c0.05 13/14 and c0.1 6/14 (the registration note recorded that the dev crash at c0.1 probably breaks this).
- (5) MET: eval obs/exp Roseburia hominis 0.42 (c0) -> 0.18 (c0.05, 4,725 of 25,900 reads) -> 0.001 (c0.1, 33 reads); Bifidobacterium adolescentis 0.22 -> 0.05 (511 of 11,100) -> 0.0001 (1 read, missed). At c0.1, 8 of 14 evaluable species are missed (E. coli, F. nucleatum, C. difficile, B. adolescentis, M. smithii, S. cerevisiae, C. albicans, A. muciniphila).
- (6) CONTRADICTED: Bray-Curtis 0.168 -> 0.336 (c0.05) -> 0.516 (c0.1).
Verdict: baseline c0 at 0.0008 stays the frozen baseline. The eval pool has now been looked at once per arm; do not re-score with other settings.

## Database-size arms (Phase 2, arm 2; REGISTERED 2026-10-16 before any run; bin/07_database_arms.py; runs on AWS)
- D17 (user chose AWS as the next experiment): does the 8 GB cap on the Kraken2 index cause the abundance bias (Roseburia hominis 0.42x, Bifidobacterium adolescentis 0.22x, E. coli 0.57x, Limosilactobacillus fermentum 1.66x at the frozen baseline)?
- Arms, all from the same Langmead index build (2026-06-26, collection PlusPF): PlusPF-8 (7.5 GB index; baseline and CONTROL), PlusPF-16 (14.9 GB), PlusPF full (109.9 GB). Everything else identical to D12: kraken2 image digest db5d772d..., --confidence 0.0, --minimum-hit-groups 2, the same dev and eval pools, the same 04_score.py and truth/alias files.
- Machine: one EC2 instance in eu-north-1 (Stockholm), r6i.4xlarge (128 GB RAM, 16 vCPU) or r5.4xlarge, Ubuntu 24.04, 200 GB gp3; terminated and its volume deleted after the run. Budget alert exists (account budget $15, alerts at 50% and 85%). Estimated cost USD 5-8.
- Control: PlusPF-8 on the instance must reproduce the local dev baseline species counts exactly (columns taxid to genus_clade_reads of reports/baseline_kraken2_pluspf8/dev_species_table.tsv); the script refuses the eval stage if it does not. This checks that machine, threads and image do not change results.
- Rule (same as D14b/D16): PlusPF-8 keeps the frozen fraction 0.0008. For PlusPF-16 and PlusPF full, score the dev pool with the block bootstrap at fractions 0.00005 ... 0.005 (the D16 grid) and choose the SMALLEST fraction whose bootstrap-mean F1 is within 0.002 of that arm's maximum.
  Selections go to reports/dbarms/dev_selection.json, committed and tagged `dbarms-frozen` BEFORE the eval pool is touched (the script refuses otherwise). Then each arm is scored ONCE on eval at its frozen fraction, and also at 0.0008 for a like-for-like column. ALL arms are reported.
- Decision rule for the question (written before the run; applied to obs/exp at fraction 0.0008 on the EVAL pool, and checked for the same verdict on dev): 'the cap is the main cause' if Roseburia hominis AND Bifidobacterium adolescentis both reach >= 0.70 with PlusPF full and rise monotonically 8 -> 16 -> full;
  'the cap is not the main cause' if both stay <= 0.50 with PlusPF full; anything else is 'partly' (reported as such). If dev and eval give different verdicts the result is 'inconclusive'.
- Expectations written before the run (low confidence):
  (1) reads classified at any rank and reads in the mock species at species rank do not fall with index size (8 -> 16 -> full), and PlusPF full classifies >= 95% of reads;
  (2) Roseburia hominis and B. adolescentis obs/exp rise monotonically with index size; PlusPF full lands in 'partly' (Roseburia 0.5-0.8; B. adolescentis stays below 0.5 because many Bifidobacterium genomes share k-mers and reads stop at genus level);
  (3) Bray-Curtis for PlusPF full is within +-0.03 of the baseline 0.168 or lower;
  (4) evaluable recall stays 14/14 in every arm at 0.0008;
  (5) human remains a false positive in every arm, so F1 stays at or near the ceiling 0.966 (all arms within 0.03); non-human false positives at 0.0008 are 0-2 for PlusPF full (a larger database can add low-level misassignments) and 0 for the others;
  (6) the dev-chosen fraction for PlusPF-16 and PlusPF full lies between 0.0002 and 0.001, not at the edge of the grid.
  Wall-clock per pool and the instance cost are recorded as results (the practical trade-off), without a prediction. The eval pool is looked at once per arm; overall this is another look at it, so the result is reported as a pre-registered comparison, not as fresh held-out data.
- Run order on the instance: install docker and git, clone the repo, download the three indexes from s3 genome-idx (kraken/k2_pluspf_08_GB_20260626, k2_pluspf_16_GB_20260626, k2_pluspf_20260626), upload both pools from WSL, run 07 --stage dev, commit and tag dbarms-frozen, run 07 --stage eval --allow-eval, copy reports/dbarms back, terminate.

## Unverified / open
- Argo paper link and basecalling model tier (see below)
- Lot-specific CoA for the sequenced tube: lot number unknown. Ask Zymo (tech@zymoresearch.com) or the Argo paper authors. Do not block the project on it; label results as against theoretical composition.
- Basecalling model tier (fast/hac/sup) not stated in the ENA record; look in the Argo paper methods or ask the authors.
- Argo paper link: reported as confirmed in the paper's data availability statement; not verified by me.
- Decisions 1 (dataset confirmation) and 2 (threshold rule, evaluable-species minimum) not yet formally settled (see Scoring section)
- LomanLab ENA accessions (if any)
- PRJNA751627, PRJEB41184, PRJNA1335824, PRJNA554461 (clinical leads, deferred)
- Disk: 647 GB free at 2026-10-15 (df, after the database download)

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
