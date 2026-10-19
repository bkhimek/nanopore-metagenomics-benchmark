# Nanopore metagenomics: benchmark-led taxonomic profiling, with a pilot resistance-gene check

An independent portfolio project: a reproducible, pre-registered benchmark of taxonomic profiling for nanopore
metagenomic reads, run on a public, characterised mock community. The aim is a defensible picture of where the
workflow works, where it fails and why, not a pipeline for its own sake.

**Research use only. Not clinically validated.** The project makes no claim about diagnostic sensitivity or
specificity, and it does not test laboratory enrichment chemistry. Results apply to the chemistry and basecaller of
the data used.

## What has been done

- **Data:** ZymoBIOMICS Gut Microbiome Standard (D6331) sequenced on a PromethION (Kit 14, Dorado v7.2.13), BioProject
  PRJNA1181840, run SRR31222153 (23.4 M reads, 128.4 Gb; data from Chen et al., *Nature Communications* 16:1744, 2025).
  The full run is not downloaded. A seeded, block-stratified sample is fetched from SRA and split into two disjoint
  pools (development and evaluation), each about 1 Gb (185,000 reads). Both pools come from one run, so the held-out
  test is within-run, not independent.
- **Ground truth:** the manufacturer's theoretical composition (17 species), not a lot-specific measurement. The lot
  number is unknown; the lot certificate of analysis has not been obtained.
- **Method:** every experiment was registered before it was run (arms, rules for choosing settings, written
  predictions), committed and tagged in Git. Settings were tuned on the development pool only, frozen by tag, and the
  evaluation pool was scored once per arm. All arms are reported, including the ones that failed. Evaluation-pool
  results are therefore pre-registered comparisons on a pool that earlier experiments had already examined, not fresh
  held-out data.

## Results in brief (evaluation pool)

| Experiment | Outcome |
|---|---|
| Baseline: Kraken2 2.17.2, PlusPF-8 index (2026-06-26) | 14 of 14 evaluable species found; human DNA the only false positive; F1 0.966 (the ceiling with human counted as wrong); Bray-Curtis 0.17 against the theoretical composition |
| Depth series (0.1, 0.25, 0.5 Gb) | F1 0.953, 0.958, 0.966; presence/absence flat down to 0.1 Gb |
| Kraken2 confidence 0 to 0.4 | Strictness hurts on ~Q18 nanopore reads; the default (0) was best |
| Database size (8 GB, 16 GB, full PlusPF) | F1 unchanged at 0.966; Bray-Curtis 0.168, 0.153, 0.137 |
| Read trace (exploratory, development pool) | Much of the departure from the theoretical abundances is in the sequenced sample itself; Kraken2 loses reads mainly for species with close relatives (E. coli, Bifidobacterium adolescentis) |
| Sylph 1.0.0 with GTDB r232 + RefSeq fungi | F1 0.963; 13 of 14 evaluable species; 0 false positives; Bray-Curtis 0.242 (Kraken2: 0.13 to 0.17). The registered B. adolescentis test is met by 0.002 (0.495 against 0.393 + 0.10), which is marginal and very likely within noise. Abundances are worse than Kraken2's |
| Reads versus DNA mass (exploratory on the development pool, then registered Addendum A2 on the evaluation pool) | Species differ in read length (Limosilactobacillus fermentum about 2 kb, most others 7 to 9 kb), so read shares are not DNA shares. Counted in bases, Kraken2's Bray-Curtis is 0.19 (0.17 on reads); Sylph keeping one genome per species is 0.18 (0.24 with the registered sum rule, which counts E. coli reads about 2.4 times). All eight written predictions were met; both pools come from one run |
| Resistance genes (pilot, registered Addenda A3 to A6) | A reference catalogue of resistance genes was annotated in the supplier genomes with AMRFinderPlus (8 high-confidence, 30 candidate; an annotation, not biological truth). Read-level detection against ResFinder: 3 of 3 evaluable gene families detected on the evaluation pool (cepA, tet(Q), tet(W)), all five registered predictions met. Only three families were evaluable, so this is a worked example, not an accuracy estimate; the family rule had to be amended after the development result (disclosed). Candidate genes are mostly absent from the ResFinder arm, and the sequence evidence alone says nothing about the carrier species or a resistant phenotype. Species attribution (Addendum A8, Kraken2 species of each supporting read against the true organism from the supplier genomes): cepA and tet(W) reads are named to the right species (dev 34/34 and 16/18, eval 22/22 and 16/16), tet(Q) reads often only to Bacteroidales (species for 12/30 and 5/16 dev, 7/25 and 5/22 eval, for B. fragilis and P. corporis); all four units attributed, but the overall correct share was 0.684 on dev and 0.593 on eval, below the registered 0.90 and 0.85. A mock-community worked example, not an accuracy estimate. Correction (Addendum A7, base-level identity instead of minimap2's approximate estimate; all predictions met): supporting reads cepA / tet(Q) / tet(W) are 36 / 47 / 18 on dev and 22 / 51 / 16 on the evaluation pool, close to the planning values; the earlier counts (29 / 15 / 14 and 17 / 22 / 14) are kept as the registered earlier arm |

A first translation table between mock-species names and GTDB names failed on the development pool (F1 0.50 for a naming
reason). It was replaced by a genome-to-genome fastANI table in a registered addendum, disclosed as written after the
development read profiles had been seen; the first version is kept unchanged and reported beside the second.

Full numbers, predictions and verdicts are in `docs/PROJECT_STATE.md`. Small result files are in `reports/`.

## Documents

- `docs/PROJECT_STATE.md`: decisions (D1 to D25), registrations, results, open items. The single source of truth.
- `docs/HANDOVER.md`: session-by-session state and next steps.
- `docs/Nanopore_Metagenomics_Project_Plan_Rev10.docx`: project plan.
- `docs/Nanopore_Metagenomics_Workflow_Guide_v9.docx`: a beginner's guide that explains each step and term
  (thresholds, precision/recall/F1, Bray-Curtis, the Sylph experiment, the resistance-gene pilot) from scratch.
- `docs/lessons_learned.md`, `reference/README.md` (composition table and checks).

## Repository layout

| Path | Content |
|---|---|
| `bin/00` to `bin/26` | numbered scripts, in the order they were written (see below) |
| `reference/` | D6331 composition, expected reads by species, taxid aliases, GTDB crosswalks (v1, v2, v2 for evaluation) |
| `reports/` | small tracked outputs: scores, species tables, false-positive audits, tuning tables, run records |
| `metadata/` | ENA run record and the sampling manifest (seed, parameters, checksums) |
| `docs/` | project documents |
| `data/`, `results/` | **not tracked** (gitignored): reads, databases, large outputs |
| `modules/`, `workflows/`, `tests/`, `nextflow.config` | placeholders; see Status |

Scripts: `00_ena_metadata.sh` (run record), `01_sra_block_sample.py` (seeded block sampling from SRA),
`02_pool_qc.py` (read QC, pool comparison), `03_classify.py` (Kraken2 in Docker), `04_score.py` (scoring),
`05_depth_series.py`, `06_confidence_arms.py`, `07_database_arms.py` (registered experiments),
`08_trace_reads.py`, `09_composition_vs_mapping.py` (exploratory read trace), `10` to `13` (Sylph: crosswalk v1,
sketch/profile/bootstrap/score arm, fastANI crosswalk v2, gap rule), `14` to `18` (diagnostics of the Sylph abundance error: shared reads, read identity, reads versus bases; Addendum A2), `19` to `22` (resistance genes: AMRFinderPlus on the reference genomes, the reference catalogue, the read-level ResFinder arm, an exploratory check of one family-name ambiguity), `23` to `25` (exploratory, dev only: where the tet(Q) read support is lost; resolved: the registered identity filter used minimap2's approximate PAF identity, not a base-level alignment). Each script documents its arguments in its header
and with `--help`.

## Reproducing

Tools run in Docker containers with pinned image tags (Kraken2, SRA tools and Sylph are pinned; digests are recorded in
the run records). Python scripts use the standard library plus a few packages; `scipy` and similar are installed with
`pip` into the base environment. The development machine is WSL2 with 12 GB RAM and 6 CPUs; every step fits within
that except the 16 GB and full PlusPF Kraken2 indexes, which were run once on a temporary 128 GB cloud machine (about
USD 2). Reads and databases must be fetched locally:

1. Run record: `bin/00_ena_metadata.sh`. Reads: `bin/01_sra_block_sample.py` (resumable; the manifest in `metadata/`
   lists the blocks and checksums).
2. Databases are large (Kraken2 PlusPF-8 about 7.5 GB; Sylph GTDB r232 about 20 GB, fungi about 0.7 GB). Their sources,
   dates and checksums are recorded in `docs/PROJECT_STATE.md`.
3. Then follow the registered order in `docs/PROJECT_STATE.md`: classify and score the development pool, tune, tag,
   then score the evaluation pool once.

Pre-registration tags (`git tag`): `freeze-d14`, `arms-registered`, `arms-frozen`, `depth-registered`,
`dbarms-registered`, `dbarms-frozen`, `sylph-registered`, `sylph-crosswalk`, `sylph-crosswalk-v2`, `sylph-frozen`,
`sylph-gaprule`, `sylph-eval-scored`, `bases-registered`, `bases-eval-scored`, `amr-registered`, `amr-a4-registered`, `amr-map-registered`, `amr-a5-registered`, `amr-frozen`, `amr-eval-scored`, `amr-a7-registered`, `amr-a7-eval-scored`, `amr-a8-registered`, `amr-a8-eval-scored`. They prove the order in which rules, predictions and results were written.

## Status

Phase 1 (baseline benchmark) and Phase 2 arms 1 to 3 (Kraken2 confidence, database size, Sylph) are done and
documented. Phase 3 (resistance genes) is done as a pilot at the level of gene families (no alleles, point mutations or species attribution). Not done: assembly-based and Argo comparisons, Nextflow wrapping
(the first modules were deliberately written as plain scripts; `modules/`, `workflows/` and `tests/` are empty
placeholders and `nextflow.config` is a skeleton), a lot-specific composition certificate, the basecalling model tier,
and any second run or chemistry. The scripts were tested on small synthetic fixtures with answers worked out by hand;
those fixtures are not yet committed.

## Known limits

- Resistance genes: the truth is an annotation of the supplier genomes with one tool, not a measurement of the sequenced sample or of phenotype; three evaluable gene families; a read-level hit shows sequence evidence, not an acquired determinant, its host organism or a resistant strain.
- One run, one chemistry, one mock community; two pools from the same run.
- Evaluation-pool results are pre-registered comparisons on a pool that was used more than once, not fresh validation.
- The theoretical composition is not a measurement of the sequenced sample. The read trace shows the sample departs
  from it for several species (for example Roseburia hominis at about 0.4 of the expected amount); whether that is the
  lot, extraction or library preparation is unresolved.
- Reads classified as human (about 0.35 to 1% depending on the database) are counted as false positives; their origin
  is unresolved.
- Abundance is judged against the manufacturer's composition, which is a share of DNA mass, while the primary Kraken2 scores count reads. Species differ in read length in this run (about 2 kb for Limosilactobacillus fermentum, 7 to 9 kb for most), so read shares and base shares differ; the base-share view was added afterwards (Addendum A2) and does not replace the registered scores. Full-index and 16 GB Kraken2 results were not recomputed in bases.
- The Sylph comparison is between two whole set-ups (reference collection, taxonomy and algorithm), not two algorithms
  alone, and it depends on a name-translation table that had to be revised.

## Licence and data

Code and documents are released under the MIT licence (see `LICENSE`). Sequencing data, reference genomes and databases
are not included and keep their own terms: the reads are public under BioProject PRJNA1181840; the Zymo reference
genomes, Kraken2 indexes, GTDB and Sylph databases have their own licences and are downloaded by the user.

## Provenance

Developed by Krzysztof Gizynski with the help of an AI assistant (Claude, Anthropic); commits carry a co-author trailer.
Dataset facts were read from the public records (ENA, the paper and its supplement) or recomputed; anything not verified
is labelled as such in `docs/PROJECT_STATE.md`.

## Citing the data

Chen X., Yin X., Xu X., Zhang T. Species-resolved profiling of antibiotic resistance genes in complex metagenomes
through long-read overlapping with Argo. *Nature Communications* 16, 1744 (2025).
https://doi.org/10.1038/s41467-025-57088-y . Data: PRJNA1181840.
