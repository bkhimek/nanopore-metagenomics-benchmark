# Lessons learned

- Verify, don't transcribe: dataset facts from AI-generated reviews were treated as unverified until
  confirmed in ENA (PRJNA1181840 confirmed 2026-10-12). Emu was wrongly suggested for shotgun data
  (it targets 16S/ITS amplicons) and corrected.

## Corrections and failed predictions (full list, 2026-10-19)

Each item is in the registered record (docs/PROJECT_STATE.md) with its original result kept beside the replacement.

| What happened | Effect | Handling |
|---|---|---|
| First GTDB name table (Sylph) | F1 0.50 on the development pool, for a naming reason | Replaced by a genome-to-genome fastANI table, disclosed as written after the development profiles had been seen; the first table is kept and reported |
| Approximate read identity in the AMR filter (minimap2 PAF estimate) | Supporting reads changed from 29 / 15 / 14 to 36 / 47 / 18 (development) and from 17 / 22 / 14 to 22 / 51 / 16 (evaluation) for cepA / tet(Q) / tet(W) | Base-level identity added in a registered addendum (A7); the earlier counts are kept as the earlier arm |
| Family rule for AMR genes | Needed amending after the development result | Amendment disclosed |
| Positional k-mer measure (A9) | Insensitive (93 to 96% of windows have no hit); all four development predictions missed | Evaluation run not done (disclosed); replaced by a registered, better measure (A10) |
| Species attribution target (A8) | Correct share 0.593 on the evaluation pool against 0.85 registered (0.684 against 0.90 on development) | Reported as missed; the mechanism was then investigated |
| minimap2 default `-p 0.8` in the genome-context check (bin/30, run 1) | Short alignments hidden beneath a long one; two hypotheses void | Rerun with `-p 0 -N 200`; both runs reported. A negative result now needs a positive control |
| Early reading of the Veillonella read (bin/31) | The 12 kb half looked only 88% identical to P. corporis (whole-read mapping), which was used against a chimera; the registered prediction (insertion) was not met | Inference withdrawn when the half, mapped alone, was 98.9% identical (bin/32); the junction check (152 spanning reads, none split) supports a chimera; the cause of the difference between the two mappings is unresolved |
| Registered Sylph test on B. adolescentis | Met by 0.002, very likely within noise | Reported as met, marginal |

Lesson: a prediction can fail because the list of explanations was too coarse (the Veillonella read was an element-to-chromosome junction, a case the first three explanations did not name).
