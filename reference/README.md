# D6331 reference composition

Source (checked 2026-10-13): ZymoBIOMICS Gut Microbiome Standard (D6331) instruction manual, Table 1
(theoretical composition) and Table 2 (strain information). PDF SHA-256:
a6fdb1b45760e81db4e86c4367e055c8018794d86ac1aaf7c0ed0ef640cc76eb

This is the THEORETICAL composition. Zymo also publishes a lot-specific certificate of analysis
(https://www.zymoresearch.com/pages/certificate-of-analysis, needs the lot number). The lot of the tube sequenced for
SRR31222153 is not known, so the lot CoA has not been obtained. Report results against the theoretical composition and say so.

Files
- d6331_composition.csv: all 21 strains, with the five manual columns (genomic DNA, 16S only, 16S and 18S, genome copy, cell number; blank = N/A in the manual), genome size, GC, rRNA copy number.
- d6331_expected_reads_by_species.csv: 17 species (E. coli strains summed). Expected reads = depth (Gb) x genomic-DNA fraction / 5,479 bp (whole-run mean read length).
  Assumes reads follow DNA mass with no extraction or sequencing bias. These are expectations, not measurements.

Checks done
- Each column sums to 100% within rounding (genomic DNA 100.011).
- Genome-copy % equals genomic-DNA % divided by genome size, renormalised (checked for all bacteria with copy % > 0.5).
- Reference genome archive (D6331.refseq.zip, SHA-256 4795b7df33ecd244986cd938778de03e1dd1df4b84c7349837f4cd3d37797d3d): 20 genome FASTAs; bacterial and archaeal lengths match Table 2 and GC matches within 0.1 percentage point.
  The two yeast files are draft assemblies: C. albicans is 13,062 contigs and 12.58 Mb (Table 2 says 14.68 Mb); S. cerevisiae is 51 contigs, 12.84 Mb (Table 2: 13.30 Mb).
  Keep the archive in data/reference/ (not committed).

Which column to use
- Read-fraction outputs (for example Kraken2 read counts): compare with genomic DNA %.
- Coverage- or copy-based abundance (for example Sylph taxonomic abundance): compare with genome copy %.

Taxonomy IDs (added 2026-10-15)
- species_taxid and ncbi_name_in_kraken_db come from the PlusPF-8 (2026-06-26) index listing (inspect.txt, rank S rows). Score by taxid, not by name:
  Zymo's "Lactobacillus fermentum" is "Limosilactobacillus fermentum" (taxid 1613) in the current NCBI taxonomy. All five E. coli strains map to species 562.

Taxonomy aliases (added 2026-10-16, D13)
- d6331_taxid_aliases.csv: extra taxids whose species-level reads count as a mock species. One entry so far: Fusobacterium animalis (76859) counts as the D6331 "Fusobacterium nucleatum" strain (851).
  The evidence is in the file (the reference genome classifies 97.1% as F. animalis with the frozen baseline). It was added AFTER the first dev scoring; the composition and expected-reads files are unchanged.
