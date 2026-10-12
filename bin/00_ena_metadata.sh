#!/usr/bin/env bash
# Phase 0: fetch ENA run-level metadata for a project accession into metadata/.
# Usage: bin/00_ena_metadata.sh PRJNA1181840
set -euo pipefail
ACC="${1:?usage: $0 <ENA/SRA project accession>}"
OUT="$(dirname "$0")/../metadata"; mkdir -p "$OUT"
F="run_accession,experiment_accession,sample_accession,instrument_model,instrument_platform,library_strategy,library_source,library_selection,library_name,study_title,experiment_title,sample_title,read_count,base_count,fastq_bytes,fastq_md5,fastq_ftp,first_public"
curl -sS "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${ACC}&result=read_run&fields=${F}&format=tsv" \
  -o "$OUT/${ACC}_read_run.tsv"
echo "wrote $OUT/${ACC}_read_run.tsv"; df -h "$HOME" | tail -1
