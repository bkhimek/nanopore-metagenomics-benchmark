#!/usr/bin/env bash
# Phase 0/1: stream a large ENA FASTQ once, never storing it, and keep two disjoint,
# seeded read pools (development and evaluation). Verifies the stream's MD5 against ENA.
#
# Usage:  bin/01_stream_subset.sh [RUN_ACCESSION] [SEED]
# Needs:  metadata/<PROJECT>_read_run.tsv from bin/00_ena_metadata.sh (columns fastq_ftp, fastq_md5, read_count)
# Env overrides: SOURCE (URL or local .fastq.gz, for testing), DEV_FRAC, EVAL_FRAC, MIN_FREE_GB
set -euo pipefail
RUN="${1:-SRR31222153}"; SEED="${2:-20261012}"
PROJECT="${PROJECT:-PRJNA1181840}"
HERE="$(cd "$(dirname "$0")/.." && pwd)"
DATA="${DATA:-$HERE/data}"
OUT="$DATA/subsets"; mkdir -p "$OUT"
META="$HERE/metadata/${PROJECT}_read_run.tsv"
DEV_FRAC="${DEV_FRAC:-0.0156}"; EVAL_FRAC="${EVAL_FRAC:-0.0156}"   # ~2 Gb each for this run
MIN_FREE_GB="${MIN_FREE_GB:-15}"

free_gb=$(df -Pk "$OUT" | awk 'NR==2{print int($4/1024/1024)}')
[ "$free_gb" -ge "$MIN_FREE_GB" ] || { echo "Only ${free_gb} GB free in $OUT (< $MIN_FREE_GB)"; exit 1; }

if [ -z "${SOURCE:-}" ]; then
  [ -f "$META" ] || { echo "Missing $META: run bin/00_ena_metadata.sh $PROJECT first"; exit 1; }
  read -r FTP MD5 < <(awk -F'\t' -v r="$RUN" 'NR==1{for(i=1;i<=NF;i++)c[$i]=i; next} $c["run_accession"]==r{print $c["fastq_ftp"], $c["fastq_md5"]}' "$META")
  [ -n "${FTP:-}" ] || { echo "Run $RUN not found in $META"; exit 1; }
  SOURCE="https://${FTP}"           # ENA lists host/path without a scheme
else
  MD5="${EXPECTED_MD5:-}"
fi
echo "Source: $SOURCE  expected md5: ${MD5:-none}  seed: $SEED  dev/eval frac: $DEV_FRAC/$EVAL_FRAC"

command -v pigz >/dev/null && GZ="pigz -dc" || GZ="gzip -dc"
MD5FILE="$(mktemp)"
fetch() { case "$SOURCE" in http*|ftp*) curl -fsS --retry 5 --retry-delay 10 "$SOURCE";; *) cat "$SOURCE";; esac; }

fetch | tee >(md5sum | cut -d' ' -f1 > "$MD5FILE") | $GZ | \
  python3 "$HERE/bin/01_stream_split_sample.py" --seed "$SEED" --dev-frac "$DEV_FRAC" --eval-frac "$EVAL_FRAC" \
    --outdir "$OUT" --prefix "$RUN" > "$OUT/${RUN}.seed${SEED}.manifest.json"

sleep 1; GOT="$(cat "$MD5FILE")"; rm -f "$MD5FILE"
if [ -n "${MD5:-}" ] && [ "$GOT" != "$MD5" ]; then
  echo "MD5 MISMATCH: got $GOT expected $MD5. Outputs are NOT trustworthy." | tee "$OUT/${RUN}.seed${SEED}.FAILED"; exit 2
fi
echo "Stream md5 OK (${GOT}). Manifest: $OUT/${RUN}.seed${SEED}.manifest.json"
( cd "$OUT" && sha256sum ${RUN}.*.seed${SEED}.fastq.gz > ${RUN}.seed${SEED}.sha256 )
