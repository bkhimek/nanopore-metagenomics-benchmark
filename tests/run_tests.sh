#!/usr/bin/env bash
# Tests for the Nextflow core path. Needs nextflow and python3 (standard library) on PATH; no Docker, no Kraken2 database.
#   1. the fixture is reproducible (tests/make_fixture.py gives the files that are committed)
#   2. nextflow run . -profile test scores the fixture: T=10, 4 called, 3 correct, 1 false positive, precision 0.75, recall 1.0, F1 0.857 (worked out by hand in make_fixture.py)
#   3. the evaluation pool is refused without --allow_eval
set -euo pipefail
cd "$(dirname "$0")/.."
python3 tests/make_fixture.py >/dev/null
git diff --exit-code --quiet -- tests/fixtures || { echo "FAIL: fixture differs from the committed one"; exit 1; }
rm -rf results_test
nextflow run . -profile test -ansi-log false >/dev/null
python3 - <<'PY'
import csv
r = list(csv.DictReader(open("results_test/arms_summary.tsv"), delimiter="\t"))
assert len(r) == 1, r
r = r[0]
exp = {"T": "10", "called": "4", "correct": "3", "fp": "1", "precision": "0.75", "recall_evaluable": "1.0"}
for k, v in exp.items():
    assert r[k] == v, (k, r[k], v)
assert abs(float(r["f1"]) - 6 / 7) < 1e-9, r["f1"]
print("PASS: fixture scored as worked out by hand")
PY
out=$(nextflow run . -profile test --pool eval -ansi-log false 2>&1 || true)
if grep -q "refusing to run the evaluation pool" <<<"$out"; then echo "PASS: evaluation pool refused without --allow_eval"; else echo "FAIL: eval guard"; exit 1; fi
