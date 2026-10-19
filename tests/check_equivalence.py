#!/usr/bin/env python3
"""Equivalence check: does the Nextflow workflow reproduce the frozen dev baseline?
Compares the workflow's species table and summary with the committed ones (reports/baseline_kraken2_pluspf8/dev_*), cell by cell.
Run after:  nextflow run . -profile docker      (PlusPF-8, confidence 0.0, fraction 0.0008 on the dev pool)
  python3 tests/check_equivalence.py
Exit 0 only if every species row and every headline number is identical. Any difference is printed; a difference means the wrapper, not the science, is wrong.
"""
import argparse, csv, json, sys
ap = argparse.ArgumentParser()
ap.add_argument("--workflow-dir", default="results/score/dev/baseline_c0_f0.0008")
ap.add_argument("--committed-table", default="reports/baseline_kraken2_pluspf8/dev_species_table.tsv")
ap.add_argument("--committed-summary", default="reports/baseline_kraken2_pluspf8/dev_score_summary.json")
a = ap.parse_args()
bad = 0
rd = lambda p: list(csv.reader(open(p), delimiter="\t"))
w, c = rd(f"{a.workflow_dir}/species_table.tsv"), rd(a.committed_table)
if w[0] != c[0]:
    print("header differs:", w[0], c[0]); bad += 1
if len(w) != len(c):
    print(f"row count differs: workflow {len(w)}, committed {len(c)}"); bad += 1
for i, (x, y) in enumerate(zip(w[1:], c[1:]), 1):
    for j, (u, v) in enumerate(zip(x, y)):
        if u != v:
            print(f"row {i} ({y[1]}), column {c[0][j]}: workflow {u} != committed {v}"); bad += 1
sw, sc = json.load(open(f"{a.workflow_dir}/score_summary.json")), json.load(open(a.committed_summary))
for k in ("total_reads", "classified", "unclassified", "reads_in_mock_species", "pct_reads_in_mock_species"):
    if sw.get(k) != sc.get(k):
        print(f"{k}: workflow {sw.get(k)} != committed {sc.get(k)}"); bad += 1
for k, v in sc["at_T"].items():
    if sw["at_T"].get(k) != v:
        print(f"at_T {k}: workflow {sw['at_T'].get(k)} != committed {v}"); bad += 1
print("IDENTICAL: the workflow reproduces the frozen dev baseline" if not bad else f"{bad} difference(s)")
sys.exit(1 if bad else 0)
