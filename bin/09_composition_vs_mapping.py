#!/usr/bin/env python3
"""09_composition_vs_mapping.py - EXPLORATORY, dev pool only.  Compare each database arm's species counts with TWO references:
  (a) the theoretical D6331 DNA composition (the registered benchmark's truth), and
  (b) the composition measured by mapping the same reads to the reference genomes (reports/trace/trace_by_species.tsv, bin/08).
(b) removes any difference between the real tube and the manufacturer's theoretical composition, so what is left is classifier error.

Method (fixed before running): for each arm, take the species table at the frozen fraction 0.0008 (status 'detected', dna_pct > 0),
observed = observed_reads, normalised over the called mock species; reference (a) = dna_pct, (b) = mapped_reads, each normalised over the SAME
species.  Bray-Curtis = 0.5 * sum |o - e|.  The registered figure in score_summary.json also counts the human false positive as a called
species with expected 0, so it is slightly larger (about +0.003); it is printed for comparison and must be within 0.01 of (a), otherwise the
script stops (the convention would then be wrong).  Also prints observed / mapped per species and arm.
Eval pool is not used (no mapping exists for it).  Standard library only."""
import argparse, json, re, sys
from pathlib import Path

ARMS = [("PlusPF-8", "db_pluspf8_f0.0008"), ("PlusPF-16", "db_pluspf16_f0.0008"), ("PlusPF full", "db_pluspf_f0.0008")]


def read_tsv(path):
    with open(path) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        return head, [dict(zip(head, line.rstrip("\n").split("\t"))) for line in fh if line.strip()]


def bc(o, e):
    return 0.5 * sum(abs(o[k] - e[k]) for k in o)


def norm(d):
    s = sum(d.values())
    return {k: v / s for k, v in d.items()}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--scoredir", default="data/aws_run/unpacked/data/score/dev")
    ap.add_argument("--trace", default="reports/trace/trace_by_species.tsv")
    ap.add_argument("--out", default="reports/trace/bc_vs_mapping.tsv")
    a = ap.parse_args()

    _, trows = read_tsv(a.trace)
    mapped = {int(r["taxid"]): int(r["mapped_reads"]) for r in trows}
    names = {int(r["taxid"]): r["species"] for r in trows}

    summary, per_sp = [], {}
    for arm, d in ARMS:
        p = Path(a.scoredir) / d
        head, rows = read_tsv(p / "species_table.tsv")
        scol = next(h for h in head if h.startswith("status_"))
        called = [r for r in rows if r[scol] == "detected" and float(r["dna_pct"]) > 0]
        obs = {int(r["taxid"]): float(r["observed_reads"]) for r in called}
        theo = {int(r["taxid"]): float(r["dna_pct"]) for r in called}
        missing = [t for t in obs if t not in mapped]
        if missing:
            sys.exit(f"ERROR: taxids {missing} are called in {d} but absent from the trace table")
        mp = {t: float(mapped[t]) for t in obs}
        o, e_th, e_mp = norm(obs), norm(theo), norm(mp)
        b_th, b_mp = bc(o, e_th), bc(o, e_mp)
        reported = float(re.search(r'"bray_curtis":\s*([0-9.]+)', (p / "score_summary.json").read_text()).group(1))
        if abs(b_th - reported) > 0.01:
            sys.exit(f"ERROR: {arm}: recomputed BC vs theoretical {b_th:.4f} differs from score_summary.json {reported:.4f} by > 0.01; convention is wrong")
        summary.append((arm, len(called), reported, b_th, b_mp))
        for t in obs:
            per_sp.setdefault(t, {})[arm] = obs[t] / mapped[t]

    lines = ["Bray-Curtis (half L1) on the dev pool, called mock species at fraction 0.0008", "",
             f"{'arm':12s} {'species':>7s} {'registered':>10s} {'vs theory':>10s} {'vs mapping':>10s}"]
    for arm, n, rep, bt, bm in summary:
        lines.append(f"{arm:12s} {n:7d} {rep:10.3f} {bt:10.3f} {bm:10.3f}")
    lines += ["", "observed Kraken2 reads / mapped reads (1.00 = Kraken2 finds as many reads as map to the genome)",
              f"{'species':32s}" + "".join(f"{arm:>13s}" for arm, _ in ARMS)]
    for t in sorted(per_sp, key=lambda x: -mapped[x]):
        lines.append(f"{names[t][:32]:32s}" + "".join(f"{per_sp[t].get(arm, float('nan')):13.2f}" for arm, _ in ARMS))
    text = "\n".join(lines)
    print(text)
    Path(a.out).write_text(text + "\n")


if __name__ == "__main__":
    main()
