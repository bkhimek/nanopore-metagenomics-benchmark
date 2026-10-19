#!/usr/bin/env python3
"""26_amr_attribution.py - AMR species attribution (ADDENDUM A8). Dev pool by default; the evaluation pool needs the guards below.

Question: a read supports a resistance-gene family (A7 rules, unchanged). Can the organism that carries the gene be named from the same read?
Method (nothing is reclassified; no tuning):
  1. Supporting reads = the reads that A7 assigns to a family (bin/21 rules, base-level identity; PAF <pool>_vs_resfinder_c.paf, ambiguous reads excluded).
  2. Attributed organism = the species of the read in the frozen baseline Kraken2 run (PlusPF-8, confidence 0; data/classify/<pool>/baseline_c0/kraken_output.txt.gz):
     the Kraken2 taxon if it is a species or below, else 'no species' (genus or higher, or unclassified).
  3. True organism (mock community only) = minimap2 (-x map-ont --secondary=no) of the supporting reads against the 21 supplier genomes, with the bin/08 rules
     (matching bases >= 0.5 x read length; ambiguous if the runner-up has >= 0.9 of the best). Mapping only the supporting reads gives the same label as mapping the
     whole pool (each read is mapped independently). A read without a label is excluded from the shares and counted.
  4. Destination of each read relative to its true organism, as bin/08: species (correct), up_genus / up_family / up_higher (right lineage, not to species),
     sister / same_family / elsewhere (a wrong organism), unclassified.
Scored reads: supporting reads of the evaluable gene families (cepA, tet(Q), tet(W); the A7 evaluable set); erm(B) reads and any other family are reported, not scored.
Units: species x gene of the evaluable high-confidence catalogue rows (cepA/B. fragilis, tet(Q)/B. fragilis, tet(Q)/P. corporis, tet(W)/F. prausnitzii). A unit is
'attributed' if >= 3 scored reads of its family have the Kraken2 species of the unit. A 'false organism call' is a (Kraken2 species, family) pair with >= 3 reads whose species
is not a catalogue carrier of that family. A 'cross-carrier swap' is a tet(Q) read whose true organism is one tet(Q) carrier and whose Kraken2 species is the other.
Predictions for DEV (P1-P6) and EVALUATION (E1-E5) are fixed in Addendum A8 of docs/PROJECT_STATE.md before any run and are scored by this script as written below.
Evaluation pool: needs --allow-eval, --suffix, the git tag amr-a8-registered, the committed dev outputs of the same suffix, existing Kraken2 output for the pool,
and refuses to overwrite eval_*<suffix>.* files in the output folder.
Standard library only; minimap2 in the Docker image used by bin/08.
"""
import argparse, csv, gzip, importlib.util, json, os, re, subprocess, sys, time
from collections import Counter, defaultdict
from pathlib import Path

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
P = lambda *x: os.path.join(ROOT, *x)
WRONG = ("sister", "same_family", "elsewhere")
ABOVE = ("up_genus", "up_family", "up_higher")
NEED = 3  # reads for a unit or a false organism call
# predictions (A8), as registered
P2_ALL, P2_FAM, P3_GAP, P4_WRONG, E2_ALL, E3_WRONG, E4_SWAP = 0.90, 0.85, 0.03, 0.05, 0.85, 0.10, 1


def load(name, fn):
    spec = importlib.util.spec_from_file_location(name, P("bin", fn))
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
    return m


def support_assignment(per, amb):
    """read -> (family, identity) for reads supporting exactly one family (A4/A7 rule: runner-up below amb x best)."""
    out = {}
    for q, d in per.items():
        ranked = sorted(d.items(), key=lambda kv: -kv[1][0])
        if len(ranked) > 1 and ranked[1][1][0] >= amb * ranked[0][1][0]:
            continue
        out[q] = (ranked[0][0], ranked[0][1][1])
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pool", choices=["dev", "eval"], default="dev")
    ap.add_argument("--allow-eval", action="store_true")
    ap.add_argument("--suffix", default="_A8")
    ap.add_argument("--merge", default="reference/amr_family_merge.tsv")
    ap.add_argument("--map", default="reports/amr/family_map_A5.tsv")
    ap.add_argument("--a7-paf", default="")
    ap.add_argument("--reads", default="")
    ap.add_argument("--kraken-run", default="")
    ap.add_argument("--db", default="data/db/pluspf_08_20260626")
    ap.add_argument("--ref-fa", default="data/trace/dev/ref_combined.fa")
    ap.add_argument("--truth-paf", default="", help="use an existing PAF of the supporting reads against the genomes (tests)")
    ap.add_argument("--workdir", default="")
    ap.add_argument("--out", default="reports/amr")
    ap.add_argument("--image", default="quay.io/biocontainers/minimap2:2.31--h118bc1c_0")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--memory", default="8g")
    a = ap.parse_args()
    a.a7_paf = a.a7_paf or f"data/amr/{a.pool}/{a.pool}_vs_resfinder_c.paf"
    a.reads = a.reads or f"data/subsets/SRR31222153.{a.pool}.seed20261012.blocks.fastq.gz"
    a.kraken_run = a.kraken_run or f"data/classify/{a.pool}/baseline_c0"
    a.workdir = a.workdir or f"data/amr/{a.pool}/attribution"
    if a.pool == "eval":
        if not a.allow_eval:
            sys.exit("evaluation pool: needs --allow-eval")
        if subprocess.run(["git", "-C", ROOT, "rev-parse", "-q", "--verify", "refs/tags/amr-a8-registered"], capture_output=True).returncode != 0:
            sys.exit("evaluation pool refused: git tag amr-a8-registered not found (commit and tag Addendum A8 first)")
        if not Path(P(a.out, f"dev_attribution{a.suffix}.txt")).exists():
            sys.exit(f"evaluation pool refused: the dev run ({a.out}/dev_attribution{a.suffix}.txt) is not there yet")
        if any(Path(P(a.out)).glob(f"eval_*{a.suffix}.*")):
            sys.exit(f"evaluation outputs are never overwritten (an eval_*{a.suffix}.* file exists in {a.out})")
        for need in (a.a7_paf, a.reads):
            if not Path(P(need)).exists():
                sys.exit(f"not found: {need}")
        if not Path(P(a.kraken_run, "kraken_output.txt.gz")).exists():
            sys.exit(f"not found: {a.kraken_run}/kraken_output.txt.gz (check with ls data/classify/eval)")
    x08, x21 = load("x08", "08_trace_reads.py"), load("x21", "21_amr_reads.py")
    for r in csv.DictReader(open(P(a.merge)), delimiter="\t"):
        x21.MERGE[r["alias"]] = r["family"]
    work = Path(P(a.workdir)); work.mkdir(parents=True, exist_ok=True)
    meta = {"script": "bin/26_amr_attribution.py", "args": vars(a), "started": time.strftime("%Y-%m-%d %H:%M:%S")}

    # catalogue: carriers per family, evaluable units
    fmap = list(csv.DictReader(open(P(a.map)), delimiter="\t"))
    carriers, units = defaultdict(set), []
    for r in fmap:
        if r["status"] == "mapped" and r["category"] == "high-confidence" and r["resfinder_family"]:
            carriers[r["resfinder_family"]].add(r["species"])
            if r["evaluable"] == "yes":
                units.append((r["species"], r["symbol"], r["resfinder_family"]))
    eval_fams = sorted({u[2] for u in units})

    # 1. supporting reads
    per, _ = x21.parse_paf(P(a.a7_paf))
    sup = support_assignment(per, x21.AMB)
    meta["supporting_reads"] = len(sup)

    # 2. extract their FASTQ records, 3. map to the genomes
    fq = work / f"{a.pool}_amr_support.fastq"
    if not a.truth_paf:
        n = 0
        with gzip.open(P(a.reads), "rt") as fh, open(fq, "w") as out:
            while True:
                h = fh.readline()
                if not h:
                    break
                rec = [h, fh.readline(), fh.readline(), fh.readline()]
                if h[1:].split()[0] in sup:
                    out.write("".join(rec)); n += 1
        if n != len(sup):
            sys.exit(f"only {n} of {len(sup)} supporting reads found in {a.reads}")
        paf = work / f"{a.pool}_support_vs_refs.paf"
        meta["minimap2"] = x08.run_minimap2(a.image, P(a.ref_fa), str(fq), paf, a.threads, a.memory)
    else:
        paf = Path(P(a.truth_paf))
    agg = x08.read_paf(str(paf))
    truth = {rid: (x08.assign(agg[rid], 0.5, 0.9) if rid in agg else "none") for rid in sup}

    # Kraken2 and taxonomy
    kraken = x08.read_kraken(Path(P(a.kraken_run)) / "kraken_output.txt.gz")
    missing = [r for r in sup if r not in kraken]
    if missing:
        sys.exit(f"{len(missing)} supporting reads are not in the Kraken2 output (example {missing[0]}); wrong pair of files?")
    parent, ranks = x08.load_nodes(x08.find_dmp(P(a.db), "nodes.dmp"))
    names, rev = x08.load_names(x08.find_dmp(P(a.db), "names.dmp"))
    lineage = x08.make_lineage(parent)
    infos = {}

    def info(n):
        if n not in infos:
            if n not in rev:
                sys.exit(f"species '{n}' not in names.dmp")
            infos[n] = x08.SpInfo(n, rev[n], lineage, ranks, rev)
        return infos[n]

    def species_of(k):
        for x in lineage(k):
            if ranks.get(x) == "species":
                return names.get(x, str(x))
        return ""

    def is_sp(k, sp):
        i = info(sp)
        lin = lineage(k) if k else ()
        return bool(k) and (i.t in lin or any(x in lin for x in i.alias))

    rows = []
    for rid, (fam, ident) in sorted(sup.items(), key=lambda kv: (kv[1][0], kv[0])):
        k = kraken[rid]; t = truth[rid]
        attributed = species_of(k) if k else ""
        if t in ("weak", "ambiguous", "none"):
            bucket = "no_truth"
        else:
            bucket = x08.destination(info(t), k, lineage, names, ranks)[0]
        rows.append(dict(read=rid, family=fam, identity=f"{ident:.3f}", true_species=t, kraken_taxid=k,
                         kraken_taxon=names.get(k, "unclassified") if k else "unclassified", kraken_rank=ranks.get(k, "") if k else "",
                         attributed_species=attributed, destination=bucket))
    scored = [r for r in rows if r["family"] in eval_fams]
    lab = [r for r in scored if r["destination"] != "no_truth"]

    def share(rs, buckets):
        return (sum(1 for r in rs if r["destination"] in buckets) / len(rs)) if rs else float("nan")

    # per family x true species table
    grp = defaultdict(Counter)
    for r in scored:
        grp[(r["family"], r["true_species"])][r["destination"]] += 1
    gl = [f"{'family':8s}{'true organism':32s}{'reads':>6s}{'correct':>8s}{'above sp':>9s}{'wrong':>6s}{'uncl.':>6s}{'no truth':>9s}"]
    for (fam, sp), c in sorted(grp.items()):
        tot = sum(c.values())
        gl.append(f"{fam:8s}{sp[:31]:32s}{tot:>6d}{c['species']:>8d}{sum(c[b] for b in ABOVE):>9d}{sum(c[b] for b in WRONG):>6d}{c['unclassified']:>6d}{c['no_truth']:>9d}")

    # units
    unit_rows, attributed_units = [], 0
    for sp, sym, fam in units:
        n = sum(1 for r in scored if r["family"] == fam and is_sp(int(r["kraken_taxid"]), sp))
        ntrue = sum(1 for r in scored if r["family"] == fam and r["true_species"] == sp)
        ok = n >= NEED; attributed_units += ok
        unit_rows.append(dict(species=sp, gene=sym, family=fam, reads_attributed_by_kraken2=n, reads_truly_from_species=ntrue, attributed="yes" if ok else "no"))
    # false organism calls, all supported reads of catalogue families
    pairs = Counter((r["attributed_species"], r["family"]) for r in rows
                    if r["attributed_species"] and r["family"] in carriers
                    and not any(is_sp(int(r["kraken_taxid"]), c) for c in carriers[r["family"]]))
    false_calls = [(s, f, n) for (s, f), n in pairs.items() if n >= NEED]
    # cross-carrier swaps among tet(Q)
    swaps = [r for r in scored if r["family"] == "tet(Q)" and r["true_species"] in carriers["tet(Q)"]
             and r["destination"] != "species" and any(c != r["true_species"] and is_sp(int(r["kraken_taxid"]), c) for c in carriers["tet(Q)"])]

    all_c = share(lab, ("species",)); wrong = share(lab, WRONG)
    famc = {f: share([r for r in lab if r["family"] == f], ("species",)) for f in eval_fams}
    bf = lambda f: share([r for r in lab if r["family"] == f and r["true_species"] == "Bacteroides fragilis"], ("species",))
    L = [f"AMR species attribution, {a.pool.upper()} pool (Addendum A8; A7 read support; Kraken2 {a.kraken_run}; truth = minimap2 to the supplier genomes for the supporting reads)", "",
         f"supporting reads {len(rows)}; of the evaluable families {len(scored)}; with a true-organism label {len(lab)} ({len(scored) - len(lab)} without: weak, ambiguous or no alignment)",
         "", *gl, "",
         f"{'unit (species x gene)':42s}{'Kraken2-attributed reads':>26s}{'truly from species':>20s}  attributed (>= {NEED})"]
    for u in unit_rows:
        L.append(f"{u['species'] + ' ' + u['gene']:42s}{u['reads_attributed_by_kraken2']:>26d}{u['reads_truly_from_species']:>20d}  {u['attributed']}")
    L += ["", f"correct species (labelled scored reads): {all_c:.3f}; by family: " + ", ".join(f"{f} {v:.3f}" for f, v in sorted(famc.items())),
          f"wrong organism (sister, same family, elsewhere): {wrong:.3f}; right lineage but not to species: {share(lab, ABOVE):.3f}; unclassified: {share(lab, ('unclassified',)):.3f}",
          f"cross-carrier swaps among tet(Q) reads: {len(swaps)}; false organism calls (>= {NEED} reads, species not a catalogue carrier): {len(false_calls)}"]
    for s, f, n in false_calls:
        L.append(f"  false organism call: {s} for {f} ({n} reads)")
    other = [r for r in rows if r["family"] not in eval_fams]
    if other:
        L.append("other supported families (reported, not scored): " + ", ".join(f"{f} {n}" for f, n in sorted(Counter(r["family"] for r in other).items())))
    L.append("")
    if a.pool == "dev":
        g = bf("cepA") - bf("tet(Q)")
        res = [("P1 all four units attributed", attributed_units == len(units), f"{attributed_units} of {len(units)}"),
               ("P2 correct species >= 0.90 overall and >= 0.85 per family", all_c >= P2_ALL and all(v >= P2_FAM for v in famc.values()), f"{all_c:.3f}; " + ", ".join(f"{f} {v:.3f}" for f, v in sorted(famc.items()))),
               (f"P3 B. fragilis tet(Q) correct share at least {P3_GAP} below cepA (mobile-element hypothesis)", g >= P3_GAP, f"cepA {bf('cepA'):.3f}, tet(Q) {bf('tet(Q)'):.3f}"),
               ("P4 wrong organism <= 0.05", wrong <= P4_WRONG, f"{wrong:.3f}"),
               ("P5 no cross-carrier swap among tet(Q) reads", len(swaps) == 0, str(len(swaps))),
               ("P6 no false organism call", len(false_calls) == 0, str(len(false_calls)))]
    else:
        res = [("E1 all four units attributed", attributed_units == len(units), f"{attributed_units} of {len(units)}"),
               ("E2 correct species >= 0.85 overall", all_c >= E2_ALL, f"{all_c:.3f}"),
               ("E3 wrong organism <= 0.10", wrong <= E3_WRONG, f"{wrong:.3f}"),
               (f"E4 cross-carrier swaps <= {E4_SWAP}", len(swaps) <= E4_SWAP, str(len(swaps))),
               ("E5 no false organism call", len(false_calls) == 0, str(len(false_calls)))]
    L.append("Registered predictions (A8):")
    for name, ok, val in res:
        L.append(f"  {name}: {'MET' if ok else 'NOT MET'}  ({val})")
    L += ["", "Mock community only: the true organism comes from the supplier genomes, which a real sample does not have. The Kraken2 database contains close relatives of every carrier, "
          "yet a mock of 21 genomes is easier than a real sample. Attribution is by the whole read, not by the gene; a gene on a mobile element can sit in sequence shared by several species. "
          "Research use only, not clinically validated."]
    txt = "\n".join(L)
    os.makedirs(P(a.out), exist_ok=True)
    open(P(a.out, f"{a.pool}_attribution{a.suffix}.txt"), "w").write(txt + "\n")
    with open(P(a.out, f"{a.pool}_attribution_reads{a.suffix}.tsv"), "w") as o:
        w = csv.DictWriter(o, fieldnames=list(rows[0]), delimiter="\t"); w.writeheader(); w.writerows(rows)
    with open(P(a.out, f"{a.pool}_attribution_units{a.suffix}.tsv"), "w") as o:
        w = csv.DictWriter(o, fieldnames=list(unit_rows[0]), delimiter="\t"); w.writeheader(); w.writerows(unit_rows)
    json.dump(meta, open(P(a.out, f"{a.pool}_attribution{a.suffix}.json"), "w"), indent=1)
    print(txt)


if __name__ == "__main__":
    main()
