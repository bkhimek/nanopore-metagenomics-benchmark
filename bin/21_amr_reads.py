#!/usr/bin/env python3
"""21_amr_reads.py - AMR read-level pilot (ADDENDUM A4; DEV pool only; the evaluation pool has no code path here).

Step 1  --map-only : family map (reports/amr/family_map.tsv). Every AMR-type catalogue row (bin/20; POINT rows are not scored) is matched to a
                     ResFinder gene family by case-insensitive exact symbol; no match = 'not in arm database' (excluded from recall, not a miss).
                     Also computes the A4 evaluable value E = carrier_coverage x (L - 0.6 g) / L (L = carrier mean read length from bin/16, g = gene
                     length from the reference coordinates); E >= 10 is 'evaluable'. No reads are used. Commit this file BEFORE step 2.
Step 2  (default)  : minimap2 -x map-ont --secondary=yes -N 100 of the dev reads against ResFinder all.fsa (per-class files duplicate it and are not used),
                     then the A4 rules: an alignment counts when matches/block >= 0.80 and aligned reference span >= 0.80 of the gene; per read and
                     family the best allele alignment (most matching bases) is kept; a read is assigned to the family with the most matching bases and is
                     'ambiguous' if the runner-up has >= 0.9 of it (not counted as support). Family = allele name before '_<digits>_<accession>'.
                     A family is 'detected' at >= 3 supporting reads; the curve at 1, 2, 3, 5, 10 reads is reported.
Outputs: reports/amr/dev_family_reads.tsv (all families with >= 1 supporting read), dev_outcomes.tsv (one row per catalogue species x family),
         dev_amr_report.txt, dev_run.json (command, image digest, ResFinder version file). Scoring unit: species x family (E. coli strains collapsed).
A5 (2026-10-19): --merge FILE (alias -> family; reference/amr_family_merge.tsv) renames ResFinder families that sequence evidence (bin/22) shows to be
the same gene under different names; --suffix _A5 writes the outputs next to the A4-rule outputs (which are kept). Without --merge the behaviour is
exactly A4 as registered.
Standard library only; minimap2 in the Docker image used by bin/08.
"""
import argparse, csv, importlib.util, json, os, re, sys
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
P = lambda *x: str(ROOT.joinpath(*x))
HDR = re.compile(r"^(.+?)_(\d+)_(.+)$")
MERGE = {}  # alias -> family (A5); empty = A4 rule
MIN_ID, MIN_SPAN, AMB, CALL, E_MIN = 0.80, 0.80, 0.9, 3, 10.0
CURVE = (1, 2, 3, 5, 10)


def family(header):
    m = HDR.match(header.lstrip(">").split()[0])
    f = m.group(1) if m else header.lstrip(">").split()[0]
    return MERGE.get(f, f)


def families_in(fsa):
    fam = Counter()
    with open(fsa) as fh:
        for line in fh:
            if line.startswith(">"):
                fam[family(line)] += 1
    return fam


def catalogue(path):
    """AMR-type rows, collapsed to species x symbol (E. coli strains merged); POINT rows listed separately."""
    keep, point = {}, []
    for r in csv.DictReader(open(path), delimiter="\t"):
        if r["subtype"] == "POINT":
            point.append(r); continue
        k = (r["species"], r["symbol"])
        g = abs(int(r["stop"]) - int(r["start"])) + 1
        o = keep.get(k)
        if o is None:
            keep[k] = dict(r, gene_len=g, strains=1, cats={r["category"]})
        else:
            o["strains"] += 1; o["cats"].add(r["category"]); o["gene_len"] = max(o["gene_len"], g)
    for o in keep.values():  # a species x symbol is high-confidence only if every strain's row is
        o["category"] = "high-confidence" if o["cats"] == {"high-confidence"} else "candidate"
    return keep, point


def build_map(a):
    fam = families_in(P(a.resfinder, "all.fsa"))
    low = {f.lower(): f for f in fam}
    cov = {r["species"]: (float(r["cov"]), float(r["len"])) for r in csv.DictReader(open(P(a.coverage)), delimiter="\t")}
    cat, point = catalogue(P(a.catalogue))
    rows = []
    for (sp, sym), o in sorted(cat.items()):
        f = low.get(sym.lower())
        c, L = cov.get(sp, (None, None))
        E = c * (L - 0.6 * o["gene_len"]) / L if c is not None else None
        rows.append({"species": sp, "symbol": sym, "category": o["category"], "class": o["class"], "strains": o["strains"], "gene_len": o["gene_len"],
                     "resfinder_family": f or "", "resfinder_alleles": fam[f] if f else 0,
                     "status": "mapped" if f else "not in arm database",
                     "carrier_cov": "" if c is None else f"{c:.2f}", "mean_read_len": "" if L is None else f"{L:.0f}",
                     "E": "" if E is None else f"{E:.1f}",
                     "evaluable": "no (carrier not in dev trace)" if E is None else ("yes" if E >= E_MIN else "low opportunity (E < 10)")})
    os.makedirs(P(a.out), exist_ok=True)
    with open(P(a.out, f"family_map{a.suffix}.tsv"), "w") as o:
        w = csv.DictWriter(o, fieldnames=list(rows[0]), delimiter="\t"); w.writeheader(); w.writerows(rows)
    L = [f"Family map (A4): {len(rows)} species x symbol units (AMR type; {len(point)} POINT rows not scored); ResFinder families {len(fam)}", ""]
    for r in rows:
        L.append(f"{r['species'][:28]:29s}{r['symbol'][:14]:15s}{r['category'][:9]:10s}{r['status'][:19]:20s}{r['resfinder_family'][:14]:15s}E {r['E'] or '-':>6s}  {r['evaluable']}")
    L.append("")
    L.append(f"high-confidence mapped: {sum(r['category']=='high-confidence' and r['status']=='mapped' for r in rows)} of {sum(r['category']=='high-confidence' for r in rows)}; "
             f"high-confidence mapped AND evaluable: {sum(r['category']=='high-confidence' and r['status']=='mapped' and r['evaluable']=='yes' for r in rows)}")
    txt = "\n".join(L); open(P(a.out, f"family_map_report{a.suffix}.txt"), "w").write(txt + "\n"); print(txt)
    return rows


def parse_paf(path):
    """read -> family -> [best matching bases, identity of that alignment]; applies the A4 alignment filter."""
    per, nal = {}, 0
    with open(path) as fh:
        for line in fh:
            f = line.split("\t", 12)
            q, qlen, tname, tlen, ts, te, mt, bl = f[0], int(f[1]), f[5], int(f[6]), int(f[7]), int(f[8]), int(f[9]), int(f[10])
            nal += 1
            if bl <= 0 or mt / bl < MIN_ID or (te - ts) / tlen < MIN_SPAN:
                continue
            d = per.setdefault(q, {}); fam = family(tname)
            if fam not in d or mt > d[fam][0]:
                d[fam] = [mt, mt / bl]
    return per, nal


def run(a, fam_rows):
    sys.path.insert(0, P("bin"))
    spec = importlib.util.spec_from_file_location("x08", P("bin", "08_trace_reads.py"))
    x08 = importlib.util.module_from_spec(spec); spec.loader.exec_module(x08)
    work = Path(P(a.workdir)); work.mkdir(parents=True, exist_ok=True)
    paf = work / "dev_vs_resfinder.paf"; meta = {"script": "bin/21_amr_reads.py", "args": vars(a)}
    if a.paf:
        paf = Path(a.paf)
    else:
        import subprocess, time
        image = a.image
        if subprocess.run(["docker", "image", "inspect", image], capture_output=True).returncode != 0:
            subprocess.run(["docker", "pull", image], check=True)
        dig = subprocess.run(["docker", "image", "inspect", image, "--format", "{{json .RepoDigests}}"], capture_output=True, text=True).stdout.strip()
        fa, rd = Path(P(a.resfinder, "all.fsa")).resolve(), Path(P(a.reads)).resolve()
        cmd = ["docker", "run", "--rm", "-u", f"{os.getuid()}:{os.getgid()}", "-m", a.memory, "--cpus", str(a.threads),
               "-v", f"{fa.parent}:/ref:ro", "-v", f"{rd.parent}:/reads:ro", "-v", f"{paf.resolve().parent}:/out",
               image, "minimap2", "-x", "map-ont", "--secondary=yes", "-N", "100", "-t", str(a.threads),
               "-o", f"/out/{paf.name}", f"/ref/{fa.name}", f"/reads/{rd.name}"]
        t0 = time.time(); subprocess.run(cmd, check=True)
        meta.update({"image": image, "repo_digests": dig, "command": " ".join(cmd), "wall_s": round(time.time() - t0, 1)})
    vf = Path(P(a.out, "resfinder_version.txt"))
    meta["resfinder_version"] = vf.read_text() if vf.exists() else "MISSING"
    per, nal = parse_paf(str(paf))
    support, ambig, idsum = Counter(), Counter(), defaultdict(float)
    for q, d in per.items():
        ranked = sorted(d.items(), key=lambda kv: -kv[1][0])
        best = ranked[0]
        if len(ranked) > 1 and ranked[1][1][0] >= AMB * best[1][0]:
            for f, _ in ranked:
                if _[0] >= AMB * best[1][0]:
                    ambig[f] += 1
            continue
        support[best[0]] += 1; idsum[best[0]] += best[1][1]
    os.makedirs(P(a.out), exist_ok=True)
    with open(P(a.out, f"dev_family_reads{a.suffix}.tsv"), "w") as o:
        o.write("family\tsupporting_reads\tambiguous_reads\tmean_identity\n")
        for f in sorted(set(support) | set(ambig), key=lambda f: -support[f]):
            o.write(f"{f}\t{support[f]}\t{ambig[f]}\t{idsum[f]/support[f]:.3f}\n" if support[f] else f"{f}\t0\t{ambig[f]}\t\n")
    outc, mapped_fams = [], set()
    for r in fam_rows:
        f = r["resfinder_family"]; mapped_fams.add(f)
        s = support[f] if f else 0
        if r["status"] != "mapped":
            res = "not in arm database"
        elif r["evaluable"] != "yes":
            res = ("detected (not scored: " if s >= CALL else "not detected (not scored: ") + r["evaluable"] + ")"
        else:
            res = "detected" if s >= CALL else "missed"
        outc.append(dict(r, support=s, ambiguous=ambig[f] if f else 0, result=res, **{f"ge{k}": int(s >= k) for k in CURVE}))
    cols = list(outc[0])
    with open(P(a.out, f"dev_outcomes{a.suffix}.tsv"), "w") as o:
        w = csv.DictWriter(o, fieldnames=cols, delimiter="\t"); w.writeheader(); w.writerows(outc)
    hc = [r for r in outc if r["category"] == "high-confidence"]
    ev = [r for r in hc if r["status"] == "mapped" and r["evaluable"] == "yes"]
    L = [f"AMR read-level pilot, DEV pool (A4 rules; {nal} alignments; {len(per)} reads pass the alignment filter; ResFinder {a.resfinder})", "",
         f"{'species':26s}{'gene':13s}{'cat':10s}{'E':>6s}{'reads':>6s}{'amb':>5s}  result"]
    for r in sorted(outc, key=lambda r: (r["category"] != "high-confidence", r["species"], r["symbol"])):
        L.append(f"{r['species'][:25]:26s}{r['symbol'][:12]:13s}{r['category'][:9]:10s}{r['E'] or '-':>6s}{r['support']:>6d}{r['ambiguous']:>5d}  {r['result']}")
    # Without read-level species attribution a family has ONE support count shared by every carrier species (tet(Q) in B. fragilis and
    # P. corporis), so recall is counted over distinct families; the species x family rows are listed above for transparency.
    evf = {}
    for r in ev:
        evf[r["resfinder_family"]] = support[r["resfinder_family"]]
    L += ["", f"Recall, evaluable high-confidence DISTINCT families ({len(ev)} species x gene rows; {len(evf)} families, since support cannot be attributed to a species): "
          f"{sum(v >= CALL for v in evf.values())} of {len(evf)} detected at >= {CALL} reads"]
    for k in CURVE:
        L.append(f"  at >= {k} reads: {sum(v >= k for v in evf.values())} of {len(evf)}")
    unexp = [(f, support[f]) for f in support if f not in mapped_fams and support[f] >= CALL]
    L += ["", f"Families with >= {CALL} supporting reads that are NOT in the reference catalogue (not automatically false positives): {len(unexp)}"]
    for f, s in sorted(unexp, key=lambda x: -x[1])[:40]:
        L.append(f"  {f:28s}{s:>6d} reads (ambiguous {ambig[f]})")
    L += ["", "Gene-family level only; no alleles, point mutations or species attribution. Reference annotation is not biological truth. Research use only."]
    txt = "\n".join(L); open(P(a.out, f"dev_amr_report{a.suffix}.txt"), "w").write(txt + "\n"); print(txt)
    json.dump(meta, open(P(a.out, f"dev_run{a.suffix}.json"), "w"), indent=1)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pool", choices=["dev"], default="dev", help="dev only; the evaluation pool needs a freeze addendum and new code")
    ap.add_argument("--map-only", action="store_true")
    ap.add_argument("--resfinder", default="data/db/resfinder_db")
    ap.add_argument("--catalogue", default="reports/amr/reference_catalogue.tsv")
    ap.add_argument("--coverage", default="reports/sylph/dev_diagnostic/base_basis_by_species.tsv")
    ap.add_argument("--reads", default="data/subsets/SRR31222153.dev.seed20261012.blocks.fastq.gz")
    ap.add_argument("--workdir", default="data/amr/dev")
    ap.add_argument("--out", default="reports/amr")
    ap.add_argument("--image", default="quay.io/biocontainers/minimap2:2.31--h118bc1c_0")
    ap.add_argument("--threads", type=int, default=6)
    ap.add_argument("--memory", default="8g")
    ap.add_argument("--paf", help="use an existing PAF (skips mapping; for tests)")
    ap.add_argument("--merge", default="", help="A5: TSV alias<TAB>family<TAB>evidence; empty = A4 rule")
    ap.add_argument("--suffix", default="", help="suffix for output names, e.g. _A5")
    a = ap.parse_args()
    if a.merge:
        for r in csv.DictReader(open(P(a.merge)), delimiter="\t"):
            MERGE[r["alias"]] = r["family"]
    rows = build_map(a)
    if not a.map_only:
        run(a, rows)


if __name__ == "__main__":
    main()
