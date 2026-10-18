#!/usr/bin/env python3
"""13_crosswalk_extend.py - gap rule of Addendum A1.2 (PROJECT_STATE.md), fixed BEFORE any eval profile exists.

Rule: any genome in a Sylph profile whose accession is NOT among the 42 candidate genomes used to build crosswalk v2
(reports/sylph/crosswalk_v2_candidates.tsv) is tested with the SAME fastANI rule against the 21 reference genomes
(ANI >= 95, aligned fraction >= 0.5). Reads never decide. The extension can only ADD (mock species, database species) pairs; it never removes
or changes a v2 row, except that a 'not_represented' row is replaced when a new genome qualifies for that mock species. A database species that
would map to two mock species aborts the script. If the profile contains no new accession the output is an exact copy of the base crosswalk.

Outputs: the extended crosswalk (default reference/d6331_gtdb_crosswalk_v2_eval.csv, same columns as v2), reports/sylph/crosswalk_v2_gap_<tag>.tsv
(every new accession with its species and best fastANI per mock species) and a JSON record. Run it on the dev profile first (expected: no new
accessions) and then on the eval profile.
"""
import argparse, csv, glob, json, os, shutil, sys, time
import importlib.util

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
P = lambda *x: os.path.join(ROOT, *x)
_s = importlib.util.spec_from_file_location("x12", P("bin", "12_crosswalk_fastani.py"))
x12 = importlib.util.module_from_spec(_s); _s.loader.exec_module(x12)
cw = x12.cw


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--profile", nargs="+", required=True, help="Sylph profile TSV(s) to check, e.g. data/sylph/eval/arm_3a.tsv")
    ap.add_argument("--tag", required=True, help="label for the report files, e.g. dev or eval")
    ap.add_argument("--base", default="reference/d6331_gtdb_crosswalk_v2.csv")
    ap.add_argument("--base-candidates", default="reports/sylph/crosswalk_v2_candidates.tsv")
    ap.add_argument("--truth", default="reference/d6331_expected_reads_by_species.csv")
    ap.add_argument("--out", default="reference/d6331_gtdb_crosswalk_v2_eval.csv")
    ap.add_argument("--refdir", default="data/reference/D6331.refseq/genomes")
    ap.add_argument("--gtdb-taxonomy", nargs="+", default=["data/db/gtdb_tax/bac120_taxonomy_r232.tsv.gz", "data/db/gtdb_tax/ar53_taxonomy_r232.tsv.gz"])
    ap.add_argument("--workdir", default="data/crosswalk_v2")
    ap.add_argument("--min-ani", type=float, default=95.0)
    ap.add_argument("--min-af", type=float, default=0.5)
    ap.add_argument("--threads", type=int, default=6)
    ap.add_argument("--fastani-tsv", help="existing fastANI output for the new genomes (tests only)")
    a = ap.parse_args()

    gtdb = cw.load_gtdb([P(x) for x in a.gtdb_taxonomy])
    base_acc = {l.split("\t")[0] for l in list(open(P(a.base_candidates)))[1:]}
    prof = x12.candidates(a.profile, gtdb)
    new = {k: v for k, v in prof.items() if k not in base_acc}
    with open(P(a.base)) as fh:
        rd = csv.DictReader(fh)
        cols, rows = rd.fieldnames, list(rd)
    tids = {int(r["mock_taxid"]): r["mock_species"] for r in rows}
    rec = {"script": "bin/13_crosswalk_extend.py", "tag": a.tag, "profiles": a.profile, "n_base_candidates": len(base_acc),
           "n_profile_accessions": len(prof), "n_new": len(new), "time": time.strftime("%Y-%m-%d %H:%M:%S")}
    os.makedirs(P("reports", "sylph"), exist_ok=True)
    added = []
    pairs = []
    if new:
        work = P(a.workdir); os.makedirs(work, exist_ok=True)
        if a.fastani_tsv:
            pf = P(a.fastani_tsv)
        else:
            gdir = os.path.join(work, "genomes")
            missing = x12.download(sorted(new), gdir)
            rec["not_downloadable"] = missing
            if missing:
                print(f"WARNING: could not download {missing}", file=sys.stderr)
            pf = os.path.join(work, f"fastani_gap_{a.tag}.tsv")
            refs = sorted(glob.glob(P(a.refdir, "*.fasta")))
            rec["fastani_version"] = x12.run_fastani(refs, [os.path.join(gdir, k + ".fna") for k in sorted(new) if os.path.exists(os.path.join(gdir, k + ".fna"))], pf, a.threads)
        names = {}
        for t in csv.DictReader(open(P(a.truth))):
            nm = (t["ncbi_name_in_kraken_db"].strip(), int(t["species_taxid"]))
            for k in (t["species"], t["ncbi_name_in_kraken_db"]):
                names[k.strip().lower()] = nm
        have = {(int(r["mock_taxid"]), r["database_species"]) for r in rows if r["status"] == "ok"}
        owner = {r["database_species"]: int(r["mock_taxid"]) for r in rows if r["status"] == "ok"}
        best = {}
        for line in open(pf):
            f = line.rstrip("\n").split("\t")
            if len(f) < 5:
                continue
            q = os.path.basename(f[0])[:-len(".fasta")]
            m = x12.ACC.search(os.path.basename(f[1]))
            if not m or m.group(0) not in new:
                sys.exit(f"fastANI reference {f[1]} is not a new accession")
            acc, ani, af = m.group(0), float(f[2]), int(f[3]) / int(f[4])
            mock = next((names[c.lower()] for c in (cw.stem_species(q), q.replace("_", " ")) if c.lower() in names), None)
            if mock is None:
                sys.exit(f"reference {q} not found in the base crosswalk")
            sp = new[acc]["species"] or "unresolved:" + acc
            pairs.append((acc, sp, mock[0], ani, af))
            if ani >= a.min_ani and af >= a.min_af:
                d = best.setdefault((mock[1], sp), {"ani": 0.0, "af": 0.0, "acc": acc, "files": set()})
                d["files"].add(q)
                if ani > d["ani"]:
                    d.update(ani=ani, af=af, acc=acc)
        for (tid, sp), d in sorted(best.items()):
            if owner.setdefault(sp, tid) != tid:
                sys.exit(f"gap rule: database species {sp} would map to two mock species; decide explicitly")
            if (tid, sp) in have:
                continue
            rows = [r for r in rows if not (int(r["mock_taxid"]) == tid and r["status"] == "not_represented")]
            rows.append({"mock_species": tids[tid], "mock_taxid": str(tid), "database_species": sp, "status": "ok", "fastani": f"{d['ani']:.2f}",
                         "aligned_fraction": f"{d['af']:.2f}", "genome_accession": d["acc"], "reference_files": ";".join(sorted(d["files"])) + " [gap rule]"})
            added.append((tids[tid], sp, d["ani"], d["af"]))
        rows.sort(key=lambda r: (r["mock_species"], -float(r["fastani"] or 0)))
    with open(P(a.out), "w", newline="") as o:
        w = csv.DictWriter(o, fieldnames=cols); w.writeheader(); w.writerows(rows)
    with open(P("reports", "sylph", f"crosswalk_v2_gap_{a.tag}.tsv"), "w") as o:
        o.write("accession\tdatabase_species\tbest_mock_species\tbest_fastani\taligned_fraction\taccepted\n")
        for acc, sp in sorted((k, v["species"]) for k, v in new.items()):
            ps = sorted((p for p in pairs if p[0] == acc), key=lambda p: -p[3])
            b = ps[0] if ps else (acc, sp, "", float("nan"), float("nan"))
            ok = any(sp == x[1] for x in added) or any(p[0] == acc and p[3] >= a.min_ani and p[4] >= a.min_af for p in ps)
            o.write(f"{acc}\t{sp}\t{b[2]}\t{b[3]:.2f}\t{b[4]:.2f}\t{'yes' if ok else 'no'}\n")
    rec["added"] = [list(x) for x in added]
    json.dump(rec, open(P("reports", "sylph", f"crosswalk_v2_gap_{a.tag}.json"), "w"), indent=2)
    print(f"gap rule ({a.tag}): {len(prof)} accessions in the profile(s), {len(new)} not among the {len(base_acc)} base candidates; {len(added)} mapping(s) added")
    for x in added:
        print(f"  added: {x[0]} <- {x[1]} (ANI {x[2]:.2f}, AF {x[3]:.2f})")
    print(f"wrote {a.out}")


if __name__ == "__main__":
    main()
