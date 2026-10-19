#!/usr/bin/env python3
"""Build the synthetic fixture for the Nextflow tests (tests/fixtures/). Deterministic; answers are worked out by hand below.
Mock: three truth species (50, 30, 20 % of the DNA) and one false positive (a fourth genus), 1,000 reads, 5 % unclassified.
  species  taxid  reads   truth share   expected reads
  Alpha one  1001   480     50 %          500
  Beta two   1002   280     30 %          300
  Gamma three 1003  150     20 %          200
  Delta four 1004    40     (not in truth: false positive)
  unclassified      50
At threshold fraction 0.01 (T = 10 reads): 4 species called, 3 correct, 1 false positive -> precision 0.75, recall (evaluable) 1.0, F1 0.857.
"""
import gzip, json, os
H = os.path.join(os.path.dirname(os.path.abspath(__file__)), "fixtures")
os.makedirs(os.path.join(H, "db"), exist_ok=True); cl = os.path.join(H, "classify", "dev", "baseline_c0"); os.makedirs(cl, exist_ok=True)
sp = [(1001, "Alpha one", 10, "Alpha", 480, 50.0), (1002, "Beta two", 11, "Beta", 280, 30.0), (1003, "Gamma three", 12, "Gamma", 150, 20.0), (1004, "Delta four", 13, "Delta", 40, 0.0)]
open(os.path.join(H, "truth.csv"), "w").write("species,species_taxid,ncbi_name_in_kraken_db,domain,n_strains,genomic_dna_pct\n" + "".join(
    f"{n},{t},{n},Bacteria,1,{d}\n" for t, n, g, gn, r, d in sp if d > 0))
open(os.path.join(H, "aliases.csv"), "w").write("species_taxid,accepted_taxid,accepted_name,evidence\n")
nodes = ["1\t|\t1\t|\tno rank\t|\n"] + [f"{g}\t|\t1\t|\tgenus\t|\n" for t, n, g, gn, r, d in sp] + [f"{t}\t|\t{g}\t|\tspecies\t|\n" for t, n, g, gn, r, d in sp]
open(os.path.join(H, "db", "nodes.dmp"), "w").write("".join(x.replace("\t|\n", "\t|\t\n") for x in nodes))
rep = ["  5.00\t50\t50\t0\t0\tU\t0\tunclassified", " 95.00\t950\t0\t0\t0\tR\t1\troot"]
for t, n, g, gn, r, d in sp:
    rep.append(f"{r / 10:6.2f}\t{r}\t0\t0\t0\tG\t{g}\t  {gn}")
    rep.append(f"{r / 10:6.2f}\t{r}\t{r}\t{r}\t{r}\tS\t{t}\t    {n}")
open(os.path.join(cl, "report.txt"), "w").write("\n".join(rep) + "\n")
with gzip.GzipFile(os.path.join(cl, "kraken_output.txt.gz"), "wb", mtime=0) as _g, __import__("io").TextIOWrapper(_g) as o:   # mtime=0: reproducible bytes
    i = 0
    for t, n, g, gn, r, d in sp:
        for _ in range(r):
            i += 1; o.write(f"C\tfix.{i}\t{t}\t1000\t1000:1\n")
    for _ in range(50):
        i += 1; o.write(f"U\tfix.{i}\t0\t1000\t0:1\n")
json.dump({"pool": "dev", "label": "baseline_c0", "database": "data/db/fixture_db", "input": "fixture_dev.fastq.gz", "input_sha256_in_manifest": "fixture",
           "image": "none (fixture)", "options": {"confidence": 0.0}}, open(os.path.join(cl, "run.json"), "w"), indent=2)
json.dump({"pools": {"dev": {"file": "fixture_dev.fastq.gz", "sha256": "fixture"}}, "blocks": [{"pool": "dev", "start": 1}, {"pool": "dev", "start": 501}]},
          open(os.path.join(H, "fixture.sra_blocks.manifest.json"), "w"), indent=2)
print("fixture written to", H)
