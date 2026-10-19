// One row per arm: the headline numbers of every score_summary.json (threshold, calls, precision, recall, F1, Bray-Curtis).
process SUMMARISE {
    container params.python_image
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path summaries

    output:
    path "arms_summary.tsv"

    script:
    """
    python3 - <<'PY'
    import glob, json
    cols = ["pool", "run", "T", "called", "correct", "fp", "precision", "recall_evaluable", "f1", "bray_curtis", "reads_in_mock_species", "pct_reads_in_mock_species"]
    rows = []
    for p in sorted(glob.glob("*.summary.json")):
        d = json.load(open(p)); a = d.get("at_T", {})
        rows.append([d.get("pool", ""), d.get("run", p), a.get("T", ""), a.get("called", ""), a.get("correct", ""), a.get("fp", ""), a.get("precision", ""),
                     a.get("recall_evaluable", ""), a.get("f1", ""), a.get("bray_curtis", ""), d.get("reads_in_mock_species", ""), d.get("pct_reads_in_mock_species", "")])
    with open("arms_summary.tsv", "w") as o:
        o.write("\\t".join(cols) + "\\n")
        for r in rows:
            o.write("\\t".join(str(x) for x in r) + "\\n")
    PY
    """

    stub:
    """
    printf 'pool\\trun\\n' > arms_summary.tsv
    """
}
