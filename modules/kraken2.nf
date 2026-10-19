// Kraken2 classification of one read pool at one --confidence value.
// Same command as bin/03_classify.py (image digest, --minimum-hit-groups 2, --report-minimizer-data, gzip of the per-read output),
// and the same run.json fields that bin/04_score.py reads (pool, label, database, input_sha256_in_manifest).
process KRAKEN2 {
    tag "${meta.pool}:${meta.label}"
    container params.kraken2_image
    cpus params.threads
    memory params.memory
    publishDir "${params.outdir}/classify/${meta.pool}/${meta.label}", mode: 'copy'

    input:
    tuple val(meta), path(reads)
    path db

    output:
    tuple val(meta), path("report.txt"), path("kraken_output.txt.gz"), path("run.json"), emit: run

    script:
    def rec = [pool: meta.pool, label: meta.label, input: meta.input, input_sha256_in_manifest: meta.input_sha256,
               image: params.kraken2_image, database: params.db_name, database_archive_md5: params.db_archive_md5,
               options: [confidence: meta.confidence as double, minimum_hit_groups: 2, threads: task.cpus, report_minimizer_data: true],
               produced_by: "nextflow modules/kraken2.nf"]
    def json = groovy.json.JsonOutput.toJson(rec)   // one line, no single quotes: safe inside printf '%s'
    """
    kraken2 --db ${db} --threads ${task.cpus} --gzip-compressed \\
        --confidence ${meta.confidence} --minimum-hit-groups 2 --report-minimizer-data \\
        --report report.txt --output kraken_output.txt ${reads}
    gzip -6 kraken_output.txt
    printf '%s\\n' '${json}' > run.json
    """

    stub:
    """
    printf '' > report.txt
    printf '' | gzip > kraken_output.txt.gz
    echo '{}' > run.json
    """
}
