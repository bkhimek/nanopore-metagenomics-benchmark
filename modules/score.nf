// Scoring of one Kraken2 run with bin/04_score.py, unchanged.
// bin/04_score.py finds everything relative to the project root, so each task builds a small project-shaped
// folder (bin/, reference/, data/classify, data/subsets, data/db) from staged inputs and runs the script there.
// The evaluation pool is scored only with --allow_eval, at one fixed fraction, with no sweep (same rule as bin/04).
process SCORE {
    tag "${meta.pool}:${meta.label}"
    container params.python_image
    publishDir "${params.outdir}/score/${meta.pool}", mode: 'copy', pattern: 'score/*', saveAs: { it.replaceFirst('^score/', '') }

    input:
    tuple val(meta), path(report), path(kout), path(runjson)
    path truth
    path aliases
    path manifest
    path nodes
    path score_script

    output:
    tuple val(meta), path("score/*"), emit: scored
    path "${meta.pool}_${meta.label}.summary.json", emit: summary

    script:
    def allow = meta.pool == 'eval' ? '--allow-eval' : ''
    """
    mkdir -p bin reference data/subsets data/classify/${meta.pool}/${meta.label} score
    ln -s \$PWD/${score_script} bin/04_score.py
    ln -s \$PWD/${truth} reference/${truth}
    ln -s \$PWD/${aliases} reference/${aliases}
    ln -s \$PWD/${manifest} data/subsets/${manifest}
    for f in ${report} ${kout} ${runjson}; do ln -s \$PWD/\$f data/classify/${meta.pool}/${meta.label}/\$f; done
    DB=\$(python3 -c "import json; print(json.load(open('${runjson}'))['database'])")
    mkdir -p "\$DB"; ln -s \$PWD/${nodes} "\$DB/nodes.dmp"
    python3 bin/04_score.py --run data/classify/${meta.pool}/${meta.label} \\
        --truth reference/${truth} --aliases reference/${aliases} \\
        --min-frac ${params.min_frac} --bootstrap ${params.bootstrap} ${allow} > score_stdout.txt
    cp -r data/score/${meta.pool}/* score/
    cp score_stdout.txt score/\$(ls score | head -1)/score_stdout.txt
    cp score/*/score_summary.json ${meta.pool}_${meta.label}.summary.json
    """

    stub:
    """
    mkdir -p score/${meta.label}_fSTUB
    echo '{"at_T": {}, "bootstrap": {}}' > score/${meta.label}_fSTUB/score_summary.json
    cp score/${meta.label}_fSTUB/score_summary.json ${meta.pool}_${meta.label}.summary.json
    """
}
