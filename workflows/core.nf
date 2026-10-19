// Core benchmark path as a workflow: classify one pool with Kraken2 at one or more --confidence values, score each run against the D6331 truth
// at a fixed threshold fraction, and collect the arms in one table. Pool sampling (bin/01) is NOT wrapped: it is a one-off, resumable SRA download
// whose seed and manifest are its record; the workflow takes the finished pool and its manifest as input.
include { KRAKEN2   } from '../modules/kraken2'
include { SCORE     } from '../modules/score'
include { SUMMARISE } from '../modules/summarise'

workflow CORE {
    if (params.pool == 'eval' && !params.allow_eval) {
        error "refusing to run the evaluation pool without --allow_eval (thresholds must be frozen on the dev pool first; bin/04 and bin/06 apply the same rule)"
    }
    manifest = file(params.manifest, checkIfExists: true)
    man = new groovy.json.JsonSlurper().parse(manifest)
    pinfo = man.pools[params.pool]
    sha = pinfo.sha256

    // one task per confidence value; label as in bin/06 (baseline_c0, c0.05, ...)
    reads = params.classified_dir ? null : (params.reads ? file(params.reads, checkIfExists: true) : file("${params.datadir}/subsets/${pinfo.file}", checkIfExists: true))
    arms = params.classified_dir ? Channel.empty() : Channel.of(params.arms.toString().tokenize(',')*.trim()).flatten().map { c ->
        def conf = c as double
        def label = conf == 0 ? 'baseline_c0' : 'c' + String.format('%s', (conf == (long) conf) ? (long) conf : conf)
        [[pool: params.pool, label: label, confidence: c, input: reads.name, input_sha256: sha], reads]
    }

    if (params.classified_dir) {
        // skip classification: use an existing run folder (report.txt, kraken_output.txt.gz, run.json); used by the tests and for rescoring
        d = file(params.classified_dir, checkIfExists: true)
        rj = new groovy.json.JsonSlurper().parse(d.resolve('run.json'))
        runs = Channel.of([[pool: rj.pool, label: rj.label], d.resolve('report.txt'), d.resolve('kraken_output.txt.gz'), d.resolve('run.json')])
    } else {
        KRAKEN2(arms, file(params.db, checkIfExists: true))
        runs = KRAKEN2.out.run
    }

    SCORE(runs, file(params.truth, checkIfExists: true), file(params.aliases, checkIfExists: true), manifest,
          file("${params.db}/nodes.dmp", checkIfExists: true), file("${projectDir}/bin/04_score.py", checkIfExists: true))
    SUMMARISE(SCORE.out.summary.collect())
}
