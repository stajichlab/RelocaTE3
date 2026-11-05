#!/usr/bin/env nextflow
/*
 * stajichlab/RelocaTE3: transposable element insertion detection.
 *
 * Usage: nextflow run . -profile test,pixi
 *        nextflow run . --help
 *        nextflow run . --help genome
 *
 * The workflow holds no detection logic: each process runs one relocaTE3
 * subcommand with the arguments `relocaTE3 run-all` passes, so a sample's
 * files match run-all's.
 */

include { paramsHelp; validateParameters; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'
include { RELOCATE3; groupLanes } from './workflows/relocate3'
include { softwareEnvironment } from './subworkflows/local/software_environment'

workflow {
    main:
    // nf-schema 3: help only through paramsHelp(), before validation so that
    // --help needs no other parameter. `--help` alone arrives as true or
    // 'true'; `--help <name>` names one parameter (nf-schema docs pattern).
    if (params.help) {
        log.info(paramsHelp(
            parameter: params.help.toString() in ['true', 'false'] ? null : params.help,
            command: 'nextflow run stajichlab/RelocaTE3 -profile <pixi|docker|apptainer|singularity> --input samples.csv --genome reference.fa --te_library te.fa --outdir results',
        ))
        exit(0)
    }
    validateParameters()
    log.info(paramsSummaryLog())

    // -profile pixi puts ${params.pixi_env}/bin first on every task's PATH
    // (nextflow.config). Stop here, not in the first task, when relocaTE3 is
    // not there. Container profiles never read pixi_env.
    if (workflow.profile.tokenize(',').contains('pixi')) {
        def exe = file("${params.pixi_env}/bin/relocaTE3")
        if (!exe.exists()) {
            error("-profile pixi: no relocaTE3 at ${exe} (--pixi_env ${params.pixi_env}). Run `pixi install -e pipeline` in the pipeline directory, or point --pixi_env at an environment that has it.")
        }
    }

    // What ran the tasks: pipeline_info/software_environment.txt.
    channel.of(softwareEnvironment())
        .collectFile(name: 'software_environment.txt', storeDir: "${params.outdir}/pipeline_info")

    // The sheet's rows are validated here, against assets/schema_input.json,
    // not by validateParameters(): the `input` param has no "schema" link,
    // because nf-schema 3.0.0 applies validation.allowParamsSubstitution
    // (on in the test profile: ${params.test_data}) only in samplesheetToList.
    // Still before any task runs.
    // Sheet rows [ meta (id), fastq_1, fastq_2, genotype_bam ] grouped into
    // one item per sample, in sheet order:
    // [ meta, [ R1 lanes ], [ R2 lanes ], genotype_bam or [] ].
    ch_samples = channel.fromList(
        groupLanes(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
    )

    RELOCATE3(ch_samples)
}
