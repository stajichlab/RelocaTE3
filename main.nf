#!/usr/bin/env nextflow
/*
 * stajichlab/RelocaTE3: transposable element insertion detection.
 *
 * Usage: nextflow run . -profile test,pixi
 *        nextflow run . --help
 *
 * The workflow holds no detection logic: each process runs one relocaTE3
 * subcommand with the arguments `relocaTE3 run-all` passes, so a sample's
 * files match run-all's.
 */

include { validateParameters; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'
include { RELOCATE3; groupLanes } from './workflows/relocate3'

workflow {
    main:
    validateParameters()
    log.info(paramsSummaryLog(workflow))

    // Sheet rows [ meta (id, genotype_bam), fastq_1, fastq_2 ] grouped into
    // one item per sample: [ meta, [ R1 lanes ], [ R2 lanes ] ], sheet order.
    ch_samples = channel.fromList(
        groupLanes(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
    )

    RELOCATE3(ch_samples)
}
