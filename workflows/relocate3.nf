/*
 * RELOCATE3: per-sample TE insertion calling.
 *
 * Skeleton only: the genome preparation, per-sample steps and lane merging are
 * wired in later (see plans/2026-11-02-nextflow-m1-plan.md, Tasks 4-5).
 */

workflow RELOCATE3 {
    take:
    ch_samples // channel: [ meta, fastq_1, fastq_2 ]

    main:
    ch_samples.count().view { n -> "RELOCATE3: ${n} sample sheet row(s) validated" }
}
