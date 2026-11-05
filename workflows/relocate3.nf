/*
 * RELOCATE3: per-sample TE insertion calling, as `relocaTE3 run-all` per sample.
 *
 *   samples (grouped by groupLanes in main.nf, before any task)
 *                     -> CAT_FASTQ (only for a sample with several lanes)
 *   genome            -> PREPARE_GENOME (once)
 *   per sample        -> ONE_SAMPLE
 *   all characterized -> RELOCATE3_COHORT (once; only if any sample was characterized)
 *   summaries         -> MULTIQC (once, at the end; unless skip_multiqc)
 *   versions.yml      -> ${params.outdir}/pipeline_info/software_versions.yml
 *
 * Published files (conf/modules.config): each sample's are run-all's, under
 * ${params.outdir}/<sample>/; the cohort's under ${params.outdir}/cohort/;
 * the report under ${params.outdir}/multiqc/.
 */

include { CAT_FASTQ      } from '../modules/local/cat_fastq/main'
include { PREPARE_GENOME } from '../subworkflows/local/prepare_genome'
include { ONE_SAMPLE     } from '../subworkflows/local/one_sample'
include { RELOCATE3_COHORT } from '../modules/local/relocate3_cohort/main'
include { MULTIQC        } from '../modules/local/multiqc/main'

workflow RELOCATE3 {
    take:
    ch_samples // channel: [ meta (id), [ R1 lanes ], [ R2 lanes ], genotype_bam or [] ], one per sample (groupLanes)

    main:
    ch_versions = channel.empty()

    // The genotype BAM travels apart from meta and the reads: only
    // CHARACTERIZE reads it, so adding, moving or touching it must not change
    // the inputs (and -resume hashes) of the sample's other tasks.
    // [ meta, BAM/CRAM, index ] for samples that bring their own alignment
    ch_genotype_bams = ch_samples
        .filter { _meta, _reads_1, _reads_2, bam -> bam }
        .map { meta, _reads_1, _reads_2, bam -> [ meta, bam, genotypeBamIndex(meta.id, bam) ] }

    ch_lanes = ch_samples
        .map { meta, reads_1, reads_2, _bam -> [ meta, reads_1, reads_2 ] }
        .branch { _meta, reads_1, _reads_2 ->
            merge: reads_1.size() > 1
            single: true
        }

    CAT_FASTQ(ch_lanes.merge)
    ch_versions = ch_versions.mix(CAT_FASTQ.out.versions)

    // [ meta, R1, R2 ], R2 is [] for single-end
    ch_reads = ch_lanes.single
        .map { meta, reads_1, reads_2 -> [ meta, reads_1[0], reads_2 ? reads_2[0] : [] ] }
        .mix(CAT_FASTQ.out.reads)

    PREPARE_GENOME(file(params.genome, checkIfExists: true))
    ch_versions = ch_versions.mix(PREPARE_GENOME.out.versions)

    ONE_SAMPLE(
        ch_reads,
        PREPARE_GENOME.out.genome_with_index,
        file(params.te_library, checkIfExists: true),
        params.repeatmasker ? file(params.repeatmasker, checkIfExists: true) : [],
        ch_genotype_bams
    )
    ch_versions = ch_versions.mix(ONE_SAMPLE.out.versions)

    // Cohort tables, once. toSortedList waits for every sample and orders the
    // tables by sample id, so the task's inputs (-resume) do not depend on
    // the order samples finish. It emits [] when no sample was characterized
    // (--genotype false, no genotype_bam): then no COHORT.
    RELOCATE3_COHORT(
        ONE_SAMPLE.out.characterized
            .map { meta, table -> [ meta.id, table ] }
            .toSortedList { a, b -> a[0] <=> b[0] }
            .filter { pairs -> pairs.size() > 0 }
            .map { pairs -> [ pairs.collect { p -> p[0] }, pairs.collect { p -> p[1] } ] }
    )
    ch_versions = ch_versions.mix(RELOCATE3_COHORT.out.versions)

    // One entry per process: identical versions.yml of several samples
    // collapse. Plus the workflow's own entry (as nf-core pipelines do).
    ch_versions_text = ch_versions
        .map { f -> f.text }
        .unique()
        .mix(channel.of(workflowVersions()))

    ch_report_versions = channel.empty()
    if (!params.skip_multiqc) {
        // MultiQC's "Software versions" section reads *_mqc_versions.yml:
        // { group: { tool: version } }. Groups are the simple process names.
        ch_mqc_versions = ch_versions_text
            .map { text -> text.replaceAll(/(?m)^"(?:[^":]*:)*([^":]+)":/, '$1:') }
            .collectFile(name: 'relocate3_software_mqc_versions.yml', sort: true)

        // Summary files only; sorted by name (unique: cohort, <id>.stats,
        // versions) so the task's inputs do not depend on completion order.
        MULTIQC(
            RELOCATE3_COHORT.out.summary
                .mix(ONE_SAMPLE.out.stats.map { _meta, stats -> stats })
                .mix(ch_mqc_versions)
                .toSortedList { a, b -> a.name <=> b.name },
            file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true),
            reportComment()
        )
        ch_report_versions = MULTIQC.out.versions.map { f -> f.text }
    }

    ch_versions_text
        .mix(ch_report_versions)
        .collectFile(name: 'software_versions.yml', storeDir: "${params.outdir}/pipeline_info", sort: true)
}

/*
 * The workflow's entry in software_versions.yml, as nf-core pipelines write it.
 */
def workflowVersions() {
    return "Workflow:\n    Nextflow: ${workflow.nextflow.version}\n    ${workflow.manifest.name}: ${workflow.manifest.version}\n"
}

/*
 * The MultiQC report comment: pipeline version and what ran the tasks
 * (pipeline_info/software_environment.txt has the full record).
 */
def reportComment() {
    def revision = workflow.commitId ? " (${workflow.revision ?: 'commit'} ${workflow.commitId.take(7)})" : ''
    def software = workflow.containerEngine
        ? "${workflow.containerEngine} image ${workflow.container}"
        : (workflow.profile.tokenize(',').contains('pixi') ? "pixi environment ${params.pixi_env}" : 'tools on the host PATH')
    return "${workflow.manifest.name} ${workflow.manifest.version}${revision}; software: ${software}. Run ${workflow.runName}, ${workflow.start}."
}

/*
 * Group sample-sheet rows ([ meta, fastq_1, fastq_2, genotype_bam ] from samplesheetToList)
 * by sample, keeping sheet order (groupBy keeps the first-seen order of keys
 * and the order of rows within a key). The rows of a sample must agree on
 * single- vs paired-end and on genotype_bam; a given genotype_bam needs its
 * index next to it. Called on the list in main.nf, so a bad sheet fails with
 * its message before any task starts (an error raised inside a channel
 * operator is reported only as InvocationTargetException).
 *
 * Returns [ [ meta, [ R1 lanes ], [ R2 lanes ], genotype_bam ], ... ]
 *   meta is [ id: sample ] only; R2 lanes [] for single-end; genotype_bam [] if none.
 */
def groupLanes(rows) {
    return rows
        .groupBy { row -> row[0].id }
        .collect { id, lanes ->
            def paired = lanes.collect { row -> row[2] ? 'paired-end' : 'single-end' }.unique()
            if (paired.size() > 1) {
                error("Sample '${id}': its sample-sheet rows mix paired-end and single-end reads (${paired.join(', ')}); all lanes of a sample must be the same")
            }
            def bams = lanes.collect { row -> row[3] ? row[3].toString() : '' }.unique()
            if (bams.size() > 1) {
                error("Sample '${id}': its sample-sheet rows give different genotype_bam values (${bams.collect { b -> b ?: '<empty>' }.join(', ')}); all lanes of a sample must give the same one")
            }
            def bam = lanes[0][3] ?: []
            if (bam) {
                genotypeBamIndex(id, bam)
            }
            [ [ id: id ], lanes.collect { row -> row[1] }, lanes.findAll { row -> row[2] }.collect { row -> row[2] }, bam ]
        }
}

/*
 * The index of a sample's genotype_bam, which must sit next to it:
 * <file>.bai or <file>.csi for BAM, <file>.crai for CRAM.
 */
def genotypeBamIndex(id, bam) {
    def suffixes = bam.name.endsWith('.cram') ? [ '.crai' ] : [ '.bai', '.csi' ]
    def candidates = suffixes.collect { s -> file("${bam}${s}") }
    def index = candidates.find { f -> f.exists() }
    if (!index) {
        error("Sample '${id}': genotype_bam ${bam} has no index; expected ${candidates.join(' or ')} (create it with `samtools index`)")
    }
    return index
}
