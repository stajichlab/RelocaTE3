/*
 * RELOCATE3: per-sample TE insertion calling, as `relocaTE3 run-all` per sample.
 *
 *   samples (grouped by groupLanes in main.nf, before any task)
 *                     -> CAT_FASTQ (only for a sample with several lanes)
 *   genome            -> PREPARE_GENOME (once)
 *   per sample        -> ONE_SAMPLE
 *   versions.yml      -> ${params.outdir}/pipeline_info/software_versions.yml
 *
 * Published files (conf/modules.config) are run-all's, under
 * ${params.outdir}/<sample>/.
 */

include { CAT_FASTQ      } from '../modules/local/cat_fastq/main'
include { PREPARE_GENOME } from '../subworkflows/local/prepare_genome'
include { ONE_SAMPLE     } from '../subworkflows/local/one_sample'

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

    // One entry per process: identical versions.yml of several samples collapse.
    ch_versions
        .map { f -> f.text }
        .unique()
        .collectFile(name: 'software_versions.yml', storeDir: "${params.outdir}/pipeline_info", sort: true)
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
