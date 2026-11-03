/*
 * ONE_SAMPLE: the per-sample steps of `relocaTE3 run-all`, one task each.
 *
 *   TE_SEARCH -> ALIGN_GENOME -> FIND_INSERTIONS -> CHARACTERIZE   (genotype, or genotype_bam given)
 *                             \-> FIND_REFERENCE                   (repeatmasker)
 *   ALIGN_READS (genotype, sample without genotype_bam), in parallel with TE_SEARCH
 *
 * Channels are joined on meta ([ id: sample ]), the same map at every step of
 * one sample. The genotype BAM is not in meta: it reaches CHARACTERIZE only,
 * so it does not change any other task's inputs (-resume).
 *
 * Take: reads              [ meta, R1, R2 ]   R2 is [] for single-end
 *       genome_with_index  [ genome, [ index files ] ]   (value channel)
 *       te_library         TE library FASTA
 *       repeatmasker       RepeatMasker .out, or [] (no FIND_REFERENCE)
 *       genotype_bams      [ meta, BAM/CRAM, index ] for samples with a
 *                          genotype_bam (index checked by the caller)
 * Emit: versions           versions.yml of every task
 */

include { RELOCATE3_TE_SEARCH       } from '../../modules/local/relocate3_te_search/main'
include { RELOCATE3_ALIGN_GENOME    } from '../../modules/local/relocate3_align_genome/main'
include { RELOCATE3_ALIGN_READS     } from '../../modules/local/relocate3_align_reads/main'
include { RELOCATE3_FIND_INSERTIONS } from '../../modules/local/relocate3_find_insertions/main'
include { RELOCATE3_FIND_REFERENCE  } from '../../modules/local/relocate3_find_reference/main'
include { RELOCATE3_CHARACTERIZE    } from '../../modules/local/relocate3_characterize/main'

workflow ONE_SAMPLE {
    take:
    reads
    genome_with_index
    te_library
    repeatmasker
    genotype_bams

    main:
    ch_versions = channel.empty()

    // run-all steps 2-3: TE-library search and trimming
    RELOCATE3_TE_SEARCH(reads, te_library)
    ch_versions = ch_versions.mix(RELOCATE3_TE_SEARCH.out.versions)

    // run-all step 4: [ meta, R1, R2, [ te_containing ], [ flanking ] ]
    RELOCATE3_ALIGN_GENOME(
        reads
            .join(RELOCATE3_TE_SEARCH.out.te_containing)
            .join(RELOCATE3_TE_SEARCH.out.flanking),
        genome_with_index
    )
    ch_versions = ch_versions.mix(RELOCATE3_ALIGN_GENOME.out.versions)

    // run-all step 5. The fullreads BAM is optional: [] when align-genome
    // wrote none, as run-all passes --fullreads-bam only if the file exists.
    ch_fullreads = RELOCATE3_ALIGN_GENOME.out.fullreads
        .map { meta, bam, bai -> [ meta, [ bam, bai ] ] }
    ch_find_insertions = RELOCATE3_ALIGN_GENOME.out.bam
        .join(RELOCATE3_TE_SEARCH.out.read_repeat)
        .join(ch_fullreads, remainder: true)
        .map { meta, bam, bai, read_repeat, fullreads -> [ meta, bam, bai, read_repeat, fullreads ?: [] ] }
    RELOCATE3_FIND_INSERTIONS(ch_find_insertions, repeatmasker)
    ch_versions = ch_versions.mix(RELOCATE3_FIND_INSERTIONS.out.versions)

    // run-all steps 0/6: reference TE calls, only with --repeatmasker
    if (params.repeatmasker) {
        RELOCATE3_FIND_REFERENCE(
            RELOCATE3_ALIGN_GENOME.out.bam.join(RELOCATE3_TE_SEARCH.out.read_repeat),
            repeatmasker
        )
        ch_versions = ch_versions.mix(RELOCATE3_FIND_REFERENCE.out.versions)
    }

    // run-all step 7: genotyping (`if args.genotype or args.genotype_bam`).
    // A sample with a genotype_bam is always characterized on it, as run-all
    // does; params.genotype only decides whether ALIGN_READS aligns the reads
    // of the samples without one.
    ch_genotype_bams = genotype_bams
    if (params.genotype) {
        RELOCATE3_ALIGN_READS(
            reads
                .join(genotype_bams.map { meta, _bam, _index -> [ meta, true ] }, remainder: true)
                .filter { _meta, _r1, _r2, has_bam -> !has_bam }
                .map { meta, r1, r2, _has_bam -> [ meta, r1, r2 ] },
            genome_with_index.map { genome, _index -> genome }
        )
        ch_versions = ch_versions.mix(RELOCATE3_ALIGN_READS.out.versions)
        ch_genotype_bams = ch_genotype_bams.mix(RELOCATE3_ALIGN_READS.out.bam)
    }

    RELOCATE3_CHARACTERIZE(
        RELOCATE3_FIND_INSERTIONS.out.sites.join(ch_genotype_bams),
        genome_with_index
    )
    ch_versions = ch_versions.mix(RELOCATE3_CHARACTERIZE.out.versions)

    emit:
    versions = ch_versions
}
