/*
 * PREPARE_GENOME: index the genome once for every sample.
 *
 * Runs `relocaTE3 index-genome` (run-all step 1): samtools .fai, minimap2 .mmi
 * and the genome aligner's index, so no per-sample task builds an index.
 *
 * Take: genome   genome FASTA (a path value)
 * Emit: genome_with_index  value channel [ genome, [ genome index files ] ]
 *       versions           versions.yml
 */

include { RELOCATE3_INDEX_GENOME } from '../../modules/local/relocate3_index_genome/main'

workflow PREPARE_GENOME {
    take:
    genome

    main:
    RELOCATE3_INDEX_GENOME(genome)

    emit:
    genome_with_index = RELOCATE3_INDEX_GENOME.out.index
    versions = RELOCATE3_INDEX_GENOME.out.versions
}
