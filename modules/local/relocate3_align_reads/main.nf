/*
 * RELOCATE3_ALIGN_READS: run-all step 7a, the genotyping alignment of all
 * original reads (only when a sample has no genotype_bam).
 *
 *   relocaTE3 align-reads -1 R1 [-2 R2] -g <genome> -n <id> -o . --threads <cpus>
 *
 * align-reads always uses minimap2 on the FASTA itself (minimap2 -x sr builds
 * its index in memory; the .mmi is not used), so only the genome is staged.
 * It is independent of the TE search and can run in parallel with it.
 *
 * Input:  [ meta, R1, R2 ]   R2 may be []
 *         genome FASTA
 * Output: bam      [ meta, genome_aln/<id>.reads.genome.bam, .bai ]
 *         versions versions.yml
 */

process RELOCATE3_ALIGN_READS {
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(reads_1), path(reads_2)
    path genome

    output:
    tuple val(meta), path("genome_aln/${meta.id}.reads.genome.bam"), path("genome_aln/${meta.id}.reads.genome.bam.bai"), emit: bam
    path "versions.yml", emit: versions

    script:
    def prefix = meta.id
    def right = reads_2 ? "-2 ${reads_2}" : ''
    """
    relocaTE3 align-reads \\
        -1 ${reads_1} \\
        ${right} \\
        -g ${genome} \\
        -n ${prefix} \\
        -o . \\
        --threads ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        relocate3: \$(relocaTE3 --version)
        minimap2: \$(minimap2 --version)
        samtools: \$(samtools --version | sed -n '1s/^samtools //p')
    END_VERSIONS
    """

    stub:
    def prefix = meta.id
    """
    mkdir -p genome_aln
    touch genome_aln/${prefix}.reads.genome.bam genome_aln/${prefix}.reads.genome.bam.bai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        relocate3: stub
    END_VERSIONS
    """
}
