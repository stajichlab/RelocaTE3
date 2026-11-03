/*
 * RELOCATE3_INDEX_GENOME: run-all step 1, once per genome.
 *
 *   relocaTE3 index-genome -g <genome> --genome-aligner <params.genome_aligner>
 *
 * Builds the samtools .fai, the minimap2 .mmi and the genome aligner's own
 * index next to the staged genome. Emits the genome together with every index
 * file, so later steps stage them side by side and never build an index
 * themselves (align-genome only builds one when it is missing).
 *
 * Input:  genome FASTA
 * Output: index    [ genome, [ genome.* index files ] ]
 *         versions versions.yml
 */

process RELOCATE3_INDEX_GENOME {
    tag "${genome.name}"
    label 'process_medium'

    input:
    path genome

    output:
    tuple path(genome), path("${genome}.*"), emit: index
    path "versions.yml", emit: versions

    script:
    """
    # Temporary files (Python tempfile) in the task directory, not /tmp.
    export TMPDIR="\${PWD}"

    relocaTE3 index-genome \\
        -g ${genome} \\
        --genome-aligner ${params.genome_aligner}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        relocate3: \$(relocaTE3 --version)
        bwa: \$(bwa 2>&1 | sed -n 's/^Version: //p')
        minimap2: \$(minimap2 --version)
        samtools: \$(samtools --version | sed -n '1s/^samtools //p')
    END_VERSIONS
    """

    stub:
    """
    touch ${genome}.fai ${genome}.mmi ${genome}.bwt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        relocate3: stub
    END_VERSIONS
    """
}
