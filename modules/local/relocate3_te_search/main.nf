/*
 * RELOCATE3_TE_SEARCH: run-all steps 2-3 (TE-library search and trimming).
 *
 *   relocaTE3 run --left R1 [--right R2] -T <te_library> -n <id> -o . --threads <cpus>
 *       --te-aligner --te-opts= --min-match --min-trimmed --mismatch   (from params)
 *
 * Arguments are those cmd_run_all passes to `run`
 * (tests/nextflow_modules_args_test.py checks every one is here).
 *
 * Input:  [ meta, R1, R2 ]   R2 may be [] (single-end)
 *         TE library FASTA
 * Output: te_bams       [ meta, [ <id>.{left,right}.bam ], [ .bai ] ]
 *         te_containing [ meta, [ te_containing/* ] ]
 *         flanking      [ meta, [ flanking/* ] ]
 *         te_portions   [ meta, [ te_portions/* ] ]
 *         read_repeat   [ meta, te_containing/<id>.read_repeat_name.txt ]  (also in te_containing)
 *         versions      versions.yml
 * The TE library's .fai, written next to the staged library, is not emitted.
 */

process RELOCATE3_TE_SEARCH {
    tag "${meta.id}"
    label 'process_high'

    input:
    tuple val(meta), path(reads_1), path(reads_2)
    path te_library

    output:
    tuple val(meta), path("${meta.id}.*.bam"), path("${meta.id}.*.bam.bai"), emit: te_bams, optional: true
    tuple val(meta), path("te_containing/*"), emit: te_containing
    tuple val(meta), path("flanking/*"), emit: flanking
    tuple val(meta), path("te_portions/*"), emit: te_portions, optional: true
    tuple val(meta), path("te_containing/${meta.id}.read_repeat_name.txt"), emit: read_repeat
    path "versions.yml", emit: versions

    script:
    def prefix = meta.id
    def right = reads_2 ? "--right ${reads_2}" : ''
    """
    relocaTE3 run \\
        --left ${reads_1} \\
        ${right} \\
        --te-library ${te_library} \\
        --name ${prefix} \\
        --outdir . \\
        --threads ${task.cpus} \\
        --te-aligner ${params.te_aligner} \\
        --te-opts='${params.te_opts}' \\
        --min-match ${params.min_match} \\
        --min-trimmed ${params.min_trimmed} \\
        --mismatch ${params.mismatch}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        relocate3: \$(relocaTE3 --version)
        blat: \$(blat 2>&1 | sed -n '1s/.*BLAT v\\. *\\([^ ]*\\).*/\\1/p')
        samtools: \$(samtools --version | sed -n '1s/^samtools //p')
    END_VERSIONS
    """

    stub:
    def prefix = meta.id
    """
    mkdir -p te_containing flanking te_portions
    touch ${prefix}.left.bam ${prefix}.left.bam.bai ${prefix}.right.bam ${prefix}.right.bam.bai
    touch te_containing/${prefix}.left.ContainingReads.fq te_containing/${prefix}.right.ContainingReads.fq
    touch te_containing/${prefix}.read_repeat_name.txt te_containing/${prefix}.te_hit_names.txt
    touch flanking/${prefix}.left.flankingReads.fq flanking/${prefix}.right.flankingReads.fq
    touch te_portions/${prefix}.five_prime.fa te_portions/${prefix}.three_prime.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        relocate3: stub
    END_VERSIONS
    """
}
