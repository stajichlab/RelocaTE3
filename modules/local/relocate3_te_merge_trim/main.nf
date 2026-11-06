/*
 * RELOCATE3_TE_MERGE_TRIM: BLAT TE search, step 3 of 3 (scatter). Merges a
 * sample's chunk BAMs in manifest order and trims: the same files as
 * RELOCATE3_TE_SEARCH (`relocaTE3 run`), published to the same paths.
 *
 *   relocaTE3 te-merge-trim --left R1 [--right R2] --te-library <te_library>
 *       --name <id> --outdir . --manifest manifest.tsv --chunk-bams <BAMs>
 *       --threads <cpus> --min-match --min-trimmed --mismatch   (from params)
 *
 * Arguments: the "run" step's that te-merge-trim accepts
 * (tests/nextflow_modules_args_test.py). Chunk BAMs may come in any order.
 *
 * Input:  [ meta, R1, R2, manifest.tsv, [ every chunk BAM of the sample ] ]
 *         R2 may be [] (single-end)
 *         TE library FASTA
 * Output: as RELOCATE3_TE_SEARCH:
 *         te_bams       [ meta, [ <id>.{left,right}.bam ], [ .bai ] ]
 *         te_containing [ meta, [ te_containing/* ] ]
 *         flanking      [ meta, [ flanking/* ] ]
 *         te_portions   [ meta, [ te_portions/* ] ]
 *         read_repeat   [ meta, te_containing/<id>.read_repeat_name.txt ]  (also in te_containing)
 *         versions      versions.yml
 * The TE library's .fai, written next to the staged library, is not emitted.
 */

process RELOCATE3_TE_MERGE_TRIM {
    tag "${meta.id}"
    label 'process_merge_trim'

    input:
    tuple val(meta), path(reads_1), path(reads_2), path(manifest), path(chunk_bams, arity: '1..*')
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
    # Temporary files (Python tempfile) in the task directory, not /tmp.
    export TMPDIR="\${PWD}"

    relocaTE3 te-merge-trim \\
        --left ${reads_1} \\
        ${right} \\
        --te-library ${te_library} \\
        --name ${prefix} \\
        --outdir . \\
        --manifest ${manifest} \\
        --chunk-bams ${chunk_bams.join(' ')} \\
        --threads ${task.cpus} \\
        --min-match ${params.min_match} \\
        --min-trimmed ${params.min_trimmed} \\
        --mismatch ${params.mismatch}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        relocate3: \$(relocaTE3 --version)
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
