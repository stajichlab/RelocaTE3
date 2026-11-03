/*
 * RELOCATE3_FIND_REFERENCE: run-all steps 0/6 (reference/shared TE calls),
 * only with params.repeatmasker.
 *
 *   relocaTE3 find-reference -b <id>.repeat.<aligner>.sorted.bam
 *       -R te_containing/<id>.read_repeat_name.txt --repeatmasker <rm.out> -n <id> -o .
 *
 * Arguments are those cmd_run_all passes to `find-reference`.
 *
 * Input:  [ meta, genome BAM, .bai, read_repeat_name.txt ]
 *         RepeatMasker .out
 * Output: existing_te [ meta, existingTE.bed ]   (run-all writes it at the sample root)
 *         results     [ meta, [ results/<id>.all_ref_insert.{gff,txt} ] ]
 *         versions    versions.yml
 */

process RELOCATE3_FIND_REFERENCE {
    tag "${meta.id}"
    label 'process_single'

    input:
    tuple val(meta), path(bam), path(bai), path(read_repeat, stageAs: 'te_containing/*')
    path repeatmasker

    output:
    tuple val(meta), path("existingTE.bed"), emit: existing_te
    tuple val(meta), path("results/${meta.id}.all_ref_insert.*"), emit: results
    path "versions.yml", emit: versions

    script:
    def prefix = meta.id
    """
    # Temporary files (Python tempfile) in the task directory, not /tmp.
    export TMPDIR="\${PWD}"

    relocaTE3 find-reference \\
        -b ${bam} \\
        -R ${read_repeat} \\
        --repeatmasker ${repeatmasker} \\
        -n ${prefix} \\
        -o .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        relocate3: \$(relocaTE3 --version)
        pysam: \$(python -c 'import pysam; print(pysam.__version__)')
    END_VERSIONS
    """

    stub:
    def prefix = meta.id
    """
    mkdir -p results
    touch existingTE.bed results/${prefix}.all_ref_insert.gff results/${prefix}.all_ref_insert.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        relocate3: stub
    END_VERSIONS
    """
}
