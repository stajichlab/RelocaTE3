/*
 * RELOCATE3_SAMPLE_STATS: one sample's reads and calls for the MultiQC report.
 *
 *   relocaTE3 sample-stats -d . -n <id> --target <target> --te-name <te_name>
 *       -o <id>.stats_mqc.tsv
 *
 * Counts lines of files the sample's tasks already wrote (src/RelocaTE3/
 * sample_stats.py lists the files and columns); changes no file and is not a
 * run-all step. The inputs are staged in run-all's layout (te_containing/,
 * flanking/, results/), so the same command works on a run-all sample dir.
 *
 * Optional inputs use an empty list: characterized [] without genotyping,
 * reference [] without params.repeatmasker; their columns are then left out.
 *
 * Input:  [ meta, read_repeat_name.txt, [ flanking FASTQs ],
 *           [ FIND_INSERTIONS results ], [ CHARACTERIZE results ] or [],
 *           [ FIND_REFERENCE results ] or [] ]
 * Output: stats    [ meta, <id>.stats_mqc.tsv ]   (MultiQC custom content)
 *         versions versions.yml
 */

process RELOCATE3_SAMPLE_STATS {
    tag "${meta.id}"
    label 'process_single'

    input:
    tuple val(meta), path(read_repeat, stageAs: 'te_containing/*'), path(flanking, stageAs: 'flanking/*'), path(nonref, stageAs: 'results/*'), path(characterized, stageAs: 'results/*', arity: '0..*'), path(reference, stageAs: 'results/*', arity: '0..*')

    output:
    tuple val(meta), path("${meta.id}.stats_mqc.tsv"), emit: stats
    path "versions.yml", emit: versions

    script:
    """
    relocaTE3 sample-stats \\
        -d . \\
        -n ${meta.id} \\
        --target ${params.target} \\
        --te-name ${params.te_name} \\
        -o ${meta.id}.stats_mqc.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        relocate3: \$(relocaTE3 --version)
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.stats_mqc.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        relocate3: stub
    END_VERSIONS
    """
}
