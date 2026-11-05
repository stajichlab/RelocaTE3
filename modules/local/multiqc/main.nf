/*
 * MULTIQC: the run's summary report, once, after every sample.
 *
 *   multiqc --no-ansi --config multiqc_config.yml [--title <multiqc_title>]
 *       --comment <comment> --filename relocate3_report.html .
 *
 * Reads only summary files ("custom content", *_mqc.tsv, and
 * *_mqc_versions.yml for the Software versions section), never calls:
 *   - cohort.summary_mqc.tsv           (RELOCATE3_COHORT)
 *   - <id>.stats_mqc.tsv               (RELOCATE3_SAMPLE_STATS, one per sample)
 *   - relocate3_software_mqc_versions.yml (the tasks' versions.yml, collected)
 * Each file is staged in its own numbered dir (?/*); the caller sorts them.
 *
 * MultiQC 1.35 comes from the pipeline env (pixi profile) or the image.
 * PYTHONNOUSERSITE=1: a host ~/.local site-packages (or a broken
 * ~/.local/bin/multiqc install) must not shadow the env's packages under
 * -profile pixi; the image sets it already.
 *
 * Input:  [ summary files ]
 *         multiqc_config.yml (assets/)
 *         comment: report comment (pipeline version, software environment)
 * Output: report   relocate3_report.html
 *         data     relocate3_report_data/
 *         versions versions.yml
 */

process MULTIQC {
    label 'process_single'

    input:
    path multiqc_files, stageAs: '?/*'
    path multiqc_config
    val comment

    output:
    path "relocate3_report.html", emit: report
    path "relocate3_report_data", emit: data
    path "versions.yml", emit: versions

    script:
    def title = params.multiqc_title ? "--title ${shellQuote(params.multiqc_title.toString())}" : ''
    """
    export PYTHONNOUSERSITE=1

    multiqc \\
        --no-ansi \\
        --force \\
        --config ${multiqc_config} \\
        ${title} \\
        --comment ${shellQuote(comment.toString())} \\
        --filename relocate3_report.html \\
        .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$(multiqc --version | sed 's/^.*version //')
    END_VERSIONS
    """

    stub:
    """
    mkdir relocate3_report_data
    touch relocate3_report.html relocate3_report_data/multiqc_data.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: stub
    END_VERSIONS
    """
}

/*
 * A single-quoted shell word: each ' becomes '\'' so the text reaches the
 * command unchanged.
 */
def shellQuote(String text) {
    return "'" + text.replace("'", "'\\''") + "'"
}
