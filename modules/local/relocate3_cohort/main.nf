/*
 * RELOCATE3_COHORT: cohort tables from every characterized sample, once.
 *
 *   relocaTE3 cohort -i <id>=<table> [<id>=<table> ...] -o .
 *
 * Not a run-all step: run-batch calls the same subcommand (cmd_cohort). The
 * logic lives in src/RelocaTE3/cohort.py; this process only stages tables.
 *
 * Every table has the same name (<target>.<te_name>.all_nonref_insert.
 * characTErized.txt), so each is staged in its own numbered dir (tables/1/,
 * tables/2/, ...). Samples are passed as NAME=PATH, so a sample with no
 * characterized calls (header-only table) keeps its column of '.' in the
 * matrix. The caller sorts the pairs by sample id, so the task's inputs (and
 * -resume hash) do not depend on the order samples finish; the output does
 * not depend on input order anyway.
 *
 * Input:  [ [ sample ids ], [ characterized tables ] ]   (same order)
 * Output: table    cohort.characTErized.txt
 *         matrix   cohort.genotype_matrix.tsv
 *         summary  cohort.summary_mqc.tsv   (MultiQC custom content)
 *         versions versions.yml
 */

process RELOCATE3_COHORT {
    label 'process_single'

    input:
    tuple val(ids), path(tables, stageAs: 'tables/?/*', arity: '1..*')

    output:
    path "cohort.characTErized.txt", emit: table
    path "cohort.genotype_matrix.tsv", emit: matrix
    path "cohort.summary_mqc.tsv", emit: summary
    path "versions.yml", emit: versions

    script:
    if (ids.size() != tables.size()) {
        error("RELOCATE3_COHORT: ${ids.size()} sample ids but ${tables.size()} tables")
    }
    def inputs = [ids, tables].transpose().collect { id, table -> "${id}=${table}" }.join(' ')
    """
    relocaTE3 cohort \\
        -i ${inputs} \\
        -o .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        relocate3: \$(relocaTE3 --version)
    END_VERSIONS
    """

    stub:
    def stub_inputs = [ids, tables].transpose().collect { id, table -> "${id}=${table}" }.join(' ')
    """
    # relocaTE3 cohort -i ${stub_inputs} -o .   (named here so tests can check the inputs)
    touch cohort.characTErized.txt cohort.genotype_matrix.tsv cohort.summary_mqc.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        relocate3: stub
    END_VERSIONS
    """
}
