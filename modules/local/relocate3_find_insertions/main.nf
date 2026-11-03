/*
 * RELOCATE3_FIND_INSERTIONS: run-all step 5 (non-reference insertion calls).
 *
 *   relocaTE3 find-insertions -b <id>.repeat.<aligner>.sorted.bam
 *       -R te_containing/<id>.read_repeat_name.txt -n <id> -o .
 *       --tsd --target --te-name --mismatch --min-mapq --size --distance
 *       --[no-]require-both-junctions                                  (from params)
 *       [--reference-ins <repeatmasker>] [--fullreads-bam genome_aln/<id>.fullreads.genome.bam]
 *
 * Arguments are those cmd_run_all passes to `find-insertions`.
 *
 * Optional inputs use an empty list, the standard Nextflow pattern for an
 * absent path:
 *   - fullreads:    [] when align-genome wrote no fullreads BAM (run-all passes
 *                   --fullreads-bam only if that file exists);
 *   - repeatmasker: [] without params.repeatmasker (run-all's
 *                   --reference-ins is then None).
 *
 * Input:  [ meta, genome BAM, .bai, read_repeat_name.txt, [ fullreads BAM, .bai ] or [] ]
 *         RepeatMasker .out, or []
 * Output: results  [ meta, [ results/<target>.<te_name>.all_nonref_* ] ]   (10 tables)
 *         sites    [ meta, results/<target>.<te_name>.all_nonref_insert.txt ] (CHARACTERIZE input)
 *         versions versions.yml
 */

process RELOCATE3_FIND_INSERTIONS {
    tag "${meta.id}"
    label 'process_single'

    input:
    tuple val(meta), path(bam), path(bai), path(read_repeat, stageAs: 'te_containing/*'), path(fullreads, stageAs: 'genome_aln/*', arity: '0..*')
    path repeatmasker

    output:
    tuple val(meta), path("results/${params.target}.${params.te_name}.all_nonref_*"), emit: results
    tuple val(meta), path("results/${params.target}.${params.te_name}.all_nonref_insert.txt"), emit: sites
    path "versions.yml", emit: versions

    script:
    def prefix = meta.id
    // Staged paths print as "genome_aln/<file>".
    def fullreads_bam = fullreads.collect { f -> f.toString() }.find { f -> f.endsWith('.bam') }
    def fullreads_arg = fullreads_bam ? "--fullreads-bam ${fullreads_bam}" : ''
    def reference_ins = repeatmasker ? "--reference-ins ${repeatmasker}" : ''
    def both_junctions = params.require_both_junctions ? '--require-both-junctions' : '--no-require-both-junctions'
    """
    relocaTE3 find-insertions \\
        -b ${bam} \\
        -R ${read_repeat} \\
        --tsd ${params.tsd} \\
        --target ${params.target} \\
        -n ${prefix} \\
        -o . \\
        --te-name ${params.te_name} \\
        ${reference_ins} \\
        ${fullreads_arg} \\
        --mismatch ${params.mismatch} \\
        --min-mapq ${params.min_mapq} \\
        ${both_junctions} \\
        --size ${params.size} \\
        --distance ${params.distance}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        relocate3: \$(relocaTE3 --version)
        pysam: \$(python -c 'import pysam; print(pysam.__version__)')
    END_VERSIONS
    """

    stub:
    def stem = "results/${params.target}.${params.te_name}"
    """
    mkdir -p results
    for s in all.gff all.txt gff high_conf.gff high_conf.txt raw.gff raw.txt txt; do
        touch ${stem}.all_nonref_insert.\${s}
    done
    touch ${stem}.all_nonref_supporting.gff ${stem}.all_nonref_supporting.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        relocate3: stub
    END_VERSIONS
    """
}
