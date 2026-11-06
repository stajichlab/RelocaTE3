/*
 * RELOCATE3_TE_SPLIT: BLAT TE search, step 1 of 3 (scatter; then
 * RELOCATE3_TE_SEARCH_CHUNK per group of chunks, RELOCATE3_TE_MERGE_TRIM).
 *
 *   relocaTE3 te-split --left R1 [--right R2] --name <id> --outdir .
 *       --chunk-reads <params.blat_chunk_reads>
 *
 * Converts each read side to FASTA chunks of blat_chunk_reads reads, as
 * `relocaTE3 run` does internally. Outputs do not depend on the chunk size.
 *
 * Input:  [ meta, R1, R2 ]   R2 may be [] (single-end)
 * Output: manifest  [ meta, <id>.te_chunks/manifest.tsv ]  (each chunk's side,
 *                   index, read count and mate)
 *         chunks    [ meta, [ <id>.te_chunks/chunk.<id>.<side>.NNNNN.fa ] ]
 *                   at least one; reads without any record fail here
 *         versions  versions.yml
 */

process RELOCATE3_TE_SPLIT {
    tag "${meta.id}"
    label 'process_single'

    input:
    tuple val(meta), path(reads_1), path(reads_2)

    output:
    tuple val(meta), path("${meta.id}.te_chunks/manifest.tsv"), emit: manifest
    tuple val(meta), path("${meta.id}.te_chunks/chunk.*.fa", arity: '1..*'), emit: chunks
    path "versions.yml", emit: versions

    script:
    def prefix = meta.id
    def right = reads_2 ? "--right ${reads_2}" : ''
    """
    # Temporary files (Python tempfile) in the task directory, not /tmp.
    export TMPDIR="\${PWD}"

    relocaTE3 te-split \\
        --left ${reads_1} \\
        ${right} \\
        --name ${prefix} \\
        --outdir . \\
        --chunk-reads ${params.blat_chunk_reads}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        relocate3: \$(relocaTE3 --version)
        seqtk: \$(seqtk 2>&1 | sed -n 's/^Version: *//p')
    END_VERSIONS
    """

    stub:
    // Three chunks per side, so tests can check the grouping.
    def prefix = meta.id
    def sides = reads_2 ? [ 'left', 'right' ] : [ 'left' ]
    def rows = sides.withIndex().collectMany { side, number ->
        (0..2).collect { i ->
            "${side}\t${i}\tchunk.${prefix}.${side}.0000${i}.fa\t1\t${reads_2 ? number + 1 : '.'}\t${number ? reads_2 : reads_1}"
        }
    }
    def chunks = sides.collectMany { side -> (0..2).collect { i -> "${prefix}.te_chunks/chunk.${prefix}.${side}.0000${i}.fa" } }
    """
    mkdir -p ${prefix}.te_chunks
    touch ${chunks.join(' ')}
    printf 'side\\tindex\\tchunk\\tn_reads\\tmate\\tread_file\\n' > ${prefix}.te_chunks/manifest.tsv
    printf '%s\\n' ${rows.collect { r -> "'${r}'" }.join(' ')} >> ${prefix}.te_chunks/manifest.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        relocate3: stub
    END_VERSIONS
    """
}
