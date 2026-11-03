/*
 * CAT_FASTQ: merge the lanes of one sample (sample-sheet rows sharing a name)
 * before the TE search.
 *
 * Gzip-safe: a gzipped lane is appended as is (gzip members concatenate into a
 * valid stream); a plain-text lane is compressed first. Lanes are written in
 * the order given, so read order follows the sample sheet. Lanes are staged in
 * numbered folders, so lanes with the same file name do not collide.
 *
 * Input:  [ meta, [ R1 lanes ], [ R2 lanes ] ]   ([] for single-end)
 * Output: reads    [ meta, <id>_1.fq.gz, [ <id>_2.fq.gz ] ]
 *                  _2 is declared with arity '0..*' (Nextflow rejects '0..1'
 *                  on outputs): a one-element list, or [] for single-end.
 *                  Both stage like the sample sheet's fastq_2 downstream.
 *         versions versions.yml
 */

process CAT_FASTQ {
    tag "${meta.id}"
    label 'process_single'

    input:
    tuple val(meta), path(reads_1, stageAs: 'lanes_1/?/*', arity: '1..*'), path(reads_2, stageAs: 'lanes_2/?/*', arity: '0..*')

    output:
    tuple val(meta), path("${meta.id}_1.fq.gz"), path("${meta.id}_2.fq.gz", arity: '0..*'), emit: reads
    path "versions.yml", emit: versions

    script:
    def prefix = meta.id
    def lanes_2 = reads_2 ? "cat_gz ${prefix}_2.fq.gz ${reads_2.join(' ')}" : ''
    """
    # Append each lane to \$1: gzip members as is, plain text compressed.
    cat_gz() {
        local out="\$1"; shift
        : > "\${out}"
        for f in "\$@"; do
            if [ "\$(head -c 2 "\${f}" | od -An -tx1 | tr -d ' \\n')" = "1f8b" ]; then
                cat "\${f}" >> "\${out}"
            else
                gzip -c -n "\${f}" >> "\${out}"
            fi
        done
    }

    cat_gz ${prefix}_1.fq.gz ${reads_1.join(' ')}
    ${lanes_2}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gzip: \$(gzip --version | sed -n '1s/^gzip //p')
    END_VERSIONS
    """

    stub:
    def prefix = meta.id
    def touch_2 = reads_2 ? "echo '' | gzip > ${prefix}_2.fq.gz" : ''
    """
    echo '' | gzip > ${prefix}_1.fq.gz
    ${touch_2}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gzip: stub
    END_VERSIONS
    """
}
