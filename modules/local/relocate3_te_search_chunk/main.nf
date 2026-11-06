/*
 * RELOCATE3_TE_SEARCH_CHUNK: BLAT TE search, step 2 of 3 (scatter). One task
 * searches a group of up to params.blat_chunks_per_task chunk FASTAs of one
 * sample, task.cpus chunks at once (BLAT is single-threaded).
 *
 *   relocaTE3 te-search-chunk --te-library <te_library> --manifest manifest.tsv
 *       --chunk <FASTAs> --outdir . --threads <cpus> --te-opts=   (from params)
 *
 * The manifest gives each chunk's mate, so read names are those of
 * `relocaTE3 run`. Arguments: the "run" step's that te-search-chunk accepts
 * (tests/nextflow_modules_args_test.py).
 *
 * Input:  [ meta, manifest.tsv, [ chunk FASTAs ] ]
 *         TE library FASTA
 * Output: bams      [ meta, [ chunk.<id>.<side>.NNNNN.bam ] ]  one per chunk FASTA
 *         versions  versions.yml
 * The chunk BAMs' .bai and the TE library's .fai (next to the staged library)
 * are not emitted: TE_MERGE_TRIM reads the BAMs whole.
 */

process RELOCATE3_TE_SEARCH_CHUNK {
    // sample + first chunk, so the trace tells a sample's chunk tasks apart
    tag "${meta.id}:${chunks[0].baseName}"
    label 'process_blat_chunk'

    input:
    tuple val(meta), path(manifest), path(chunks, arity: '1..*')
    path te_library

    output:
    tuple val(meta), path("chunk.*.bam", arity: '1..*'), emit: bams
    path "versions.yml", emit: versions

    script:
    // Single-quoted in the shell: escape each ' as '\'' so the value reaches relocaTE3 unchanged.
    def te_opts = params.te_opts.toString().replace("'", "'\\''")
    """
    # Temporary files (Python tempfile) in the task directory, not /tmp.
    export TMPDIR="\${PWD}"

    relocaTE3 te-search-chunk \\
        --te-library ${te_library} \\
        --manifest ${manifest} \\
        --chunk ${chunks.join(' ')} \\
        --outdir . \\
        --threads ${task.cpus} \\
        --te-opts='${te_opts}'

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        relocate3: \$(relocaTE3 --version)
        blat: \$(blat 2>&1 | sed -n '1s/.*BLAT v\\. *\\([^ ]*\\).*/\\1/p')
        samtools: \$(samtools --version | sed -n '1s/^samtools //p')
    END_VERSIONS
    """

    stub:
    """
    for fa in ${chunks.join(' ')}; do
        touch "\${fa%.fa}.bam" "\${fa%.fa}.bam.bai"
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        relocate3: stub
    END_VERSIONS
    """
}
