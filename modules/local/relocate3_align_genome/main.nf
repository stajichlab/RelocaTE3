/*
 * RELOCATE3_ALIGN_GENOME: run-all step 4 (place flanks on the genome).
 *
 *   relocaTE3 align-genome -g <genome> -f flanking/<id>.*.flankingReads.fq -n <id> -o .
 *       -1 R1 [-2 R2] --genome-aligner --genome-opts= --threads <cpus>
 *
 * Arguments are those cmd_run_all passes to `align-genome` (paired=False is
 * the flag's default, so --paired is never given). The flank files are passed
 * sorted by path, as run-all does (Python sorted() on the same strings).
 * align-genome reads te_containing/<id>.te_hit_names.txt from -o, so
 * te_containing/ and flanking/ are staged under those names. The genome is
 * staged with its index files (RELOCATE3_INDEX_GENOME), so no index is built.
 *
 * Input:  [ meta, R1, R2, [ te_containing files ], [ flanking files ] ]   R2 may be []
 *         [ genome, [ genome index files ] ]
 * Output: bam       [ meta, <id>.repeat.<aligner>.sorted.bam, .bai ]
 *         fullreads [ meta, genome_aln/<id>.fullreads.genome.bam, .bai ]
 *                   optional: align-genome skips it when no junction read aligns
 *         versions  versions.yml
 */

process RELOCATE3_ALIGN_GENOME {
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(reads_1), path(reads_2), path(te_containing, stageAs: 'te_containing/*'), path(flanking, stageAs: 'flanking/*', arity: '1..*')
    tuple path(genome), path(genome_index)

    output:
    tuple val(meta), path("${meta.id}.repeat.*.sorted.bam"), path("${meta.id}.repeat.*.sorted.bam.bai"), emit: bam
    tuple val(meta), path("genome_aln/${meta.id}.fullreads.genome.bam"), path("genome_aln/${meta.id}.fullreads.genome.bam.bai"), emit: fullreads, optional: true
    path "versions.yml", emit: versions

    script:
    // Single-quoted in the shell: escape each ' as '\'' so the value reaches relocaTE3 unchanged.
    def genome_opts = params.genome_opts.toString().replace("'", "'\\''")
    def prefix = meta.id
    // A staged path prints as its stageAs location ("flanking/<file>"), the
    // string run-all sorts; match run-all's glob <id>.*.flankingReads.fq on the base name.
    def flanks = flanking
        .collect { f -> f.toString() }
        .findAll { f -> def base = f.tokenize('/')[-1]; base.startsWith("${prefix}.") && base.endsWith('.flankingReads.fq') }
        .sort()
    def right = reads_2 ? "-2 ${reads_2}" : ''
    """
    # Temporary files (Python tempfile) in the task directory, not /tmp.
    export TMPDIR="\${PWD}"

    relocaTE3 align-genome \\
        -g ${genome} \\
        -f ${flanks.join(' ')} \\
        -n ${prefix} \\
        -o . \\
        -1 ${reads_1} \\
        ${right} \\
        --genome-aligner ${params.genome_aligner} \\
        --genome-opts='${genome_opts}' \\
        --threads ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        relocate3: \$(relocaTE3 --version)
        bwa: \$(bwa 2>&1 | sed -n 's/^Version: //p')
        samtools: \$(samtools --version | sed -n '1s/^samtools //p')
    END_VERSIONS
    """

    stub:
    def prefix = meta.id
    """
    mkdir -p genome_aln
    touch ${prefix}.repeat.${params.genome_aligner}.sorted.bam ${prefix}.repeat.${params.genome_aligner}.sorted.bam.bai
    touch genome_aln/${prefix}.fullreads.genome.bam genome_aln/${prefix}.fullreads.genome.bam.bai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        relocate3: stub
    END_VERSIONS
    """
}
