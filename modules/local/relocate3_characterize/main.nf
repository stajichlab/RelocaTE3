/*
 * RELOCATE3_CHARACTERIZE: run-all step 7 (zygosity of each insertion).
 *
 *   relocaTE3 characterize -s results/<target>.<te_name>.all_nonref_insert.txt
 *       -b <genotype BAM/CRAM> -g <genome> -o results [--excision]
 *       --samtools samtools --bcftools bcftools
 *
 * Arguments are those cmd_run_all passes to `characterize`. samtools and
 * bcftools are run-all's defaults: the tools come from the task environment.
 *
 * The BAM is an explicit path input: a sample's genotype_bam lives in meta
 * (nf-schema), and files in meta are not staged. Task wiring passes either
 * that BAM/CRAM with its index or RELOCATE3_ALIGN_READS' BAM.
 *
 * The sites table is staged as results/<name> and is not emitted; only the
 * new results/*.characTErized.{txt,gff} are.
 *
 * Input:  [ meta, sites table, BAM or CRAM, .bai/.crai ]
 *         [ genome, [ genome index files ] ]   (.fai needed for CRAM)
 * Output: results  [ meta, [ results/*.characTErized.{gff,txt} ] ]
 *         versions versions.yml
 */

process RELOCATE3_CHARACTERIZE {
    tag "${meta.id}"
    label 'process_single'

    input:
    tuple val(meta), path(sites, stageAs: 'results/*'), path(bam), path(bam_index)
    tuple path(genome), path(genome_index)

    output:
    tuple val(meta), path("results/*.characTErized.*"), emit: results
    path "versions.yml", emit: versions

    script:
    def excision = params.excision ? '--excision' : ''
    """
    relocaTE3 characterize \\
        -s ${sites} \\
        -b ${bam} \\
        -g ${genome} \\
        -o results \\
        ${excision} \\
        --samtools samtools \\
        --bcftools bcftools

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        relocate3: \$(relocaTE3 --version)
        samtools: \$(samtools --version | sed -n '1s/^samtools //p')
    END_VERSIONS
    """

    stub:
    def stem = sites.toString().replaceFirst(/\.txt$/, '')
    """
    mkdir -p results
    touch ${stem}.characTErized.gff ${stem}.characTErized.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        relocate3: stub
    END_VERSIONS
    """
}
