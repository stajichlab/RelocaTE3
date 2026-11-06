/*
 * TE_SEARCH_SCATTER: the BLAT TE search (run-all steps 2-3, `relocaTE3 run`)
 * spread over tasks; its outputs are RELOCATE3_TE_SEARCH's, file for file.
 *
 *   TE_SPLIT         one task per sample: chunk FASTAs + manifest.tsv
 *   TE_SEARCH_CHUNK  one task per group of params.blat_chunks_per_task chunks
 *   TE_MERGE_TRIM    one task per sample, once all its chunk BAMs exist
 *
 * Grouping: a sample's chunk FASTAs (both sides) sorted by file name
 * (chunk.<id>.<side>.NNNNN.fa: left before right, then by index) and cut into
 * consecutive groups. The groups depend only on the split, never on task
 * timing, so -resume reuses every TE_SEARCH_CHUNK task. Each group's task
 * stages the manifest, which gives each chunk's mate.
 *
 * Gathering: the number of groups is known when the sample is split, so each
 * sample's chunk BAMs are collected with groupKey(meta, n_groups): its
 * TE_MERGE_TRIM starts when its own last group finishes, not when every
 * sample's do. The BAMs are sorted by name (te-merge-trim orders them by the
 * manifest anyway; sorting keeps the task's inputs stable for -resume).
 *
 * Take: reads       [ meta, R1, R2 ]   R2 is [] for single-end
 *       te_library  TE library FASTA
 * Emit: te_containing, flanking, read_repeat, versions: as RELOCATE3_TE_SEARCH
 */

include { RELOCATE3_TE_SPLIT        } from '../../modules/local/relocate3_te_split/main'
include { RELOCATE3_TE_SEARCH_CHUNK } from '../../modules/local/relocate3_te_search_chunk/main'
include { RELOCATE3_TE_MERGE_TRIM   } from '../../modules/local/relocate3_te_merge_trim/main'

workflow TE_SEARCH_SCATTER {
    take:
    reads
    te_library

    main:
    ch_versions = channel.empty()

    RELOCATE3_TE_SPLIT(reads)
    ch_versions = ch_versions.mix(RELOCATE3_TE_SPLIT.out.versions)

    // [ meta, manifest, [ [ chunk FASTAs of group 1 ], [ group 2 ], ... ] ]
    ch_groups = RELOCATE3_TE_SPLIT.out.manifest
        .join(RELOCATE3_TE_SPLIT.out.chunks)
        .map { meta, manifest, chunks ->
            [ meta, manifest, chunks.toSorted { a, b -> a.name <=> b.name }.collate(params.blat_chunks_per_task as int) ]
        }

    RELOCATE3_TE_SEARCH_CHUNK(
        ch_groups.flatMap { meta, manifest, groups -> groups.collect { group -> [ meta, manifest, group ] } },
        te_library
    )
    ch_versions = ch_versions.mix(RELOCATE3_TE_SEARCH_CHUNK.out.versions)

    // [ meta, [ every chunk BAM of the sample ] ], per sample as soon as its
    // n_groups tasks are done
    ch_chunk_bams = RELOCATE3_TE_SEARCH_CHUNK.out.bams
        .combine(ch_groups.map { meta, _manifest, groups -> [ meta, groups.size() ] }, by: 0)
        .map { meta, bams, n_groups -> [ groupKey(meta, n_groups), bams ] }
        .groupTuple()
        .map { key, bams -> [ key.getGroupTarget(), bams.flatten().toSorted { a, b -> a.name <=> b.name } ] }

    RELOCATE3_TE_MERGE_TRIM(
        reads
            .join(RELOCATE3_TE_SPLIT.out.manifest)
            .join(ch_chunk_bams),
        te_library
    )
    ch_versions = ch_versions.mix(RELOCATE3_TE_MERGE_TRIM.out.versions)

    emit:
    te_containing = RELOCATE3_TE_MERGE_TRIM.out.te_containing
    flanking      = RELOCATE3_TE_MERGE_TRIM.out.flanking
    read_repeat   = RELOCATE3_TE_MERGE_TRIM.out.read_repeat
    versions      = ch_versions
}
