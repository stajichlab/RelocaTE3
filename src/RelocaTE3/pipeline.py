"""Single-sample orchestration for RelocaTE3.

This is the zero-workflow-engine path: run the pipeline for one sample on a
laptop or a single HPC node. A workflow engine (Nextflow) can instead call the
individual CLI subcommands to scatter work across chromosomes/samples.

Currently implemented: Step 3 (TE read identification + trimming). Genome
re-alignment and insertion calling will be added as those modules land.
"""

from __future__ import annotations

from pathlib import Path

from RelocaTE3 import logger
from RelocaTE3.align import Aligner
from RelocaTE3.characterize import characterize_insertions, write_characterized
from RelocaTE3.genome_align import align_to_genome, read_read_repeat
from RelocaTE3.insertions import (
    find_insertions,
    write_insertions_gff,
    write_insertions_txt,
)
from RelocaTE3.librelocate import RelocaTE
from RelocaTE3.ReadLibrary import ReadLibrary
from RelocaTE3.reference_te import (
    filter_reference_overlaps,
    find_reference_insertions,
    write_existing_te_bed_from_rm,
)


def run_sample(
    reads: ReadLibrary,
    te_library: str,
    genome: str,
    outdir: str | Path,
    repeatmasker: str | None = None,
    genotype: bool = False,
    threads: int = 1,
    len_cut_match: int = 10,
    len_cut_trim: int = 10,
    mismatch_allowance: int = 0,
    required_junction_reads: int = 1,
    verbose: int = 0,
) -> Path:
    """Run the available pipeline steps for a single sample.

    Steps: (3) identify and trim TE-containing reads, (4) re-align flanking reads
    and supporting mates to the genome, (5) cluster junction/supporting reads and
    call non-reference insertions, (0/6) — when a RepeatMasker ``.out`` is supplied
    — build the reference TE annotation, call reference/shared insertions, and drop
    non-reference calls overlapping a known reference TE, and (7) — when
    ``genotype`` is set — align the original reads to the genome and classify each
    insertion as homozygous/heterozygous/somatic. Returns the path to the
    non-reference GFF.
    """
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    relocate = RelocaTE(TElib=te_library, threads=threads, verbose=verbose)
    logger.info(
        "Step 3: identifying and trimming TE-containing reads for %s", reads.name
    )
    relocate.identify_TE_reads(
        reads,
        outdir,
        TE_library=te_library,
        len_cut_match=len_cut_match,
        len_cut_trim=len_cut_trim,
        mismatch_allowance=mismatch_allowance,
    )

    logger.info("Step 4: re-aligning flanking reads to the genome for %s", reads.name)
    genome_bam, fullreads_bam = align_to_genome(reads, genome, outdir, threads=threads)

    logger.info("Step 5: clustering reads and calling insertions for %s", reads.name)
    rr_path = outdir / "te_containing" / f"{reads.name}.read_repeat_name.txt"
    read_repeat = read_read_repeat(rr_path) if rr_path.exists() else {}
    insertions = find_insertions(
        str(genome_bam),
        read_repeat,
        genome,
        fullreads_bam=str(fullreads_bam) if fullreads_bam else None,
        required_junction_reads=required_junction_reads,
    )
    results_dir = outdir / "results"
    results_dir.mkdir(parents=True, exist_ok=True)

    if repeatmasker:
        logger.info(
            "Steps 0/6: reference TE annotation and reference insertions for %s",
            reads.name,
        )
        bed_path = outdir / "existingTE.bed"
        reference_tes = write_existing_te_bed_from_rm(repeatmasker, bed_path)
        insertions = filter_reference_overlaps(insertions, reference_tes)
        ref_insertions = find_reference_insertions(
            str(genome_bam), read_repeat, reference_tes
        )
        write_insertions_gff(
            ref_insertions, results_dir / f"{reads.name}.all_ref_insert.gff", reads.name
        )
        write_insertions_txt(
            ref_insertions, results_dir / f"{reads.name}.all_ref_insert.txt"
        )

    gff_path = results_dir / f"{reads.name}.all_nonref_insert.gff"
    write_insertions_gff(insertions, gff_path, reads.name)
    write_insertions_txt(
        insertions, results_dir / f"{reads.name}.all_nonref_insert.txt"
    )

    if genotype:
        logger.info("Step 7: genotyping insertions for %s", reads.name)
        reads_bam = outdir / "genome_aln" / f"{reads.name}.reads.genome.bam"
        Aligner(threads).map_library_to_genome(
            genome, reads, str(reads_bam), cpu_threads=threads
        )
        characterize_insertions(insertions, str(reads_bam))
        write_characterized(
            insertions,
            results_dir / f"{reads.name}.all_nonref_insert.characTErized.gff",
            results_dir / f"{reads.name}.all_nonref_insert.characTErized.txt",
            reads.name,
        )

    return gff_path
