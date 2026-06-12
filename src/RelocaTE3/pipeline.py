"""Single-sample orchestration for RelocaTE3.

This is the zero-workflow-engine path: run the pipeline for one sample on a
laptop or a single HPC node. A workflow engine (Nextflow) can instead call the
individual CLI subcommands to scatter work across chromosomes/samples.

Use :func:`run_sample` for one sample or :func:`run_samples` to process
multiple samples in parallel with a thread pool.
"""

from __future__ import annotations

import subprocess
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from RelocaTE3 import __version__, logger
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


def _tool_version(tool: str, *extra_flags: str) -> str:
    """Return the first output line of ``tool --version`` (or first stderr line).

    Falls back to ``"unknown"`` if the binary is not found or fails.
    """
    flags = list(extra_flags) if extra_flags else ["--version"]
    try:
        result = subprocess.run(
            [tool, *flags],
            capture_output=True,
            text=True,
            timeout=10,
        )
        output = result.stdout or result.stderr
        first = next((ln for ln in output.splitlines() if ln.strip()), "")
        return first.strip() or "unknown"
    except Exception:
        return "unknown"


def log_provenance(te_library: str, genome: str, threads: int) -> None:
    """Log RelocaTE3 version, external tool versions, and key run parameters."""
    logger.info("RelocaTE3 %s", __version__)
    logger.info("  minimap2 : %s", _tool_version("minimap2", "--version"))
    logger.info("  samtools : %s", _tool_version("samtools", "--version"))
    logger.info("  bedtools : %s", _tool_version("bedtools", "--version"))
    logger.info("  TE library : %s", te_library)
    logger.info("  genome     : %s", genome)
    logger.info("  threads    : %d", threads)


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
    """Run the full pipeline for a single sample.

    Steps executed in order:

    - **Step 3** — identify TE-containing reads, trim the TE portion, emit
      flanking reads and ``read_repeat_name.txt``.
    - **Step 4** — re-align flanking reads and supporting mates to the genome.
    - **Step 5** — cluster junction/supporting reads and call non-reference
      insertions.
    - **Steps 0/6** *(optional, requires ``repeatmasker``)* — build the
      reference TE annotation, call reference/shared insertions, and drop
      non-reference calls overlapping a known reference TE of the same family.
    - **Step 7** *(optional, requires ``genotype=True``)* — align the original
      reads to the genome and classify each insertion as
      homozygous/heterozygous/somatic.

    Returns the path to the non-reference insertion GFF.
    """
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    t0 = time.monotonic()
    log_provenance(te_library, genome, threads)

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

    elapsed = time.monotonic() - t0
    logger.info(
        "Pipeline complete for %s in %.1f s; results at %s",
        reads.name,
        elapsed,
        gff_path,
    )
    return gff_path


def run_samples(
    samples: list[ReadLibrary],
    te_library: str,
    genome: str,
    outdir: str | Path,
    repeatmasker: str | None = None,
    genotype: bool = False,
    sample_threads: int = 1,
    step_threads: int = 1,
    len_cut_match: int = 10,
    len_cut_trim: int = 10,
    mismatch_allowance: int = 0,
    required_junction_reads: int = 1,
    verbose: int = 0,
) -> list[Path | None]:
    """Run the pipeline for multiple samples, up to ``sample_threads`` at a time.

    Each sample runs :func:`run_sample` in its own subdirectory
    ``<outdir>/<sample.name>/``.  ``step_threads`` controls the CPU threads
    given to each individual alignment step; the total CPU load is at most
    ``sample_threads × step_threads``.

    Returns a list of GFF paths in the same order as ``samples``.  Samples
    that fail are logged as errors and their position in the list is ``None``.
    """
    outdir = Path(outdir)
    results: dict[str, Path | None] = {}

    def _run(reads: ReadLibrary) -> tuple[str, Path]:
        gff = run_sample(
            reads,
            te_library,
            genome,
            outdir / reads.name,
            repeatmasker=repeatmasker,
            genotype=genotype,
            threads=step_threads,
            len_cut_match=len_cut_match,
            len_cut_trim=len_cut_trim,
            mismatch_allowance=mismatch_allowance,
            required_junction_reads=required_junction_reads,
            verbose=verbose,
        )
        return reads.name, gff

    with ThreadPoolExecutor(max_workers=max(1, sample_threads)) as pool:
        futures = {pool.submit(_run, reads): reads.name for reads in samples}
        for future in as_completed(futures):
            name = futures[future]
            try:
                _, gff = future.result()
                results[name] = gff
                logger.info("Sample %s finished → %s", name, gff)
            except Exception as exc:
                logger.error("Sample %s failed: %s", name, exc)
                results[name] = None

    return [results[s.name] for s in samples]
