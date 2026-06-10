"""Step 7: genotype non-reference insertions as homo/heterozygous/somatic.

A pysam reimplementation of RelocaTE2's ``characterizer.pl``. For each insertion
with junction support on both sides, it counts *spanner* reads — original reads
that map cleanly (full-match, no mismatches) across the insertion site, i.e. the
reference allele lacking the TE — and compares them to the *flanker* (junction)
reads supporting the TE allele to assign a zygosity/status.

This needs a whole-genome alignment of the *original untrimmed* reads, produced by
:meth:`RelocaTE3.align.Aligner.map_library_to_genome` (or supplied by the user).
The optional excision-with-footprint VCF analysis from the Perl tool is not ported.
"""

from __future__ import annotations

from pathlib import Path

import pysam

from RelocaTE3 import logger
from RelocaTE3.models import Insertion

# a spanner must extend at least this far past the site on both sides
SPAN_MARGIN = 5


def count_spanners(
    bam: pysam.AlignmentFile, chrom: str, pos: int, margin: int = SPAN_MARGIN
) -> int:
    """Count reference-allele reads spanning ``pos`` on ``chrom``.

    A spanner starts at least ``margin`` bp before and ends at least ``margin`` bp
    after the site, with a pure-match CIGAR (no soft-clips or indels) and no
    mismatches (NM == 0) — i.e. the locus reads as reference, with no TE present.
    """
    spanners = 0
    for rec in bam.fetch(chrom, max(pos - 1, 0), pos):
        if rec.is_unmapped:
            continue
        start = rec.reference_start + 1  # 1-based
        end = rec.reference_end  # 1-based inclusive
        if end < pos + margin or start > pos - margin:
            continue
        cigar = rec.cigartuples
        if not cigar or len(cigar) != 1 or cigar[0][0] != 0:  # single M operation
            continue
        nm = rec.get_tag("NM") if rec.has_tag("NM") else 0
        if nm == 0:
            spanners += 1
    return spanners


def classify_status(avg_flankers: float, spanners: int) -> str:
    """Assign a genotype/status from flanker vs spanner counts (RelocaTE2 ladder)."""
    if spanners == 0:
        return "homozygous"
    if avg_flankers >= 5 and spanners < 5:
        return "homozygous/excision_no_footprint"
    if spanners < avg_flankers * 0.2 and spanners <= 10:
        return "homozygous/excision_no_footprint"
    if avg_flankers <= 2 and spanners > 10:
        return "somatic_insertion"
    if abs(avg_flankers - spanners) <= 5:
        return "heterozygous"
    if abs(avg_flankers - spanners) - ((avg_flankers + spanners) / 2) <= 10:
        return "heterozygous"
    if avg_flankers > 10 and spanners > 10:
        return "heterozygous"
    if (spanners - avg_flankers) > (spanners + avg_flankers) / 2 and avg_flankers <= 10:
        return "somatic_insertion"
    return "other"


def characterize_insertions(
    insertions: list[Insertion], genome_reads_bam: str
) -> list[Insertion]:
    """Annotate insertions with spanner counts and genotype status in place.

    Only insertions with junction reads on both sides are genotyped (matching
    RelocaTE2); others keep an empty status. Returns the same list for chaining.
    """
    with pysam.AlignmentFile(genome_reads_bam, "rb") as bam:
        for ins in insertions:
            if ins.left_junction_reads >= 1 and ins.right_junction_reads >= 1:
                ins.spanners = count_spanners(bam, ins.chrom, ins.end)
                ins.status = classify_status(ins.avg_flankers, ins.spanners)
    n = sum(1 for i in insertions if i.status)
    logger.info("Genotyped %d insertions with two-sided junction support", n)
    return insertions


def write_characterized(
    insertions: list[Insertion], gff_path: str | Path, txt_path: str | Path, sample: str
) -> None:
    """Write the genotyped insertions as GFF3 and a tab-delimited table."""
    with open(txt_path, "w") as txt:
        txt.write(
            "strain\tTE\tTSD\tchromosome.pos\tstrand\tavg_flankers\tspanners\tstatus\n"
        )
        for ins in insertions:
            if not ins.status:
                continue
            txt.write(
                f"{sample}\t{ins.te_name}\t{ins.tsd}\t{ins.chrom}:{ins.start}..{ins.end}\t"
                f"{ins.strand}\t{ins.avg_flankers}\t{ins.spanners}\t{ins.status}\n"
            )

    with open(gff_path, "w") as gff:
        gff.write("##gff-version 3\n")
        for ins in insertions:
            if not ins.status:
                continue
            attrs = (
                f"ID={ins.chrom}.{ins.end}.spanners;avg_flankers={ins.avg_flankers};"
                f"spanners={ins.spanners};type={ins.status};TE={ins.te_name};TSD={ins.tsd}"
            )
            gff.write(
                f"{ins.chrom}\t{sample}\ttransposable_element_attribute\t"
                f"{ins.start}\t{ins.end}\t.\t{ins.strand}\t.\t{attrs}\n"
            )
