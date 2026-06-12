"""Characterize transposon insertion sites as homozygous, heterozygous, somatic or excised.

This is a Python port of ``characterizer.pl`` from RelocaTE2. It reads a
RelocaTE non-reference insertion table together with the BAM file(s) of the
original (pre-trimming) reads aligned to the reference genome. For every
candidate insertion site it counts the reads that *span* the site without a
clipped transposon junction ("spanners") and compares this to the number of
reads that flank the site ("flankers") to decide the zygosity / origin of the
insertion. Optionally it inspects spanning reads with indels for excision
footprints.
"""

from __future__ import annotations

import re
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path

import pysam

from RelocaTE3 import logger

# A read with no clipping or indels: a single fully-matched block.
_ALL_MATCH_CIGAR = re.compile(r"^\d+M$")
# Parse "Chr1 1446..1448" style coordinates -> the insertion (end) position.
_COORD_END = re.compile(r"\d+\.\.(\d+)")
_COORD_START = re.compile(r"(\d+)\.\.\d+")


def _format_number(value: float) -> str:
    """Render a number the way the original Perl did (no trailing ``.0``)."""
    return f"{value:g}"


class Characterizer:
    """Characterize RelocaTE insertion sites using read support from BAM files."""

    def __init__(
        self,
        samtools: str = "samtools",
        bcftools: str = "bcftools",
        verbose: int = 0,
    ):
        """Initialize the Characterizer.

        Args:
            samtools: path to the ``samtools`` executable (used for excision mode).
            bcftools: path to the ``bcftools`` executable (used for excision mode).
            verbose: verbosity level.
        """
        self.samtools = samtools
        self.bcftools = bcftools
        self.verbose = verbose

    # ------------------------------------------------------------------
    # public API
    # ------------------------------------------------------------------
    def characterize(
        self,
        sites_file: Path,
        bam_files: list[Path],
        genome_fasta: Path | None = None,
        outdir: Path | None = None,
        excision: bool = False,
    ) -> tuple[Path, Path]:
        """Characterize insertion sites and write ``.characTErized.txt``/``.gff``.

        Args:
            sites_file: RelocaTE non-reference insertion table (e.g.
                ``SAMPLE.mping.all_nonref.txt``).
            bam_files: BAM file(s) of the original reads aligned to the reference
                genome (before the TE was trimmed).
            genome_fasta: reference genome FASTA (required for ``excision``).
            outdir: directory for the output files (defaults to the directory of
                ``sites_file``).
            excision: also search for excision events that leave a footprint.

        Returns:
            Tuple of ``(txt_path, gff_path)`` for the written output files.
        """
        sites_file = Path(sites_file)
        if outdir is None:
            outdir = sites_file.parent
        outdir = Path(outdir)
        outdir.mkdir(parents=True, exist_ok=True)

        stem = sites_file.name
        if stem.endswith(".txt"):
            stem = stem[: -len(".txt")]
        txt_path = outdir / f"{stem}.characTErized.txt"
        gff_path = outdir / f"{stem}.characTErized.gff"

        alignments = [pysam.AlignmentFile(str(b), "rb") for b in bam_files]
        try:
            tsds: dict[str, dict[int, str]] = defaultdict(dict)
            # indel-containing spanning reads kept for excision analysis,
            # keyed by "chrom.pos"; the inner dict de-duplicates SAM lines.
            indel_reads: dict[str, dict[str, int]] = defaultdict(dict)
            statuses: dict[str, str] = {}
            to_print: dict = defaultdict(lambda: defaultdict(dict))

            self._score_sites(
                sites_file, alignments, tsds, indel_reads, statuses, to_print
            )

            if excision:
                if genome_fasta is None:
                    raise ValueError(
                        "excision analysis requires a reference genome FASTA"
                    )
                self._find_excision_footprints(
                    Path(genome_fasta), outdir, tsds, indel_reads, to_print
                )

            self._write_outputs(txt_path, gff_path, to_print)
        finally:
            for aln in alignments:
                aln.close()

        logger.info("Wrote %s and %s", txt_path, gff_path)
        return txt_path, gff_path

    # ------------------------------------------------------------------
    # scoring
    # ------------------------------------------------------------------
    def _score_sites(
        self, sites_file, alignments, tsds, indel_reads, statuses, to_print
    ):
        """Read the sites table and compute spanner/flanker support per site."""
        with open(sites_file) as handle:
            for line in handle:
                if "TE\tTSD\tExper\tchromosome\tinsertion_site" in line or re.search(
                    r"TE.TSD.Exper.chromosome.insertion_site", line
                ):
                    continue
                if not line.strip():
                    continue
                line = line.rstrip("\n")

                # mping TTA A119 Chr1 1446..1448 + T:1 R:0 L:1
                fields = line.split("\t")
                if len(fields) < 9:
                    continue
                (
                    te,
                    tsd,
                    exp,
                    chromosome,
                    coor,
                    te_orient,
                    total_string,
                    right_string,
                    left_string,
                ) = fields[:9]

                end_match = _COORD_END.search(coor)
                if not end_match:
                    continue
                pos = int(end_match.group(1))
                tsds[chromosome][pos] = tsd

                total_count = self._extract_count(total_string, "T")
                left_count = self._extract_count(left_string, "L")
                right_count = self._extract_count(right_string, "R")

                if not (
                    (left_count >= 1 and right_count >= 1)
                    or tsd == "supporting_junction"
                ):
                    continue

                site = f"{chromosome}.{pos}"
                spanners = self._count_spanners(
                    alignments, chromosome, pos, site, indel_reads
                )
                average_flankers = total_count / 2
                status = self._classify(average_flankers, spanners)

                statuses[site] = status
                record = to_print[chromosome][pos]
                record[tsd] = {
                    "TE": te,
                    "flank": average_flankers,
                    "span": spanners,
                    "status": status,
                    "strain": exp,
                    "coor": coor,
                    "TE_orient": te_orient,
                }

    @staticmethod
    def _extract_count(text: str, prefix: str) -> int:
        """Parse a "T:5"/"L:2"/"R:0" style count, defaulting to 0."""
        match = re.search(rf"{prefix}:(\d+)", text)
        return int(match.group(1)) if match else 0

    def _count_spanners(self, alignments, chromosome, pos, site, indel_reads) -> int:
        """Count reads that fully span the insertion site without clipping.

        Reads carrying indels (CIGAR I/N/D) over the site are stashed for the
        optional excision-footprint analysis rather than counted as spanners.
        """
        # samtools view chr:pos-pos (1-based inclusive) == fetch(chr, pos-1, pos)
        sam_lines: dict[str, int] = {}
        for aln in alignments:
            try:
                fetched = aln.fetch(chromosome, pos - 1, pos)
            except (ValueError, KeyError):
                # contig not present in this BAM
                continue
            for record in fetched:
                sam_lines.setdefault(record.to_string(), len(sam_lines))

        spanners = 0
        for sam_line in sam_lines:
            cols = sam_line.split("\t")
            cigar = cols[5]
            seq_len = len(cols[9]) if cols[9] != "*" else 0
            start = int(cols[3])
            end = start + seq_len - 1
            if end < pos + 5 or start > pos - 5:
                continue
            if _ALL_MATCH_CIGAR.match(cigar):
                xm = self._sam_tag(cols, "XM")
                nm = self._sam_tag(cols, "NM")
                if xm is not None:
                    if xm == 0:
                        spanners += 1
                elif nm is not None:
                    if nm == 0:
                        spanners += 1
                else:
                    spanners += 1
            elif re.search(r"[IND]", cigar):
                indel_reads[site][sam_line] = 1
        return spanners

    @staticmethod
    def _sam_tag(cols: list[str], tag: str) -> int | None:
        """Return an integer SAM tag value (e.g. ``NM:i:0``) or ``None``."""
        for col in cols[11:]:
            if col.startswith(f"{tag}:i:"):
                return int(col.split(":", 2)[2])
        return None

    @staticmethod
    def _classify(average_flankers: float, spanners: int) -> str:
        """Decide insertion status from flanker vs spanner support.

        Mirrors the cascading thresholds from ``characterizer.pl`` exactly.
        """
        if spanners == 0:
            return "homozygous"
        if average_flankers >= 5 and spanners < 5:
            return "homozygous/excision_no_footprint"
        if spanners < (average_flankers * 0.2) and spanners <= 10:
            return "homozygous/excision_no_footprint"
        if average_flankers <= 2 and spanners > 10:
            return "somatic_insertion"
        if abs(average_flankers - spanners) <= 5:
            return "heterozygous"
        if abs(average_flankers - spanners) - ((average_flankers + spanners) / 2) <= 10:
            return "heterozygous"
        if average_flankers > 10 and spanners > 10:
            return "heterozygous"
        if (spanners - average_flankers) > (
            spanners + average_flankers
        ) / 2 and average_flankers <= 10:
            return "somatic_insertion"
        return "other"

    # ------------------------------------------------------------------
    # excision footprint analysis
    # ------------------------------------------------------------------
    def _find_excision_footprints(
        self, genome_fasta, outdir, tsds, indel_reads, to_print
    ):
        """Call variants among indel spanning reads and flag excision footprints."""
        fai = Path(f"{genome_fasta}.fai")
        if not fai.exists():
            pysam.faidx(str(genome_fasta))
        header = self._fasta_header(genome_fasta)

        excision_info = outdir / "excisions_with_footprint.vcfinfo"
        with tempfile.TemporaryDirectory() as workdir, open(
            excision_info, "a"
        ) as info_out:
            for site, sam_dict in indel_reads.items():
                chromosome, loc = site.rsplit(".", 1)
                loc = int(loc)
                if chromosome not in to_print or loc not in to_print[chromosome]:
                    continue
                sam_lines = list(sam_dict.keys())
                if len(sam_lines) <= 1:
                    continue

                vcf_path = self._call_variants(
                    site, sam_lines, header, genome_fasta, Path(workdir)
                )
                if vcf_path is None:
                    continue
                self._parse_excision_vcf(
                    vcf_path, chromosome, loc, tsds, to_print, info_out
                )

    @staticmethod
    def _fasta_header(genome_fasta) -> pysam.AlignmentHeader:
        """Build a BAM header from the genome FASTA index (mirrors ``view -bT``)."""
        sq = []
        with open(f"{genome_fasta}.fai") as handle:
            for line in handle:
                name, length = line.split("\t")[:2]
                sq.append({"SN": name, "LN": int(length)})
        return pysam.AlignmentHeader.from_dict(
            {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": sq}
        )

    def _call_variants(
        self, site, sam_lines, header, genome_fasta, workdir
    ) -> Path | None:
        """Write per-site reads to a sorted BAM and call variants with bcftools."""
        bam_path = workdir / f"{site}.bam"
        sorted_bam = workdir / f"{site}.sorted.bam"
        vcf_path = workdir / f"{site}.var.flt.vcf"

        with pysam.AlignmentFile(str(bam_path), "wb", header=header) as out:
            for sam_line in sam_lines:
                out.write(pysam.AlignedSegment.fromstring(sam_line, header))
        pysam.sort("-o", str(sorted_bam), str(bam_path))
        pysam.index(str(sorted_bam))

        try:
            mpileup = subprocess.run(
                [
                    self.bcftools,
                    "mpileup",
                    "-C50",
                    "-f",
                    str(genome_fasta),
                    str(sorted_bam),
                ],
                capture_output=True,
                check=True,
            )
            with open(vcf_path, "wb") as vcf_out:
                subprocess.run(
                    [self.bcftools, "call", "-mv", "-Ov"],
                    input=mpileup.stdout,
                    stdout=vcf_out,
                    check=True,
                )
        except (subprocess.CalledProcessError, FileNotFoundError) as err:
            logger.warning("Variant calling failed for %s: %s", site, err)
            return None
        return vcf_path

    def _parse_excision_vcf(
        self, vcf_path, insert_ref, insert_pos, tsds, to_print, info_out
    ):
        """Inspect a per-site VCF for indels that indicate an excision footprint."""
        tsd = tsds.get(insert_ref, {}).get(insert_pos, "")
        tsd_len = len(tsd)
        with open(vcf_path) as handle:
            for line in handle:
                if line.startswith("#"):
                    continue
                line = line.rstrip("\n")
                cols = line.split("\t")
                if len(cols) < 5:
                    continue
                first_base = int(cols[1])
                ref_seq = cols[3]
                strain_seq = cols[4]

                aln_start = first_base - insert_pos - 1
                aln_end_ref = first_base - len(ref_seq) - insert_pos - 1
                aln_end_strain = first_base - len(strain_seq) - insert_pos - 1

                aln_start_near = abs(aln_start) <= tsd_len + 1
                aln_end_ref_near = abs(aln_end_ref) <= tsd_len + 1
                aln_end_strain_near = abs(aln_end_strain) <= tsd_len + 1

                end_ref = first_base + len(ref_seq) - 1
                end_strain = first_base + len(strain_seq) + 1
                insert_bwt_ends = (end_ref < insert_pos < end_strain) or (
                    end_strain < insert_pos < end_ref
                )
                all_after_insertion = (
                    (first_base - tsd_len + 1) > insert_pos
                    and end_ref > insert_pos
                    and end_strain > insert_pos
                )
                # ``aln_start_near`` is computed to mirror the original script,
                # which leaves it unused in the footprint decision.
                del aln_start_near

                if (
                    aln_end_ref_near or aln_end_strain_near or insert_bwt_ends
                ) and not all_after_insertion:
                    info_out.write(f"{insert_ref}.{insert_pos}\t{line}\n")
                    status = to_print[insert_ref][insert_pos][tsd]["status"]
                    if "/excision_with_footprint" not in status:
                        to_print[insert_ref][insert_pos][tsd]["status"] = (
                            status + "/excision_with_footprint"
                        )

    # ------------------------------------------------------------------
    # output
    # ------------------------------------------------------------------
    def _write_outputs(self, txt_path, gff_path, to_print):
        """Write the tabular and GFF3 characterization outputs."""
        with open(txt_path, "w") as txt_out, open(gff_path, "w") as gff_out:
            txt_out.write(
                "strain\tTE\tTSD\tchromosome.pos\tstrand\tavg_flankers\tspanners\tstatus\n"
            )
            gff_out.write("##gff-version 3\n")

            for chrom in sorted(to_print):
                for pos in sorted(to_print[chrom]):
                    for tsd in sorted(to_print[chrom][pos]):
                        rec = to_print[chrom][pos][tsd]
                        flankers = _format_number(rec["flank"])
                        spanners = rec["span"]
                        coor = rec["coor"]
                        start_match = _COORD_START.search(coor)
                        start = start_match.group(1) if start_match else str(pos)
                        txt_out.write(
                            f"{rec['strain']}\t{rec['TE']}\t{tsd}\t{chrom}:{coor}\t"
                            f"{rec['TE_orient']}\t{flankers}\t{spanners}\t{rec['status']}\n"
                        )
                        gff_out.write(
                            f"{chrom}\t{rec['strain']}\ttransposable_element_attribute\t"
                            f"{start}\t{pos}\t{rec['TE_orient']}\t.\t.\t"
                            f"ID={chrom}.{pos}.spanners;avg_flankers={flankers};"
                            f"spanners={spanners};type={rec['status']};TE={rec['TE']};TSD={tsd}\n"
                        )
