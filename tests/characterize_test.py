"""Tests for the insertion-site Characterizer."""

import os
import tempfile
import unittest
from pathlib import Path

import pysam

from RelocaTE3.characterize import Characterizer


def _write_bam(bam_path: str, contig: str, length: int, reads: list[dict]):
    """Build a sorted, indexed BAM from a list of simple read specs."""
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": contig, "LN": length}],
    }
    reads = sorted(reads, key=lambda r: r["start"])
    with pysam.AlignmentFile(bam_path, "wb", header=header) as out:
        for i, spec in enumerate(reads):
            seg = pysam.AlignedSegment(out.header)
            seg.query_name = spec.get("name", f"read{i}")
            seg.query_sequence = "A" * spec["len"]
            seg.flag = 0
            seg.reference_id = 0
            seg.reference_start = spec["start"] - 1  # spec start is 1-based
            seg.mapping_quality = 60
            seg.cigarstring = spec["cigar"]
            seg.query_qualities = pysam.qualitystring_to_array("I" * spec["len"])
            seg.set_tag("NM", spec.get("nm", 0), value_type="i")
            out.write(seg)
    pysam.index(bam_path)


def _bam_to_cram(bam_path: str, cram_path: str, genome_fasta: str):
    """Convert a test BAM to an indexed CRAM using the supplied reference."""
    with pysam.AlignmentFile(bam_path, "rb") as source:
        with pysam.AlignmentFile(
            cram_path,
            "wc",
            template=source,
            reference_filename=genome_fasta,
        ) as out:
            for read in source:
                out.write(read)
    pysam.index(cram_path)


class TestCharacterizer(unittest.TestCase):
    """Exercise the core spanner/flanker classification."""

    def test_heterozygous_call(self):
        """Three clean spanning reads against 2 avg flankers -> heterozygous."""
        with tempfile.TemporaryDirectory() as workdir:
            bam_path = os.path.join(workdir, "reads.bam")
            # three reads fully spanning pos 1000 with no mismatches (spanners)
            spanning = [
                {"start": 985, "len": 40, "cigar": "40M", "nm": 0, "name": f"span{i}"}
                for i in range(3)
            ]
            # a clipped read at the junction should NOT count as a spanner
            clipped = [
                {"start": 990, "len": 40, "cigar": "10S30M", "nm": 0, "name": "clip"}
            ]
            _write_bam(bam_path, "Chr1", 2000, spanning + clipped)

            sites_file = os.path.join(workdir, "HEG4.mping.all_nonref.txt")
            with open(sites_file, "w") as fh:
                fh.write(
                    "TE\tTSD\tExper\tchromosome\tinsertion_site\tstrand\tT\tR\tL\n"
                )
                fh.write("mping\tTTA\tHEG4\tChr1\t998..1000\t+\tT:4\tR:2\tL:2\n")

            characterizer = Characterizer()
            txt_path, gff_path = characterizer.characterize(
                sites_file=Path(sites_file),
                bam_files=[Path(bam_path)],
                outdir=Path(workdir),
            )

            self.assertTrue(Path(txt_path).exists())
            self.assertTrue(Path(gff_path).exists())
            lines = Path(txt_path).read_text().splitlines()
            # header + one data row
            self.assertEqual(len(lines), 2)
            row = lines[1].split("\t")
            # strain, TE, TSD, chromosome.pos, strand, avg_flankers, spanners, status
            self.assertEqual(row[0], "HEG4")
            self.assertEqual(row[1], "mping")
            self.assertEqual(row[2], "TTA")
            self.assertEqual(row[3], "Chr1:998..1000")
            self.assertEqual(row[5], "2")  # avg_flankers = 4/2
            self.assertEqual(row[6], "3")  # three spanners
            self.assertEqual(row[7], "heterozygous")

    def test_cram_input_uses_reference_fasta(self):
        """A CRAM input is detected and produces the same spanner call as BAM."""
        with tempfile.TemporaryDirectory() as workdir:
            genome_fasta = os.path.join(workdir, "genome.fa")
            with open(genome_fasta, "w") as fh:
                fh.write(">Chr1\n")
                fh.write("A" * 2000 + "\n")
            pysam.faidx(genome_fasta)

            bam_path = os.path.join(workdir, "reads.bam")
            cram_path = os.path.join(workdir, "reads.cram")
            spanning = [
                {
                    "start": 985,
                    "len": 40,
                    "cigar": "40M",
                    "nm": 0,
                    "name": f"span{i}",
                }
                for i in range(3)
            ]
            _write_bam(bam_path, "Chr1", 2000, spanning)
            _bam_to_cram(bam_path, cram_path, genome_fasta)

            sites_file = os.path.join(workdir, "HEG4.mping.all_nonref.txt")
            with open(sites_file, "w") as fh:
                fh.write("mping\tTTA\tHEG4\tChr1\t998..1000\t+\tT:4\tR:2\tL:2\n")

            txt_path, _ = Characterizer().characterize(
                sites_file=Path(sites_file),
                bam_files=[Path(cram_path)],
                genome_fasta=Path(genome_fasta),
                outdir=Path(workdir),
            )

            row = Path(txt_path).read_text().splitlines()[1].split("\t")
            self.assertEqual(row[6], "3")
            self.assertEqual(row[7], "heterozygous")

    def test_homozygous_no_spanners(self):
        """No spanning reads over the site -> homozygous insertion."""
        with tempfile.TemporaryDirectory() as workdir:
            bam_path = os.path.join(workdir, "reads.bam")
            # read that does not span the +/-5bp window around pos 1000
            _write_bam(
                bam_path,
                "Chr1",
                2000,
                [{"start": 800, "len": 40, "cigar": "40M", "nm": 0}],
            )

            sites_file = os.path.join(workdir, "HEG4.mping.all_nonref.txt")
            with open(sites_file, "w") as fh:
                fh.write("mping\tTTA\tHEG4\tChr1\t998..1000\t+\tT:8\tR:4\tL:4\n")

            characterizer = Characterizer()
            txt_path, _ = characterizer.characterize(
                sites_file=Path(sites_file),
                bam_files=[Path(bam_path)],
                outdir=Path(workdir),
            )
            row = Path(txt_path).read_text().splitlines()[1].split("\t")
            self.assertEqual(row[6], "0")  # no spanners
            self.assertEqual(row[7], "homozygous")

    def test_skips_single_sided_multi_read_junction(self):
        """A one-sided call is dropped even when it carries a real TSD.

        RelocaTE2's gate is ``characterizer.pl:91``::

            if ( ($left_count >= 1 and $right_count >= 1)
                 or $TSD eq 'supporting_junction' )

        so junction reads on one side only never reach the characterized
        output, whatever the TSD column holds. This gate previously used ``or``
        plus a ``total_count > 1`` guard, which admitted every one-sided
        multi-read cluster and produced the riceTElib precision collapse.
        """
        with tempfile.TemporaryDirectory() as workdir:
            bam_path = os.path.join(workdir, "reads.bam")
            _write_bam(
                bam_path,
                "Chr1",
                2000,
                [{"start": 985, "len": 40, "cigar": "40M", "nm": 0}],
            )

            sites_file = os.path.join(workdir, "HEG4.mping.all_nonref.txt")
            with open(sites_file, "w") as fh:
                # single-sided (R:0) but multi-read (T:2 L:2) with a real TSD
                fh.write("mping\tTTA\tHEG4\tChr1\t998..1000\t+\tT:2\tR:0\tL:2\n")

            characterizer = Characterizer()
            txt_path, _ = characterizer.characterize(
                sites_file=Path(sites_file),
                bam_files=[Path(bam_path)],
                outdir=Path(workdir),
            )
            # header only -- the one-sided site was dropped
            self.assertEqual(len(Path(txt_path).read_text().splitlines()), 1)

    def test_keeps_single_sided_supporting_junction(self):
        """The ``supporting_junction`` sentinel is RelocaTE2's one-sided exception.

        It is the only one-sided class ``characterizer.pl:91`` admits, and the
        only one ``clean_false_positive.py:99,107`` does not grep out.
        """
        with tempfile.TemporaryDirectory() as workdir:
            bam_path = os.path.join(workdir, "reads.bam")
            _write_bam(
                bam_path,
                "Chr1",
                2000,
                [{"start": 985, "len": 40, "cigar": "40M", "nm": 0}],
            )

            sites_file = os.path.join(workdir, "HEG4.mping.all_nonref.txt")
            with open(sites_file, "w") as fh:
                fh.write(
                    "mping\tsupporting_junction\tHEG4\tChr1\t998..1000\t+"
                    "\tT:2\tR:0\tL:2\n"
                )

            characterizer = Characterizer()
            txt_path, _ = characterizer.characterize(
                sites_file=Path(sites_file),
                bam_files=[Path(bam_path)],
                outdir=Path(workdir),
            )
            # header + one data row
            self.assertEqual(len(Path(txt_path).read_text().splitlines()), 2)

    def test_skips_singleton(self):
        """A singleton (T:1) has no real junction support and must be skipped."""
        with tempfile.TemporaryDirectory() as workdir:
            bam_path = os.path.join(workdir, "reads.bam")
            _write_bam(
                bam_path,
                "Chr1",
                2000,
                [{"start": 985, "len": 40, "cigar": "40M", "nm": 0}],
            )

            sites_file = os.path.join(workdir, "HEG4.mping.all_nonref.txt")
            with open(sites_file, "w") as fh:
                fh.write("mping\tsingleton\tHEG4\tChr1\t998..1000\t+\tT:1\tR:0\tL:1\n")

            characterizer = Characterizer()
            txt_path, _ = characterizer.characterize(
                sites_file=Path(sites_file),
                bam_files=[Path(bam_path)],
                outdir=Path(workdir),
            )
            # only the header should be present
            self.assertEqual(len(Path(txt_path).read_text().splitlines()), 1)

    def test_preserves_te_family_evidence_metadata(self):
        """Step 7 appends family evidence without changing its first eight columns."""
        with tempfile.TemporaryDirectory() as workdir:
            bam_path = os.path.join(workdir, "reads.bam")
            _write_bam(bam_path, "Chr1", 2000, [])

            sites_file = os.path.join(workdir, "HEG4.mping.all_nonref.txt")
            with open(sites_file, "w") as fh:
                fh.write(
                    "mPing\tTTA\tHEG4\tChr1\t998..1000\t+\tT:3\tR:2\tL:1\t"
                    "ST:0\tSR:0\tSL:0\t"
                    "TE_family_support:mPing=2,RIRE3=1\t"
                    "TE_family_confidence:0.666667\t"
                    "TE_family_status:dominant\t"
                    "TE_supporting_family_support:RIRE3=2\t"
                    "TE_supporting_family_confidence:1.000000\t"
                    "TE_supporting_family_status:unique\t"
                    "TE_family_concordance:discordant\n"
                )

            txt_path, gff_path = Characterizer().characterize(
                sites_file=Path(sites_file),
                bam_files=[Path(bam_path)],
                outdir=Path(workdir),
            )

            header, row = Path(txt_path).read_text().splitlines()
            assert header.split("\t")[-7:] == [
                "TE_family_support",
                "TE_family_confidence",
                "TE_family_status",
                "TE_supporting_family_support",
                "TE_supporting_family_confidence",
                "TE_supporting_family_status",
                "TE_family_concordance",
            ]
            assert row.split("\t")[-7:] == [
                "mPing=2,RIRE3=1",
                "0.666667",
                "dominant",
                "RIRE3=2",
                "1.000000",
                "unique",
                "discordant",
            ]
            gff = Path(gff_path).read_text()
            assert "TE_family_support=mPing=2,RIRE3=1;" in gff
            assert "TE_family_confidence=0.666667;" in gff
            assert "TE_family_status=dominant;" in gff
            assert "TE_supporting_family_support=RIRE3=2;" in gff
            assert "TE_supporting_family_confidence=1.000000;" in gff
            assert "TE_supporting_family_status=unique;" in gff
            assert "TE_family_concordance=discordant" in gff


if __name__ == "__main__":
    unittest.main()
