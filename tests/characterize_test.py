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

    def test_skips_single_sided_support(self):
        """Sites lacking both left and right flanker support are skipped."""
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
                # R:0 -> no right support, not a supporting_junction -> skipped
                fh.write("mping\tTTA\tHEG4\tChr1\t998..1000\t+\tT:2\tR:0\tL:2\n")

            characterizer = Characterizer()
            txt_path, _ = characterizer.characterize(
                sites_file=Path(sites_file),
                bam_files=[Path(bam_path)],
                outdir=Path(workdir),
            )
            # only the header should be present
            self.assertEqual(len(Path(txt_path).read_text().splitlines()), 1)


if __name__ == "__main__":
    unittest.main()
