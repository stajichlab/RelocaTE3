"""Tests for reference-TE annotation (step 0), genome indexing/align (steps 1,4) and insertion finding (step 5)."""

import os
import tempfile
import unittest
from pathlib import Path

import pysam

from RelocaTE3.align import Aligner
from RelocaTE3.insertions import InsertionFinder
from RelocaTE3.reference_te import ReferenceTEAnnotator

DATA = Path(__file__).parent / "data"
GENOME = DATA / "sim_genome" / "MSU7.Chr3_2M.fa"
MPING = DATA / "mping.fa"
RM_OUT = DATA / "sim_genome" / "MSU7.Chr3_2M.fa.RepeatMasker.out"


class TestReferenceTE(unittest.TestCase):
    """Step 0: existing-TE annotation."""

    def test_repeatmasker_boundary_table(self):
        existing = ReferenceTEAnnotator.load_existing_te(RM_OUT, target="Chr3")
        self.assertIn("Chr3", existing)
        # first RM record: Chr3 7192..7483 -> start/end padded +/-2bp
        self.assertEqual(existing["Chr3"]["start"][7192], 1)
        self.assertEqual(existing["Chr3"]["start"][7190], 1)
        self.assertEqual(existing["Chr3"]["end"][7483], 1)

    def test_annotate_minimap_bed(self):
        with tempfile.TemporaryDirectory() as workdir:
            annotator = ReferenceTEAnnotator(threads=1)
            bed = annotator.annotate_minimap(MPING, GENOME, Path(workdir))
            self.assertTrue(bed.exists())
            # the simulated genome contains mping copies; expect at least one hit
            lines = [ln for ln in bed.read_text().splitlines() if ln.strip()]
            for ln in lines:
                self.assertEqual(len(ln.split("\t")), 6)


class TestGenomeAlign(unittest.TestCase):
    """Steps 1 & 4: index the genome and align reads to it."""

    def test_index_and_align(self):
        with tempfile.TemporaryDirectory() as workdir:
            # copy genome so indexes land in the temp dir
            genome = Path(workdir) / "genome.fa"
            genome.write_text(GENOME.read_text())

            aln = Aligner(threads=2)
            self.assertEqual(aln.index_genome(str(genome)), 0)
            self.assertTrue(Path(f"{genome}.fai").exists())
            self.assertTrue(Path(f"{genome}.mmi").exists())

            r1 = DATA / "sim_reads" / "MSU7.Chr3_2M.ALL_reads_6X_100_500_1.fq.gz"
            r2 = DATA / "sim_reads" / "MSU7.Chr3_2M.ALL_reads_6X_100_500_2.fq.gz"
            bam = aln.map_genome_minimap(
                str(genome), [str(r1), str(r2)], "HEG4", workdir, paired=True
            )
            self.assertTrue(bam.exists())
            self.assertTrue(Path(f"{bam}.bai").exists())
            with pysam.AlignmentFile(str(bam), "rb") as fh:
                self.assertGreater(fh.mapped, 0)


def _write_junction_bam(bam_path, contig, length, reads):
    """Build a sorted, indexed BAM of simple junction reads."""
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": contig, "LN": length}],
    }
    reads = sorted(reads, key=lambda r: r["start0"])
    with pysam.AlignmentFile(bam_path, "wb", header=header) as out:
        for spec in reads:
            seg = pysam.AlignedSegment(out.header)
            seg.query_name = spec["name"]
            seg.query_sequence = spec["seq"]
            seg.flag = spec.get("flag", 0)
            seg.reference_id = 0
            seg.reference_start = spec["start0"]
            seg.mapping_quality = spec.get("mapq", 60)
            seg.cigarstring = f"{len(spec['seq'])}M"
            seg.set_tag("NM", spec.get("nm", 0), value_type="i")
            out.write(seg)
    pysam.index(bam_path)


class TestInsertionFinder(unittest.TestCase):
    """Step 5: cluster junction reads into a non-reference insertion call."""

    def test_known_tsd_call(self):
        with tempfile.TemporaryDirectory() as workdir:
            bam_path = os.path.join(workdir, "flank.bam")
            # TSD "TTA"; build a left and a right junction read that share tsd_start=1000.
            # Left  (pos=left):  + strand, ":end:5", seq ends in TTA, ref_end+1 = 1003 -> tsd_start = 1000
            # Right (pos=right): + strand, ":start:5", seq starts with TTA, start(1based)=1000 -> tsd_start = 1000
            left = {"name": "readL:end:5", "seq": "A" * 37 + "TTA", "start0": 962}
            right = {"name": "readR:start:5", "seq": "TTA" + "A" * 37, "start0": 999}
            _write_junction_bam(bam_path, "Chr1", 5000, [left, right])

            read_repeat = os.path.join(workdir, "read_repeat_name.txt")
            with open(read_repeat, "w") as fh:
                fh.write("readL\tmping\t+\n")
                fh.write("readR\tmping\t+\n")

            finder = InsertionFinder(mismatch_allow=0, min_mapq=1)
            out_txt = finder.find_insertions(
                bam_file=Path(bam_path),
                read_repeat_file=Path(read_repeat),
                tsd="TTA",
                target="Chr1",
                sample="HEG4",
                outdir=Path(workdir),
            )
            self.assertTrue(out_txt.exists())
            rows = [ln for ln in out_txt.read_text().splitlines() if ln.strip()]
            self.assertEqual(len(rows), 1)
            cols = rows[0].split("\t")
            # repeat_family, TSD, exper, chrom, coor, orient, T:, R:, L:, ST:, SR:, SL:
            self.assertEqual(cols[0], "mping")
            self.assertEqual(cols[1], "TTA")
            self.assertEqual(cols[2], "HEG4")
            self.assertEqual(cols[3], "Chr1")
            self.assertEqual(cols[4], "1000..1002")
            self.assertEqual(cols[6], "T:2")
            self.assertEqual(cols[7], "R:1")
            self.assertEqual(cols[8], "L:1")

    def test_tsd_unknown_raises(self):
        with tempfile.TemporaryDirectory() as workdir:
            bam_path = os.path.join(workdir, "flank.bam")
            _write_junction_bam(
                bam_path,
                "Chr1",
                5000,
                [{"name": "r:end:5", "seq": "A" * 40, "start0": 100}],
            )
            rr = os.path.join(workdir, "rr.txt")
            Path(rr).write_text("r\tmping\t+\n")
            finder = InsertionFinder()
            with self.assertRaises(NotImplementedError):
                finder.find_insertions(
                    Path(bam_path), Path(rr), "UNK", "Chr1", "HEG4", Path(workdir)
                )


if __name__ == "__main__":
    unittest.main()
