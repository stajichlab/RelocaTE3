"""Contract tests for the pluggable aligner backends (Step 3/4 alignment)."""

import random
import shutil
import tempfile
import unittest
from pathlib import Path

import pysam

from RelocaTE3.aligners import get_aligner, psl_to_sam
from RelocaTE3.ReadLibrary import ReadLibrary

_BINARIES = {
    "minimap2": ["minimap2"],
    "bwa": ["bwa"],
    "bwamem2": ["bwa-mem2"],
    "bowtie2": ["bowtie2", "bowtie2-build"],
    "blat": ["blat"],
}


def _available(name):
    return shutil.which("samtools") is not None and all(
        shutil.which(b) is not None for b in _BINARIES[name]
    )


def _refseq(n=2000, seed=7):
    r = random.Random(seed)
    return "".join(r.choice("ACGT") for _ in range(n))


REF = _refseq()


def _write_fa(path, seqs):
    with open(path, "w") as fh:
        for name, seq in seqs.items():
            fh.write(f">{name}\n{seq}\n")


def _write_fq(path, reads):
    with open(path, "w") as fh:
        for name, seq in reads:
            fh.write(f"@{name}\n{seq}\n+\n{'I' * len(seq)}\n")


def _assert_contract(test, bam):
    test.assertTrue(Path(bam).exists())
    test.assertTrue(Path(f"{bam}.bai").exists())
    with pysam.AlignmentFile(str(bam), "rb") as fh:
        test.assertEqual(fh.header["HD"].get("SO"), "coordinate")
        mapped = 0
        for rec in fh.fetch(until_eof=True):
            test.assertFalse(rec.is_unmapped)  # mapped-only
            test.assertTrue(rec.query_name)  # names preserved
            test.assertTrue(rec.has_tag("NM"))  # NM present
            mapped += 1
        test.assertGreater(mapped, 0)


class TestBackendContract(unittest.TestCase):
    """Every installed backend must satisfy the downstream BAM contract."""

    def _run_backend(self, name):
        with tempfile.TemporaryDirectory() as d:
            d = Path(d)
            ref = d / "ref.fa"
            _write_fa(ref, {"chr1": REF})
            fq = d / "reads.fq"
            _write_fq(fq, [(f"r{i}", REF[i * 100 : i * 100 + 80]) for i in range(1, 6)])

            aln = get_aligner(name, threads=1)
            aln.index(ref)
            bam = aln.map_genome(ref, [fq], d / "g.bam", paired=False, threads=1)
            _assert_contract(self, bam)

            # TE-library stage (single-end ReadLibrary -> one "left" BAM)
            telib = d / "te.fa"
            _write_fa(telib, {"teA": REF[500:1000]})
            tefq = d / "te_reads.fq"
            _write_fq(
                tefq,
                [(f"t{i}", REF[500 + i * 20 : 500 + i * 20 + 80]) for i in range(1, 6)],
            )
            rl = ReadLibrary([str(tefq)], "S1")
            bams = aln.map_te_library(rl, telib, d, threads=1)
            self.assertEqual([b.name for b in bams], ["S1.left.bam"])
            _assert_contract(self, bams[0])

    @unittest.skipUnless(_available("minimap2"), "minimap2 not available")
    def test_minimap2(self):
        self._run_backend("minimap2")

    @unittest.skipUnless(_available("bwa"), "bwa not available")
    def test_bwa(self):
        self._run_backend("bwa")

    @unittest.skipUnless(_available("bwamem2"), "bwa-mem2 not available")
    def test_bwamem2(self):
        self._run_backend("bwamem2")

    @unittest.skipUnless(_available("bowtie2"), "bowtie2 not available")
    def test_bowtie2(self):
        self._run_backend("bowtie2")

    @unittest.skipUnless(_available("blat"), "blat not available")
    def test_blat(self):
        with tempfile.TemporaryDirectory() as d:
            d = Path(d)
            telib = d / "te.fa"
            _write_fa(telib, {"teA": REF[500:1000]})
            tefq = d / "te.fq"
            _write_fq(
                tefq,
                [(f"t{i}", REF[500 + i * 20 : 500 + i * 20 + 80]) for i in range(1, 6)],
            )
            rl = ReadLibrary([str(tefq)], "S1")
            aln = get_aligner("blat")
            bams = aln.map_te_library(rl, telib, d)
            _assert_contract(self, bams[0])
            # SEQ must be populated (not "*"): the trim step reads it to emit
            # flanking reads. blat's PSL carries no sequence, so the backend must
            # fill it from the query.
            with pysam.AlignmentFile(str(bams[0]), "rb") as fh:
                recs = list(fh.fetch(until_eof=True))
            self.assertTrue(recs and all(r.query_sequence for r in recs),
                            "blat BAM records have no SEQ")
            with self.assertRaises(NotImplementedError):
                aln.map_genome(telib, [tefq], d / "x.bam")


class TestBwaMateSuffix(unittest.TestCase):
    """bwa mem strips a trailing /1,/2 from read names. map_te_library maps each
    side separately (left=R1, right=R2), so it must restore the mate suffix
    (left -> /1, right -> /2); downstream flank-pairing needs it to find each
    junction flank's genomic mate. minimap2/bowtie2 keep the suffix natively."""

    @unittest.skipUnless(_available("bwa"), "bwa not available")
    def test_map_te_library_restores_mate_suffix(self):
        with tempfile.TemporaryDirectory() as d:
            d = Path(d)
            telib = d / "te.fa"
            _write_fa(telib, {"teA": REF[500:1000]})
            r1, r2 = d / "r1.fq", d / "r2.fq"
            _write_fq(r1, [(f"t{i}/1", REF[500 + i * 20 : 500 + i * 20 + 80]) for i in range(1, 6)])
            _write_fq(r2, [(f"t{i}/2", REF[600 + i * 20 : 600 + i * 20 + 80]) for i in range(1, 6)])
            rl = ReadLibrary([str(r1), str(r2)], "S1")
            aln = get_aligner("bwa", threads=1)
            bams = {b.name: b for b in aln.map_te_library(rl, telib, d, threads=1)}
            self.assertIn("S1.left.bam", bams)
            self.assertIn("S1.right.bam", bams)
            with pysam.AlignmentFile(str(bams["S1.left.bam"]), "rb") as fh:
                left = [r.query_name for r in fh.fetch(until_eof=True)]
            with pysam.AlignmentFile(str(bams["S1.right.bam"]), "rb") as fh:
                right = [r.query_name for r in fh.fetch(until_eof=True)]
            self.assertTrue(left and all(n.endswith("/1") for n in left), left)
            self.assertTrue(right and all(n.endswith("/2") for n in right), right)


def _revcomp(seq):
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]


class TestMinimapPairedGenome(unittest.TestCase):
    """minimap2 must honor ``paired=True`` in ``map_genome`` so junction flanks are
    anchored by their mate (short-read paired mode ``-ax sr``). Previously the
    backend concatenated R1+R2 and mapped single-end, so no read was paired."""

    @unittest.skipUnless(_available("minimap2"), "minimap2 not available")
    def test_map_genome_paired_produces_pairs(self):
        with tempfile.TemporaryDirectory() as d:
            d = Path(d)
            ref = d / "ref.fa"
            _write_fa(ref, {"chr1": REF})
            r1 = d / "r1.fq"
            r2 = d / "r2.fq"
            # two distinct pairs so a naive concat (all R1s then all R2s) cannot
            # accidentally interleave mates; real usage looks like this.
            _write_fq(r1, [("p", REF[500:600]), ("q", REF[900:1000])])
            _write_fq(
                r2,
                [("p", _revcomp(REF[700:800])), ("q", _revcomp(REF[1100:1200]))],
            )
            aln = get_aligner("minimap2", threads=1)
            aln.index(ref)
            bam = aln.map_genome(ref, [str(r1), str(r2)], d / "g.bam", paired=True, threads=1)
            with pysam.AlignmentFile(str(bam), "rb") as bf:
                recs = [r for r in bf.fetch(until_eof=True) if not r.is_unmapped]
            self.assertTrue(recs, "no reads mapped")
            self.assertTrue(
                any(r.is_paired for r in recs),
                "minimap2 map_genome(paired=True) did not pair the reads",
            )


class TestBowtie2JunctionSoftclip(unittest.TestCase):
    """Regression: bowtie2 TE-library mapping must soft-clip junction reads.

    A TE-junction read is part TE, part genomic flank -- exactly the reads that
    carry the non-reference insertion signal. In bowtie2's default end-to-end
    mode such a read must align over its full length, so the foreign flank half
    makes it align 0 times and no soft-clipped flank is recovered (find-insertions
    then returns nothing). ``Bowtie2Backend.map_te_library`` passes ``--local`` so
    bowtie2 soft-clips the flank and keeps the read, matching the bwa-mem and
    minimap2 backends. This test fails if ``--local`` is dropped: the junction
    reads stop aligning and no soft-clipped alignment is produced.

    (The ``TestBackendContract`` bowtie2 case uses reads wholly inside the TE,
    which align end-to-end, so it does not exercise this path.)
    """

    @unittest.skipUnless(_available("bowtie2"), "bowtie2 not available")
    def test_te_mapping_softclips_junction_reads(self):
        te = REF[500:1000]  # 500 bp TE
        flank = _refseq(n=500, seed=99)  # independent foreign flank
        # Junction reads: 40 bp TE + 40 bp foreign flank, in both orientations.
        reads = []
        for i in range(1, 6):
            te_chunk = te[100 + i * 20 : 100 + i * 20 + 40]
            fl_chunk = flank[i * 20 : i * 20 + 40]
            reads.append((f"j{i}_5", te_chunk + fl_chunk))  # TE then flank
            reads.append((f"j{i}_3", fl_chunk + te_chunk))  # flank then TE
        with tempfile.TemporaryDirectory() as d:
            d = Path(d)
            telib = d / "te.fa"
            _write_fa(telib, {"teA": te})
            tefq = d / "junction.fq"
            _write_fq(tefq, reads)
            rl = ReadLibrary([str(tefq)], "S1")
            aln = get_aligner("bowtie2", threads=1)
            bams = aln.map_te_library(rl, telib, d, threads=1)
            mapped = 0
            softclipped = 0
            with pysam.AlignmentFile(str(bams[0]), "rb") as fh:
                for rec in fh.fetch(until_eof=True):
                    if rec.is_unmapped:
                        continue
                    mapped += 1
                    if any(op == 4 for op, _ in (rec.cigartuples or [])):
                        softclipped += 1
            self.assertGreater(
                mapped, 0, "no junction reads aligned (end-to-end regression?)"
            )
            self.assertGreater(
                softclipped,
                0,
                "no soft-clipped junction alignments (bowtie2 needs --local)",
            )


class TestPslToSam(unittest.TestCase):
    """The PSL->SAM converter (BLAT backend) tested without the blat binary."""

    def test_single_block_exact(self):
        psl = "80\t0\t0\t0\t0\t0\t0\t0\t+\tt1\t80\t0\t80\tteA\t500\t100\t180\t1\t80,\t0,\t100,"
        sam = psl_to_sam([psl])
        self.assertEqual(len(sam), 1)
        f = sam[0].split("\t")
        self.assertEqual(f[0], "t1")
        self.assertEqual(f[1], "0")  # + strand
        self.assertEqual(f[2], "teA")
        self.assertEqual(f[3], "101")  # 1-based POS
        self.assertEqual(f[5], "80M")
        self.assertIn("NM:i:0", sam[0])

    def test_softclip_and_mismatch(self):
        psl = "68\t2\t0\t0\t0\t0\t0\t0\t+\tt2\t80\t5\t75\tteA\t500\t200\t270\t1\t70,\t5,\t200,"
        sam = psl_to_sam([psl])
        f = sam[0].split("\t")
        self.assertEqual(f[5], "5S70M5S")
        self.assertEqual(f[3], "201")
        self.assertIn("NM:i:2", sam[0])

    def test_skips_header_lines(self):
        self.assertEqual(psl_to_sam(["psLayout version 3", "", "match\tmis"]), [])


if __name__ == "__main__":
    unittest.main()
