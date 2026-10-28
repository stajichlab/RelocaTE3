"""Tests for reference-TE annotation (step 0), genome indexing/align (steps 1,4) and insertion finding (step 5)."""

import os
import tempfile
import unittest
from pathlib import Path

import pysam

from RelocaTE3.align import Aligner
from RelocaTE3.insertions import InsertionFinder, _pair_breakpoints
from RelocaTE3.reference_te import ReferenceTEAnnotator


def _bp(counts):
    """Build a {position: [reads]} breakpoint map with the given per-position support."""
    return {pos: [object()] * n for pos, n in counts.items()}


class TestPairBreakpoints(unittest.TestCase):
    """RelocaTE2's nearest-right sub-clustering (TSD_from_read_depth:603-770)."""

    def test_each_left_pairs_with_its_nearest_right(self):
        # RelocaTE2 emits one sub-cluster per left position; a minor left beside
        # a dominant one is NOT absorbed, it gets its own pair (:610-676).
        pairs = _pair_breakpoints(_bp({1003: 5, 1005: 1}), _bp({1000: 5}))
        self.assertEqual(pairs, [(1003, 1000), (1005, 1000)])

    def test_lone_breakpoints_within_the_window_still_pair(self):
        # 60 bp apart is inside RelocaTE2's 100 bp pairing window, so it pairs.
        # The call is suppressed later, at emission, when TSD_len_calculate
        # returns a non-positive length (relocaTE_insertionFinder.py:818).
        pairs = _pair_breakpoints(_bp({100: 3}), _bp({40: 3}))
        self.assertEqual(pairs, [(100, 40)])

    def test_lone_breakpoints_beyond_the_window_are_split(self):
        pairs = _pair_breakpoints(_bp({1000: 3}), _bp({500: 3}))
        self.assertEqual(sorted(pairs, key=lambda p: (p[0] or 0, p[1] or 0)),
                         [(None, 500), (1000, None)])

    def test_distinct_sites_stay_separate(self):
        # two real insertions ~50 bp apart: each left takes its own nearest right.
        pairs = _pair_breakpoints(_bp({1003: 3, 1053: 3}), _bp({1000: 3, 1050: 3}))
        self.assertEqual(pairs, [(1003, 1000), (1053, 1050)])

DATA = Path(__file__).parent / "data"
GENOME = DATA / "sim_genome" / "MSU7.Chr3_2M.fa"
MPING = DATA / "TE_lib" / "mping.fa"
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
    """Build a sorted, indexed BAM of simple junction reads on one contig."""
    return _write_multi_contig_bam(bam_path, [(contig, length)], reads)


def _write_multi_contig_bam(bam_path, contigs, reads):
    """Build a sorted, indexed BAM across one or more contigs.

    Args:
        bam_path: output BAM path.
        contigs: list of (name, length) tuples; index in this list is the
            reference_id used by each read spec via the ``ref`` key (default 0).
        reads: list of read specs. Each must include ``name``, ``seq``,
            ``start0``; may include ``ref`` (contig index), ``flag``, ``mapq``,
            ``nm``.
    """
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": name, "LN": length} for name, length in contigs],
    }
    reads = sorted(reads, key=lambda r: (r.get("ref", 0), r["start0"]))
    with pysam.AlignmentFile(bam_path, "wb", header=header) as out:
        for spec in reads:
            seg = pysam.AlignedSegment(out.header)
            seg.query_name = spec["name"]
            seg.query_sequence = spec["seq"]
            seg.flag = spec.get("flag", 0)
            seg.reference_id = spec.get("ref", 0)
            seg.reference_start = spec["start0"]
            seg.mapping_quality = spec.get("mapq", 60)
            seg.cigarstring = f"{len(spec['seq'])}M"
            seg.set_tag("NM", spec.get("nm", 0), value_type="i")
            out.write(seg)
    pysam.index(bam_path)


class TestInsertionFinder(unittest.TestCase):
    """Step 5: cluster junction reads into a non-reference insertion call."""

    def test_target_all_preserves_per_cluster_chrom(self):
        """With --target ALL, the chrom column must be the cluster's real chrom, not 'ALL'."""
        with tempfile.TemporaryDirectory() as workdir:
            bam_path = os.path.join(workdir, "flank.bam")
            # Two independent insertion sites on two different chromosomes,
            # each with a left + right junction sharing tsd_start.
            reads = [
                # Chr1 site at tsd_start=1000
                {
                    "name": "rL1:end:5",
                    "seq": "A" * 37 + "TTA",
                    "start0": 962,
                    "ref": 0,
                },
                {
                    "name": "rR1:start:5",
                    "seq": "TTA" + "A" * 37,
                    "start0": 999,
                    "ref": 0,
                },
                # Chr2 site at tsd_start=2000
                {
                    "name": "rL2:end:5",
                    "seq": "A" * 37 + "TTA",
                    "start0": 1962,
                    "ref": 1,
                },
                {
                    "name": "rR2:start:5",
                    "seq": "TTA" + "A" * 37,
                    "start0": 1999,
                    "ref": 1,
                },
            ]
            _write_multi_contig_bam(bam_path, [("Chr1", 5000), ("Chr2", 5000)], reads)

            read_repeat = os.path.join(workdir, "read_repeat_name.txt")
            with open(read_repeat, "w") as fh:
                for n in ("rL1", "rR1", "rL2", "rR2"):
                    fh.write(f"{n}\tmping\t+\n")

            finder = InsertionFinder(mismatch_allow=0, min_mapq=1)
            out_txt = finder.find_insertions(
                bam_file=Path(bam_path),
                read_repeat_file=Path(read_repeat),
                tsd="UNK",
                target="ALL",
                sample="HEG4",
                outdir=Path(workdir),
            )
            rows = [ln for ln in out_txt.read_text().splitlines() if ln.strip()]
            self.assertEqual(len(rows), 2)
            chroms = sorted(r.split("\t")[3] for r in rows)
            self.assertEqual(chroms, ["Chr1", "Chr2"])
            for r in rows:
                self.assertNotEqual(r.split("\t")[3], "ALL")

    def test_tsd_unknown_infers_variable_length_tsd(self):
        with tempfile.TemporaryDirectory() as workdir:
            bam_path = os.path.join(workdir, "flank.bam")
            # A 5-bp TSD at 101..105: right start=101; left end=105.
            _write_junction_bam(
                bam_path,
                "Chr1",
                5000,
                [
                    {"name": "rR:start:5", "seq": "GATCA" + "A" * 35, "start0": 100},
                    {"name": "rL:end:5", "seq": "A" * 35 + "GATCA", "start0": 65},
                ],
            )
            rr = os.path.join(workdir, "rr.txt")
            Path(rr).write_text("rR\tmping\t+\nrL\tmping\t+\n")
            finder = InsertionFinder()
            output = finder.find_insertions(
                Path(bam_path), Path(rr), "UNK", "Chr1", "HEG4", Path(workdir)
            )
            row = output.read_text().strip().split("\t")
            self.assertEqual(row[0], "mping")
            self.assertEqual(row[1], "GATCA")
            self.assertEqual(row[4], "101..105")

    def test_tsd_unknown_pools_same_start_subcandidates(self):
        """Nearest-right pairs sharing a TSD start become one RelocaTE2-style call."""
        with tempfile.TemporaryDirectory() as workdir:
            bam_path = os.path.join(workdir, "flank.bam")
            reads = [
                # Site A: two left breakpoints share one right breakpoint. Both
                # pairs have two reads, so the shorter resolved TSD wins.
                {
                    "name": "aLeftShort:end:5",
                    "seq": "A" * 38 + "TA",
                    "start0": 61,
                },
                {
                    "name": "aLeftLong:end:5",
                    "seq": "A" * 35 + "TACTC",
                    "start0": 64,
                },
                {
                    "name": "aRight:start:5",
                    "seq": "TACTC" + "A" * 35,
                    "start0": 99,
                },
                # Site B: the longer-TSD pair has an extra left read, so its
                # three reads beat the two-read short-TSD pair.
                {
                    "name": "bLeftShort:end:5",
                    "seq": "A" * 38 + "TA",
                    "start0": 2961,
                },
                {
                    "name": "bLeftLong1:end:5",
                    "seq": "A" * 35 + "TACTC",
                    "start0": 2964,
                },
                {
                    "name": "bLeftLong2:end:5",
                    "seq": "A" * 35 + "TACTC",
                    "start0": 2964,
                },
                {
                    "name": "bRight:start:5",
                    "seq": "TACTC" + "A" * 35,
                    "start0": 2999,
                },
            ]
            _write_junction_bam(bam_path, "Chr1", 5000, reads)

            rr = Path(workdir) / "rr.txt"
            rr.write_text(
                "aLeftShort\tfamilyAlpha\t+\n"
                "aLeftLong\tfamilyBeta\t+\n"
                "aRight\tfamilyShared\t+\n"
                "bLeftShort\tfamilyShared\t+\n"
                "bLeftLong1\tfamilyShared\t+\n"
                "bLeftLong2\tfamilyShared\t+\n"
                "bRight\tfamilyShared\t+\n"
            )

            output = InsertionFinder().find_insertions(
                Path(bam_path), rr, "UNK", "Chr1", "HEG4", Path(workdir)
            )
            rows = [line.split("\t") for line in output.read_text().splitlines()]

            self.assertEqual(len(rows), 2)
            self.assertEqual(rows[0][0], "familyShared")
            self.assertEqual(rows[0][1], "TA")
            self.assertEqual(rows[0][4], "100..101")
            self.assertEqual(rows[0][6:9], ["T:4", "R:2", "L:2"])
            self.assertEqual(rows[1][1], "TACTC")
            self.assertEqual(rows[1][4], "3000..3004")
            self.assertEqual(rows[1][6:9], ["T:5", "R:2", "L:3"])

    def test_require_both_junctions_drops_single_sided(self):
        """require_both_junctions drops one-sided calls (RelocaTE2 parity) while
        keeping two-sided insertions; default keeps both."""
        with tempfile.TemporaryDirectory() as workdir:
            bam_path = os.path.join(workdir, "flank.bam")
            # Site A: two-sided (left + right) at tsd_start~1000.
            # Site C: left-only junction ~1000 bp away, with three reads.
            #
            # Three, not one: RelocaTE2's cluster arbitration
            # (write_output:271-283) discards a one-sided candidate backed by
            # fewer than MIN_ONE_SIDED_JUNCTIONS reads whenever the same cluster
            # also holds a two-sided candidate. A single-read site C would be
            # removed there and never reach the require_both_junctions test.
            reads = [
                {"name": "aR:start:5", "seq": "TTA" + "A" * 37, "start0": 999},
                {"name": "aL:end:5", "seq": "A" * 37 + "TTA", "start0": 962},
                {"name": "cL1:end:5", "seq": "A" * 37 + "TTA", "start0": 1962},
                {"name": "cL2:end:5", "seq": "A" * 37 + "TTA", "start0": 1962},
                {"name": "cL3:end:5", "seq": "A" * 37 + "TTA", "start0": 1962},
            ]
            _write_junction_bam(bam_path, "Chr1", 5000, reads)
            read_repeat = os.path.join(workdir, "rr.txt")
            with open(read_repeat, "w") as fh:
                for n in ("aR", "aL", "cL1", "cL2", "cL3"):
                    fh.write(f"{n}\tmping\t+\n")

            def _run(finder):
                out = finder.find_insertions(
                    bam_file=Path(bam_path), read_repeat_file=Path(read_repeat),
                    tsd="UNK", target="Chr1", sample="HEG4", outdir=Path(workdir),
                )
                return [ln for ln in out.read_text().splitlines() if ln.strip()]

            # Permissive policy must be asked for explicitly: requiring both
            # junctions is the default (see the flag's docstring for why).
            permissive_rows = _run(InsertionFinder(
                mismatch_allow=0, min_mapq=1, require_both_junctions=False))
            self.assertEqual(len(permissive_rows), 2)  # two-sided + single-sided

            default_rows = _run(InsertionFinder(mismatch_allow=0, min_mapq=1))
            self.assertEqual(len(default_rows), 1)  # single-sided dropped by default

            both_rows = _run(InsertionFinder(
                mismatch_allow=0, min_mapq=1, require_both_junctions=True))
            self.assertEqual(both_rows, default_rows)  # explicit == default
            cols = both_rows[0].split("\t")
            self.assertNotEqual(cols[7], "R:0")  # kept call has a right junction
            self.assertNotEqual(cols[8], "L:0")  # ...and a left junction

    def test_supporting_junction_keeps_one_sided_call_with_opposite_support(self):
        """The strict policy retains RelocaTE2's narrow one-sided exception."""
        with tempfile.TemporaryDirectory() as workdir:
            bam_path = os.path.join(workdir, "flank.bam")
            reads = [
                # Right junction at the left edge of the legacy 3-base sentinel.
                {"name": "junction:start:5", "seq": "A" * 40, "start0": 999},
                # Unpaired, forward read ending to its left: opposite-side support.
                {"name": "support", "seq": "A" * 40, "start0": 900},
                # Left junction at 3000; its sentinel begins one base beyond the
                # aligned flank, matching RelocaTE2's reference_end + 1 convention.
                {"name": "left-junction:end:5", "seq": "A" * 40, "start0": 2960},
                # Reverse read beginning to its right: opposite-side support.
                {
                    "name": "right-support",
                    "seq": "A" * 40,
                    "start0": 3100,
                    "flag": 16,
                },
            ]
            _write_junction_bam(bam_path, "Chr1", 5000, reads)
            repeat = Path(workdir) / "read_repeat_name.txt"
            repeat.write_text("junction\tmping\t+\nleft-junction\tmping\t+\n")

            out = InsertionFinder(mismatch_allow=0, min_mapq=1).find_insertions(
                bam_file=Path(bam_path),
                read_repeat_file=repeat,
                tsd="UNK",
                target="Chr1",
                sample="HEG4",
                outdir=Path(workdir),
            )
            rows = [line.split("\t") for line in out.read_text().splitlines()]
            self.assertEqual(len(rows), 2)
            self.assertEqual(rows[0][1], "supporting_junction")
            self.assertEqual(rows[0][4], "1000..1002")
            self.assertEqual(rows[0][7:12], ["R:1", "L:0", "ST:1", "SR:0", "SL:1"])
            self.assertEqual(rows[1][1], "supporting_junction")
            self.assertEqual(rows[1][4], "3001..3003")
            self.assertEqual(rows[1][7:12], ["R:0", "L:1", "ST:1", "SR:1", "SL:0"])

    def test_paired_nonjunction_record_is_not_supporting_evidence(self):
        """RelocaTE2 excludes paired, untagged BAM records from support counts."""
        with tempfile.TemporaryDirectory() as workdir:
            bam_path = os.path.join(workdir, "flank.bam")
            reads = [
                {"name": "junction:start:5", "seq": "A" * 40, "start0": 999},
                {"name": "paired-mate", "seq": "A" * 40, "start0": 900, "flag": 1},
            ]
            _write_junction_bam(bam_path, "Chr1", 5000, reads)
            repeat = Path(workdir) / "read_repeat_name.txt"
            repeat.write_text("junction\tmping\t+\n")

            out = InsertionFinder(mismatch_allow=0, min_mapq=1).find_insertions(
                bam_file=Path(bam_path),
                read_repeat_file=repeat,
                tsd="UNK",
                target="Chr1",
                sample="HEG4",
                outdir=Path(workdir),
            )
            self.assertEqual(out.read_text(), "")


if __name__ == "__main__":
    unittest.main()
