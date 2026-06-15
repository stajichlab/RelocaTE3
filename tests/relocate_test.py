"""Run Relocate tests."""

import os
import unittest
from pathlib import Path

from RelocaTE3.align import Aligner
from RelocaTE3.librelocate import RelocaTE
from RelocaTE3.ReadLibrary import ReadLibrary


class TestRelocaTE(unittest.TestCase):
    """Test the RelocaTE class."""

    def test_trim_TE_reads(self):
        """Map reads to the TE library then trim, checking coordinates are produced."""
        rl = ReadLibrary(
            [
                os.path.join(
                    Path(__file__).parent,
                    "data",
                    "sim_reads",
                    "MSU7.Chr3_2M.ALL_reads_6X_100_500_1.fq.gz",
                ),
                os.path.join(
                    Path(__file__).parent,
                    "data",
                    "sim_reads",
                    "MSU7.Chr3_2M.ALL_reads_6X_100_500_2.fq.gz",
                ),
            ],
            "HEG4",
        )
        TElib = os.path.join(Path(__file__).parent, "data", "mping.fa")

        outdir = os.path.join(Path(__file__).parent, "results")
        if not os.path.exists(outdir):
            os.mkdir(outdir)

        aligner = Aligner(threads=2)
        bamfiles = aligner.map_minimap_library(rl, outdir, TElib)
        self.assertEqual(len(bamfiles), 2)

        relocate = RelocaTE(TElib=TElib, threads=2)
        trimlib = relocate.trim_TE_reads(rl, bamfiles)
        self.assertEqual(trimlib.name, "HEG4")
        # reads matching the mping TE library should produce trimmed coordinates
        self.assertGreater(len(trimlib.trimmed_coordinates), 0)
