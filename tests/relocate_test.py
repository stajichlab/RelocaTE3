"""Run Relocate tests."""

import unittest
from RelocaTE3.ReadLibrary import ReadLibrary
from RelocaTE3.librelocate import RelocaTE
from pathlib import Path
import os


class TestRelocaTE(unittest.TestCase):
    """Test the RelocaTE class."""

    def test_relocate(self):
        rl = ReadLibrary([
            os.path.join(Path(__file__).parent, "data", "sim_reads", "MSU7.Chr3_2M.ALL_reads_6X_100_500_1.fq.gz"),
            os.path.join(Path(__file__).parent, "data", "sim_reads", "MSU7.Chr3_2M.ALL_reads_6X_100_500_2.fq.gz"),
            ], 'HEG4')
        TElib = os.path.join(Path(__file__).parent, "data", "mping.fa")
        # with tempfile.TemporaryDirectory() as outdir:
        relocate = RelocaTE(threads=2)
        print(os.path.join(Path(__file__).parent))
        for outdir in [os.path.join(Path(__file__).parent, 'results')]:
            if not os.path.exists(outdir):
                os.mkdir(outdir)
            relocate.identify_TE_reads(rl, outdir, TElib, "minimap2")
