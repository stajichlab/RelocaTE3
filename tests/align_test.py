"""Run aligner tests."""

import unittest
from RelocaTE3.ReadLibrary import ReadLibrary
from RelocaTE3.align import Aligner
from pathlib import Path
import os
# import tempfile


class TestAligner(unittest.TestCase):
    """Test the Alignment class"""
    def test_minimap(self):

        rl = ReadLibrary([
            os.path.join(Path(__file__).parent, "data", "sim_reads", "MSU7.Chr3_2M.ALL_reads_6X_100_500_1.fq.gz"),
            os.path.join(Path(__file__).parent, "data", "sim_reads", "MSU7.Chr3_2M.ALL_reads_6X_100_500_2.fq.gz"),
            ], 'HEG4')
        self.assertEqual("HEG4", rl.name)
        mm = Aligner()
        self.assertEqual("minimap2", mm.minimap)

        TElib = os.path.join(Path(__file__).parent, "data", "mping.fa")
        self.assertTrue(os.path.exists(TElib))
        # with tempfile.TemporaryDirectory() as outdir:
        for outdir in [os.path.join(Path(__file__).parent, 'results')]:
            if not os.path.exists(outdir):
                os.mkdir(outdir)
            self.assertTrue(mm.index_minimap(TElib, TElib+".mmi", True) >= 0)
            # now run with debug flags as On to see STDERR
            # mm.verbose = True
            # self.assertTrue(mm._index_minimap(TElib, TElib+".mmi", True) >= 0)
            mm.verbose = False
            bamfiles = mm.map_minimap_library(rl, outdir, TElib)
            self.assertEqual(len(bamfiles), 2)

    def test_relocate_TE_read(self):
        rl = ReadLibrary([
            os.path.join(Path(__file__).parent, "data", "sim_reads", "MSU7.Chr3_2M.ALL_reads_6X_100_500_1.fq.gz"),
            os.path.join(Path(__file__).parent, "data", "sim_reads", "MSU7.Chr3_2M.ALL_reads_6X_100_500_2.fq.gz"),
            ], 'HEG4')
        mm = Aligner()

        TElib = os.path.join(Path(__file__).parent, "data", "mping.fa")
        # with tempfile.TemporaryDirectory() as outdir:
        for outdir in [os.path.join(Path(__file__).parent, 'results')]:
            if not os.path.exists(outdir):
                os.mkdir(outdir)

            bamfiles = mm.map_minimap_library(rl, outdir, TElib)
            self.assertEquals(2, len(bamfiles))


if __name__ == '__main__':
    unittest.main()
