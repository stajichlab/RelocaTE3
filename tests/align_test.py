"""Run aligner tests."""

import unittest
import RelocaTE3.ReadLibrary as ReadLibrary
from RelocaTE3.align import Aligner
from pathlib import Path
import os
import tempfile


class TestAligner(unittest.TestCase):
    """Test the Alignment class"""
    def test_build(self):
        rl = ReadLibrary.ReadLibrary(['A123_R1.fq.gz', 'A123_R2.fq.gz'], 'A123')
        self.assertEqual("A123", rl.name)
        mm = Aligner()
        self.assertEqual("minimap2", mm.minimap)

        TElib = os.path.join(Path(__file__).parent, "data", "mping.fa")
        self.assertTrue(os.path.exists(TElib))
        with tempfile.TemporaryDirectory() as outdir:
            print(outdir)
            self.assertTrue(mm._index_minimap(TElib, TElib+".mmi", True) >= 0)
            # now run with debug flags as On to see STDERR
            mm.verbose = True
            self.assertTrue(mm._index_minimap(TElib, TElib+".mmi", True) >= 0)
            mm.verbose = False

            mm._map_minimap_library(TElib+".mmi", rl, outdir)


if __name__ == '__main__':
    unittest.main()
