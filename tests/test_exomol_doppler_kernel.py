"""Verify LTE line-strength normalization and Doppler bin fractions."""

import math
from pathlib import Path
import shutil
import struct
import subprocess
from tempfile import TemporaryDirectory
import unittest


@unittest.skipUnless(shutil.which("cc"), "a C compiler is required for the offline opacity test")
class DopplerKernelTests(unittest.TestCase):
    def test_single_line_at_shared_bin_edge(self):
        source = Path(__file__).resolve().parents[1] / "tools" / "exomol_doppler.c"
        with TemporaryDirectory() as directory:
            root = Path(directory)
            binary = root / "exomol_doppler"
            subprocess.run(["cc", "-std=c99", str(source), "-o", str(binary), "-lm"], check=True)
            config = root / "config"
            config.write_text("2 1 2 28 200000\n1999 2000 2001\n3500 1200\n1 200 3\n2 2200 5\n")
            result = subprocess.run([binary, config], input=b"2 1 10.0\n",
                                    capture_output=True, check=True)
        bins = struct.unpack("=2d", result.stdout)
        c2_over_t = 1.438776877 / 3500
        strength = (
            10 * 5 * math.exp(-c2_over_t*200) * (1-math.exp(-c2_over_t*2000))
            / (8*math.pi*2.99792458e10 * 2000**2 * 1200)
        )
        for cross_section in bins:
            self.assertTrue(math.isclose(cross_section, strength/2, rel_tol=1e-12))
