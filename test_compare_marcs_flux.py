"""Checks for matching and integrating native MARCS-grid samples."""

from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

from compare_marcs_flux import sampled_bands


class SampledBandTests(unittest.TestCase):
    def test_bounding_samples_integrate_and_grid_must_match(self):
        with TemporaryDirectory() as directory:
            path = Path(directory) / "native.csv"
            path.write_text(
                "wavelength_vacuum_A,reference_flux_per_A\n1000,3\n2000,5\n3000,7\n",
                encoding="utf-8",
            )
            self.assertEqual(sampled_bands(path, [1000., 2000., 3000.], [(150, 250)]),
                             [(150, 250, 5000.)])
            with self.assertRaisesRegex(ValueError, "complete, ordered MARCS"):
                sampled_bands(path, [1000., 3000.], [(150, 250)])
