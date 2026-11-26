"""Regression checks for molecular cross-section diagnostics."""

import json
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest

from giantlayers.molecular_cross_section import MolecularCrossSection


class MolecularCrossSectionTests(unittest.TestCase):
    def test_hot_omission_only_above_the_upper_table_boundary(self):
        data = {
            "species": "TiO",
            "temperature_K": [1800, 3400],
            "pressure_dyn_cm2": [10, 100],
            "wavelength_vacuum_A": [5000, 6000],
            "cross_section_cm2": [[[1, 1], [1, 1]], [[2, 2], [2, 2]]],
        }
        with TemporaryDirectory() as directory:
            path = Path(directory) / "tio.json"
            path.write_text(json.dumps(data), encoding="utf-8")
            baseline = MolecularCrossSection(path)
            diagnostic = MolecularCrossSection(path, omit_hot=True)
        grid = baseline.on_grid([5500, 6500])
        for temperature, expected in [(1700, 1), (3400, 2), (3401, 0)]:
            with self.subTest(temperature=temperature):
                gas = SimpleNamespace(temperature=temperature, mass_density=1,
                                      number_densities={"TiO": 1})
                self.assertEqual(diagnostic.mass_opacity(gas, 10, grid), [expected, 0])
        self.assertEqual(baseline.mass_opacity(gas, 10, grid), [2, 0])
