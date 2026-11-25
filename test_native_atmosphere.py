"""Independent initialization, equilibrium, and export checks."""

import math
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

from giantlayers import initialize
from lte_chemistry import Thermochemistry
from marcs_writer import read_document, write_document
from native_marcs import build_document
from spherical_hydrostatic import reconstruct_for_transfer


class NativeAtmosphereChecks(unittest.TestCase):
    def test_warm_equilibrium_matches_cold_solution(self):
        chemistry = Thermochemistry()
        abundance = {"H": 12, "He": 10.93, "C": 9.0, "O": 8.83, "Fe": 7.50}
        seed = chemistry.solve(3499, 250.1, abundance)
        warm = chemistry.solve(3500, 250, abundance, initial=seed)
        cold = chemistry.solve(3500, 250, abundance)
        for a, b in zip(warm.number_densities.values(), cold.number_densities.values()):
            self.assertTrue(math.isclose(a, b, rel_tol=2e-9))
        self.assertTrue(math.isclose(warm.mass_density, cold.mass_density, rel_tol=2e-9))

    def test_native_structure_exports_calculated_depth_and_abundances(self):
        seed, boundary = initialize(3500, 1, 1, {"H": 12, "He": 10.93, "C": 9.0}, layers=5)
        model, _ = reconstruct_for_transfer(seed, boundary=boundary)
        # Native transfer normally provides these Rosseland and convection fields.
        convection = [[0, 0, 0, 0, 0.1, 0, 0, 0, 0, 0, 0, 0, 0]] * 5
        document = build_document(model, (None,) * 7 + (convection,), (), (), boundary, 2)
        with TemporaryDirectory() as directory:
            path = Path(directory) / "test.mod.gz"
            write_document(document, path)
            reread = read_document(path)
        self.assertEqual(reread.abundances[5], 9.0)
        self.assertEqual(len(reread.structure), 5)
        self.assertEqual(len(reread.partial_pressures), 3)
        self.assertEqual(reread.header[4][:5], "  2.0")
        self.assertEqual(reread.header[6][:12], "-99.00-99.00")
        self.assertGreater(float(reread.header[7][:12]), 0)
        self.assertTrue(all(a[0] < b[0] for a, b in zip(reread.structure, reread.structure[1:])))
        self.assertTrue(all(row[5] > 0 for row in reread.structure))


if __name__ == "__main__":
    unittest.main()
