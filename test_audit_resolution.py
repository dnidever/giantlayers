"""Regression checks for comparisons between different radial grids."""

import unittest

from audit_resolution import audit


class ResolutionAuditTests(unittest.TestCase):
    def test_detects_an_interior_difference_even_when_surface_flux_agrees(self):
        coarse = [(0, 3000, 1, .001), (2, 4000, 100, .002)]
        fine = [(0, 3000, 1, .001), (1, 3650, 10, -.05), (2, 4000, 100, .002)]
        result = audit(coarse, fine, temperature_tolerance=100, pressure_tolerance=.01)
        self.assertEqual(result["max_temperature_difference_K"], 150)
        self.assertEqual(result["fine"]["worst_layer"], 2)
        self.assertEqual(result["max_gas_pressure_difference_dex"], 0)
        self.assertEqual(result["checks"], dict(same_mass_extent=True, both_in_luminosity_balance=False,
                                               temperature_change_within_tolerance=False,
                                               gas_pressure_change_within_tolerance=True))

    def test_rejects_disjoint_grids(self):
        with self.assertRaisesRegex(ValueError, "share a nonzero"):
            audit([(0, 3000, 1, 0), (1, 4000, 2, 0)], [(2, 4000, 3, 0), (3, 5000, 4, 0)])
