"""Regression checks for the bounded deepest-layer luminosity correction."""

import unittest

from spherical_heating import bracket_bottom


class BottomFluxBracketTests(unittest.TestCase):
    def test_stiff_flux_root_is_found_inside_the_bounds(self):
        visited = []

        def transfer(x):
            visited.append(x)
            return None, None, (0.12, 140 * (x - 0.0021) + 3000 * (x - 0.0021)**2)

        value, result = bracket_bottom(transfer, 0, transfer(0), -0.01, 0.01, 0.001)
        self.assertTrue(all(-0.01 <= x <= 0.01 for x in visited))
        self.assertAlmostEqual(value, 0.0021, delta=0.00001)
        self.assertLess(abs(result[2][-1]), 0.001)
        self.assertEqual(result[2][0], 0.12)

    def test_unbracketed_worse_endpoints_preserve_original_trial(self):
        def transfer(x):
            return None, None, (0.02, 0.1 + x*x)

        original = transfer(0)
        self.assertEqual(bracket_bottom(transfer, 0, original, -0.01, 0.01, 0.001), (0, original))
