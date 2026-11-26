"""Checks for selecting exact opacity table samples."""

import unittest

from giantlayers.exomol_import import select


class ExomolAxisTests(unittest.TestCase):
    def test_pressure_points_use_cgs_and_keep_source_order(self):
        self.assertEqual(select((1e-5, 2e-5, 1e-4), None, 1e6, "Pressure", (100, 10)),
                         ([0, 2], [10, 100]))

    def test_unavailable_points_are_rejected(self):
        with self.assertRaisesRegex(ValueError, "does not match"):
            select((1800, 2000, 2200), None, 1, "Temperature", (1900, 2200))
