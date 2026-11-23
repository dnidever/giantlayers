"""Analytic checks for the spherical formal transfer and sampled line kernels."""

import math
import unittest

from lte_opacity import C_ANGSTROM_S, line_opacity, sample_line_opacity
from spherical_transfer import disk_rays, surface_fluxes


class TransferChecks(unittest.TestCase):
    def test_transparent_shell_preserves_core_luminosity(self):
        radii = (3.0, 2.8, 2.5, 2.0)
        rays = disk_rays(radii, core_rays=12)
        sphere, _ = surface_fluxes((0, 0, 0, 1), (0, 0, 0), rays)
        self.assertAlmostEqual(sum(ray[-1] for ray in rays), math.pi)
        self.assertAlmostEqual(sphere * radii[0]**2, math.pi * radii[-1]**2)

    def test_thin_shell_approaches_plane_parallel_limit(self):
        rays = disk_rays((1 + 2e-10, 1 + 1e-10, 1), core_rays=12)
        sphere, plane = surface_fluxes((1, 1, 1), (10, 10), rays)
        self.assertTrue(math.isclose(sphere, plane, rel_tol=1e-8))

    def test_indexed_line_sampling_matches_direct_sum(self):
        wavelengths = [4000 + i * 0.01 for i in range(101)]
        transitions = (
            (C_ANGSTROM_S / 4000.2, 2e9, 0.0, 10),
            (C_ANGSTROM_S / 4000.4, 1e9, 0.3, 15),
        )
        expected = [line_opacity(C_ANGSTROM_S / w, transitions) for w in wavelengths]
        actual = sample_line_opacity(wavelengths, transitions)
        for sampled, direct in zip(actual, expected):
            self.assertTrue(math.isclose(sampled, direct, rel_tol=1e-14))


if __name__ == "__main__":
    unittest.main()
