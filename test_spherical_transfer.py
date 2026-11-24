"""Analytic checks for the spherical formal transfer and sampled line kernels."""

import math
import unittest
from unittest.mock import patch

from lte_opacity import C_ANGSTROM_S, line_opacity, sample_line_opacity
from marcs_writer import MarcsDocument
from spherical_heating import (FormalTransfer, angular_geometry, correct_temperatures,
                               shell_energy_budget, solve_wavelength)
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

    def test_tangent_ray_diagonal_matches_source_perturbations(self):
        rays = angular_geometry((3.0, 2.8, 2.5, 2.0), core_rays=12)
        thermal = (1.0, 1.3, 2.0, 3.0)
        field = FormalTransfer(rays, (0.2, 0.5, 0.6), thermal)
        mean, _ = field.moments(thermal)
        for i, original in enumerate(mean):
            changed = list(thermal)
            changed[i] += 1
            perturbed, _ = field.moments(changed)
            self.assertAlmostEqual(perturbed[i] - original, field.diagonal[i], places=12)

    def test_pure_scattering_produces_no_thermal_heating(self):
        rays = angular_geometry((3.0, 2.8, 2.5, 2.0), core_rays=12)
        _, _, heating, iterations = solve_wavelength(
            rays, (0.2, 0.3, 0.5), (1.0,) * 4, (1.0,) * 4, (0.0,) * 4, (1.0,) * 4
        )
        self.assertEqual(heating, [0.0] * 4)
        self.assertGreater(iterations, 0)

    def test_shell_heating_accounts_for_spherical_area(self):
        radii = (3.0, 2.5, 2.0)
        flux = tuple(10 / (4 * math.pi * r**2) for r in radii)
        _, shells, mismatch = shell_energy_budget(radii, (1.0,) * 3, flux, (0.0,) * 3)
        self.assertAlmostEqual(sum(shell[1] for shell in shells), 0.0)
        self.assertAlmostEqual(mismatch, 0.0)

    def test_temperature_step_reduces_shell_luminosity_residual(self):
        header = ("", "", "1", "", "", "", "", "1")
        structure = tuple((0, 0, 0, t) for t in (3000, 3100, 3200))
        model = MarcsDocument(header, (), structure, ((0,) * 7,) * 3, ())
        targets = (3050, 3150, 3200)

        def linear_transfer(trial, *args):
            luminosity = [
                4 * math.pi * (1 + (row[3] - expected) * 1e-4)
                for row, expected in zip(trial.structure, targets)
            ]
            return None, None, luminosity

        with patch("spherical_heating.evaluate", side_effect=linear_transfer):
            corrected, _, history = correct_temperatures(model, (), (), {}, 2)
        self.assertLess(history[-1][0], 1e-12)
        for row, expected in zip(corrected.structure, targets):
            self.assertAlmostEqual(row[3], expected)
        self.assertEqual(corrected.structure[-1], model.structure[-1])

    def test_absolute_luminosity_uses_independent_photospheric_mode(self):
        header = ("", "3500", "1", "", "", "", "", "1")
        structure = ((-4, 0, 0, 3000), (0, 0, 0, 3500), (2, 0, 0, 5000))
        model = MarcsDocument(header, (), structure, ((0,) * 7,) * 3, ())
        desired = (3040, 3450, 5000)

        def linear_transfer(trial, *args):
            return None, None, [
                4 * math.pi * (1 + 1e-4 * (row[3] - target))
                for row, target in zip(trial.structure, desired)
            ]

        with patch("spherical_heating.evaluate", side_effect=linear_transfer):
            updated, _, history = correct_temperatures(model, (), (), {}, 1, luminosity_target=True)
        self.assertLess(history[-1][0], 1e-12)
        self.assertLess(abs(history[-1][1]), 1e-12)
        for row, expected in zip(updated.structure, desired):
            self.assertAlmostEqual(row[3], expected)

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
