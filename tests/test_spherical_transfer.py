from importlib.resources import files
"""Analytic checks for the spherical formal transfer and sampled line kernels."""

from dataclasses import replace
import math
from pathlib import Path
import shutil
import subprocess
from tempfile import TemporaryDirectory
import unittest
from types import SimpleNamespace
from unittest.mock import patch

from giantlayers.lte_opacity import C_ANGSTROM_S, LineProfiles, line_opacity, sample_line_opacity
from giantlayers.marcs_writer import MarcsDocument
from giantlayers.spherical_heating import (FormalTransfer, angular_geometry, balance_all_layers, calibrate_heating, correct_temperatures,
                               merit_descent, shell_energy_budget, solve_wavelength, temperature_direction)
from giantlayers.spherical_kernel_loader import CompiledKernel
from giantlayers.spherical_transfer import disk_rays, surface_fluxes


class TransferChecks(unittest.TestCase):
    def test_heating_calibration_preserves_absorbed_to_emitted_power(self):
        sampled = (-.2, .1, -1.)
        full = (-.4, .3, -.95)
        scales, offsets = calibrate_heating(sampled, full)
        for a, b, scale, offset in zip(sampled, full, scales, offsets):
            self.assertAlmostEqual(scale*a+offset, b)
        self.assertAlmostEqual(scales[0], .75)
        self.assertEqual(scales[-1], 1.)

    def test_large_temperature_correction_reaches_step_limit(self):
        self.assertAlmostEqual(temperature_direction([[1.0]], [0.4])[0], -0.03)
        strong, weak = temperature_direction([[10.0, 0.0], [0.0, 0.001]], [-0.2, 0.2])
        self.assertAlmostEqual(strong, 0.02, places=6)
        self.assertAlmostEqual(weak, -0.03, places=6)

    def test_fallback_descends_the_line_search_merit(self):
        residual = [.2, -.3, .1]
        columns = [[1., .4, 0.], [0., 1., -2.]]
        weight = .8
        direction = merit_descent(columns, residual, 2, weight)
        shifted = [r+sum(d*c[i] for d, c in zip(direction, columns))
                   for i, r in enumerate(residual)]
        def merit(values):
            return math.sqrt(sum(x*x for x in values)/2)+weight*max(map(abs, values))
        self.assertLess(merit(shifted), merit(residual))
        self.assertLessEqual(max(map(abs, direction)), .01)

    def test_direct_scattering_matches_converged_iteration(self):
        rays = angular_geometry((10.0, 9.9, 9.8, 9.7), 12)
        inputs = (rays, (0.1, 0.1, 0.1), (1, 1, 1, 1), (1, 1.1, 1.3, 1.7),
                  (0.2, 0.2, 0.2, 0.2), (50, 50, 50, 50))
        direct = solve_wavelength(*inputs, max_iterations=0)
        iterative = solve_wavelength(*inputs, tolerance=1e-11, max_iterations=500)
        self.assertLess(iterative[-1], 500)
        for x, y in zip(direct[:3], iterative[:3]):
            for a, b in zip(x, y):
                self.assertAlmostEqual(a, b, places=8)

    @unittest.skipUnless(shutil.which("cc"), "C compiler unavailable")
    def test_compiled_direct_scattering_matches_python(self):
        with TemporaryDirectory() as directory:
            library = Path(directory)/"giantlayers-transfer.so"
            source = files("giantlayers").joinpath("spherical_kernel.c")
            subprocess.run(["cc", "-std=c99", "-O2", "-shared", "-fPIC", str(source),
                            "-o", str(library), "-lm"], check=True, capture_output=True)
            rays = angular_geometry((10.0, 9.95, 9.9, 9.85, 9.8), 12)
            inputs = (rays, (0.05,)*4, (1,)*5, (1, 1.1, 1.3, 1.6, 1.7),
                      (1e-5,)*5, (500,)*5)
            expected = solve_wavelength(*inputs, max_iterations=0)
            actual = CompiledKernel(library, rays)(*inputs, max_iterations=0)
            for a, b in zip(expected[:3], actual[:3]):
                for x, y in zip(a, b):
                    self.assertTrue(math.isclose(x, y, rel_tol=1e-10, abs_tol=1e-12))

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

        with patch("giantlayers.spherical_heating.evaluate", side_effect=linear_transfer):
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

        with patch("giantlayers.spherical_heating.evaluate", side_effect=linear_transfer):
            updated, _, history = correct_temperatures(model, (), (), {}, 1, luminosity_target=True)
        self.assertLess(history[-1][0], 1e-12)
        self.assertLess(abs(history[-1][1]), 1e-12)
        for row, expected in zip(updated.structure, desired):
            self.assertAlmostEqual(row[3], expected)

    def test_temperature_trials_reuse_the_original_hydrostatic_boundary(self):
        header = ("", "3500", "1", "", "", "", "", "1")
        structure = ((-4, 0, 0, 3000, 0, 1), (0, 0, 0, 3500, 0, 1))
        model = MarcsDocument(header, (), structure, ((0,) * 7,) * 2, ())
        boundary = SimpleNamespace(reference_radius=1)
        pressures = []

        def reconstruct(trial, *, boundary):
            self.assertIs(boundary, original_boundary)
            rows = tuple((*row[:5], row[3] / 100) for row in trial.structure)
            return replace(trial, structure=rows), None

        def linear_transfer(trial, *args):
            pressure = trial.structure[0][5]
            pressures.append(pressure)
            return None, None, [4 * math.pi * pressure / 30.5, 4 * math.pi]

        original_boundary = boundary
        with patch("giantlayers.spherical_heating.boundary_from_marcs", return_value=(boundary, None)) as grid:
            with patch("giantlayers.spherical_heating.reconstruct_for_transfer", side_effect=reconstruct):
                with patch("giantlayers.spherical_heating.evaluate", side_effect=linear_transfer):
                    updated, _, history = correct_temperatures(model, (), (), {}, 1, hydrostatic=True)
        grid.assert_called_once()
        self.assertGreater(len(set(pressures)), 2)
        self.assertAlmostEqual(updated.structure[0][3], 3050)
        self.assertAlmostEqual(updated.structure[0][5], 30.5)
        self.assertLess(history[-1][0], 1e-12)

    def test_native_convection_controls_deep_luminosity_residual(self):
        header = ("", "3500", "1", "", "", "", "", "1")
        rows = ((-4, 0, 0, 3000), (0, 0, 0, 3500), (2, 0, 0, 4300), (3, 0, 0, 4500))
        model = MarcsDocument(header, (), rows, ((0,) * 7,) * 4, ())
        convection_temperatures = []

        def transfer(trial, *args):
            outer, photo, _, deep = [row[3] for row in trial.structure]
            convection_temperatures.append(deep)
            surface = 1 + (photo - 3450) * 1e-4
            radiation = [4 * math.pi * x for x in
                         (surface, surface + (outer - 3040) * 1e-4, 1.0, 0.6)]
            convection = (0.0, 0.0, 0.0, 0.4 + (deep - 4600) * 1e-4)
            return None, None, radiation, None, None, None, convection

        with patch("giantlayers.spherical_heating.evaluate", side_effect=transfer):
            corrected, _, history = correct_temperatures(
                model, (), (), {}, 1, luminosity_target=True, native_convection=True, deep_layers=1
            )
        self.assertGreater(len(set(convection_temperatures)), 2)
        self.assertLess(max(map(abs, history[-1])), 1e-12)
        for row, expected in zip(corrected.structure, (3040, 3450, 4300, 4600)):
            self.assertAlmostEqual(row[3], expected)

    def test_global_balance_updates_radiative_and_convective_layers(self):
        header = ("", "3500", "1", "", "", "", "", "1")
        rows = tuple((i, 0, 0, 3000 + 100*i) for i in range(8))
        model = MarcsDocument(header, (), rows, ((0,) * 7,) * 8, ())
        convection = (0.0, 0.0, 0.001, 0.0, 0.0, 0.01, 0.1, 0.4)

        def transfer(trial, *args):
            radiative = [
                4*math.pi*(1 - conv + 1e-4*(row[3] - 1.01*initial[3]))
                for row, initial, conv in zip(trial.structure, rows, convection)
            ]
            return None, None, radiative, None, None, None, convection

        boundary = SimpleNamespace(reference_radius=1)
        progress = []
        with (
            patch("giantlayers.spherical_heating.boundary_from_marcs", return_value=(boundary, None)),
            patch("giantlayers.spherical_heating.reconstruct_for_transfer", side_effect=lambda profile, **_: (profile, None)),
            patch("giantlayers.spherical_heating.evaluate", side_effect=transfer),
        ):
            corrected, _, history, knots, start = balance_all_layers(
                model, (), (), {}, upper_knots=2, progress=lambda step, metrics: progress.append((step, metrics))
            )
            with patch("giantlayers.spherical_heating.temperature_direction", return_value=[-.01]*5):
                _, _, fallback_history, _, _ = balance_all_layers(model, (), (), {}, upper_knots=2)

        self.assertEqual(start, 5)
        self.assertEqual(knots, (0, 4, 5, 6, 7))
        self.assertEqual(len(history), 2)
        self.assertEqual(progress, list(enumerate(history)))
        self.assertLess(history[-1][1], 1e-8)
        self.assertEqual(len(fallback_history), 2)
        self.assertLess(fallback_history[-1][0], fallback_history[0][0])
        for row, initial in zip(corrected.structure, rows):
            self.assertAlmostEqual(row[3], 1.01*initial[3], places=5)

    def test_full_grid_heating_offset_changes_surface_temperature(self):
        header = ("", "3500", "1", "", "", "", "", "1")
        rows = tuple((i, 0, 0, 3000.) for i in range(4))
        model = MarcsDocument(header, (), rows, ((0,) * 7,) * 4, ())
        boundary = SimpleNamespace(reference_radius=1, column_mass=(0, 1, 2, 3))

        def transfer(trial, *args):
            heating = ((trial.structure[0][3]-3000)/100, 0., 0., 0.)
            return (None, None, (4*math.pi,)*4, None, None, None, (0.,)*4,
                    ((0., 0., 0., 0., 1.),)*4, heating)

        with (
            patch("giantlayers.spherical_heating.boundary_from_marcs", return_value=(boundary, None)),
            patch("giantlayers.spherical_heating.reconstruct_for_transfer", side_effect=lambda profile, **_: (profile, None)),
            patch("giantlayers.spherical_heating.evaluate", side_effect=transfer),
        ):
            corrected, _, _, _, _ = balance_all_layers(
                model, (), (), {}, outer_heating_weight=1, heating_offset=(.03, 0., 0., 0.)
            )
        self.assertAlmostEqual(corrected.structure[0][3], 2997., places=5)

    def test_indexed_line_sampling_matches_direct_sum(self):
        wavelengths = [4000 + i * 0.01 for i in range(101)]
        transitions = (
            (C_ANGSTROM_S / 4000.2, 2e9, 0.0, 10),
            (C_ANGSTROM_S / 4000.4, 1e9, 0.3, 15),
            (C_ANGSTROM_S / 1000, 2e9, 0.0, 12),
            (C_ANGSTROM_S / 11000, 1e9, 0.8, 16),
        )
        expected = [line_opacity(C_ANGSTROM_S / w, transitions) for w in wavelengths]
        actual = sample_line_opacity(wavelengths, LineProfiles(transitions))
        for sampled, direct in zip(actual, expected):
            self.assertTrue(math.isclose(sampled, direct, rel_tol=1e-14))


if __name__ == "__main__":
    unittest.main()
