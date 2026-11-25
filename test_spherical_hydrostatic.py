"""Analytic boundary and refinement checks for spherical hydrostatic balance."""

import unittest
from types import SimpleNamespace

from spherical_hydrostatic import G, solve


class ConstantDensity:
    def solve(self, temperature, pressure, abundances, initial=None):
        return SimpleNamespace(mass_density=1e-8, electron_density=0.0)


class HydrostaticChecks(unittest.TestCase):
    def test_radiative_acceleration_reduces_gas_pressure_gradient(self):
        radius = 1e9
        gm = 1e20
        mass = (0, 0.5, 1.0)
        arguments = (mass, (3500,)*3, {}, gm/G, gm/radius**2, 5.0, 1, ConstantDensity())
        gravity_only = solve(*arguments)
        supported = solve(*arguments, radiative_acceleration=(20,)*3)
        for m, p, q in zip(mass, gravity_only.pressure, supported.pressure):
            self.assertAlmostEqual(p-q, 20*m)
        with self.assertRaisesRegex(RuntimeError, "acceleration reaches gravity"):
            solve(*arguments, radiative_acceleration=(100,)*3)

    def test_fractional_photosphere_anchor(self):
        radius = 1e9
        gm = 1e20
        profile = solve((0, 0.5, 1.0), (3500,)*3, {}, gm/G, gm/radius**2,
                        5.0, 1.25, ConstantDensity())
        self.assertAlmostEqual(0.75*profile.radius[1] + 0.25*profile.radius[2], radius)

    def test_constant_density_sphere(self):
        reference_radius = 1e9
        gm = 1e20
        errors = []
        for n in (21, 41):
            mass = [i / (n - 1) for i in range(n)]
            profile = solve(mass, (3500,) * n, {}, gm / G, gm / reference_radius**2,
                            5.0, n // 2, ConstantDensity())
            self.assertEqual(profile.pressure[0], 5.0)
            self.assertEqual(profile.radius[n // 2], reference_radius)
            expected_radius = [reference_radius + (0.5 - m) / 1e-8 for m in mass]
            self.assertLess(max(abs(a - b) for a, b in zip(profile.radius, expected_radius)), 1e-6)
            expected_pressure = [5 + gm * 1e-8 * (1 / r - 1 / expected_radius[0])
                                 for r in expected_radius]
            errors.append(max(abs(a - b) for a, b in zip(profile.pressure, expected_pressure)))
        self.assertLess(errors[1], errors[0] / 3.9)


if __name__ == "__main__":
    unittest.main()
