"""Check mixing-length stability and the Henyey algebraic closure."""

import math
import unittest

from giantlayers.gray_atmosphere import SIGMA
from giantlayers.mixing_length import henyey_flux, slopes


class MixingLengthChecks(unittest.TestCase):
    def test_gradient_recovers_quadratic_at_interior_and_boundaries(self):
        axis = [0, 1, 3, 5]
        gradient = slopes(axis, [x*x for x in axis], quadratic_endpoints=True)
        for derivative, x in zip(gradient, axis):
            self.assertAlmostEqual(derivative, 2*x)

    def test_stable_layer_carries_no_convective_flux(self):
        self.assertEqual(henyey_flux(3500, 100, 1e-8, 10, 0.1, 0.2, 2e8, 1, 0.4), (0, 0))

    def test_convective_velocity_and_flux_satisfy_henyey_equations(self):
        t, p, rho, g, kappa, gradient, cp, q, adiabatic = (
            7500, 6000, 1e-8, 10, 2, 0.8, 2e8, 1.2, 0.2
        )
        alpha, y, nu = 1.5, 0.076, 8
        flux, v = henyey_flux(t, p, rho, g, kappa, gradient, cp, q, adiabatic, alpha, y, nu)
        hp = p / (g*rho)
        tau = kappa*rho*alpha*hp
        gamma = v*rho*cp*(1+y*tau**2)/(8*SIGMA*t**3*tau)
        contrast = gamma/(1+gamma)*(gradient-adiabatic)
        self.assertTrue(math.isclose(v, alpha*math.sqrt(g*hp*q*contrast/nu), rel_tol=1e-14))
        self.assertTrue(math.isclose(flux, .5*rho*cp*t*v*alpha*contrast, rel_tol=1e-14))


if __name__ == "__main__":
    unittest.main()
