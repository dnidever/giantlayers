"""Check the residual-change metric used to refresh the Jacobian."""

import unittest

from giantlayers.spherical_heating import jacobian_prediction_error


class JacobianPredictionTests(unittest.TestCase):
    def test_predicts_a_linear_step(self):
        self.assertAlmostEqual(jacobian_prediction_error(
            [[2, 0], [0, 1]], [0.25, 0.5], [1, 1], [1.5, 1.5]
        ), 0)

    def test_detects_an_inaccurate_or_reversed_prediction(self):
        self.assertGreater(jacobian_prediction_error(
            [[2, 0], [0, 1]], [0.25, 0.5], [1, 1], [0.5, 0.5]
        ), 0.5)
