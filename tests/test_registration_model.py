import unittest
import numpy as np

from src.mip_record_conditioned_actualisation import (
    RegistrationParameters,
    bounded_actualisation_rate,
    gaussian_spatial_kernel,
    integrate_registration_load,
    stochastic_actualisation_step,
    trace_distance,
)


class RegistrationModelTests(unittest.TestCase):
    def test_trace_distance_extremes(self):
        rho0 = np.diag([1.0, 0.0]).astype(complex)
        rho1 = np.diag([0.0, 1.0]).astype(complex)
        self.assertAlmostEqual(trace_distance(rho0, rho0), 0.0)
        self.assertAlmostEqual(trace_distance(rho0, rho1), 1.0)

    def test_load_recovers_after_source_stops(self):
        times = np.linspace(0.0, 10.0, 101)
        source = np.zeros_like(times)
        source[:11] = 1.0
        load = integrate_registration_load(times, source, tau_r=1.0)
        self.assertGreater(load[10], load[-1])
        self.assertTrue(np.all(load >= 0.0))

    def test_rate_is_bounded(self):
        params = RegistrationParameters(1.0, 0.1, 0.9, 2.0, 2.0)
        rate = bounded_actualisation_rate(np.array([0.0, 2.0, 100.0]), params)
        self.assertAlmostEqual(rate[0], 0.1)
        self.assertTrue(np.all(np.diff(rate) > 0.0))
        self.assertLess(rate[-1], 1.0)

    def test_kernel_and_state_normalisation(self):
        kernel = gaussian_spatial_kernel(np.array([0.0, 1.0, 2.0]), xi=1.0)
        self.assertAlmostEqual(float(kernel.sum()), 1.0)
        psi = np.array([np.sqrt(0.3), np.sqrt(0.7)], dtype=complex)
        updated = stochastic_actualisation_step(
            psi,
            np.zeros((2, 2), dtype=complex),
            np.diag([1.0, -1.0]).astype(complex),
            gamma=0.2,
            dt=1e-3,
            d_wiener=0.01,
        )
        self.assertAlmostEqual(float(np.linalg.norm(updated)), 1.0, places=12)


if __name__ == "__main__":
    unittest.main()
