import pytest
import numpy as np
from src.dq_plus import DQPlus

@pytest.mark.parametrize("K, r, q, vol, tau_nodes, expected_first_value", [
    (100, 0.05, 0.02, 0.2, np.linspace(0, 1, 10), 100 * min(1, 0.05 / 0.02)),
    (120, 0.03, 0.03, 0.25, np.linspace(0, 2, 20), 120 * min(1, 0.03 / 0.03)),
    (80, 0.07, 0.01, 0.3, np.linspace(0, 1.5, 15), 80 * min(1, 0.07 / 0.01)),
])
def test_initialize_boundary(K, r, q, vol, tau_nodes, expected_first_value):
    """
    Test the initialization of the exercise boundary.
    """
    dqplus = DQPlus(K, r, q, vol, tau_nodes)
    dqplus.initialize_boundary()

    initial_boundary = dqplus.get_boundary_values()
    assert initial_boundary[0] == pytest.approx(expected_first_value, rel=1e-6)

@pytest.mark.parametrize("K, r, q, vol, tau_nodes", [
    (100, 0.05, 0.02, 0.2, np.linspace(0, 1, 10)),
    (120, 0.03, 0.03, 0.25, np.linspace(0, 2, 20)),
    (80, 0.07, 0.01, 0.3, np.linspace(0, 1.5, 15)),
])
def test_fixed_point_iteration_convergence(K, r, q, vol, tau_nodes):
    """
    Test the fixed-point iteration for convergence.
    """
    dqplus = DQPlus(K, r, q, vol, tau_nodes)
    dqplus.initialize_boundary()
    dqplus.fixed_point_iteration()

    refined_boundary = dqplus.get_boundary_values()
    assert all(refined_boundary >= 0)
    assert all(refined_boundary[i] >= refined_boundary[i + 1] for i in range(len(refined_boundary) - 1))

@pytest.mark.parametrize("K, r, q, vol, tau_nodes", [
    (100, 0.05, 0.02, 0.2, np.linspace(0, 1, 10)),
    (120, 0.03, 0.03, 0.25, np.linspace(0, 2, 20)),
    (80, 0.07, 0.01, 0.3, np.linspace(0, 1.5, 15)),
])
def test_compute_H(K, r, q, vol, tau_nodes):
    """
    Test the compute_H() method to ensure H(sqrt(tau)) values are correctly calculated.
    """
    dqplus = DQPlus(K, r, q, vol, tau_nodes)
    dqplus.initialize_boundary()
    H_values = dqplus.compute_H()

    # Check if H_values are non-negative
    assert all(H >= 0 for H in H_values)

@pytest.mark.parametrize("K, r, q, vol, tau_nodes", [
    (100, 0.05, 0.02, 0.2, np.linspace(0, 1, 10)),
    (120, 0.03, 0.03, 0.25, np.linspace(0, 2, 20)),
    (80, 0.07, 0.01, 0.3, np.linspace(0, 1.5, 15)),
])
def test_initialize_chebyshev_interpolation(K, r, q, vol, tau_nodes):
    """
    Test the initialize_chebyshev_interpolation() method for correct Chebyshev coefficients.
    """
    dqplus = DQPlus(K, r, q, vol, tau_nodes)
    dqplus.initialize_boundary()
    H_values = dqplus.compute_H()
    a_coefficients = dqplus.initialize_chebyshev_interpolation(H_values)

    # Check if the Chebyshev coefficients are finite numbers
    assert all(np.isfinite(a) for a in a_coefficients)

    # Ensure the coefficients are not all zero (a basic check for validity)
    assert np.any(a_coefficients != 0)

if __name__ == "__main__":
    pytest.main(["-v"])
