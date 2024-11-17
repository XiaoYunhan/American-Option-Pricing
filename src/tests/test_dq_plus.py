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
    assert all(np.isfinite(a) for a in a_coefficients)
    assert np.any(a_coefficients != 0)

@pytest.mark.parametrize("coefficients, z, expected_value", [
    ([1, 0.5, 0.25], 0.5, 1.125),
    ([2, -1, 0.5], -0.5, 2.25),
    ([0, 1, 1], 0, -1.0),
])
def test_clenshaw_algorithm(coefficients, z, expected_value):
    """
    Test the clenshaw_algorithm() method to ensure correct evaluation of Chebyshev polynomials.
    """
    dqplus = DQPlus(0, 0, 0, 0, np.array([]))
    dqplus.chebyshev_coefficients = np.array(coefficients)
    result = dqplus.clenshaw_algorithm(z)
    assert result == pytest.approx(expected_value, rel=1e-6)

@pytest.mark.parametrize("K, r, q, vol, tau_nodes, y_nodes", [
    (100, 0.05, 0.02, 0.2, np.linspace(0, 1, 10), np.linspace(-1, 1, 5)),
    (120, 0.03, 0.03, 0.25, np.linspace(0, 2, 20), np.linspace(-1, 1, 7)),
    (80, 0.07, 0.01, 0.3, np.linspace(0, 1.5, 15), np.linspace(-1, 1, 4)),
    (100, 0.05, 0.02, 0.2, np.linspace(0, 10, 10), np.linspace(-1, 1, 5)),
])
def test_evaluate_boundary(K, r, q, vol, tau_nodes, y_nodes):
    """
    Test the evaluate_boundary() method for correct boundary evaluation using Clenshaw's algorithm.
    """
    dqplus = DQPlus(K, r, q, vol, tau_nodes)
    dqplus.initialize_boundary()
    H_values = dqplus.compute_H()
    dqplus.initialize_chebyshev_interpolation(H_values)
    tau_test = tau_nodes[len(tau_nodes) // 2]
    B_values = dqplus.evaluate_boundary(tau_test, y_nodes)
    assert all(np.isfinite(B) for B in B_values)
    assert all(B >= 0 for B in B_values)

if __name__ == "__main__":
    pytest.main(["-v"])

