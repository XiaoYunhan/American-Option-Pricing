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

    # Check the first value of the boundary initialization
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

    # Check that the refined boundary values are non-negative and decreasing
    assert all(refined_boundary >= 0)
    assert all(refined_boundary[i] >= refined_boundary[i + 1] for i in range(len(refined_boundary) - 1))

if __name__ == "__main__":
    pytest.main(["-v", "src/tests/test_dq_plus.py"])
