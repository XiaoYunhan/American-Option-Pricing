import pytest
import numpy as np
from src.dq_plus import DQPlus

@pytest.mark.parametrize(
    "K, r, q, vol, tau_nodes, m_iterations",
    [
        (100, 0.05, 0.02, 0.2, np.linspace(0, 1, 10), 10),
        (120, 0.03, 0.03, 0.25, np.linspace(0, 2, 20), 15),
        (80, 0.07, 0.01, 0.3, np.linspace(0, 1.5, 15), 20),
    ]
)
def test_run_full_algorithm(K, r, q, vol, tau_nodes, m_iterations):
    """
    Test the full DQPlus algorithm to ensure the final boundary values are valid.
    """
    dqplus = DQPlus(K, r, q, vol, tau_nodes)

    # Initialize the boundary
    dqplus.initialize_boundary()

    # Run the full algorithm with specified iterations
    final_B_values = dqplus.run_full_algorithm(m_iterations)

    # Validate the final boundary values
    assert all(final_B_values >= 0), "Negative boundary value found in the final result."
    assert all(np.isfinite(final_B_values)), "Non-finite boundary value found in the final result."

    # Optionally, you could check the length of the final_B_values to match tau_nodes
    assert len(final_B_values) == len(tau_nodes), "Mismatch in length of final boundary values and tau nodes."

if __name__ == "__main__":
    pytest.main(["-v"])