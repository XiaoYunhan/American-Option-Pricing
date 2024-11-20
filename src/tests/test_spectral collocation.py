import pytest
import numpy as np
from src.dq_plus import DQPlus

@pytest.mark.parametrize(
    "K, r, q, vol, n, l, m_iterations",
    [
        (100, 0.05, 0.02, 0.2, 10, 5, 10),
        (120, 0.03, 0.03, 0.25, 20, 7, 15),
        (80, 0.07, 0.01, 0.3, 15, 4, 50),
    ]
)
def test_run_full_algorithm(K, r, q, vol, n, l, m_iterations):
    """
    Test the full DQPlus algorithm to ensure the final boundary values are valid.
    """
    dqplus = DQPlus(K, r, q, vol, n, l)

    # Initialize the boundary
    dqplus.initialize_boundary()

    # Initialize Chebyshev interpolation to ensure coefficients are available
    H_values = dqplus.compute_H()
    dqplus.initialize_chebyshev_interpolation(H_values)

    # Run the full algorithm with specified iterations
    final_B_values = dqplus.run_full_algorithm(m_iterations)

    # Validate the final boundary values
    assert all(final_B_values >= 0), "Negative boundary value found in the final result."
    assert all(np.isfinite(final_B_values)), "Non-finite boundary value found in the final result."

    # Validate the length of the final_B_values matches the number of collocation nodes (n)
    assert len(final_B_values) == n, "Mismatch in length of final boundary values and collocation points."

if __name__ == "__main__":
    pytest.main(["-v"])
