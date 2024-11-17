# File: src/tests/test_chebyshev_interpolator.py

import pytest
import numpy as np
from src.chebyshev_interpolator import ChebyshevInterpolator

@pytest.mark.parametrize("n, tau_max, expected_first_node, expected_last_node", [
    (9, 1.0, 0.0, np.sqrt(2.0)),  # Expected sqrt(2.0) based on the formula
    (5, 2.0, 0.0, np.sqrt(4.0)),  # Expected sqrt(4.0) based on the formula
    (7, 0.5, 0.0, np.sqrt(1.0)),  # Expected sqrt(1.0) based on the formula
])
def test_chebyshev_interpolator_nodes(n, tau_max, expected_first_node, expected_last_node):
    """
    Test the ChebyshevInterpolator class to ensure correct node computation.
    This test checks:
    - Correct computation of Chebyshev interpolation nodes.
    - Correct calculation of collocation times.
    """
    interpolator = ChebyshevInterpolator(n, tau_max)
    interpolator.compute_nodes()

    x_nodes, tau_nodes = interpolator.get_nodes()

    # Check the first and last nodes for correctness
    assert x_nodes[0] == pytest.approx(expected_first_node, rel=1e-6)
    assert tau_nodes[0] == pytest.approx(expected_first_node ** 2, rel=1e-6)
    assert x_nodes[-1] == pytest.approx(expected_last_node, rel=1e-6)
    assert tau_nodes[-1] == pytest.approx(expected_last_node ** 2, rel=1e-6)

if __name__ == "__main__":
    pytest.main(["-v", "src/tests/test_chebyshev_interpolator.py"])

