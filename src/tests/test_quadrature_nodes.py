import pytest
import numpy as np
from src.quadrature_nodes import QuadratureNodes

@pytest.mark.parametrize("l, expected_first_node, expected_last_node", [
    (5, -0.90618, 0.90618),  # Basic test case with l = 5
    (3, -0.774597, 0.774597),  # Test case with fewer quadrature points
    (7, -0.949108, 0.949108),  # Test case with more quadrature points
])
def test_quadrature_nodes(l, expected_first_node, expected_last_node):
    """
    Test the QuadratureNodes class to ensure correct node and weight computation.
    This test checks:
    - Correct computation of quadrature nodes.
    - Correct calculation of quadrature weights.
    """
    quadrature = QuadratureNodes(l)
    quadrature.compute_legendre_nodes()

    y_nodes, w_weights = quadrature.get_nodes_and_weights()

    # Check the first and last nodes for correctness
    assert y_nodes[0] == pytest.approx(expected_first_node, rel=1e-5)
    assert y_nodes[-1] == pytest.approx(expected_last_node, rel=1e-5)

    # Check that the sum of weights is approximately 2 (integral over [-1, 1])
    assert np.sum(w_weights) == pytest.approx(2.0, rel=1e-5)

if __name__ == "__main__":
    pytest.main(["-v", "src/tests/test_quadrature_nodes.py"])

