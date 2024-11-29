# import pytest
# import numpy as np
# from src.chebyshev_interpolator import ChebyshevInterpolator

# @pytest.mark.parametrize("n, tau_max, a, b, expected_first_node, expected_last_node", [
#     (9, 1.0, 0.0, 1.0, 0.0, 1.0),  # Default interval [0, 1]
#     (5, 2.0, -1.0, 1.0, -1.0, 1.0),  # Interval [-1, 1]
#     (7, 0.5, 2.0, 4.0, 2.0, 4.0),  # Custom interval [2, 4]
# ])
# def test_chebyshev_interpolator_nodes(n, tau_max, a, b, expected_first_node, expected_last_node):
#     """
#     Test the ChebyshevInterpolator class to ensure correct node computation.
#     This test checks:
#     - Correct computation of Chebyshev interpolation nodes.
#     - Correct calculation of collocation times.
#     """
#     interpolator = ChebyshevInterpolator(n, tau_max)
#     interpolator.compute_nodes()

#     x_nodes, tau_nodes = interpolator.get_nodes()

#     # Check the first and last nodes for correctness
#     assert x_nodes[0] == pytest.approx(expected_first_node, rel=1e-6)
#     assert tau_nodes[0] == pytest.approx(expected_first_node ** 2, rel=1e-6)
#     assert x_nodes[-1] == pytest.approx(expected_last_node, rel=1e-6)
#     assert tau_nodes[-1] == pytest.approx(expected_last_node ** 2, rel=1e-6)

# if __name__ == "__main__":
#     pytest.main(["-v", "src/tests/test_chebyshev_interpolator.py"])
