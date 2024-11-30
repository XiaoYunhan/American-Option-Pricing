import pytest
import numpy as np
from src.chebyshev_interpolator import ChebyshevInterpolator
import re

@pytest.mark.parametrize(
    "n, tau_max, expected_z_first, expected_z_last",
    [
        (5, 1.0, 1.0, -1.0),  # Small number of nodes
        (10, 2.0, 1.0, -1.0),  # Medium number of nodes
        (20, 0.5, 1.0, -1.0),  # Larger number of nodes
    ]
)
def test_compute_nodes(n, tau_max, expected_z_first, expected_z_last):
    """
    Test the computation of Chebyshev nodes and collocation times.
    """
    interpolator = ChebyshevInterpolator(n, tau_max)
    interpolator.compute_nodes()

    # Check the z_nodes values
    z_nodes = interpolator.z_nodes
    assert len(z_nodes) == n + 1, f"Expected {n + 1} z_nodes, but got {len(z_nodes)}."
    assert np.isclose(z_nodes[0], expected_z_first), "First z_node does not match expected value."
    assert np.isclose(z_nodes[-1], expected_z_last), "Last z_node does not match expected value."

    # Check the x_nodes values
    x_nodes = interpolator.x_nodes
    assert len(x_nodes) == n + 1, f"Expected {n + 1} x_nodes, but got {len(x_nodes)}."
    assert np.all(x_nodes >= 0), "x_nodes should all be non-negative."

    # Check the tau_nodes values
    tau_nodes = interpolator.tau_nodes
    assert len(tau_nodes) == n + 1, f"Expected {n + 1} tau_nodes, but got {len(tau_nodes)}."
    assert np.all(tau_nodes >= 0), "tau_nodes should all be non-negative."
    assert np.allclose(tau_nodes, x_nodes**2), "tau_nodes should be the square of x_nodes."

def test_get_nodes_before_computation():
    """
    Test that an error is raised if get_nodes() is called before compute_nodes().
    """
    interpolator = ChebyshevInterpolator(5, 1.0)
    with pytest.raises(ValueError, match=re.escape("Nodes have not been computed yet. Call compute_nodes() first.")):
        interpolator.get_nodes()

@pytest.mark.parametrize(
    "n, tau_max",
    [
        (5, 1.0),
        (10, 2.0),
        (20, 0.5),
    ]
)
def test_get_nodes_after_computation(n, tau_max):
    """
    Test that get_nodes() returns the correct x_nodes and tau_nodes after computation.
    """
    interpolator = ChebyshevInterpolator(n, tau_max)
    interpolator.compute_nodes()
    x_nodes, tau_nodes = interpolator.get_nodes()

    # Validate that the returned nodes match the computed ones
    assert np.array_equal(x_nodes, interpolator.x_nodes), "x_nodes from get_nodes() do not match computed x_nodes."
    assert np.array_equal(tau_nodes, interpolator.tau_nodes), "tau_nodes from get_nodes() do not match computed tau_nodes."

if __name__ == "__main__":
    pytest.main(["-v"])
