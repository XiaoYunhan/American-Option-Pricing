import numpy as np

class QuadratureNodes:
    """
    This class handles the computation of Gauss-Legendre quadrature nodes
    and weights for numerical integration.
    
    Attributes:
    - l (int): Number of quadrature points.
    - y_nodes (numpy array): Computed quadrature nodes.
    - w_weights (numpy array): Computed quadrature weights.
    """

    def __init__(self, l):
        """
        Initialize the QuadratureNodes with the number of points (l).
        
        Parameters:
        - l (int): Number of quadrature points.
        """
        self.l = l
        self.y_nodes = None
        self.w_weights = None

    def compute_legendre_nodes(self):
        """
        Compute the Gauss-Legendre quadrature nodes and weights.
        This method uses the roots of the Legendre polynomial.
        """
        # Compute nodes and weights using numpy's leggauss function
        self.y_nodes, self.w_weights = np.polynomial.legendre.leggauss(self.l)

    def get_nodes_and_weights(self):
        """
        Retrieve the computed quadrature nodes and weights.
        
        Returns:
        - y_nodes (numpy array): Quadrature nodes.
        - w_weights (numpy array): Quadrature weights.
        """
        if self.y_nodes is None or self.w_weights is None:
            raise ValueError("Nodes and weights have not been computed yet. Call compute_legendre_nodes() first.")
        
        return self.y_nodes, self.w_weights
    
    def compute_tanh_sinh_nodes(self):
        """
        Compute the Tanh_sinh quadrature nodes and weights.
        """
        n = self.l
        nodes = np.zeros(n)
        weights = np.zeros(n)

        # Compute the nodes and weights
        for k in range(n):
            # Compute the k-th node
            nodes[k] = np.tanh(np.pi * (k + 0.5) / n)

            # Compute the weight for the k-th node
            weights[k] = (np.pi / n) / (np.cosh(np.pi * (k + 0.5) / n) ** 2)

        self.y_nodes, self.w_weights = nodes, weights
