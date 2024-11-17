import numpy as np
from scipy.interpolate import BarycentricInterpolator
from src.chebyshev_interpolator import ChebyshevInterpolator

class DQPlus:
    """
    The core class for American option pricing using DQ+ method.
    This class handles initialization of exercise boundary, fixed-point iteration, and interpolation.
    """

    def __init__(self, K, r, q, vol, tau_nodes):
        """
        Initialize the DQPlus engine with option parameters and collocation nodes.

        Parameters:
        - K (float): Strike price.
        - r (float): Risk-free interest rate.
        - q (float): Dividend yield.
        - vol (float): Volatility of the underlying asset.
        - tau_nodes (numpy array): Collocation time points (tau_i).
        """
        self.K = K
        self.r = r
        self.q = q
        self.vol = vol
        self.tau_nodes = tau_nodes
        self.n = len(tau_nodes)
        self.initial_boundary = np.zeros(self.n)
        self.chebyshev_interpolator = None

    def initialize_boundary(self):
        """
        Initialize the early exercise boundary using QD+ approximation.
        """
        # Initialize B^(0)(tau_0)
        B_tau_0 = self.K * min(1, self.r / self.q) if self.q > 0 else self.K
        self.initial_boundary[0] = B_tau_0

        # Estimate B^(0)(tau_i) for i = 1, ..., n using exponential decay model
        for i in range(1, self.n):
            tau_i = self.tau_nodes[i]
            decay_factor = np.exp(-self.r * tau_i)
            self.initial_boundary[i] = B_tau_0 * decay_factor
    
    def compute_H(self):
        """
        Compute H(sqrt(tau)) for each collocation point based on the previous boundary values.

        Returns:
        - H_values (numpy array): The computed H(sqrt(tau)) values.
        """
        H_values = np.zeros(self.n)
        for i in range(self.n):
            tau = self.tau_nodes[i]
            B_tau = self.initial_boundary[i]
            H_values[i] = (np.log(B_tau / self.K))**2

        return H_values

    def fixed_point_iteration(self, max_iter=10, tol=1e-8):
        """
        Perform fixed-point iteration using Newton's method to refine the exercise boundary.

        Parameters:
        - max_iter (int): Maximum number of iterations.
        - tol (float): Convergence tolerance.
        """
        boundary_values = self.initial_boundary.copy()

        for iteration in range(max_iter):
            prev_boundary = boundary_values.copy()

            # Create a new interpolator for the current boundary values
            interpolator = BarycentricInterpolator(self.tau_nodes, boundary_values)

            for i in range(1, self.n):
                tau = self.tau_nodes[i]
                boundary_estimate = interpolator(tau)

                # Newton's update step
                boundary_values[i] = boundary_estimate - (boundary_estimate - prev_boundary[i]) / 2

            # Check for convergence
            max_diff = np.max(np.abs(boundary_values - prev_boundary))
            if max_diff < tol:
                break

        self.initial_boundary = boundary_values

    def get_boundary_values(self):
        """
        Retrieve the computed exercise boundary values.

        Returns:
        - initial_boundary (numpy array): Refined exercise boundary values.
        """
        return self.initial_boundary
    
    def initialize_chebyshev_interpolation(self, H_values):
        """
        Initialize the Chebyshev interpolation by computing the coefficients a_k.
        """
        self.chebyshev_interpolator = ChebyshevInterpolator(self.n - 1, np.max(self.tau_nodes))
        self.chebyshev_interpolator.compute_nodes()

        x_nodes, tau_nodes = self.chebyshev_interpolator.get_nodes()
        z_nodes = -np.cos(np.pi * np.arange(self.n) / (self.n - 1))

        a_coefficients = np.zeros(self.n)
        for k in range(self.n):
            sum_value = 0.0
            for i in range(1, self.n - 1):
                cos_term = np.cos(i * k * np.pi / self.n)
                sum_value += H_values[i] * cos_term

            a_coefficients[k] = (1 / (2 * self.n)) * (H_values[0] + (-1)**self.n * H_values[-1]) + (2 / self.n) * sum_value

        return a_coefficients
