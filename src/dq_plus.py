import numpy as np
from scipy.interpolate import BarycentricInterpolator
from scipy.integrate import quad
from src.chebyshev_interpolator import ChebyshevInterpolator

class DQPlus:
    """
    The core class for American option pricing using DQ+ method.
    This class handles initialization of exercise boundary, fixed-point iteration, and interpolation.
    """

    def __init__(self, K, r, q, vol, tau_nodes, eta=0.5):
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
        self.eta = eta
        self.initial_boundary = np.zeros(self.n)
        self.chebyshev_interpolator = None
        self.chebyshev_coefficients = None

    def initialize_boundary(self):
        """
        Initialize the early exercise boundary using QD+ approximation.
        """
        # Initialize B^(0)(tau_0)
        B_tau_0 = self.K * min(1, self.r / self.q) if self.q > 0 else self.K
        self.initial_boundary[0] = B_tau_0

        # Estimate B^(0)(tau_i) for i = 1, ..., n using exponential decay model
        for i in range(1, self.n):
            decay_factor = np.exp(-self.r * self.tau_nodes[i])
            self.initial_boundary[i] = B_tau_0 * decay_factor
    
    def compute_H(self):
        """
        Compute H(sqrt(tau)) for each collocation point based on the previous boundary values.

        Returns:
        - H_values (numpy array): The computed H(sqrt(tau)) values.
        """
        H_values = np.zeros(self.n)
        for i in range(self.n):
            # tau = self.tau_nodes[i]
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
        self.chebyshev_coefficients = a_coefficients
        return a_coefficients
    
    def clenshaw_algorithm(self, z):
        """
        Evaluate the Chebyshev polynomial using Clenshaw's algorithm.

        Parameters:
        - z (float): The input value to evaluate the polynomial.

        Returns:
        - (float): The evaluated value of the Chebyshev polynomial.
        """
        if self.chebyshev_coefficients is None:
            raise ValueError("Chebyshev coefficients are not initialized.")
        
        n = len(self.chebyshev_coefficients)
        b_next = 0.0
        b_curr = 0.0

        for j in range(n - 1, 0, -1):
            b_temp = b_curr
            b_curr = 2 * z * b_curr - b_next + self.chebyshev_coefficients[j]
            b_next = b_temp

        return z * b_curr - b_next + self.chebyshev_coefficients[0]

    def evaluate_boundary(self, tau, y_nodes):
        """
        For each tau_i, evaluate B(tau) using Clenshaw algorithm at the adjusted points.

        Parameters:
        - tau (float): The current time node tau_i.
        - y_nodes (numpy array): The quadrature nodes y_k.

        Returns:
        - B_values (numpy array): Evaluated boundary values at adjusted points.
        """
        if self.chebyshev_coefficients is None:
            H_values = self.compute_H()
            self.initialize_chebyshev_interpolation(H_values)

        B_values = np.zeros(len(y_nodes))

        for k, y_k in enumerate(y_nodes):
            adjusted_tau = tau - tau * (1 + y_k)**2 / 4
            z = 2 * np.sqrt(adjusted_tau / np.max(self.tau_nodes)) - 1
            B_values[k] = max(self.clenshaw_algorithm(z), 1e-10)  # Ensure B_values are positive
            
        return B_values
    
    def compute_ND_values(self, tau, B_values):
        """
        Compute N(tau, B) and D(tau, B) using numerical quadrature.
        """
        vol_sqrt_tau = max(self.vol * np.sqrt(tau), 1e-8) # Increase minimum value protection to avoid division by zero

        def integrand_N(u, B):
            if B <= 1e-10:
                return 0.0
            exp_factor = np.exp((self.q - self.r) * u)
            m = (np.log(max(B / self.K, 1e-10)) + (self.r - self.q) * u) / vol_sqrt_tau + 0.5 * vol_sqrt_tau
            if not np.isfinite(m):
                return 0.0
            return exp_factor * np.exp(-0.5 * m**2) / (np.sqrt(2 * np.pi))

        def integrand_D(u, B):
            if B <= 1e-10:
                return 0.0
            exp_factor = np.exp(self.q * u)
            m = (np.log(max(B / self.K, 1e-10)) + (self.r - self.q) * u) / vol_sqrt_tau - 0.5 * vol_sqrt_tau
            if not np.isfinite(m):
                return 0.0
            return exp_factor * np.exp(-0.5 * m**2) / (np.sqrt(2 * np.pi))

        # Calculate N and D for each B value
        try:
            N_values = np.array([quad(integrand_N, 0, tau, args=(B,), limit=100)[0] for B in B_values])
            D_values = np.array([quad(integrand_D, 0, tau, args=(B,), limit=100)[0] for B in B_values])
        except Exception as e:
            print(f"Integration error: {e}")
            N_values = np.zeros(len(B_values))
            D_values = np.ones(len(B_values)) * 1e-10

        # Use np.clip to protect the lower bound of D_values to avoid division by zero
        D_values = np.clip(D_values, 1e-10, None)

        return N_values, D_values
    
    def compute_f_values(self, tau, B_values):
        """
        Compute f(tau, B) values based on N and D.
        """
        N_values, D_values = self.compute_ND_values(tau, B_values)
        f_values = N_values / D_values
        return f_values
    
    def compute_f_derivative(self, tau, B_values, h=1e-5):
        """
        Compute f'(tau, B) using numerical differentiation (central difference method).

        Parameters:
        - tau (float): The current time node tau_i.
        - B_values (numpy array): Evaluated boundary values at adjusted points.
        - h (float): Step size for numerical differentiation.

        Returns:
        - f_derivative (numpy array): Computed derivative values of f with respect to B.
        """
        f_values = self.compute_f_values(tau, B_values)
        f_derivative = np.zeros(len(B_values))

        for i, B in enumerate(B_values):
            B_forward = B + h
            B_backward = B - h

            # Compute f values at B + h and B - h
            f_forward = self.compute_f_values(tau, np.array([B_forward]))[0]
            f_backward = self.compute_f_values(tau, np.array([B_backward]))[0]

            # Central difference approximation
            f_derivative[i] = (f_forward - f_backward) / (2 * h)

        return f_derivative
    
    def update_boundary(self, B_values, max_iter=5):
        """
        Update the boundary values using the Jacobi-Newton scheme.

        Parameters:
        - B_values (numpy array): Current boundary values B^{(j)}(\tau_i).

        Returns:
        - B_next (numpy array): Updated boundary values B^{(j+1)}(\tau_i).
        """
        B_next = np.zeros(self.n)

        for i in range(self.n):
            tau = self.tau_nodes[i]
            B_current = B_values[i]

            for iter_count in range(max_iter):

                # Calculate f(tau, B) and f'(tau, B)
                f_value = self.compute_f_values(tau, np.array([B_current]))[0]
                f_derivative = self.compute_f_derivative(tau, np.array([B_current]))[0]

                # Avoid division by zero
                denominator = np.clip(f_derivative - 1.0, 1e-8, None)

                # Jacobi-Newton update formula
                numerator = B_current - f_value
                delta_B = self.eta * (numerator / denominator)

                # Update step limit
                max_step = 0.1 * B_current
                delta_B = np.clip(delta_B, -max_step, max_step)

                # Update boundary value and check non-negativity
                B_updated = B_current + delta_B
                if B_updated >= 0:
                    B_next[i] = B_updated
                    break
                else:
                    # If the result is negative, decrease the step size \(\eta\)
                    self.eta *= 0.5
                    B_current = max(B_current, 1e-10)
                    
            # Ensure the boundary value is non-negative
            B_next[i] = max(B_next[i], 1e-10)

        return B_next
    
    def run_full_algorithm(self, m=10):
        """
        Run the complete Jacobi-Newton iterative scheme for pricing American options.

        Parameters:
        - m (int): Number of iterations for the Jacobi-Newton scheme.

        Returns:
        - B_values (numpy array): The final boundary values after m iterations.
        """
        # Initialize boundary values using the approximate method (B^{(0)})
        B_values = self.initial_boundary.copy()

        for j in range(1, m + 1):
            print(f"Starting iteration {j}/{m}")

            # Step 5: Compute H(sqrt(tau)) and initialize Chebyshev interpolation
            H_values = self.compute_H()
            self.initialize_chebyshev_interpolation(H_values)

            # Step 6: Evaluate boundary using Clenshaw algorithm at adjusted points
            y_nodes = np.linspace(-1, 1, len(self.tau_nodes))
            evaluated_B_values = np.zeros((self.n, len(y_nodes)))

            for i in range(self.n):
                tau = self.tau_nodes[i]
                evaluated_B_values[i] = self.evaluate_boundary(tau, y_nodes)

            # Step 7: Compute N(tau_i, B) and D(tau_i, B), then compute f(tau_i, B)
            f_values = np.zeros(self.n)
            for i in range(self.n):
                tau = self.tau_nodes[i]
                N_values, D_values = self.compute_ND_values(tau, evaluated_B_values[i])
                f_values[i] = np.mean(N_values / D_values)

            # Step 8: Compute the derivative f'(tau_i, B) for each tau_i
            f_derivative_values = np.zeros(self.n)
            for i in range(self.n):
                tau = self.tau_nodes[i]
                f_derivative_values[i] = self.compute_f_derivative(tau, np.array([B_values[i]]))[0]

            # Step 9: Update B_values using the Jacobi-Newton scheme
            B_values = self.update_boundary(B_values)

            print(f"Iteration {j}/{m} completed.")

        print("Jacobi-Newton iterations completed.")
        return B_values
