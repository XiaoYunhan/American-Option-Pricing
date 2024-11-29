# American Option Pricing

## Project Requirement
1. Implement Spectral Collocation Method for pricing American options
2. Replicate Table 2 in the paper [High Performance American Option Pricing](/docs/High_Performance_American_Option_Pricing.pdf)
<!-- ![Table2](./plot/table2.jpg) -->
| (l, m, n)  | p   | American Premium  | Relative Error | CPU Seconds |
|------------|------|-------------------|----------------|-------------|
| (5, 1, 4)  | 15   | 0.106783919132    | 1.5E-03        | 9.9E-06     |
| (7, 2, 5)  | 20   | 0.106934846948    | 1.7E-04        | 2.3E-05     |
| (11, 2, 5) | 31   | 0.106939863585    | 1.2E-04        | 3.1E-05     |
| (15, 2, 6) | 41   | 0.106954833468    | 2.0E-05        | 4.5E-05     |
| (15, 3, 7) | 41   | 0.106952928838    | 2.1E-06        | 7.1E-05     |
| (25, 4, 9) | 51   | 0.106952731254    | 2.7E-07        | 1.7E-04     |
| (25, 5, 12)| 61   | 0.106952704598    | 1.7E-08        | 2.9E-04     |
| (25, 6, 15)| 61   | 0.106952703049    | 2.8E-09        | 4.6E-04     |
| (35, 8, 16)| 81   | 0.106952702764    | 1.6E-10        | 8.4E-04     |
| (51, 8, 24)| 101  | 0.106952702748    | 1.5E-11        | 2.1E-03     |
| (65, 8, 32)| 101  | 0.106952702747    | 8.5E-13        | 3.0E-03     |

*Table 2*: Estimated 1-year American premium for \( K = S = 100 \). Model settings were \( r = q = 5\% \) and \( \sigma = 0.25 \). All numbers were computed using fixed point system A, with \( (l, m, n) \) and \( p \) as given in the table. Relative errors are measured against the American premium computed with \( (l, m, n) = (201, 16, 64) \) and \( p = 201 \). Results for \( (l, m, n) = (5, 1, 4) \) and \( (l, m, n) = (7, 2, 5) \) were computed with Gauss-Legendre quadrature; all other results were computed with tanh-sinh quadrature.

1. Analyze the results in terms of accuracy, numerical stability and convergence speed for the difference choices of model parameters l, m, n and spot price S, interest rate r, dividend q and time to maturity τ .
2. (optional) Implement Crank-Nicolson method and compare it with Spectral Collocation Method.

## Spectral Collocation Method for Pricing American Options

In this project, we implement the **Spectral Collocation Method** for pricing American options using the Jacobi-Newton iterative scheme combined with Chebyshev interpolation. This numerical method is efficient and accurate for solving free-boundary problems such as pricing American options, where the early exercise boundary must be determined.

### Overview

The Spectral Collocation Method solves for the early exercise boundary of American options by:
1. Establishing a collocation grid using Chebyshev nodes.
2. Using numerical quadrature and Chebyshev interpolation to compute boundary values iteratively.
3. Refining the boundary values using a Jacobi-Newton scheme to achieve convergence.

### Steps of the Algorithm

The complete algorithm involves the following steps:

1. **Compute Chebyshev Nodes**:
   - We compute the Chebyshev nodes $x_i$ based on the extremum points of the Chebyshev polynomial.
   - The collocation time points $\tau_i$ are determined using $\tau_i = x_i^2$, creating a grid of points for the numerical solution.

2. **Compute Quadrature Nodes and Weights**:
   - Using the Gauss-Legendre quadrature method, we compute the quadrature nodes $y_k$ and weights $w_k$ for accurate numerical integration.

3. **Initialize the Early Exercise Boundary**:
   - The initial guess $B^{(0)}(\tau_i)$ is set using an approximate method (QD+ approach).
   - The initial boundary value at $\tau_0$ is $B^{(0)}(\tau_0) = K \min\left(1, \frac{r}{q}\right)$, and subsequent values are estimated using an exponential decay model.

4. **Compute $H(\sqrt{\tau})$ and Initialize Chebyshev Interpolation**:
   - We compute $H(\sqrt{\tau}) = \left(\ln\left(\frac{B^{(j-1)}(\tau)}{K}\right)\right)^2$ at each collocation point.
   - The Chebyshev interpolation coefficients $a_k$ are determined using the discrete Chebyshev transform.

5. **Evaluate the Boundary Using Clenshaw Algorithm**:
   - For each $\tau_i$, we use the Clenshaw algorithm to evaluate the boundary values at adjusted points $\tau_i - \tau_i(1 + y_k)^2 / 4$.
   - This provides an efficient way to compute the polynomial evaluation using the Chebyshev coefficients.

6. **Compute $N(\tau_i, B)$ and $D(\tau_i, B)$ Using Numerical Quadrature**:
   - We use numerical integration to compute the numerator $N(\tau_i, B)$ and denominator $D(\tau_i, B)$ based on the integrands involving normal distributions.
   - These values are used to determine the function $f(\tau, B)$, which is central to the fixed-point iteration scheme.

7. **Compute the Derivative of $f(\tau_i, B)$**:
   - We approximate the derivative $f'(\tau, B)$ using finite differences.
   - This derivative is critical for the Jacobi-Newton update step.

8. **Update Boundary Values Using Jacobi-Newton Scheme**:
   - The boundary values are updated iteratively using the Jacobi-Newton formula:
     $B^{(j+1)}(\tau) = B^{(j)}(\tau) + \eta \frac{B^{(j)}(\tau) - f(\tau, B^{(j)})}{f'(\tau, B^{(j)}) - 1}$
   - The hyper-parameter $\eta$ controls the step size, ensuring stability of the iteration.

9. **Iterate Until Convergence**:
   - The algorithm iterates through steps 5 to 8 for a fixed number of iterations (typically $m$ iterations).
   - The process stops when the boundary values converge within a specified tolerance.

### Code Structure

The core implementation of the Spectral Collocation Method is encapsulated in the `DQPlus` class. Here is a brief overview of the key components:

- **`DQPlus` Class**:
  - **Initialization** (`__init__`): Sets up the option parameters and collocation grid.
  - **Boundary Initialization** (`initialize_boundary`): Initializes the early exercise boundary using the QD+ approximation.
  - **Chebyshev Interpolation** (`initialize_chebyshev_interpolation`): Computes the Chebyshev coefficients for interpolation.
  - **Clenshaw Algorithm** (`clenshaw_algorithm`): Efficiently evaluates the Chebyshev polynomial.
  - **Boundary Evaluation** (`evaluate_boundary`): Evaluates the boundary values at adjusted points using the Clenshaw algorithm.
  - **Numerical Quadrature** (`compute_ND_values`): Computes $N$ and $D$ using numerical integration.
  - **Fixed-Point Iteration** (`fixed_point_iteration`): Refines the boundary values iteratively using Newton's method.
  - **Jacobi-Newton Update** (`update_boundary`): Updates the boundary values using the Jacobi-Newton formula.
  - **Full Algorithm Execution** (`run_full_algorithm`): Executes the complete algorithm for a fixed number of iterations.

### Example Usage

The following is an example of how to use the `DQPlus` class to compute the early exercise boundary for an American option:

```python
import numpy as np
from src.dq_plus import DQPlus

# Define option parameters
K = 100       # Strike price
r = 0.05      # Risk-free rate
q = 0.02      # Dividend yield
vol = 0.2     # Volatility
tau_nodes = np.linspace(0, 1, 10)  # Collocation grid

# Initialize the DQPlus engine
dqplus = DQPlus(K, r, q, vol, tau_nodes)

# Initialize the boundary and run the full algorithm
dqplus.initialize_boundary()
final_boundary = dqplus.run_full_algorithm(m=10)

# Output the final computed boundary values
print("Final early exercise boundary values:", final_boundary)
```

### Performance and Convergence

The Spectral Collocation Method combined with the Jacobi-Newton scheme provides:
- **High Accuracy**: The use of Chebyshev interpolation minimizes interpolation errors.
- **Fast Convergence**: The Jacobi-Newton update efficiently refines the boundary values.
- **Numerical Stability**: The method incorporates safeguards against non-finite values and ensures positivity of the boundary values.

This approach is well-suited for pricing American options where accurate determination of the early exercise boundary is critical.

# Crank-Nicolson Method for American Option Pricing

## Overview
This project implements the **Crank-Nicolson method** to price **American options** using finite difference techniques. The method solves the Black-Scholes Partial Differential Equation (PDE) with early exercise constraints applied using **Projected Successive Over-Relaxation (PSOR)**.

The Crank-Nicolson method is a second-order, implicit time-stepping scheme that is stable and efficient for pricing options. The solver is flexible and supports both **call** and **put** options with user-defined parameters for grid size, time step, and boundary conditions.

---

## Steps of the Algorithm
1. **Discretize the PDE**:
   Transform the Black-Scholes PDE:
   $\frac{\partial V}{\partial t} + \frac{1}{2} \sigma^2 S^2 \frac{\partial^2 V}{\partial S^2} + (r - q) S \frac{\partial V}{\partial S} - rV = 0$
   into a tridiagonal system of equations using the transformation $\( x = \ln(S) \)$.

2. **Set Initial and Boundary Conditions**:
   - Initial condition: Payoff at \( t = T \):
     $V(S, T) = \max(K - S, 0) \quad \text{(Put)}, \quad V(S, T) = \max(S - K, 0) \quad \text{(Call)}.$
   - Boundary conditions at $\( S \to 0 \)$ and $\( S \to \infty \)$ are derived based on option type and time decay.

3. **Iterate Backwards in Time**:
   - Solve the discretized equations for each time step $\( t \)$ using the Crank-Nicolson scheme:
     $A X^{n+1} = B X^n$,
     where $\( A \)$ and $\( B \)$ are tridiagonal matrices.

4. **Apply Early Exercise Constraints**:
   Ensure the option price satisfies the early exercise condition:
   $V(S, t) \geq \max(K - S, 0) \quad \text{(Put)}, \quad V(S, t) \geq \max(S - K, 0) \quad \text{(Call)}.$

5. **Output the Results**:
   Interpolate the final option prices for any input stock price $\( S_0 \)$.

---

## Code Structure
1. **Class: `CrankNicolsonSolver`**
   - Encapsulates the entire Crank-Nicolson method, including grid setup, time-stepping, and PSOR.

2. **Key Methods**:
   - `__init__`: Initializes parameters and grids.
   - `solve(S0)`: Computes the option price for a given stock price $\( S_0 \)$.
   - `solvePDE()`: Iterates backward in time to solve the PDE.
   - `setInitialCondition()`: Sets the payoff at maturity.
   - `setCoeff(dt, t)`: Sets up the tridiagonal matrix coefficients.
   - `solveLinearSystem()`: Solves the linear system using a banded matrix solver.
   - `solvePSOR()`: Applies PSOR to enforce early exercise constraints.

3. **Dependencies**:
   - `numpy` for numerical operations.
   - `scipy.linalg` for solving tridiagonal systems.

---

## Example Usage

```python
from crank_nicolson import CrankNicolsonSolver

# Initialize solver
option_type = "put"  # Change to "call" for call options
solver = CrankNicolsonSolver(
    riskfree=0.05,         # Risk-free rate (5%)
    dividend=0.02,         # Dividend yield (2%)
    volatility=0.2,        # Volatility (20%)
    strike=50,             # Strike price
    maturity=1,            # Time to maturity (1 year)
    option_type=option_type
)
solver.max_dt = 0.01       # Maximum time step
solver.USE_PSOR = True     # Enable PSOR for early exercise

# Solve for option price
S0 = 50
price = solver.solve(S0)
print(f"The {option_type} option price for S0 = {S0} is: {price:.2f}")
```

## Reference
https://github.com/antdvid/FastAmericanOptionPricing

https://github.com/lballabio/QuantLib/blob/master/ql/pricingengines/vanilla/qdfpamericanengine.cpp