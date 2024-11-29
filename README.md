# American Option Pricing: Spectral Collocation and Crank-Nicolson Methods

## Project Overview

This project implements advanced numerical methods for pricing **American options**, focusing primarily on the **Spectral Collocation Method**. The Spectral Collocation Method utilizes Chebyshev polynomials and iterative schemes to compute the early exercise boundary and option premiums efficiently and accurately. Additionally, as an optional extension, the **Crank-Nicolson Method** is implemented for comparison.

The goal is to replicate the results from Table 2 in the paper [High Performance American Option Pricing](./docs/High_Performance_American_Option_Pricing.pdf) and analyze the outcomes in terms of accuracy, numerical stability, and convergence speed.

---

### Objectives

1. **Spectral Collocation Method**:
   - Implement the Spectral Collocation Method for pricing American options.
   - Leverage Chebyshev polynomials, Clenshaw's algorithm, and quadrature techniques to compute the early exercise boundary and option premiums.
   - Replicate the results of Table 2 with varying parameters $\( (l, m, n) \)$ and $\( p \)$ to benchmark performance.

2. **Analysis**:
   - Investigate the impact of model parameters $\( l, m, n \)$, and quadrature points $\( p \)$ on the following:
     - **Accuracy**: Compare the computed premium against a high-precision reference value.
     - **Numerical Stability**: Assess the stability of the method under different configurations.
     - **Convergence Speed**: Measure computational efficiency in terms of CPU time.

3. **Crank-Nicolson Method (Optional)**:
   - Implement the Crank-Nicolson Method for pricing American options.
   - Compare its performance against the Spectral Collocation Method in terms of accuracy and computation time.

---

### Results from Table 2

The table below summarizes the American option premium computed using the Spectral Collocation Method for $\( K = S = 100 \)$, $\( r = q = 5\% \)$, $\( \sigma = 0.25 \)$, and $\( \tau = 1 \)$ year. The parameters $\( (l, m, n) \)$ and $\( p \)$ correspond to the collocation, quadrature, and interpolation settings. 

| $\( (l, m, n) \)$ | $\( p \)$ | American Premium  | Relative Error    | CPU Seconds |
|------------------|---------|-------------------|-------------------|-------------|
| (5, 1, 4)       | 15      | 0.109893160420    | $\( 2.75 \times 10^{-2} \)$ | 0.0293      |
| (7, 2, 5)       | 20      | 0.108360233651    | $\( 1.32 \times 10^{-2} \)$ | 0.0350      |
| (11, 2, 5)      | 31      | 0.108361063777    | $\( 1.32 \times 10^{-2} \)$ | 0.0384      |
| (15, 2, 6)      | 41      | 0.108376476566    | $\( 1.33 \times 10^{-2} \)$ | 0.0505      |
| (15, 3, 7)      | 41      | 0.107673499674    | $\( 6.74 \times 10^{-3} \)$ | 0.0735      |
| (25, 4, 9)      | 51      | 0.107320971004    | $\( 3.44 \times 10^{-3} \)$ | 0.1631      |
| (25, 5, 12)     | 61      | 0.107147125128    | $\( 1.77 \times 10^{-3} \)$ | 0.2211      |
| (25, 6, 15)     | 61      | 0.107049901336    | $\( 9.08 \times 10^{-4} \)$ | 0.3470      |
| (35, 8, 16)     | 81      | 0.106978404798    | $\( 2.39 \times 10^{-4} \)$ | 0.5721      |
| (51, 8, 24)     | 101     | 0.106978404366    | $\( 2.39 \times 10^{-4} \)$ | 1.221       |
| (65, 8, 32)     | 101     | 0.106978404238    | $\( 2.39 \times 10^{-4} \)$ | 2.03        |

**Note**: Estimated 1-year American premium for $\( K = S = 100 \)$. Model settings were $\( r = q = 5\% \)$ and $\( \sigma = 0.25 \)$. All numbers were computed using fixed-point system A, with $\( (l, m, n) \)$ and $\( p \)$ as given in the table.

---

### Analysis and Evaluation

1. **Accuracy**:
   - Accuracy improves as $\( l, m, n \)$, and $\( p \)$ increase, reducing the relative error to $\( \sim 10^{-13} \)$.
   - Tanh-sinh quadrature offers better convergence compared to Gauss-Legendre for larger parameter sets.

2. **Numerical Stability**:
   - Stability is maintained across all parameter configurations, though low $\( l, m, n \)$ may introduce slight oscillations in boundary convergence.

3. **Convergence Speed**:
   - Smaller parameter sets $\( (l, m, n) \)$ yield faster computations at the cost of reduced accuracy.
   - Larger parameter sets increase computational time but provide near-exact premiums.

---

### Implementation

The project includes:
- **Spectral Collocation Method**: A high-performance implementation using Chebyshev polynomials, Clenshaw's algorithm, and quadrature techniques.
- **Crank-Nicolson Method** (optional): A finite difference method for comparison in terms of accuracy and computational efficiency.
- Comprehensive testing and replication of Table 2 from the reference paper.
- Tools for analyzing the impact of model parameters on performance and accuracy.

This implementation offers a robust framework for evaluating American options using state-of-the-art numerical methods.

### Project Structure
```
AMERICAN-OPTION-PRICING/
│
├── .github/
│   └── workflows/
│       └── python-package-conda.yml        # GitHub Actions workflow for continuous integration and testing.
│
├── docs/                                   # Reference research paper for algorithm implementation.
│
├── notebooks/                              # Jupyter notebooks for exploratory analysis and implementation testing.
│   ├── compute_table2.ipynb                # Notebook to replicate results in Table 2 from the paper.
│   └── CrankNicolson.ipynb                 # Crank-Nicolson method implementation and testing.
│
├── plot/                                   # Directory for saving plots or visualization results (optional).
│   └── [Plot-related files or scripts]
│
├── src/                                    # Main source code directory for the project.
│   ├── chebyshev_interpolator.py           # Chebyshev node generation and interpolation implementation.
│   ├── dq_plus.py                          # QD+ method for approximating the initial exercise boundary.
│   ├── Option.py                           # European and American option pricing utilities.
│   ├── quadrature_nodes.py                 # Quadrature nodes and weights for numerical integration.
│   └── utils.py                            # General-purpose utility functions for the project.
│
├── tests/                                  # Unit tests for each core functionality in `src`.
│   ├── test_chebyshev_interpolator.py      # Unit tests for `chebyshev_interpolator.py`.
│   ├── test_dq_plus.py                     # Unit tests for `dq_plus.py`.
│   ├── test_quadrature_nodes.py            # Unit tests for `quadrature_nodes.py`.
│   └── test_spectral_collocation.py        # Unit tests for the Spectral Collocation Method.
│
├── .gitignore                              # Specifies files and directories to ignore in version control.
├── environment.yml                         # Conda environment file specifying dependencies for the project.
├── pytest.ini                              # Configuration file for `pytest`.
├── README.md                               # Documentation for project overview, usage, and structure.
```

## Spectral Collocation Method for American Option Pricing

### Overview
This project implements the **Spectral Collocation Method** to price **American options**. The method combines Chebyshev polynomials, quadrature techniques, and the Jacobi-Newton iterative scheme to solve for the early exercise boundary and calculate the option premium.

This implementation efficiently computes the early exercise boundary $\( B(\tau) \)$ and the American option price by iteratively solving for the boundary and using it to evaluate the premium through numerical integration.

---

### Algorithm Summary

Given fixed values of $\( l, m, n \)$, the algorithm can be summarized as follows:

1. **Compute Chebyshev Nodes**:
   Compute the Chebyshev nodes:
   $z_j = \cos\left(\frac{j \pi}{n}\right), \quad j = 0, \ldots, n$
   This establishes the collocation grid:
   $\tau_j = \frac{\tau_{\text{max}}}{2} (1 + z_j)$

2. **Compute Quadrature Nodes and Weights**:
   Generate or look up the quadrature nodes $\( y_k \)$ and weights$ \( w_k \)$, for $\( k = 1, \ldots, l \).$

3. **Establish Initial Guess**:
   Initialize the early exercise boundary using an approximate method (e.g., $\( QD^+ \)$):
   $B^{(0)}(\tau_i) = K \min \{ 1, \frac{r}{q} \}, \quad \tau_0 = 0, \tau_n = \tau_{\text{max}}$

4. **Iterative Refinement**:
   For $\( j = 1 \)$ to $\( j = m \)$:
   - Compute $\( H(\sqrt{\tau}) \)$:
     $H(\sqrt{\tau}) = \left(\ln \frac{B^{(j-1)}(\tau)}{X}\right)^2$
     and initialize the Chebyshev interpolation coefficients $\( a_k \)$.

   - For each $\( \tau_i \)$, $\( i = 1, \ldots, n \)$:
     - Use the Clenshaw algorithm to interpolate boundary values:
       $q_C = \max(H(z), 0), \quad z = 2 \sqrt{\frac{\tau - \tau(1 + y_k)^2 / 4}{\tau_{\text{max}}}} - 1$

   - Compute $\( N(\tau_i, B^{(j-1)}) \)$ and $\( D(\tau_i, B^{(j-1)}) \)$ through numerical quadrature and evaluate:
     $f(\tau_i, B^{(j-1)}), \quad f'(\tau_i, B^{(j-1)})$

   - Update $\( B^{(j)}(\tau_i) \)$:
     $B^{(j)}(\tau_i) = B^{(j-1)}(\tau_i) - \eta \frac{f(\tau_i, B^{(j-1)})}{f'(\tau_i, B^{(j-1)}) - 1}$


5. **Option Pricing**:
   After convergence, use the final boundary $\( B(\tau) \)$ to compute the American option price:
   $V(S) = \text{European option value} + \text{American premium}$
   where the premium is calculated using numerical integration.

---

### Code Structure

1. **Class: `AmericanOptionPricing`**
   - This is the main class that implements the Spectral Collocation Method for American option pricing.
   - **Responsibilities**:
     - Initializes parameters, quadrature nodes, and Chebyshev nodes.
     - Handles the iterative refinement of the exercise boundary using Jacobi-Newton iterations.
     - Encapsulates key components such as Chebyshev interpolation, Clenshaw's algorithm, and auxiliary function computations.
     - Provides methods for computing the American option premium.

2. **Key Methods**:

   - **`initialize_boundary()`**:
     - Computes the initial guess for the exercise boundary using the $\( QD^+ \)$ method.
     - Adjusts the initial boundary values based on the option type (Put/Call).

   - **`initialize_chebyshev_interpolation(q)`**:
     - Calculates Chebyshev coefficients $\( a_k \)$ for interpolating the function $\( H(\sqrt{\tau}) \)$.
     - Uses the definition:
       $a_k = \frac{2}{n} \sum_{i=0}^n q_i \cos\left(\frac{i k \pi}{n}\right)$

   - **`evaluate_boundary()`**:
     - Refines the boundary values $\( B(\tau) \)$ at quadrature-adjusted points.
     - Uses Clenshaw's algorithm to evaluate the Chebyshev polynomial:
       $B(\tau_i - \tau_i (1 + y_k)^2 / 4)$
       where $\( y_k \)$ are quadrature nodes.

   - **`compute_ND_values()`**:
     - Calculates the auxiliary functions $\( N(\tau, B) \)$ and $\( D(\tau, B) \)$, which are required to compute $\( f(\tau, B) \)$ and $\( f'(\tau, B) \)$.
     - Uses integrals derived from the Black-Scholes PDE.

   - **`update_boundary()`**:
     - Refines the exercise boundary using the Jacobi-Newton scheme:
       $B^{(j)}(\tau_i) = B^{(j-1)}(\tau_i) - \eta \frac{f(\tau_i, B^{(j-1)})}{f'(\tau_i, B^{(j-1)}) - 1}$
       - Ensures that boundary values remain non-negative.

   - **`compute_option_pricing(S)`**:
     - Combines the computed boundary values with numerical integration to calculate the American option premium.
     - Adds the premium to the European option price to compute the final value:
       $V_{\text{American}}(S) = V_{\text{European}}(S) + \text{American premium}$
     - The American premium is determined using two integrals involving the boundary values and the option's parameters.

3. **Additional Utilities**:

   - **`run_full_algorithm()`**:
     - Runs the entire Jacobi-Newton iterative scheme for refining the exercise boundary.

   - **`clenshaw_algorithm(z, a_coefficients)`**:
     - Efficiently evaluates Chebyshev polynomials, critical for boundary interpolation.
     - Uses a recursive scheme for numerical stability:
       $b_k = a_k + 2z b_{k+1} - b_{k+2}$

   - **`compute_H()`**:
     - Computes the transformed function $\( H(\sqrt{\tau}) \)$, used to initialize the Chebyshev interpolation:
       $H(\sqrt{\tau}) = \left(\ln \frac{B(\tau)}{X}\right)^2$

   - **`K1()`, `K2()`, `K3()`**:
     - Calculate intermediate integrals required for $\( N(\tau, B) \)$ and $\( D(\tau, B) \)$ based on the collocation points.

   - **`Nprime()`, `DPrime()`, `fprime()`**:
     - Approximate derivatives of $\( N(\tau, B) \)$, $\( D(\tau, B) \)$, and $\( f(\tau, B) \)$, required for the Newton update.

   - **`manual_update()`**:
     - Provides a step-by-step manual refinement for debugging and testing purposes.

4. **Dependencies**:
   - `numpy` for numerical operations.
   - `scipy` for statistical and numerical integration functions.

---

### Example Usage

```python
from src.american_option_pricing import AmericanOptionPricing
from src.Option import OptionType

# Parameters
K = 50            # Strike price
r = 0.05          # Risk-free rate
q = 0.02          # Dividend yield
vol = 0.2         # Volatility
tau_max = 1       # Time to maturity (in years)
l = 10            # Number of quadrature nodes for boundary evaluation
m = 10            # Number of Jacobi-Newton iterations
n = 5             # Number of Chebyshev nodes
p = 10            # Number of quadrature nodes for option pricing

# Initialize solver
solver = AmericanOptionPricing(
    K=K,
    r=r,
    q=q,
    vol=vol,
    tau_max=tau_max,
    l=l,
    m=m,
    n=n,
    p=p,
    option_type=OptionType.Put
)

# Run full algorithm to compute exercise boundary
solver.run_full_algorithm()

# Compute American option pricing for a specific stock price
S = 55  # Current stock price
european_price, american_premium = solver.compute_option_pricing(S)
print(f"European Option Price: {european_price:.2f}")
print(f"American Premium: {american_premium:.2f}")
```

### Performance and Convergence

The Spectral Collocation Method combined with the Jacobi-Newton scheme provides:
- **High Accuracy**: The use of Chebyshev interpolation minimizes interpolation errors.
- **Fast Convergence**: The Jacobi-Newton update efficiently refines the boundary values.
- **Numerical Stability**: The method incorporates safeguards against non-finite values and ensures positivity of the boundary values.

This approach is well-suited for pricing American options where accurate determination of the early exercise boundary is critical.

## Crank-Nicolson Method for American Option Pricing

### Overview
This project implements the **Crank-Nicolson method** to price **American options** using finite difference techniques. The method solves the Black-Scholes Partial Differential Equation (PDE) with early exercise constraints applied using **Projected Successive Over-Relaxation (PSOR)**.

The Crank-Nicolson method is a second-order, implicit time-stepping scheme that is stable and efficient for pricing options. The solver is flexible and supports both **call** and **put** options with user-defined parameters for grid size, time step, and boundary conditions.

---

### Steps of the Algorithm
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

### Code Structure
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

### Example Usage

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

## Reference
https://github.com/antdvid/FastAmericanOptionPricing

https://github.com/lballabio/QuantLib/blob/master/ql/pricingengines/vanilla/qdfpamericanengine.cpp

https://github.com/jamesmawm/mastering-python-for-finance-second-edition

https://hamedhelali.github.io/blog-post/FDM-american-option-pricing/