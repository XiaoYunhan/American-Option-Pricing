import pytest
import numpy as np
from crank_nicolson import CrankNicolsonSolver

@pytest.mark.parametrize(
    "S0, K, r, q, sigma, T, N, M, option_type",
    [
        (100, 100, 0.05, 0.02, 0.2, 1, 50, 50, "put"),
        (120, 100, 0.03, 0.03, 0.25, 0.5, 100, 100, "call"),
        (80, 100, 0.07, 0.01, 0.3, 2, 150, 150, "put"),
    ]
)
def test_crank_nicolson_solver(S0, K, r, q, sigma, T, N, M, option_type):
    """
    Test the Crank-Nicolson solver for American options to ensure the final option price is valid.
    """
    # Initialize the Crank-Nicolson solver
    solver = CrankNicolsonSolver(
        riskfree=r,
        dividend=q,
        volatility=sigma,
        strike=K,
        maturity=T,
        option_type=option_type,
        N=N,
        M=M
    )

    # Solve the PDE for the given initial stock price S0
    option_price = solver.solve(S0)

    # Validate the option price
    assert option_price >= 0, "Negative option price computed."
    assert np.isfinite(option_price), "Non-finite option price computed."

    # Validate the option price for deep in-the-money and out-of-the-money cases
    if option_type == "put":
        if S0 == 0:
            assert np.isclose(option_price, K * np.exp(-r * T), atol=1e-2), "Put option price does not match boundary condition for S -> 0."
        if S0 > 2 * K:
            assert np.isclose(option_price, 0, atol=1e-2), "Put option price does not match boundary condition for S -> ∞."
    elif option_type == "call":
        if S0 == 0:
            assert np.isclose(option_price, 0, atol=1e-2), "Call option price does not match boundary condition for S -> 0."
        if S0 > 2 * K:
            assert np.isclose(option_price, S0 - K * np.exp(-r * T), atol=1e-2), "Call option price does not match boundary condition for S -> ∞."

if __name__ == "__main__":
    pytest.main(["-v"])
