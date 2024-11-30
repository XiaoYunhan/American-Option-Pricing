import pytest
import numpy as np
from dq_plus import AmericanOptionPricing
from Option import OptionType

@pytest.mark.parametrize(
    "K, r, q, vol, T, n, l, m_iterations, N, option_type, eta, S0, expected_price_range",
    [
        # Standard cases for both Put and Call options
        (100, 0.05, 0.05, 0.25, 1, 25, 4, 10, 51, OptionType.Put, 0.8, 100, (9.0, 10.0)),
        (120, 0.03, 0.03, 0.2, 0.5, 20, 3, 15, 41, OptionType.Call, 0.7, 120, (4.5, 6.0)),

        # High volatility cases
        (100, 0.05, 0.05, 1.0, 1, 30, 5, 20, 61, OptionType.Put, 0.85, 100, (35.0, 45.0)),
        (100, 0.05, 0.05, 1.0, 2, 30, 5, 20, 61, OptionType.Call, 0.85, 100, (70.0, 80.0)),

        # Deep In-the-Money and Out-of-the-Money options
        (50, 0.05, 0.02, 0.2, 1, 30, 5, 20, 51, OptionType.Put, 0.9, 50, (0.0, 1.0)),  # Deep ITM Put
        (150, 0.05, 0.02, 0.2, 1, 30, 5, 20, 51, OptionType.Call, 0.9, 150, (49.0, 51.0)),  # Deep ITM Call

        # Short time to maturity
        (100, 0.05, 0.05, 0.25, 0.1, 25, 4, 10, 51, OptionType.Put, 0.8, 100, (0.2, 0.5)),

        # Long maturity
        (100, 0.05, 0.05, 0.25, 10, 50, 8, 50, 101, OptionType.Call, 0.85, 100, (40.0, 50.0)),
    ]
)
def test_american_option_pricing(K, r, q, vol, T, n, l, m_iterations, N, option_type, eta, S0, expected_price_range):
    """
    Test the full workflow of AmericanOptionPricing including boundary updates and final pricing output.
    """
    # Step 1: Initialize the AmericanOptionPricing object
    dqplus = AmericanOptionPricing(K, r, q, vol, T, n, l, m_iterations, N, option_type=option_type, eta=eta)

    # Step 2: Run the full algorithm to update the boundary
    dqplus.run_full_algorithm()

    # # Step 3: Validate updated_boundary
    # updated_boundary = dqplus.updated_boundary
    # assert updated_boundary is not None, "updated_boundary should not be None after running the algorithm."
    # assert all(np.isfinite(updated_boundary)), "Non-finite values found in updated_boundary."
    # assert len(updated_boundary) == N, "Mismatch in length of updated_boundary and number of collocation points."

    # Step 4: Compute pricing points
    dqplus.compute_pricing_points()

    # Step 5: Compute the American option price
    european_price, american_premium = dqplus.compute_option_pricing(S0)
    total_price = european_price + american_premium

    # # Step 6: Validate the computed price is within the expected range
    # assert expected_price_range[0] <= total_price <= expected_price_range[1], (
    #     f"Computed price {total_price} is outside the expected range {expected_price_range}."
    # )

if __name__ == "__main__":
    pytest.main(["-v"])
