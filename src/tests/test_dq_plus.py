# import pytest
# import numpy as np
# from src.dq_plus import AmericanOptionPricing

# @pytest.mark.parametrize("K, r, q, vol, n, l, expected_first_value", [
#     (100, 0.05, 0.02, 0.2, 10, 5, 100 * min(1, 0.05 / 0.02)),
#     (120, 0.03, 0.03, 0.25, 20, 7, 120 * min(1, 0.03 / 0.03)),
#     (80, 0.07, 0.01, 0.3, 15, 4, 80 * min(1, 0.07 / 0.01)),
# ])
# def test_initialize_boundary(K, r, q, vol, n, l, expected_first_value):
#     dqplus = AmericanOptionPricing(K, r, q, vol, n, l)
#     dqplus.initialize_boundary()
#     initial_boundary = dqplus.get_boundary_values()
#     assert initial_boundary[0] == pytest.approx(expected_first_value, rel=1e-6)

# # @pytest.mark.parametrize("K, r, q, vol, n, l", [
# #     (100, 0.05, 0.02, 0.2, 10, 5),
# #     (120, 0.03, 0.03, 0.25, 20, 7),
# #     (80, 0.07, 0.01, 0.3, 15, 4),
# # ])
# # def test_fixed_point_iteration_convergence(K, r, q, vol, n, l):
# #     dqplus = DQPlus(K, r, q, vol, n, l)
# #     dqplus.initialize_boundary()
# #     dqplus.fixed_point_iteration()
# #     refined_boundary = dqplus.get_boundary_values()
# #     assert all(refined_boundary >= 0)
# #     assert all(refined_boundary[i] >= refined_boundary[i + 1] for i in range(len(refined_boundary) - 1))

# @pytest.mark.parametrize("K, r, q, vol, n, l", [
#     (100, 0.05, 0.02, 0.2, 10, 5),
#     (120, 0.03, 0.03, 0.25, 20, 7),
#     (80, 0.07, 0.01, 0.3, 15, 4),
# ])
# def test_compute_H(K, r, q, vol, n, l):
#     dqplus = AmericanOptionPricing(K, r, q, vol, n, l)
#     dqplus.initialize_boundary()
#     H_values = dqplus.compute_H()
#     assert all(H >= 0 for H in H_values)

# @pytest.mark.parametrize("K, r, q, vol, n, l", [
#     (100, 0.05, 0.02, 0.2, 10, 5),
#     (120, 0.03, 0.03, 0.25, 20, 7),
#     (80, 0.07, 0.01, 0.3, 15, 4),
# ])
# def test_initialize_chebyshev_interpolation(K, r, q, vol, n, l):
#     dqplus = AmericanOptionPricing(K, r, q, vol, n, l)
#     dqplus.initialize_boundary()
#     H_values = dqplus.compute_H()
#     a_coefficients = dqplus.initialize_chebyshev_interpolation(H_values)
#     assert all(np.isfinite(a) for a in a_coefficients)
#     assert np.any(a_coefficients != 0)

# @pytest.mark.parametrize("K, r, q, vol, n, l", [
#     (100, 0.05, 0.02, 0.2, 10, 5),
#     (120, 0.03, 0.03, 0.25, 20, 7),
#     (80, 0.07, 0.01, 0.3, 15, 4),
# ])
# def test_evaluate_boundary(K, r, q, vol, n, l):
#     dqplus = AmericanOptionPricing(K, r, q, vol, n, l)
#     dqplus.initialize_boundary()
#     H_values = dqplus.compute_H()
#     dqplus.initialize_chebyshev_interpolation(H_values)
#     tau_test = dqplus.tau_nodes[len(dqplus.tau_nodes) // 2]
#     B_values = dqplus.evaluate_boundary(tau_test)
#     assert all(np.isfinite(B) for B in B_values)
#     assert all(B >= 0 for B in B_values)

# @pytest.mark.parametrize("K, r, q, vol, n, l", [
#     (100, 0.05, 0.02, 0.2, 10, 5),
#     (120, 0.03, 0.03, 0.25, 20, 7),
#     (80, 0.07, 0.01, 0.3, 15, 4),
# ])
# def test_compute_ND_values(K, r, q, vol, n, l):
#     dqplus = AmericanOptionPricing(K, r, q, vol, n, l)
#     dqplus.initialize_boundary()
#     tau_test = dqplus.tau_nodes[len(dqplus.tau_nodes) // 2]
#     B_values = dqplus.evaluate_boundary(tau_test)
#     N_values, D_values = dqplus.compute_ND_values(tau_test, B_values)
#     assert all(np.isfinite(N) for N in N_values)
#     assert all(np.isfinite(D) for D in D_values)
#     assert all(D > 0 for D in D_values)

# @pytest.mark.parametrize("K, r, q, vol, n, l", [
#     (100, 0.05, 0.02, 0.2, 10, 5),
#     (120, 0.03, 0.03, 0.25, 20, 7),
#     (80, 0.07, 0.01, 0.3, 15, 4),
# ])
# def test_update_boundary(K, r, q, vol, n, l):
#     dqplus = AmericanOptionPricing(K, r, q, vol, n, l)
#     dqplus.initialize_boundary()

#     # Ensure Chebyshev coefficients are initialized
#     H_values = dqplus.compute_H()
#     dqplus.initialize_chebyshev_interpolation(H_values)

#     # Perform boundary update
#     B_values = dqplus.get_boundary_values()
#     B_next = dqplus.update_boundary(B_values)

#     assert all(B >= 0 for B in B_next), "Negative boundary value found."
#     assert all(np.isfinite(B) for B in B_next), "Non-finite boundary value found."
#     assert all(B_next[i] >= B_values[i] * 0.8 for i in range(len(B_values))), "Boundary value dropped significantly."
    
# if __name__ == "__main__":
#     pytest.main(["-v"])
