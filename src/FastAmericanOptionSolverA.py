from src.FastAmericanOptionSolverBase import AmericanOptionSolver
from src.Option import OptionType

import numpy as np

class FastAmericanOptionSolverA(AmericanOptionSolver):
    # all self.quadrature_num here: self.p?
    def N_func(self, tau, B):
        K3 = self.K3(tau, B)
        if self.option_type == OptionType.Put:
            return self.PDF_dminus(tau, B/self.K)/(self.vol * np.sqrt(tau)) + self.r * K3
        else:
            return self.PDF_neg_dminus(tau, B / self.K) / (self.vol * np.sqrt(tau)) + self.r * K3
            #return self.PDF_dminus(tau, B/self.K)/ (self.vol * np.sqrt(tau)) + self.r * K3

    def D_func(self, tau, B):
        K12 = self.K12(tau, B)
        if self.option_type == OptionType.Put:
            return self.PDF_dplus(tau, B/self.K)/(self.vol * np.sqrt(tau)) + self.CDF_pos_dplus(tau, B/self.K) + self.q * (K12)
        else:
            return self.PDF_neg_dplus(tau, B/self.K)/(self.vol * np.sqrt(tau)) - self.CDF_neg_dplus(tau, B/self.K) + self.q * (K12)
            #return self.PDF_dplus(tau, B/self.K)/(self.vol * np.sqrt(tau)) + self.CDF_pos_dplus(tau, B/self.K) + self.q * (K12) - np.exp(self.q * tau)

    def K12(self, tau, B):
        return self.quadrature_sum(self.K12_integrand, tau, B, self.p) #?

    def K12_integrand(self, tau, B_tau, u, B_u):
        if self.option_type == OptionType.Put:
            return np.exp(self.q * u) * self.CDF_pos_dplus(tau - u, B_tau/B_u) \
                + np.exp(self.q * u)/(self.vol * np.sqrt(tau - u)) * self.PDF_dplus(tau-u, B_tau/B_u)
        else:
            return -np.exp(self.q * u) *  self.CDF_neg_dplus(tau - u, B_tau/B_u) \
                   + np.exp(self.q * u) / (self.vol * np.sqrt(tau - u)) * self.PDF_neg_dplus(tau-u, B_tau/B_u)
            #return np.exp(self.q * u) * self.CDF_pos_dplus(tau - u, B_tau/B_u) \
             #   + np.exp(self.q * u)/(self.vol * np.sqrt(tau - u)) * self.PDF_dplus(tau-u, B_tau/B_u)

    def K3(self, tau, B):
            return self.quadrature_sum(self.K3_integrand, tau, B, self.p)

    def K3_integrand(self, tau, B_tau, u, B_u):
        ans = 0
        if self.option_type == OptionType.Put:
            ans = np.exp(self.r * u)/(self.vol * np.sqrt(tau - u)) * self.PDF_dminus(tau-u, B_tau/B_u)
        else:
            ans = np.exp(self.r * u) / (self.vol * np.sqrt(tau - u)) * self.PDF_neg_dminus(tau-u, B_tau/B_u)
            #return np.exp(self.r * u)/(self.vol * np.sqrt(tau - u)) * self.PDF_dminus(tau-u, B_tau/B_u)
        return ans

    def Nprime_func(self, tau, Q):
        tau = max(1e-10, tau)
        if self.option_type == OptionType.Put:
            return -self.dminus(tau, Q/self.K) * self.PDF_dminus(tau, Q/self.K) / (Q * self.vol * self.vol * tau) \
                    - self.r * self.quadrature_sum(self.Nprime_integrand, tau, Q, self.p)
        else:
            return -self.dminus(tau, Q/self.K) * self.PDF_neg_dminus(tau, Q/self.K) / (Q * self.vol * self.vol * tau) \
                 - self.r * self.quadrature_sum(self.Nprime_integrand, tau, Q, self.p)

    def Nprime_integrand(self, tau, B_tau, u, B_u):
        tau = max(1e-10, tau)
        if self.option_type == OptionType.Put:
            return np.exp(self.r * u) * self.dminus(tau-u, B_tau/B_u) / (B_tau * self.vol * self.vol * (tau - u))\
               * self.PDF_dminus(tau-u, B_tau/B_u)
        else:
            return np.exp(self.r * u) * self.dminus(tau - u, B_tau / B_u) / (
                        B_tau * self.vol * self.vol * (tau - u)) \
                   * self.PDF_neg_dminus(tau - u, B_tau / B_u)

    def Dprime_func(self, tau, Q):
        tau = max(1e-10, tau)
        if self.option_type == OptionType.Put:
            return self.PDF_dplus(tau, Q/self.K)/(self.vol * np.sqrt(tau) * Q) \
               * (1 - self.dplus(tau, Q/self.K)/(self.vol * np.sqrt(tau))) \
               + self.q * self.quadrature_sum(self.Dprime_integrand, tau, Q, self.p)
        else:
            return self.PDF_neg_dplus(tau, Q / self.K) / (self.vol * np.sqrt(tau) * Q) \
                   * (1 - self.dplus(tau, Q / self.K) / (self.vol * np.sqrt(tau))) \
                   + self.q * self.quadrature_sum(self.Dprime_integrand, tau, Q, self.p)

    def Dprime_integrand(self, tau, B_tau, u, B_u):
        tau = max(1e-10, tau)
        if self.option_type == OptionType.Put:
            return np.exp(self.q * u) * self.PDF_dplus(tau-u, B_tau/B_u)/(self.vol * np.sqrt(tau - u) * B_tau) \
            * (1 - self.dplus(tau - u, B_tau/B_u)/(self.vol * np.sqrt(tau - u)))
        else:
            return np.exp(self.q * u) * self.PDF_neg_dplus(tau - u, B_tau / B_u) / (self.vol * np.sqrt(tau - u) * B_tau) \
                   * (1 - self.dplus(tau - u, B_tau / B_u) / (self.vol * np.sqrt(tau - u)))
