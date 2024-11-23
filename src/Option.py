from scipy.stats import norm
import numpy as np
from enum import Enum

from src.quadrature_nodes import QuadratureNodes

class OptionType(Enum):
    Call = 'call'
    Put = 'put'

class EuropeanOption:
    def __init__(self, 
                 S0, K, tau, r, q, sigma
                 ):
        """
        Initializes a European option.

        Parameters:
        S0 (float): Current stock price
        K (float): Strike price
        T (float): Time to maturity (in years)
        r (float): Risk-free interest rate
        q (float): Dividend yield
        sigma (float): Volatility of the underlying asset
        """
        self.S0 = S0
        self.K = K
        self.tau = tau
        self.r = r
        self.q = q
        self.sigma = sigma

    def _d1(self):
        return (np.log(self.S0 / self.K) + (self.r - self.q + 0.5 * self.sigma ** 2) * self.tau
                ) / (self.sigma * np.sqrt(self.tau))
    
    def _d2(self):
        return self._d1() - self.sigma * np.sqrt(self.tau)

    def put_price(self):
        d1 = self._d1()
        d2 = self._d2()
        price = (
            self.K * np.exp(-self.r * self.tau) * norm.cdf(-d2) -
            self.S0 * np.exp(-self.q * self.tau) * norm.cdf(-d1)
        )
        return price
    def call_value(self):
        d1 = self._d1()
        d2 = self._d2()
        price = (
            self.S0 * np.exp(-self.q * self.tau) * norm.cdf(d1) -
            self.K * np.exp(-self.r * self.tau) * norm.cdf(d2)
        )
        return price