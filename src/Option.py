from scipy.stats import norm
import numpy as np
from enum import Enum

class OptionType(Enum):
    Call = 'call'
    Put = 'put'

class EuropeanOption:
    def __init__(self, 
                 S0, K, T, r, q, sigma
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
        self.T = T
        self.r = r
        self.q = q
        self.sigma = sigma

    def _d1(self):
        return (np.log(self.S0 / self.K) + (self.r - self.q + 0.5 * self.sigma ** 2) * self.T
                ) / (self.sigma * np.sqrt(self.T))
    
    def _d2(self):
        return self._d1() - self.sigma * np.sqrt(self.T)

    def put_price(self):
        d1 = self._d1()
        d2 = self._d2()
        price = (
            self.K * np.exp(-self.r * self.T) * norm.cdf(-d2) - 
            self.S0 * np.exp(-self.q * self.T) * norm.cdf(-d1)
        )
        return price
    
