# American Option Pricing

## Project Requirement
1. Implement Spectral Collocation Method for pricing American options
2. Replicate Table 2 in the paper [High Performance American Option Pricing](/docs/High_Performance_American_Option_Pricing.pdf)
![Table2](./plot/table2.jpg)
1. Analyze the results in terms of accuracy, numerical stability and convergence speed for the difference choices of model parameters l, m, n and spot price S, interest rate r, dividend q and time to maturity τ .
2. (optional) Implement Crank-Nicolson method and compare it with Spectral Collocation Method.

## Reference
https://github.com/antdvid/FastAmericanOptionPricing

https://github.com/lballabio/QuantLib/blob/master/ql/pricingengines/vanilla/qdfpamericanengine.cpp