#%%
"""
Black-Scholes Implied Volatility Calculation using Newton-Raphson Method.

This module implements the Newton-Raphson method to compute implied volatility
from market option prices. Implied volatility is the volatility parameter that,
when input into the Black-Scholes formula, yields the observed market price.

Created on Wed Oct 24 20:37:32 2018
@author: Lech A. Grzelak
"""
import numpy as np
import scipy.stats as st


def ImpliedVolatility(CP, S_0, K, sigma, tau, r, V_market):
    """
    Calculate implied volatility using Newton-Raphson method.
    
    This function iteratively solves for the volatility parameter σ that makes
    the Black-Scholes option price equal to the observed market price. The
    Newton-Raphson method uses the vega (∂V/∂σ) to rapidly converge to the solution.
    
    Args:
        CP (str): Option type - 'c' or '1' for call, 'p' or '-1' for put.
        S_0 (float): Current stock price.
        K (float): Strike price.
        sigma (float): Initial guess for implied volatility.
        tau (float): Time to maturity (in years).
        r (float): Risk-free interest rate (annualized).
        V_market (float): Observed market price of the option.
    
    Returns:
        float: Implied volatility that matches the market price.
    
    Note:
        The algorithm uses the Newton-Raphson update:
            σ_new = σ_old - (V_BS(σ) - V_market) / vega(σ)
        where vega = ∂V/∂σ
    """
    error = 1e10  # Initial error
    # Handy lambda expressions
    optPrice = lambda sigma: BS_Call_Option_Price(CP, S_0, K, sigma, tau, r)
    vega = lambda sigma: dV_dsigma(S_0, K, sigma, tau, r)
    
    # While the difference between the model and the market price is large
    # follow the iteration
    n = 1
    while error > 10e-10:
        g = optPrice(sigma) - V_market
        g_prim = vega(sigma)
        sigma_new = sigma - g / g_prim
    
        error = abs(g)
        sigma = sigma_new
        
        print(f'Iteration {n} with error = {error:.10f}')
        
        n += 1
    return sigma


def dV_dsigma(S_0, K, sigma, tau, r):
    """
    Calculate vega: the derivative of option price with respect to volatility.
    
    Vega measures the sensitivity of the option price to changes in volatility.
    For both calls and puts, vega is the same and always positive.
    
    Args:
        S_0 (float): Current stock price.
        K (float): Strike price.
        sigma (float): Volatility parameter.
        tau (float): Time to maturity (in years).
        r (float): Risk-free interest rate (annualized).
    
    Returns:
        float: Vega value (∂V/∂σ).
    
    Note:
        Vega = K * exp(-r*τ) * φ(d2) * √τ
        where φ is the standard normal PDF and d2 is from the BS formula.
    """
    d2 = (np.log(S_0 / float(K)) + (r - 0.5 * np.power(sigma, 2.0)) * tau) / \
         float(sigma * np.sqrt(tau))
    value = K * np.exp(-r * tau) * st.norm.pdf(d2) * np.sqrt(tau)
    return value


def BS_Call_Option_Price(CP, S_0, K, sigma, tau, r):
    """
    Calculate Black-Scholes option price for European call or put.
    
    This is the classic Black-Scholes formula for pricing European options
    under the assumption of geometric Brownian motion for the stock price.
    
    Args:
        CP (str): Option type - 'c' or '1' for call, 'p' or '-1' for put.
        S_0 (float): Current stock price.
        K (float): Strike price.
        sigma (float): Volatility parameter (annualized).
        tau (float): Time to maturity (in years).
        r (float): Risk-free interest rate (annualized).
    
    Returns:
        float: Black-Scholes option price.
    
    Note:
        Call: V = S_0*Φ(d1) - K*exp(-r*τ)*Φ(d2)
        Put:  V = K*exp(-r*τ)*Φ(-d2) - S_0*Φ(-d1)
        where Φ is the standard normal CDF.
    """
    d1 = (np.log(S_0 / float(K)) + (r + 0.5 * np.power(sigma, 2.0)) * tau) / \
         float(sigma * np.sqrt(tau))
    d2 = d1 - sigma * np.sqrt(tau)
    if str(CP).lower() == "c" or str(CP).lower() == "1":
        value = st.norm.cdf(d1) * S_0 - st.norm.cdf(d2) * K * np.exp(-r * tau)
    elif str(CP).lower() == "p" or str(CP).lower() == "-1":
        value = st.norm.cdf(-d2) * K * np.exp(-r * tau) - st.norm.cdf(-d1) * S_0
    return value


def main():
    """
    Main calculation routine for implied volatility example.
    
    Demonstrates the calculation of implied volatility from a given market price
    and verifies the result by computing the option price using the implied vol.
    """
    # Initial parameters and market quotes
    V_market = 2        # Market call option price
    K = 120             # Strike price
    tau = 1             # Time-to-maturity
    r = 0.05            # Interest rate
    S_0 = 100           # Today's stock price
    sigmaInit = 0.25    # Initial guess for implied volatility
    CP = "c"            # 'c' for call, 'p' for put
    
    sigma_imp = ImpliedVolatility(CP, S_0, K, sigmaInit, tau, r, V_market)
    message = f'''Implied volatility for CallPrice= {V_market}, strike K={K}, 
      maturity T= {tau}, interest rate r= {r} and initial stock S_0={S_0} 
      equals to sigma_imp = {sigma_imp:.7f}'''
                
    print(message)
    
    # Verify: compute option price using implied volatility
    val = BS_Call_Option_Price(CP, S_0, K, sigma_imp, tau, r)
    print(f'Option Price for implied volatility of {sigma_imp:.7f} is equal to {val:.7f}')
    print(f'Market price was {V_market}, difference: {abs(val - V_market):.10f}')


if __name__ == "__main__":
    main()
