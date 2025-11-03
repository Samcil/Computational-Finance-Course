#%%
"""
Asian Option Pricing using Monte Carlo Simulation.

This module demonstrates the pricing of Asian options, which are path-dependent
derivatives where the payoff depends on the average price of the underlying
asset over a period. Asian options have lower variance in payoff compared to
vanilla options, making them valuable for hedging and speculation.

Created on Thu Jan 16 2019
@author: Lech A. Grzelak
"""
import numpy as np


def PayoffValuation(S, T, r, payoff):
    """
    Calculate the present value of a payoff using Monte Carlo samples.
    
    Args:
        S (np.ndarray): Monte Carlo samples at time T.
        T (float): Time to maturity.
        r (float): Risk-free rate.
        payoff (callable): Payoff function.
    
    Returns:
        float: Discounted expected payoff.
    """
    return np.exp(-r * T) * np.mean(payoff(S))


def GeneratePathsGBMEuler(NoOfPaths, NoOfSteps, T, r, sigma, S_0):
    """
    Generate GBM paths using Euler discretization.
    
    Args:
        NoOfPaths (int): Number of Monte Carlo paths.
        NoOfSteps (int): Number of time steps.
        T (float): Time horizon.
        r (float): Risk-free rate.
        sigma (float): Volatility.
        S_0 (float): Initial stock price.
    
    Returns:
        dict: Dictionary with 'time' and 'S' arrays.
    """
    Z = np.random.normal(0.0, 1.0, [NoOfPaths, NoOfSteps])
    W = np.zeros([NoOfPaths, NoOfSteps+1])
   
    # Euler Approximation
    S1 = np.zeros([NoOfPaths, NoOfSteps+1])
    S1[:, 0] = S_0
    
    time = np.zeros([NoOfSteps+1])
        
    dt = T / float(NoOfSteps)
    for i in range(0, NoOfSteps):
        # Making sure that samples from normal have mean 0 and variance 1
        if NoOfPaths > 1:
            Z[:, i] = (Z[:, i] - np.mean(Z[:, i])) / np.std(Z[:, i])
        W[:, i+1] = W[:, i] + np.power(dt, 0.5) * Z[:, i]
        
        S1[:, i+1] = S1[:, i] + r * S1[:, i] * dt + sigma * S1[:, i] * (W[:, i+1] - W[:, i])
        time[i+1] = time[i] + dt
        
    paths = {"time": time, "S": S1}
    return paths


def mainCalculation():
    """
    Main calculation comparing vanilla European vs Asian option pricing.
    
    Asian options (average price options) have payoffs based on the arithmetic
    average of the asset price over the option's life:
        Payoff = max(A(T) - K, 0)
    where A(T) = (1/T) ∫_0^T S(t) dt
    
    Asian options typically have lower prices than vanilla options because
    averaging reduces volatility.
    """
    NoOfPaths = 5000
    NoOfSteps = 250
   
    S0 = 100.0
    r = 0.05
    T = 5
    sigma = 0.2
     
    paths = GeneratePathsGBMEuler(NoOfPaths, NoOfSteps, T, r, sigma, S0)
    S_paths = paths["S"]
    S_T = S_paths[:, -1]  # Terminal stock prices
    
    # Strike price
    K = 100.0
    
    # Call option payoff
    payoff = lambda S: np.maximum(S - K, 0.0)  
        
    # Vanilla European call option
    val_t0 = PayoffValuation(S_T, T, r, payoff)
    print(f"\n=== Option Pricing Results ===")
    print(f"European Call Option value: {val_t0:.4f}")
    
    # Asian option: payoff based on arithmetic average
    A_T = np.mean(S_paths, axis=1)  # Arithmetic average along each path
    valAsian_t0 = PayoffValuation(A_T, T, r, payoff)
    print(f"Asian Call Option value: {valAsian_t0:.4f}")
    
    # Variance comparison
    var_ST = np.var(S_T)
    var_AT = np.var(A_T)
    print(f"\n=== Variance Analysis ===")
    print(f"Variance of S(T): {var_ST:.4f}")
    print(f"Variance of A(T): {var_AT:.4f}")
    print(f"Variance reduction: {(1 - var_AT/var_ST) * 100:.2f}%")
    
    print(f"\n=== Key Insight ===")
    print(f"Asian option ({valAsian_t0:.4f}) < European option ({val_t0:.4f})")
    print(f"Price reduction: {((val_t0 - valAsian_t0) / val_t0 * 100):.2f}%")
    print("This is because averaging reduces the volatility of the payoff.")


if __name__ == "__main__":
    mainCalculation()