#%%
"""
Delta Hedging with the Black-Scholes Model.

This module demonstrates dynamic delta hedging strategies for European options
using the Black-Scholes framework. It simulates the hedging process by continuously
rebalancing a portfolio to maintain delta neutrality.

Created on Thu Dec 12 2018
@author: Lech A. Grzelak
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as st
import enum 
from mpl_toolkits import mplot3d
from scipy.interpolate import RegularGridInterpolator


class OptionType(enum.Enum):
    """Enumeration for option types."""
    CALL = 1.0
    PUT = -1.0

def GeneratePathsGBM(NoOfPaths, NoOfSteps, T, r, sigma, S_0):
    """
    Generate Geometric Brownian Motion paths for stock prices.
    
    Args:
        NoOfPaths (int): Number of Monte Carlo paths.
        NoOfSteps (int): Number of time steps.
        T (float): Time horizon.
        r (float): Risk-free rate.
        sigma (float): Volatility.
        S_0 (float): Initial stock price.
    
    Returns:
        dict: Dictionary with 'time' and 'S' (stock price paths).
    """
    Z = np.random.normal(0.0, 1.0, [NoOfPaths, NoOfSteps])
    X = np.zeros([NoOfPaths, NoOfSteps+1])
    W = np.zeros([NoOfPaths, NoOfSteps+1])
    time = np.zeros([NoOfSteps+1])
        
    X[:, 0] = np.log(S_0)
    
    dt = T / float(NoOfSteps)
    for i in range(0, NoOfSteps):
        # Making sure that samples from normal have mean 0 and variance 1
        if NoOfPaths > 1:
            Z[:, i] = (Z[:, i] - np.mean(Z[:, i])) / np.std(Z[:, i])
        W[:, i+1] = W[:, i] + np.power(dt, 0.5) * Z[:, i]
        X[:, i+1] = X[:, i] + (r - 0.5 * sigma * sigma) * dt + sigma * (W[:, i+1] - W[:, i])
        time[i+1] = time[i] + dt
        
    S = np.exp(X)
    paths = {"time": time, "S": S}
    return paths


def BS_Call_Put_Option_Price(CP, S_0, K, sigma, t, T, r):
    """
    Black-Scholes option pricing formula.
    
    Args:
        CP (OptionType): CALL or PUT.
        S_0 (float): Current stock price.
        K (list/array): Strike prices.
        sigma (float): Volatility.
        t (float): Current time.
        T (float): Maturity time.
        r (float): Risk-free rate.
    
    Returns:
        np.ndarray: Option prices.
    """
    K = np.array(K).reshape([len(K), 1])
    d1 = (np.log(S_0 / K) + (r + 0.5 * np.power(sigma, 2.0)) * 
          (T - t)) / (sigma * np.sqrt(T - t))
    d2 = d1 - sigma * np.sqrt(T - t)
    if CP == OptionType.CALL:
        value = st.norm.cdf(d1) * S_0 - st.norm.cdf(d2) * K * np.exp(-r * (T - t))
    elif CP == OptionType.PUT:
        value = st.norm.cdf(-d2) * K * np.exp(-r * (T - t)) - st.norm.cdf(-d1) * S_0
    return value


def BS_Delta(CP, S_0, K, sigma, t, T, r):
    """
    Calculate Black-Scholes delta (∂V/∂S).
    
    Delta measures the rate of change of option price with respect to
    the underlying asset price. It's used for hedging to create a
    delta-neutral portfolio.
    
    Args:
        CP (OptionType): CALL or PUT.
        S_0 (float): Current stock price.
        K (list/array): Strike prices.
        sigma (float): Volatility.
        t (float): Current time.
        T (float): Maturity time.
        r (float): Risk-free rate.
    
    Returns:
        np.ndarray: Delta values.
    """
    # Handle edge case where time grid slightly exceeds maturity
    if t - T > 10e-20 and T - t < 10e-7:
        t = T
    K = np.array(K).reshape([len(K), 1])
    d1 = (np.log(S_0 / K) + (r + 0.5 * np.power(sigma, 2.0)) * 
          (T - t)) / (sigma * np.sqrt(T - t))
    if CP == OptionType.CALL:
        value = st.norm.cdf(d1)
    elif CP == OptionType.PUT:
        value = st.norm.cdf(d1) - 1.0
    return value


def mainCalculation():
    """
    Main calculation for delta hedging simulation.
    
    This function simulates a delta hedging strategy where we:
    1. Sell a call option
    2. Maintain a delta-neutral portfolio by holding Δ units of stock
    3. Rebalance continuously
    4. At maturity, deliver the payoff and unwind the hedge
    
    In a perfect Black-Scholes world with continuous hedging, the final
    P&L should be zero (self-financing strategy).
    """
    NoOfPaths = 5000
    NoOfSteps = 1000

    T = 1.0
    r = 0.1
    sigma = 0.2
    s0 = 1.0
    K = [0.95]
    CP = OptionType.CALL
    
    np.random.seed(1)
    Paths = GeneratePathsGBM(NoOfPaths, NoOfSteps, T, r, sigma, s0)
    time = Paths["time"]
    S = Paths["S"]
    
    # Lambda functions for option price and delta
    C = lambda t, K, S0: BS_Call_Put_Option_Price(CP, S0, K, sigma, t, T, r)
    Delta = lambda t, K, S0: BS_Delta(CP, S0, K, sigma, t, T, r)
    
    # Initialize portfolio: sell option, buy delta shares, invest/borrow difference
    PnL = np.zeros([NoOfPaths, NoOfSteps+1])
    delta_init = Delta(0.0, K, s0)
    PnL[:, 0] = C(0.0, K, s0) - delta_init * s0  # Initial cash position
            
    CallM = np.zeros([NoOfPaths, NoOfSteps+1])
    CallM[:, 0] = C(0.0, K, s0)
    DeltaM = np.zeros([NoOfPaths, NoOfSteps+1])
    DeltaM[:, 0] = Delta(0, K, s0)
    
    # Dynamic hedging: rebalance at each time step
    for i in range(1, NoOfSteps+1):
        dt = time[i] - time[i-1]
        delta_old = Delta(time[i-1], K, S[:, i-1])
        delta_curr = Delta(time[i], K, S[:, i])
        
        # Update P&L: accrue interest and subtract cost of rebalancing
        PnL[:, i] = PnL[:, i-1] * np.exp(r * dt) - (delta_curr - delta_old) * S[:, i]
        CallM[:, i] = C(time[i], K, S[:, i])
        DeltaM[:, i] = delta_curr
    
    # Final settlement: pay off option if ITM and sell remaining hedge
    PnL[:, -1] = PnL[:, -1] - np.maximum(S[:, -1] - K, 0) + DeltaM[:, -1] * S[:, -1]
    
    # Visualization of one path
    path_id = 13
    plt.figure(1, figsize=(12, 6))
    plt.plot(time, S[path_id, :], label='Stock Price', linewidth=2)
    plt.plot(time, CallM[path_id, :], label='Call Price', linewidth=2)
    plt.plot(time, DeltaM[path_id, :], label='Delta', linewidth=2)
    plt.plot(time, PnL[path_id, :], label='P&L', linewidth=2)
    plt.legend()
    plt.grid()
    plt.xlabel('Time')
    plt.ylabel('Value')
    plt.title(f'Delta Hedging Simulation (Path {path_id})')
    
    # Histogram of final P&L
    plt.figure(2, figsize=(10, 6))
    plt.hist(PnL[:, -1], 50, edgecolor='black')
    plt.grid(True, alpha=0.3)
    plt.xlim([-0.1, 0.1])
    plt.xlabel('Final P&L')
    plt.ylabel('Frequency')
    plt.title('Distribution of Hedging P&L at Maturity')
    plt.axvline(x=0, color='r', linestyle='--', linewidth=2, label='Perfect hedge')
    plt.legend()
    
    # Summary statistics
    print(f"\nHedging Performance Summary:")
    print(f"Mean P&L: {np.mean(PnL[:, -1]):.6f}")
    print(f"Std Dev P&L: {np.std(PnL[:, -1]):.6f}")
    print(f"Min P&L: {np.min(PnL[:, -1]):.6f}")
    print(f"Max P&L: {np.max(PnL[:, -1]):.6f}")
    print(f"\nNote: With continuous hedging in BS world, P&L should be ~0")


if __name__ == "__main__":
    mainCalculation()