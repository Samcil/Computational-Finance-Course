#%%
"""
Barrier Option Pricing and Variance Reduction.

This module demonstrates the pricing of up-and-out barrier options using Monte
Carlo simulation. Barrier options are path-dependent derivatives that are
knocked out (become worthless) if the underlying asset crosses a barrier level.

Created on Thu Jan 16 2019
@author: Lech A. Grzelak
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as st
import enum 


def DigitalPayoffValuation(S, T, r, payoff):
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


def UpAndOutBarrier(S, T, r, payoff, Su):
    """
    Price an up-and-out barrier option using Monte Carlo simulation.
    
    An up-and-out barrier option is knocked out (becomes worthless) if the
    underlying asset price ever exceeds the barrier level Su during the
    option's life.
    
    Args:
        S (np.ndarray): Stock price paths of shape (NoOfPaths, NoOfSteps+1).
        T (float): Time to maturity.
        r (float): Risk-free rate.
        payoff (callable): Payoff function.
        Su (float): Upper barrier level.
    
    Returns:
        float: Option value at t=0.
    
    Note:
        The option pays off only if S(t) ≤ Su for all t ∈ [0,T].
        Barrier options are cheaper than vanilla options due to the
        knock-out feature.
    """
    # Check which paths hit the barrier
    n1, n2 = S.shape
    barrier = np.zeros([n1, n2]) + Su
    
    hitM = S > barrier
    hitVec = np.sum(hitM, 1)  # Count barrier hits per path
    hitVec = (hitVec == 0.0).astype(int)  # 1 if no hit, 0 if hit
    
    # Payoff is zero if barrier was hit
    V_0 = np.exp(-r * T) * np.mean(payoff(S[:, -1] * hitVec))
    
    return V_0


def mainCalculation():
    """
    Main calculation comparing vanilla and barrier option prices.
    
    Demonstrates:
    1. Vanilla call option pricing
    2. Up-and-out barrier call option pricing
    3. Variance reduction from barrier feature
    """
    NoOfPaths = 10000
    NoOfSteps = 250
   
    S0 = 100.0
    r = 0.05
    T = 5
    sigma = 0.2
    Su = 150  # Barrier level
    
    paths = GeneratePathsGBMEuler(NoOfPaths, NoOfSteps, T, r, sigma, S0)
    S_paths = paths["S"]
    S_T = S_paths[:, -1]
    
    # Strike prices
    K = 100.0
    K2 = 140.0
    
    # Call option payoff
    payoff = lambda S: np.maximum(S - K, 0.0)
    
    # Visualize payoff
    S_T_grid = np.linspace(50, S0 * 1.5, 200)
    
    plt.figure(1, figsize=(10, 6))
    plt.plot(S_T_grid, payoff(S_T_grid), 'b-', linewidth=2)
    plt.axvline(x=Su, color='r', linestyle='--', linewidth=2, label=f'Barrier Su={Su}')
    plt.axvline(x=K, color='g', linestyle=':', linewidth=1.5, label=f'Strike K={K}')
    plt.grid(True, alpha=0.3)
    plt.xlabel('Stock Price S(T)')
    plt.ylabel('Payoff')
    plt.title('Call Option Payoff')
    plt.legend()
    
    # Vanilla option valuation
    val_t0 = DigitalPayoffValuation(S_T, T, r, payoff)
    print(f"\n=== Option Pricing Results ===")
    print(f"Vanilla Call Option: {val_t0:.4f}")
    
    # Barrier option valuation
    barrier_price = UpAndOutBarrier(S_paths, T, r, payoff, Su)
    print(f"Up-and-Out Barrier Call Option: {barrier_price:.4f}")
    
    # Analysis
    discount = ((val_t0 - barrier_price) / val_t0) * 100
    print(f"\n=== Barrier Impact ===")
    print(f"Price reduction from barrier: {discount:.2f}%")
    print(f"Barrier level: Su = {Su} (S0 = {S0})")
    
    # Count paths that hit barrier
    hit_count = np.sum(np.any(S_paths > Su, axis=1))
    hit_pct = (hit_count / NoOfPaths) * 100
    print(f"Paths hitting barrier: {hit_count}/{NoOfPaths} ({hit_pct:.2f}%)")


if __name__ == "__main__":
    mainCalculation()
