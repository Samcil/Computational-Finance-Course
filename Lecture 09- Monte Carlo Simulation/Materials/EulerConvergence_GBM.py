#%%
"""
Euler Scheme Convergence Analysis for Geometric Brownian Motion.

This module analyzes the convergence properties of the Euler discretization
scheme for GBM by comparing it with the exact solution. Both weak convergence
(convergence of expectations) and strong convergence (pathwise convergence)
are examined.

Created on Jan 20 2019
@author: Lech A. Grzelak
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as st


def GeneratePathsGBMEuler(NoOfPaths, NoOfSteps, T, r, sigma, S_0):
    """
    Generate GBM paths using Euler scheme and exact solution.
    
    The Euler scheme approximates the GBM SDE dS = rS dt + σS dW as:
        S_{i+1} = S_i + r*S_i*Δt + σ*S_i*ΔW_i
    
    The exact solution is:
        S(t) = S_0 * exp((r - 0.5σ²)t + σW(t))
    
    Args:
        NoOfPaths (int): Number of Monte Carlo paths to simulate.
        NoOfSteps (int): Number of time steps for discretization.
        T (float): Time horizon.
        r (float): Drift parameter (risk-free rate).
        sigma (float): Volatility parameter.
        S_0 (float): Initial stock price.
    
    Returns:
        dict: A dictionary containing:
            - 'time' (np.ndarray): Time grid.
            - 'S1' (np.ndarray): Euler approximation paths.
            - 'S2' (np.ndarray): Exact solution paths.
    
    Note:
        Euler scheme has:
        - Weak convergence order: O(Δt)
        - Strong convergence order: O(√Δt)
    """
    Z = np.random.normal(0.0, 1.0, [NoOfPaths, NoOfSteps])
    W = np.zeros([NoOfPaths, NoOfSteps+1])
   
    # Euler approximation
    S1 = np.zeros([NoOfPaths, NoOfSteps+1])
    S1[:, 0] = S_0
    
    # Exact solution
    S2 = np.zeros([NoOfPaths, NoOfSteps+1])
    S2[:, 0] = S_0
    
    time = np.zeros([NoOfSteps+1])
        
    dt = T / float(NoOfSteps)
    for i in range(0, NoOfSteps):
        # Making sure that samples from normal have mean 0 and variance 1
        if NoOfPaths > 1:
            Z[:, i] = (Z[:, i] - np.mean(Z[:, i])) / np.std(Z[:, i])
        W[:, i+1] = W[:, i] + np.power(dt, 0.5) * Z[:, i]
        
        # Euler scheme
        S1[:, i+1] = S1[:, i] + r * S1[:, i] * dt + sigma * S1[:, i] * (W[:, i+1] - W[:, i])
        # Exact solution
        S2[:, i+1] = S2[:, i] * np.exp((r - 0.5 * sigma**2.0) * dt + sigma * (W[:, i+1] - W[:, i]))
        time[i+1] = time[i] + dt
        
    paths = {"time": time, "S1": S1, "S2": S2}
    return paths

def mainCalculation():
    """
    Main calculation routine for analyzing Euler scheme convergence.
    
    This function:
    1. Visualizes sample paths from Euler and exact solutions
    2. Analyzes weak and strong convergence as Δt → 0
    
    Weak convergence: |E[S_Euler(T)] - E[S_Exact(T)]| → 0 as Δt → 0
    Strong convergence: E[|S_Euler(T) - S_Exact(T)|] → 0 as Δt → 0
    """
    NoOfPaths = 25
    NoOfSteps = 25
    T = 1
    r = 0.06
    sigma = 0.3
    S_0 = 50
    
    # Simulated paths - visualization
    Paths = GeneratePathsGBMEuler(NoOfPaths, NoOfSteps, T, r, sigma, S_0)
    timeGrid = Paths["time"]
    S1 = Paths["S1"]
    S2 = Paths["S2"]
    
    plt.figure(1, figsize=(10, 6))
    plt.plot(timeGrid, np.transpose(S1), 'k-', alpha=0.5, linewidth=0.8, label='Euler approximation')   
    plt.plot(timeGrid, np.transpose(S2), '--r', alpha=0.7, linewidth=1.5, label='Exact solution')   
    plt.grid()
    plt.xlabel("Time")
    plt.ylabel("S(t)")
    plt.title("Euler vs Exact GBM Paths")
    plt.legend()
    
    # Convergence analysis
    NoOfStepsV = range(1, 500, 1)
    NoOfPaths = 250
    errorWeak = np.zeros([len(NoOfStepsV), 1])
    errorStrong = np.zeros([len(NoOfStepsV), 1])
    dtV = np.zeros([len(NoOfStepsV), 1])
    
    print("Computing convergence rates...")
    for idx, NoOfSteps in enumerate(NoOfStepsV):
        Paths = GeneratePathsGBMEuler(NoOfPaths, NoOfSteps, T, r, sigma, S_0)
        # Get the paths at maturity T
        S1_atT = Paths["S1"][:, -1]
        S2_atT = Paths["S2"][:, -1]
        
        # Weak error: difference in expectations
        errorWeak[idx] = np.abs(np.mean(S1_atT) - np.mean(S2_atT))
        # Strong error: expectation of absolute difference
        errorStrong[idx] = np.mean(np.abs(S1_atT - S2_atT))
        dtV[idx] = T / NoOfSteps
        
    plt.figure(2, figsize=(10, 6))
    plt.loglog(dtV, errorWeak, 'b-', linewidth=2, label='Weak error ~ O(Δt)')
    plt.loglog(dtV, errorStrong, '--r', linewidth=2, label='Strong error ~ O(√Δt)')
    # Add reference lines
    plt.loglog(dtV, dtV, ':', color='gray', label='O(Δt) reference')
    plt.loglog(dtV, np.sqrt(dtV), ':', color='lightcoral', label='O(√Δt) reference')
    plt.grid(True, which="both", ls="-", alpha=0.3)
    plt.xlabel("Time step Δt")
    plt.ylabel("Error")
    plt.title("Euler Scheme Convergence Analysis")
    plt.legend()
    
    print(f"Final weak error (smallest Δt): {errorWeak[-1][0]:.6e}")
    print(f"Final strong error (smallest Δt): {errorStrong[-1][0]:.6e}")


if __name__ == "__main__":
    mainCalculation()