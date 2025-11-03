#%%
"""
Correlated Brownian Motion Path Generation.

This module demonstrates the simulation of correlated Brownian motion paths
using the Cholesky decomposition method. Correlated Brownian motions are
essential for modeling multi-asset derivatives and stochastic volatility models.

Created on Feb 09 2019
@author: Lech A. Grzelak
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as st


def GeneratePathsCorrelatedBM(NoOfPaths, NoOfSteps, T, rho):
    """
    Generate correlated Brownian motion paths.
    
    This function generates two correlated Brownian motions W1(t) and W2(t)
    with correlation coefficient ρ using the Cholesky decomposition:
        W1(t) = ∫_0^t Z1(s) ds
        W2(t) = ∫_0^t (ρ*Z1(s) + √(1-ρ²)*Z2(s)) ds
    where Z1 and Z2 are independent standard normal variables.
    
    Args:
        NoOfPaths (int): Number of Monte Carlo paths to simulate.
        NoOfSteps (int): Number of time steps for discretization.
        T (float): Time horizon.
        rho (float): Correlation coefficient between W1 and W2, must be in [-1, 1].
    
    Returns:
        dict: A dictionary containing:
            - 'time' (np.ndarray): Time grid of shape (NoOfSteps+1,).
            - 'W1' (np.ndarray): First Brownian motion paths of shape (NoOfPaths, NoOfSteps+1).
            - 'W2' (np.ndarray): Second Brownian motion paths of shape (NoOfPaths, NoOfSteps+1).
    
    Note:
        - Corr(W1(t), W2(t)) = ρ * t
        - The correlation is induced using: Z2_corr = ρ * Z1 + √(1-ρ²) * Z2
    """
    Z1 = np.random.normal(0.0, 1.0, [NoOfPaths, NoOfSteps])
    Z2 = np.random.normal(0.0, 1.0, [NoOfPaths, NoOfSteps])
    W1 = np.zeros([NoOfPaths, NoOfSteps+1])
    W2 = np.zeros([NoOfPaths, NoOfSteps+1])    
    
    dt = T / float(NoOfSteps)
    time = np.zeros([NoOfSteps+1])
    for i in range(0, NoOfSteps):
        # Making sure that samples from normal have mean 0 and variance 1
        if NoOfPaths > 1:
            Z1[:, i] = (Z1[:, i] - np.mean(Z1[:, i])) / np.std(Z1[:, i])
            Z2[:, i] = (Z2[:, i] - np.mean(Z2[:, i])) / np.std(Z2[:, i])
        
        # Correlate noises using Cholesky decomposition
        Z2[:, i] = rho * Z1[:, i] + np.sqrt(1.0 - rho**2) * Z2[:, i]
        
        W1[:, i+1] = W1[:, i] + np.power(dt, 0.5) * Z1[:, i]
        W2[:, i+1] = W2[:, i] + np.power(dt, 0.5) * Z2[:, i]
        
        time[i+1] = time[i] + dt
        
    # Store the results
    paths = {"time": time, "W1": W1, "W2": W2}
    return paths

def mainCalculation():
    """
    Main calculation routine demonstrating correlated Brownian motions.
    
    This function generates and plots correlated Brownian motion pairs
    with three different correlation coefficients:
    - Negative correlation (ρ = -0.9): paths move in opposite directions
    - Positive correlation (ρ = 0.9): paths move together
    - Zero correlation (ρ = 0): paths are independent
    """
    NoOfPaths = 1
    NoOfSteps = 500
    T = 1.0
    
    # Negative correlation: paths tend to move in opposite directions
    rho = -0.9
    Paths = GeneratePathsCorrelatedBM(NoOfPaths, NoOfSteps, T, rho)
    timeGrid = Paths["time"]
    W1 = Paths["W1"]
    W2 = Paths["W2"]
    
    plt.figure(1)
    plt.plot(timeGrid, np.transpose(W1), label='W1(t)')   
    plt.plot(timeGrid, np.transpose(W2), label='W2(t)')   
    plt.grid()
    plt.xlabel("Time")
    plt.ylabel("W(t)")
    plt.title(f"Correlated Brownian Motions (ρ = {rho})")
    plt.legend()
    
    # Positive correlation: paths tend to move together
    rho = 0.9
    Paths = GeneratePathsCorrelatedBM(NoOfPaths, NoOfSteps, T, rho)
    timeGrid = Paths["time"]
    W1 = Paths["W1"]
    W2 = Paths["W2"]
    
    plt.figure(2)
    plt.plot(timeGrid, np.transpose(W1), label='W1(t)')   
    plt.plot(timeGrid, np.transpose(W2), label='W2(t)')   
    plt.grid()
    plt.xlabel("Time")
    plt.ylabel("W(t)")
    plt.title(f"Correlated Brownian Motions (ρ = {rho})")
    plt.legend()
    
    # Zero correlation: paths are independent
    rho = 0.0
    Paths = GeneratePathsCorrelatedBM(NoOfPaths, NoOfSteps, T, rho)
    timeGrid = Paths["time"]
    W1 = Paths["W1"]
    W2 = Paths["W2"]
    
    plt.figure(3)
    plt.plot(timeGrid, np.transpose(W1), label='W1(t)')   
    plt.plot(timeGrid, np.transpose(W2), label='W2(t)')   
    plt.grid()
    plt.xlabel("Time")
    plt.ylabel("W(t)")
    plt.title(f"Correlated Brownian Motions (ρ = {rho})")
    plt.legend()


if __name__ == "__main__":
    mainCalculation()