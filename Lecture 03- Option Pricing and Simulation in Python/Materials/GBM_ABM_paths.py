#%%
"""
Geometric Brownian Motion (GBM) and Arithmetic Brownian Motion (ABM) Path Generation.

This module demonstrates the simulation of GBM and ABM paths using Monte Carlo methods.
GBM is commonly used to model stock prices in quantitative finance.

Created on Thu Nov 27 2018
@author: Lech A. Grzelak
"""
import numpy as np
import matplotlib.pyplot as plt


def GeneratePathsGBMABM(NoOfPaths, NoOfSteps, T, r, sigma, S_0):
    """
    Generate sample paths for Geometric Brownian Motion (GBM) and Arithmetic Brownian Motion (ABM).
    
    This function simulates stock price paths under the risk-neutral measure using the
    Euler-Maruyama discretization scheme. The ABM process X(t) follows:
        dX(t) = (r - 0.5*sigma^2)*dt + sigma*dW(t)
    And the GBM process S(t) is obtained as:
        S(t) = exp(X(t))
    
    Args:
        NoOfPaths (int): Number of Monte Carlo paths to simulate.
        NoOfSteps (int): Number of time steps for discretization.
        T (float): Time to maturity (in years).
        r (float): Risk-free interest rate (annualized).
        sigma (float): Volatility parameter (annualized).
        S_0 (float): Initial stock price at time t=0.
    
    Returns:
        dict: A dictionary containing:
            - 'time' (np.ndarray): Time grid of shape (NoOfSteps+1,).
            - 'X' (np.ndarray): ABM paths of shape (NoOfPaths, NoOfSteps+1).
            - 'S' (np.ndarray): GBM paths of shape (NoOfPaths, NoOfSteps+1).
    
    Note:
        - Random seed is fixed to 1 for reproducibility.
        - Normal samples are standardized to have mean 0 and variance 1.
    """
    # Fixing random seed for reproducibility
    np.random.seed(1)
        
    Z = np.random.normal(0.0, 1.0, [NoOfPaths, NoOfSteps])
    X = np.zeros([NoOfPaths, NoOfSteps+1])
    S = np.zeros([NoOfPaths, NoOfSteps+1])
    time = np.zeros([NoOfSteps+1])
        
    X[:, 0] = np.log(S_0)
    
    dt = T / float(NoOfSteps)
    for i in range(0, NoOfSteps):
        # Making sure that samples from normal have mean 0 and variance 1
        if NoOfPaths > 1:
            Z[:, i] = (Z[:, i] - np.mean(Z[:, i])) / np.std(Z[:, i])
            
        X[:, i+1] = X[:, i] + (r - 0.5 * sigma**2) * dt + sigma * np.power(dt, 0.5) * Z[:, i]
        time[i+1] = time[i] + dt
        
    # Compute exponent of ABM to get GBM
    S = np.exp(X)
    paths = {"time": time, "X": X, "S": S}
    return paths

def mainCalculation():
    """
    Main calculation and visualization routine for GBM and ABM paths.
    
    Simulates paths using predefined parameters and generates plots showing
    both the logarithmic (ABM) and exponential (GBM) processes.
    """
    NoOfPaths = 25
    NoOfSteps = 500
    T = 1
    r = 0.05
    sigma = 0.4
    S_0 = 100
    
    Paths = GeneratePathsGBMABM(NoOfPaths, NoOfSteps, T, r, sigma, S_0)
    timeGrid = Paths["time"]
    X = Paths["X"]
    S = Paths["S"]
    
    plt.figure(1)
    plt.plot(timeGrid, np.transpose(X))   
    plt.grid()
    plt.xlabel("time")
    plt.ylabel("X(t)")
    plt.title("Arithmetic Brownian Motion (ABM) Paths")
    
    plt.figure(2)
    plt.plot(timeGrid, np.transpose(S))   
    plt.grid()
    plt.xlabel("time")
    plt.ylabel("S(t)")
    plt.title("Geometric Brownian Motion (GBM) Paths")


if __name__ == "__main__":
    mainCalculation()