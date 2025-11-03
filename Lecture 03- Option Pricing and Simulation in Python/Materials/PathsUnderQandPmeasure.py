#%%
"""
Monte Carlo Path Generation under P (Physical) and Q (Risk-Neutral) Measures.

This module demonstrates the difference between paths generated under the
physical measure (P-measure) with drift μ and the risk-neutral measure
(Q-measure) with drift r. It visualizes the martingale property under Q.

Created on Thu Nov 28 2018
@author: Lech A. Grzelak
"""
import numpy as np
import matplotlib.pyplot as plt


def GeneratePathsGBM(NoOfPaths, NoOfSteps, T, r, sigma, S_0):
    """
    Generate Geometric Brownian Motion paths with specified drift.
    
    This function simulates stock price paths using GBM with a given drift parameter.
    Under the risk-neutral measure (Q), the drift should be r (risk-free rate).
    Under the physical measure (P), the drift should be μ (expected return).
    
    Args:
        NoOfPaths (int): Number of Monte Carlo paths to simulate.
        NoOfSteps (int): Number of time steps for discretization.
        T (float): Time to maturity (in years).
        r (float): Drift parameter (r for Q-measure, μ for P-measure).
        sigma (float): Volatility parameter (annualized).
        S_0 (float): Initial stock price at time t=0.
    
    Returns:
        dict: A dictionary containing:
            - 'time' (np.ndarray): Time grid of shape (NoOfSteps+1,).
            - 'S' (np.ndarray): GBM paths of shape (NoOfPaths, NoOfSteps+1).
    
    Note:
        The SDE solved is: dS(t) = r*S(t)*dt + sigma*S(t)*dW(t)
    """
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
     
        X[:, i+1] = X[:, i] + (r - 0.5 * sigma * sigma) * dt + sigma * \
            np.power(dt, 0.5) * Z[:, i]
        time[i+1] = time[i] + dt
        
    # Compute exponent of ABM to get GBM
    S = np.exp(X)
    paths = {"time": time, "S": S}
    return paths
    
def MainCode():
    """
    Main calculation demonstrating paths under Q-measure and P-measure.
    
    This function generates stock price paths under two probability measures:
    - Q-measure (risk-neutral): drift = r (risk-free rate)
    - P-measure (physical): drift = μ (expected return)
    
    Under Q-measure, the discounted stock price S(t)/M(t) is a martingale,
    meaning E^Q[S(t)/M(t)] = S_0 (constant over time).
    
    Under P-measure, the discounted stock price grows/declines depending on
    the difference between μ and r.
    """
    NoOfPaths = 8
    NoOfSteps = 1000
    S_0       = 1
    r         = 0.05
    mu        = 0.15
    sigma     = 0.1
    T         = 10
    
    # Money market account (numeraire)
    M = lambda t: np.exp(r * t)
    
    # Generate Monte Carlo Paths under both measures
    pathsQ = GeneratePathsGBM(NoOfPaths, NoOfSteps, T, r, sigma, S_0)
    S_Q = pathsQ["S"]
    pathsP = GeneratePathsGBM(NoOfPaths, NoOfSteps, T, mu, sigma, S_0)
    S_P = pathsP["S"]
    time = pathsQ["time"]    
    
    # Compute discounted stock paths
    S_Qdisc = np.zeros([NoOfPaths, NoOfSteps+1])
    S_Pdisc = np.zeros([NoOfPaths, NoOfSteps+1])
    for i, ti in enumerate(time):
        S_Qdisc[:, i] = S_Q[:, i] / M(ti) 
        S_Pdisc[:, i] = S_P[:, i] / M(ti) 
    
    # Plot Q-measure: S(T)/M(T) with stock growing at rate r
    plt.figure(1)
    plt.grid()
    plt.xlabel("Time")
    plt.ylabel("Discounted Stock Price S(t)/M(t)")
    plt.title("Stock Paths under Q-measure (Risk-Neutral)")
    eSM_Q = lambda t: S_0 * np.exp(r * t) / M(t)
    plt.plot(time, eSM_Q(time), 'r--', linewidth=2)
    plt.plot(time, np.transpose(S_Qdisc), 'blue', alpha=0.6)   
    plt.legend(['E^Q[S(t)/M(t)] = S_0', 'Sample paths S(t)/M(t)'])
    
    # Plot P-measure: S(T)/M(T) with stock growing at rate μ
    plt.figure(2)
    plt.grid()
    plt.xlabel("Time")
    plt.ylabel("Discounted Stock Price S(t)/M(t)")
    plt.title("Stock Paths under P-measure (Physical)")
    eSM_P = lambda t: S_0 * np.exp(mu * t) / M(t)
    plt.plot(time, eSM_P(time), 'r--', linewidth=2)
    plt.plot(time, np.transpose(S_Pdisc), 'blue', alpha=0.6)   
    plt.legend(['E^P[S(t)/M(t)]', 'Sample paths S(t)/M(t)'])


if __name__ == "__main__":
    MainCode()