#%%
"""
Stochastic Integral Computation - Integrated Brownian Motion.

This module demonstrates the Monte Carlo simulation of stochastic integrals
of the form: ∫_0^T g(W(t)) dW(t), where W(t) is a Brownian motion.

Created on Thu Nov 27 2018
@author: Lech A. Grzelak
"""
import numpy as np
import matplotlib.pyplot as plt


def ComputeIntegrals(NoOfPaths, NoOfSteps, T, g):
    """
    Compute stochastic integrals using Monte Carlo simulation.
    
    This function simulates Brownian motion paths and computes the stochastic
    integral: I_1(T) = ∫_0^T g(W(t)) dW(t) for each path.
    
    Args:
        NoOfPaths (int): Number of Monte Carlo paths to simulate.
        NoOfSteps (int): Number of time steps for discretization.
        T (float): Time horizon for integration.
        g (callable): Function to integrate against the Brownian motion.
                     Should accept a scalar or array and return same shape.
    
    Returns:
        dict: A dictionary containing:
            - 'time' (np.ndarray): Time grid of shape (NoOfSteps+1,).
            - 'W' (np.ndarray): Brownian motion paths of shape (NoOfPaths, NoOfSteps+1).
            - 'I1' (np.ndarray): Stochastic integral paths of shape (NoOfPaths, NoOfSteps+1).
    
    Note:
        The discretization uses the Euler scheme: I_1(t_{i+1}) = I_1(t_i) + g(W(t_i)) * ΔW(t_i)
    """
    Z = np.random.normal(0.0, 1.0, [NoOfPaths, NoOfSteps])
    W = np.zeros([NoOfPaths, NoOfSteps+1])
    I1 = np.zeros([NoOfPaths, NoOfSteps+1])
    time = np.zeros([NoOfSteps+1])
    
    dt = T / float(NoOfSteps)
    for i in range(0, NoOfSteps):
        # Making sure that samples from normal have mean 0 and variance 1
        if NoOfPaths > 1:
            Z[:, i] = (Z[:, i] - np.mean(Z[:, i])) / np.std(Z[:, i])
        W[:, i+1] = W[:, i] + np.power(dt, 0.5) * Z[:, i]
        
        I1[:, i+1] = I1[:, i] + g(W[:, i]) * (W[:, i+1] - W[:, i])
        time[i+1] = time[i] + dt
        
    paths = {"time": time, "W": W, "I1": I1}
    return paths


def main():
    """
    Main calculation routine for stochastic integral computation.
    
    Computes the stochastic integral ∫_0^T W(t) dW(t) and displays
    its distribution and statistical properties.
    """
    NoOfPaths = 100000
    NoOfSteps = 1000
    T = 2
    
    # Define integrand function g(t) = t
    g = lambda t: t
    
    output = ComputeIntegrals(NoOfPaths, NoOfSteps, T, g)
    timeGrid = output["time"]
    G_T = output["I1"]
    
    plt.figure(1)
    plt.grid()
    plt.hist(G_T[:, -1], 50)
    plt.xlabel("Integral value")
    plt.ylabel("Frequency")
    plt.title(f"Distribution of Stochastic Integral at T={T}")
    
    EX = np.mean(G_T[:, -1])
    Var = np.var(G_T[:, -1])
    print(f'Mean = {EX:.6f} and variance = {Var:.6f}')
    print(f'For the integral ∫_0^T W(t) dW(t):')
    print(f'  Theoretical mean = 0')
    print(f'  Theoretical variance = T^3/3 = {T**3/3:.6f}')


if __name__ == "__main__":
    main()
