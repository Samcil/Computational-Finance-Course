#%%
"""
Poisson Process and Compensated Poisson Process Path Generation.

This module demonstrates the simulation of Poisson process paths and their
compensated versions. The compensated Poisson process is a martingale obtained
by subtracting the deterministic drift λt from the Poisson process N(t).

Created on Thu Jan 11 2019
@author: Lech A. Grzelak
"""
import numpy as np
import matplotlib.pyplot as plt


def GeneratePathsPoisson(NoOfPaths, NoOfSteps, T, xiP):
    """
    Generate paths for Poisson process and compensated Poisson process.
    
    The Poisson process N(t) counts the number of jumps up to time t with
    intensity λ (xiP). The compensated process is defined as:
        Ñ(t) = N(t) - λt
    which is a martingale.
    
    Args:
        NoOfPaths (int): Number of Monte Carlo paths to simulate.
        NoOfSteps (int): Number of time steps for discretization.
        T (float): Time horizon.
        xiP (float): Intensity parameter λ of the Poisson process (jumps per unit time).
    
    Returns:
        dict: A dictionary containing:
            - 'time' (np.ndarray): Time grid of shape (NoOfSteps+1,).
            - 'X' (np.ndarray): Poisson process paths of shape (NoOfPaths, NoOfSteps+1).
            - 'Xcomp' (np.ndarray): Compensated Poisson process paths of shape (NoOfPaths, NoOfSteps+1).
    
    Note:
        - E[N(t)] = λt
        - Var[N(t)] = λt
        - E[Ñ(t)] = 0 (compensated process is a martingale)
    """
    # Create empty matrices for Poisson process and compensated Poisson process
    X = np.zeros([NoOfPaths, NoOfSteps+1])
    Xc = np.zeros([NoOfPaths, NoOfSteps+1])
    time = np.zeros([NoOfSteps+1])
                
    dt = T / float(NoOfSteps)
    
    # Generate Poisson increments
    Z = np.random.poisson(xiP * dt, [NoOfPaths, NoOfSteps])
    
    for i in range(0, NoOfSteps):
        # Standard Poisson process: cumulative sum of jumps
        X[:, i+1] = X[:, i] + Z[:, i]
        # Compensated Poisson process: subtract drift λ*dt
        Xc[:, i+1] = Xc[:, i] - xiP * dt + Z[:, i]
        time[i+1] = time[i] + dt
        
    paths = {"time": time, "X": X, "Xcomp": Xc}
    return paths

def mainCalculation():
    """
    Main calculation routine for Poisson process visualization.
    
    Generates and plots sample paths of both the standard Poisson process
    and the compensated Poisson process. The compensated process fluctuates
    around zero (martingale property), while the standard process grows
    linearly on average.
    """
    NoOfPaths = 25
    NoOfSteps = 500
    T = 30
    xiP = 1  # Intensity: 1 jump per unit time on average
        
    Paths = GeneratePathsPoisson(NoOfPaths, NoOfSteps, T, xiP)
    timeGrid = Paths["time"]
    X = Paths["X"]
    Xc = Paths["Xcomp"]
    
    # Plot standard Poisson process
    plt.figure(1)
    plt.plot(timeGrid, np.transpose(X), '-b', alpha=0.6)
    plt.plot(timeGrid, xiP * timeGrid, 'r--', linewidth=2)
    plt.grid()
    plt.xlabel("Time")
    plt.ylabel("N(t)")
    plt.title(f"Poisson Process (λ = {xiP})")
    plt.legend(['Sample paths', f'Expected value λt = {xiP}t'])
    
    # Plot compensated Poisson process
    plt.figure(2)
    plt.plot(timeGrid, np.transpose(Xc), '-b', alpha=0.6)
    plt.axhline(y=0, color='r', linestyle='--', linewidth=2)
    plt.grid()
    plt.xlabel("Time")
    plt.ylabel("Ñ(t)")
    plt.title(f"Compensated Poisson Process Ñ(t) = N(t) - λt")
    plt.legend(['Sample paths', 'Expected value = 0'])


if __name__ == "__main__":
    mainCalculation()