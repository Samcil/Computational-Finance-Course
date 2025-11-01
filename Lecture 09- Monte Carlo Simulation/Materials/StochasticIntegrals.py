#%%
"""
Stochastic and Riemann Integrals of Brownian Motion.

This module demonstrates two types of integrals involving Brownian motion:
1. Riemann integral: ∫_0^t W(s) ds (deterministic integral of a stochastic process)
2. Stochastic integral: ∫_0^t W(s) dW(s) (Itô integral)

Created on Thu Nov 27 2018
@author: Lech A. Grzelak
"""
import numpy as np
import matplotlib.pyplot as plt


def ComputeIntegrals(NoOfPaths, NoOfSteps, T):
    """
    Compute Brownian motion and its Riemann and stochastic integrals.
    
    This function simulates Brownian motion W(t) and computes:
    - I1(t) = ∫_0^t W(s) ds  (Riemann integral)
    - I2(t) = ∫_0^t W(s) dW(s)  (Itô stochastic integral)
    
    Args:
        NoOfPaths (int): Number of Monte Carlo paths to simulate.
        NoOfSteps (int): Number of time steps for discretization.
        T (float): Time horizon.
    
    Returns:
        dict: A dictionary containing:
            - 'time' (np.ndarray): Time grid of shape (NoOfSteps+1,).
            - 'W' (np.ndarray): Brownian motion paths of shape (NoOfPaths, NoOfSteps+1).
            - 'I1' (np.ndarray): Riemann integral ∫W(s)ds of shape (NoOfPaths, NoOfSteps+1).
            - 'I2' (np.ndarray): Stochastic integral ∫W(s)dW(s) of shape (NoOfPaths, NoOfSteps+1).
    
    Note:
        - I1(t) is normally distributed with E[I1(t)] = 0 and Var[I1(t)] = t³/3
        - I2(t) = (W(t)² - t)/2 (Itô's lemma)
    """
    Z = np.random.normal(0.0, 1.0, [NoOfPaths, NoOfSteps])
    W = np.zeros([NoOfPaths, NoOfSteps+1])
    I1 = np.zeros([NoOfPaths, NoOfSteps+1])
    I2 = np.zeros([NoOfPaths, NoOfSteps+1])
    time = np.zeros([NoOfSteps+1])
    
    dt = T / float(NoOfSteps)
    for i in range(0, NoOfSteps):
        # Making sure that samples from normal have mean 0 and variance 1
        if NoOfPaths > 1:
            Z[:, i] = (Z[:, i] - np.mean(Z[:, i])) / np.std(Z[:, i])
        W[:, i+1] = W[:, i] + np.power(dt, 0.5) * Z[:, i]
        I1[:, i+1] = I1[:, i] + W[:, i] * dt  # Riemann integral
        I2[:, i+1] = I2[:, i] + W[:, i] * (W[:, i+1] - W[:, i])  # Itô integral
        time[i+1] = time[i] + dt
        
    paths = {"time": time, "W": W, "I1": I1, "I2": I2}
    return paths


def main():
    """
    Main calculation routine for visualizing stochastic integrals.
    
    Displays three processes on the same plot:
    - W(t): Brownian motion
    - ∫W(s)ds: Riemann integral of Brownian motion
    - ∫W(s)dW(s): Stochastic integral (should equal (W² - t)/2)
    """
    NoOfPaths = 1
    NoOfSteps = 1000
    T = 1
    
    W_t = ComputeIntegrals(NoOfPaths, NoOfSteps, T)
    timeGrid = W_t["time"]
    Ws = W_t["W"]
    intWsds = W_t["I1"]
    intWsdWs = W_t["I2"]
    
    plt.figure(1, figsize=(10, 6))
    plt.plot(timeGrid, np.transpose(Ws), 'b-', linewidth=2, label='W(t)')
    plt.plot(timeGrid, np.transpose(intWsds), 'r-', linewidth=2, label='∫₀ᵗ W(s)ds')
    plt.plot(timeGrid, np.transpose(intWsdWs), 'k-', linewidth=2, label='∫₀ᵗ W(s)dW(s)')
    plt.grid()
    plt.xlabel("Time")
    plt.ylabel("Value")
    plt.title("Brownian Motion and Its Integrals")
    plt.legend()
    
    # Verify Itô's lemma: ∫W(s)dW(s) = (W(T)² - T)/2
    W_T = Ws[0, -1]
    integral_value = intWsdWs[0, -1]
    expected_value = (W_T**2 - T) / 2
    print(f"Stochastic integral ∫₀ᵀ W(s)dW(s) = {integral_value:.6f}")
    print(f"Expected from Itô's lemma: (W(T)² - T)/2 = {expected_value:.6f}")
    print(f"Difference: {abs(integral_value - expected_value):.6f}")


if __name__ == "__main__":
    main()
