#%%
"""
Cox-Ingersoll-Ross (CIR) Process - Exact Simulation Method.

This module implements exact simulation of the CIR process using the noncentral
chi-squared distribution. The CIR process is commonly used to model stochastic
volatility and interest rates.

Created on Jan 20 2019
@author: Lech A. Grzelak
"""
import numpy as np
import matplotlib.pyplot as plt


def CIR_Sample(NoOfPaths, kappa, gamma, vbar, s, t, v_s):
    """
    Generate exact samples from the CIR process using noncentral chi-squared distribution.
    
    The CIR process satisfies: dV(t) = κ(v̄ - V(t))dt + γ√V(t)dW(t)
    
    The exact distribution of V(t) | V(s) is known and can be sampled using the
    noncentral chi-squared distribution with appropriate parameters.
    
    Args:
        NoOfPaths (int): Number of samples to generate.
        kappa (float): Mean reversion speed (κ).
        gamma (float): Volatility of volatility (γ).
        vbar (float): Long-term mean (v̄).
        s (float): Starting time.
        t (float): Ending time.
        v_s (np.ndarray or float): Value(s) at time s.
    
    Returns:
        np.ndarray: Samples of V(t) given V(s).
    
    Note:
        The transformation uses:
        - δ = 4κv̄/γ² (degrees of freedom)
        - c = γ²(1-exp(-κ(t-s)))/(4κ)
        - λ̄ = 4κV(s)exp(-κ(t-s))/(γ²(1-exp(-κ(t-s)))) (non-centrality parameter)
    """
    delta = 4.0 * kappa * vbar / gamma / gamma
    c = 1.0 / (4.0 * kappa) * gamma * gamma * (1.0 - np.exp(-kappa * (t - s)))
    kappaBar = 4.0 * kappa * v_s * np.exp(-kappa * (t - s)) / (gamma * gamma * (1.0 - np.exp(-kappa * (t - s))))
    sample = c * np.random.noncentral_chisquare(delta, kappaBar, NoOfPaths)
    return sample


def GeneratePathsCIRExact(NoOfPaths, NoOfSteps, T, kappa, v0, vbar, gamma):
    """
    Generate CIR process paths using exact simulation method.
    
    This function uses the exact distributional properties of the CIR process
    to generate sample paths without discretization bias.
    
    Args:
        NoOfPaths (int): Number of Monte Carlo paths.
        NoOfSteps (int): Number of time steps.
        T (float): Time horizon.
        kappa (float): Mean reversion speed.
        v0 (float): Initial variance.
        vbar (float): Long-term mean variance.
        gamma (float): Volatility of volatility.
    
    Returns:
        dict: Dictionary with 'time' and 'VExact' arrays.
    
    Note:
        The Feller condition 2κv̄ ≥ γ² ensures the process stays positive.
    """
    Z = np.random.normal(0.0, 1.0, [NoOfPaths, NoOfSteps])
    W = np.zeros([NoOfPaths, NoOfSteps+1])
    V = np.zeros([NoOfPaths, NoOfSteps+1])
    V[:, 0] = v0
    
    time = np.zeros([NoOfSteps+1])
        
    dt = T / float(NoOfSteps)
    for i in range(0, NoOfSteps):
        # Making sure that samples from normal have mean 0 and variance 1
        if NoOfPaths > 1:
            Z[:, i] = (Z[:, i] - np.mean(Z[:, i])) / np.std(Z[:, i])
        W[:, i+1] = W[:, i] + np.power(dt, 0.5) * Z[:, i]
        
        # Exact simulation using noncentral chi-squared distribution
        V[:, i+1] = CIR_Sample(NoOfPaths, kappa, gamma, vbar, 0, dt, V[:, i])
                       
        time[i+1] = time[i] + dt
        
    paths = {"time": time, "VExact": V}
    return paths


def mainCalculation():
    """
    Main calculation routine for CIR exact simulation demonstration.
    
    Generates and visualizes a sample path of the CIR process using exact
    simulation. This method provides unbiased samples at each time point.
    """
    NoOfPaths = 1
    NoOfSteps = 2000
    T = 1
    kappa = 0.7
    v0 = 0.1
    vbar = 0.1
    gamma = 0.7

    # Check Feller condition
    feller = 2 * kappa * vbar
    print(f"Feller condition check: 2κv̄ = {feller:.4f}, γ² = {gamma**2:.4f}")
    if feller >= gamma**2:
        print("✓ Feller condition satisfied: process stays positive")
    else:
        print("⚠ Feller condition not satisfied: process may hit zero")

    np.random.seed(10)
    Paths = GeneratePathsCIRExact(NoOfPaths, NoOfSteps, T, kappa, v0, vbar, gamma)
    timeGrid = Paths["time"]
    V_exact = Paths["VExact"]
    
    plt.figure(1, figsize=(10, 6))
    plt.plot(timeGrid, np.transpose(V_exact), 'b-', linewidth=2)
    plt.axhline(y=vbar, color='r', linestyle='--', label=f'Long-term mean v̄={vbar}')
    plt.xlabel("Time")
    plt.ylabel("V(t)")
    plt.title(f"CIR Process Exact Simulation (κ={kappa}, γ={gamma}, v̄={vbar})")
    plt.legend(['V(t) - exact simulation', f'Long-term mean v̄={vbar}'])
    plt.grid(True, alpha=0.3)


if __name__ == "__main__":
    mainCalculation()
