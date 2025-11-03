#%%
"""
CIR Process with Different Boundary Conditions - Euler Discretization.

This module compares two different boundary condition treatments for the CIR
process when using Euler discretization. The boundary conditions prevent the
variance process from becoming negative.

Created on Jan 20 2019
@author: Lech A. Grzelak
"""
import numpy as np
import matplotlib.pyplot as plt


def GeneratePathsCIREuler2Schemes(NoOfPaths, NoOfSteps, T, kappa, v0, vbar, gamma):
    """
    Generate CIR paths using Euler scheme with two boundary condition treatments.
    
    The CIR process: dV(t) = κ(v̄ - V(t))dt + γ√V(t)dW(t)
    
    When the Feller condition (2κv̄ ≥ γ²) is not satisfied, the Euler scheme
    may produce negative values. Two boundary conditions are implemented:
    
    1. Truncation: V(t) = max(V(t), 0) - set negative values to zero
    2. Reflection: V(t) = |V(t)| - take absolute value
    
    Args:
        NoOfPaths (int): Number of Monte Carlo paths.
        NoOfSteps (int): Number of time steps.
        T (float): Time horizon.
        kappa (float): Mean reversion speed.
        v0 (float): Initial variance.
        vbar (float): Long-term mean variance.
        gamma (float): Volatility of volatility.
    
    Returns:
        dict: Dictionary with 'time', 'Vtruncated', and 'Vreflected' arrays.
    
    Note:
        - Truncation scheme: Simple but may introduce bias
        - Reflection scheme: Preserves symmetry but may not be theoretically justified
        - For accurate simulation when Feller condition fails, use exact methods
    """
    Z = np.random.normal(0.0, 1.0, [NoOfPaths, NoOfSteps])
    W = np.zeros([NoOfPaths, NoOfSteps+1])
    V1 = np.zeros([NoOfPaths, NoOfSteps+1])
    V2 = np.zeros([NoOfPaths, NoOfSteps+1])
    V1[:, 0] = v0
    V2[:, 0] = v0
    time = np.zeros([NoOfSteps+1])
        
    dt = T / float(NoOfSteps)
    for i in range(0, NoOfSteps):
        # Making sure that samples from normal have mean 0 and variance 1
        if NoOfPaths > 1:
            Z[:, i] = (Z[:, i] - np.mean(Z[:, i])) / np.std(Z[:, i])
        W[:, i+1] = W[:, i] + np.power(dt, 0.5) * Z[:, i]
        
        # Truncated boundary condition: max(V, 0)
        V1[:, i+1] = V1[:, i] + kappa * (vbar - V1[:, i]) * dt + gamma * np.sqrt(V1[:, i]) * (W[:, i+1] - W[:, i])
        V1[:, i+1] = np.maximum(V1[:, i+1], 0.0)
        
        # Reflecting boundary condition: |V|
        V2[:, i+1] = V2[:, i] + kappa * (vbar - V2[:, i]) * dt + gamma * np.sqrt(V2[:, i]) * (W[:, i+1] - W[:, i])
        V2[:, i+1] = np.absolute(V2[:, i+1])
        time[i+1] = time[i] + dt
        
    paths = {"time": time, "Vtruncated": V1, "Vreflected": V2}
    return paths


def mainCalculation():
    """
    Main calculation demonstrating different boundary condition schemes.
    
    Uses deliberately small number of steps to show the effect of boundary
    conditions when the discretized process attempts to go negative.
    """
    NoOfPaths = 1
    NoOfSteps = 20
    T = 1
    kappa = 0.5
    v0 = 0.1
    vbar = 0.1
    gamma = 0.8
    
    # Check Feller condition
    feller = 2 * kappa * vbar
    print(f"=== CIR Boundary Condition Analysis ===")
    print(f"Feller condition: 2κv̄ = {feller:.4f} vs γ² = {gamma**2:.4f}")
    if feller < gamma**2:
        print("⚠ Feller condition NOT satisfied - boundary issues expected")
        print("  Process may hit zero with Euler discretization")
    else:
        print("✓ Feller condition satisfied")
    
    np.random.seed(210)
    Paths = GeneratePathsCIREuler2Schemes(NoOfPaths, NoOfSteps, T, kappa, v0, vbar, gamma)
    timeGrid = Paths["time"]
    V_truncated = Paths["Vtruncated"]
    V_reflected = Paths["Vreflected"]
    
    plt.figure(1, figsize=(10, 6))
    plt.plot(timeGrid, np.transpose(V_truncated), 'b-o', linewidth=2, markersize=6, label='Truncated: max(V, 0)')
    plt.plot(timeGrid, np.transpose(V_reflected), '--r^', linewidth=2, markersize=6, label='Reflected: |V|')
    plt.axhline(y=vbar, color='g', linestyle=':', linewidth=1.5, label=f'Long-term mean v̄={vbar}')
    plt.axhline(y=0, color='k', linestyle='-', linewidth=0.5, alpha=0.5)
    plt.grid(True, alpha=0.3)
    plt.xlabel("Time")
    plt.ylabel("V(t)")
    plt.title(f"CIR Euler Schemes with Different Boundary Conditions\n(κ={kappa}, γ={gamma}, v̄={vbar}, steps={NoOfSteps})")
    plt.legend()
    
    print(f"\nMinimum values before boundary treatment would show negative values")
    print("Both schemes ensure V(t) ≥ 0 but use different approaches")


if __name__ == "__main__":
    mainCalculation()
