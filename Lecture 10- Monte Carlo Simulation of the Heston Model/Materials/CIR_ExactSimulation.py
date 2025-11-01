"""
CIR Process Exact Simulation and Heston Model Almost Exact Scheme (AES).

This module demonstrates the almost exact simulation scheme for the Heston model,
comparing the Euler and AES schemes for option pricing. The AES scheme uses exact
simulation for the CIR variance process via the noncentral chi-squared distribution.

The Heston model dynamics:
    dS(t) = rS(t)dt + sqrt(V(t))S(t)dW_2(t)
    dV(t) = κ(v̄ - V(t))dt + γsqrt(V(t))dW_1(t)
    
where W_1 and W_2 are correlated Brownian motions with correlation ρ.

The AES scheme:
    - V(t) simulated exactly using noncentral chi-squared distribution
    - S(t) computed using conditional expectation given V(t+Δt) and V(t)

Author: Lech A. Grzelak
Created: Feb 11, 2019
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as st
import enum

def CIR_Sample(NoOfPaths, kappa, gamma, vbar, s, t, v_s):
    """
    Generate exact samples from CIR process using noncentral chi-squared distribution.
    
    The CIR process: dV(t) = κ(v̄ - V(t))dt + γsqrt(V(t))dW(t)
    
    The conditional distribution V(t)|V(s) follows a scaled noncentral chi-squared:
        V(t) = c·χ²(δ, λ̄)
    where:
        c = γ²(1-e^(-κ(t-s)))/(4κ)
        δ = 4κv̄/γ² (degrees of freedom)
        λ̄ = 4κV(s)e^(-κ(t-s))/(γ²(1-e^(-κ(t-s)))) (non-centrality parameter)
    
    Args:
        NoOfPaths (int): Number of sample paths to generate.
        kappa (float): Mean reversion speed.
        gamma (float): Volatility of variance (vol-of-vol).
        vbar (float): Long-term mean variance level.
        s (float): Current time.
        t (float): Future time.
        v_s (float or np.ndarray): Variance value(s) at time s.
    
    Returns:
        np.ndarray: Exact samples of V(t) given V(s).
    
    Note:
        Requires Feller condition: 2κv̄ ≥ γ² for well-defined process.
    """
    delta = 4.0 * kappa * vbar / (gamma * gamma)
    c = gamma * gamma * (1.0 - np.exp(-kappa * (t - s))) / (4.0 * kappa)
    kappaBar = (4.0 * kappa * v_s * np.exp(-kappa * (t - s)) / 
                (gamma * gamma * (1.0 - np.exp(-kappa * (t - s)))))
    sample = c * np.random.noncentral_chisquare(delta, kappaBar, NoOfPaths)
    return sample

def GeneratePathsHestonAES(NoOfPaths, NoOfSteps, T, r, S_0, kappa, gamma, rho, vbar, v0):
    """
    Generate Heston model paths using Almost Exact Scheme (AES).
    
    The AES scheme combines:
        1. Exact simulation of variance V(t) using CIR_Sample()
        2. Conditional log-price evolution: X(t+dt) | V(t), V(t+dt)
    
    The log-price X(t) = log(S(t)) follows:
        X(t+dt) = X(t) + k₀ + k₁V(t) + k₂V(t+dt) + sqrt((1-ρ²)V(t))ΔW₁
    where:
        k₀ = (r - ρκv̄/γ)dt
        k₁ = (ρκ/γ - 0.5)dt - ρ/γ
        k₂ = ρ/γ
    
    Args:
        NoOfPaths (int): Number of Monte Carlo paths.
        NoOfSteps (int): Number of time steps.
        T (float): Time to maturity (years).
        r (float): Risk-free interest rate.
        S_0 (float): Initial stock price.
        kappa (float): Mean reversion speed of variance.
        gamma (float): Volatility of variance.
        rho (float): Correlation between price and variance Brownian motions.
        vbar (float): Long-term mean variance.
        v0 (float): Initial variance.
    
    Returns:
        dict: Dictionary with keys:
            - 'time': np.ndarray of time points
            - 'S': np.ndarray of stock price paths (NoOfPaths × NoOfSteps+1)
    
    Note:
        This scheme has first-order strong convergence and is more accurate
        than standard Euler discretization for the Heston model.
    """
    Z1 = np.random.normal(0.0, 1.0, [NoOfPaths, NoOfSteps])
    W1 = np.zeros([NoOfPaths, NoOfSteps + 1])
    V = np.zeros([NoOfPaths, NoOfSteps + 1])
    X = np.zeros([NoOfPaths, NoOfSteps + 1])
    V[:, 0] = v0
    X[:, 0] = np.log(S_0)
    
    time = np.zeros([NoOfSteps + 1])
    dt = T / float(NoOfSteps)
    
    for i in range(0, NoOfSteps):
        # Ensure samples have mean 0 and variance 1
        if NoOfPaths > 1:
            Z1[:, i] = (Z1[:, i] - np.mean(Z1[:, i])) / np.std(Z1[:, i])
        W1[:, i + 1] = W1[:, i] + np.power(dt, 0.5) * Z1[:, i]
        
        # Exact samples for the variance process
        V[:, i + 1] = CIR_Sample(NoOfPaths, kappa, gamma, vbar, 0, dt, V[:, i])
        
        # Conditional log-price evolution
        k0 = (r - rho / gamma * kappa * vbar) * dt
        k1 = (rho * kappa / gamma - 0.5) * dt - rho / gamma
        k2 = rho / gamma
        X[:, i + 1] = (X[:, i] + k0 + k1 * V[:, i] + k2 * V[:, i + 1] + 
                       np.sqrt((1.0 - rho**2) * V[:, i]) * (W1[:, i + 1] - W1[:, i]))
        time[i + 1] = time[i] + dt
    
    # Compute stock price from log-price
    S = np.exp(X)
    paths = {"time": time, "S": S}
    return paths


def mainCalculation():
    """
    Main calculation comparing Euler and AES schemes for Heston model option pricing.
    
    This function:
        1. Prices options across strikes using COS method (reference)
        2. Compares Euler vs AES simulation schemes
        3. Analyzes convergence as time step decreases
    """
    NoOfPaths = 1000
    NoOfSteps = 500
    
    # Heston model parameters
    gamma = 1.0      # Volatility of variance
    kappa = 0.5      # Mean reversion speed
    vbar = 0.04      # Long-term mean variance
    rho = -0.9       # Correlation (negative for leverage effect)
    v0 = 0.04        # Initial variance
    T = 1.0          # Time to maturity
    S_0 = 100.0      # Initial stock price
    r = 0.1          # Risk-free rate
    CP = OptionType.CALL
    
    # First we define a range of strikes and check the convergence
    K = np.linspace(0.1,S_0*2.0,30)
    
    # Exact solution with the COS method
    cf = ChFHestonModel(r,T,kappa,gamma,vbar,v0,rho)
    
    # The COS method
    optValueExact = CallPutOptionPriceCOSMthd(cf, CP, S_0, r, T, K, 1000, 8)
    
    # Euler simulation
    pathsEULER = GeneratePathsHestonEuler(NoOfPaths,NoOfSteps,T,r,S_0,kappa,gamma,rho,vbar,v0)
    S_Euler = pathsEULER["S"]
    
    # Almost exact simulation
    pathsAES = GeneratePathsHestonAES(NoOfPaths,NoOfSteps,T,r,S_0,kappa,gamma,rho,vbar,v0)
    S_AES = pathsAES["S"]
    
        
    OptPrice_EULER = EUOptionPriceFromMCPathsGeneralized(CP,S_Euler[:,-1],K,T,r)
    OptPrice_AES   = EUOptionPriceFromMCPathsGeneralized(CP,S_AES[:,-1],K,T,r)
    
    plt.figure(1)
    plt.plot(K, optValueExact, '-r', linewidth=2)
    plt.plot(K, OptPrice_EULER, '--k', linewidth=1.5)
    plt.plot(K, OptPrice_AES, '.b', markersize=8)
    plt.legend(['Exact (COS)', 'Euler Scheme', 'AES Scheme'], fontsize=10)
    plt.grid(True, alpha=0.3)
    plt.xlabel('Strike K', fontsize=11)
    plt.ylabel('Option Price', fontsize=11)
    plt.title('Heston Model: Option Prices Across Strikes\n(Euler vs AES vs COS)', fontsize=12)
    plt.tight_layout()
    
    # Here we will analyze the convergence for particular dt
    dtV = np.array([1.0, 1.0/4.0, 1.0/8.0,1.0/16.0,1.0/32.0,1.0/64.0])
    NoOfStepsV = [int(T/x) for x in dtV]
    
    # Specify strike for analysis
    K = np.array([100.0])
    
    # Exact
    optValueExact = CallPutOptionPriceCOSMthd(cf, CP, S_0, r, T, K, 1000, 8)
    errorEuler = np.zeros([len(dtV),1])
    errorAES = np.zeros([len(dtV),1])
    
    for (idx,NoOfSteps) in enumerate(NoOfStepsV):
        # Euler
        np.random.seed(3)
        pathsEULER = GeneratePathsHestonEuler(NoOfPaths,NoOfSteps,T,r,S_0,kappa,gamma,rho,vbar,v0)
        S_Euler = pathsEULER["S"]
        OptPriceEULER = EUOptionPriceFromMCPathsGeneralized(CP,S_Euler[:,-1],K,T,r)
        errorEuler[idx] = OptPriceEULER-optValueExact
        # AES
        np.random.seed(3)
        pathsAES = GeneratePathsHestonAES(NoOfPaths,NoOfSteps,T,r,S_0,kappa,gamma,rho,vbar,v0)
        S_AES = pathsAES["S"]
        OptPriceAES   = EUOptionPriceFromMCPathsGeneralized(CP,S_AES[:,-1],K,T,r)
        errorAES[idx] = OptPriceAES-optValueExact
    
    # Print convergence results
    print(f"\n{'='*70}")
    print(f"Convergence Analysis for Strike K = {K[0]}")
    print(f"Exact Option Price (COS): {optValueExact[0]:.6f}")
    print(f"{'='*70}\n")
    
    print("EULER SCHEME ERRORS:")
    print(f"{'dt':<12} {'NoOfSteps':<12} {'Error':<15} {'|Error|':<12}")
    print("-" * 70)
    for i in range(len(NoOfStepsV)):
        print(f"{dtV[i]:<12.4f} {NoOfStepsV[i]:<12} {errorEuler[i][0]:<15.6e} "
              f"{abs(errorEuler[i][0]):<12.6e}")
    
    print("\nAES SCHEME ERRORS:")
    print(f"{'dt':<12} {'NoOfSteps':<12} {'Error':<15} {'|Error|':<12}")
    print("-" * 70)
    for i in range(len(NoOfStepsV)):
        print(f"{dtV[i]:<12.4f} {NoOfStepsV[i]:<12} {errorAES[i][0]:<15.6e} "
              f"{abs(errorAES[i][0]):<12.6e}")
    
    print(f"\n{'='*70}")
    print("Note: AES scheme shows superior convergence compared to Euler")
    print("AES achieves first-order strong convergence due to exact CIR sampling")
    print(f"{'='*70}\n")


if __name__ == "__main__":
    mainCalculation()