#%%
"""
Pathwise Method for Computing Option Greeks (Delta and Vega).

This module demonstrates the pathwise differentiation method for computing
option sensitivities (Greeks) using Monte Carlo simulation. The pathwise
method provides unbiased estimators with lower variance than finite difference
methods.

Created on 08 Mar 2019
@author: Lech A. Grzelak
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as st
import enum 


class OptionType(enum.Enum):
    """Enumeration for option types."""
    CALL = 1.0
    PUT = -1.0


def BS_Call_Put_Option_Price(CP, S_0, K, sigma, t, T, r):
    """
    Compute Black-Scholes option price.
    
    Args:
        CP (OptionType): CALL or PUT.
        S_0 (float): Current stock price.
        K (array): Strike prices.
        sigma (float): Volatility.
        t (float): Current time.
        T (float): Maturity time.
        r (float): Risk-free rate.
    
    Returns:
        np.ndarray: Option prices.
    """
    K = np.array(K).reshape([len(K), 1])
    d1 = (np.log(S_0 / K) + (r + 0.5 * np.power(sigma, 2.0)) * (T - t)) / \
         (sigma * np.sqrt(T - t))
    d2 = d1 - sigma * np.sqrt(T - t)
    if CP == OptionType.CALL:
        value = st.norm.cdf(d1) * S_0 - st.norm.cdf(d2) * K * np.exp(-r * (T - t))
    elif CP == OptionType.PUT:
        value = st.norm.cdf(-d2) * K * np.exp(-r * (T - t)) - st.norm.cdf(-d1) * S_0
    return value


def BS_Delta(CP, S_0, K, sigma, t, T, r):
    """
    Compute Black-Scholes delta (∂V/∂S).
    
    Args:
        CP (OptionType): CALL or PUT.
        S_0 (float): Current stock price.
        K (array): Strike prices.
        sigma (float): Volatility.
        t (float): Current time.
        T (float): Maturity time.
        r (float): Risk-free rate.
    
    Returns:
        np.ndarray: Delta values.
    """
    K = np.array(K).reshape([len(K), 1])
    d1 = (np.log(S_0 / K) + (r + 0.5 * np.power(sigma, 2.0)) * (T - t)) / \
         (sigma * np.sqrt(T - t))
    if CP == OptionType.CALL:
        value = st.norm.cdf(d1)
    elif CP == OptionType.PUT:
        value = st.norm.cdf(d1) - 1
    return value


def BS_Gamma(S_0, K, sigma, t, T, r):
    """
    Compute Black-Scholes gamma (∂²V/∂S²).
    
    Args:
        S_0 (float): Current stock price.
        K (array): Strike prices.
        sigma (float): Volatility.
        t (float): Current time.
        T (float): Maturity time.
        r (float): Risk-free rate.
    
    Returns:
        np.ndarray: Gamma values.
    """
    K = np.array(K).reshape([len(K), 1])
    d1 = (np.log(S_0 / K) + (r + 0.5 * np.power(sigma, 2.0)) * (T - t)) / \
         (sigma * np.sqrt(T - t))
    return st.norm.pdf(d1) / (S_0 * sigma * np.sqrt(T - t))


def BS_Vega(S_0, K, sigma, t, T, r):
    """
    Compute Black-Scholes vega (∂V/∂σ).
    
    Args:
        S_0 (float): Current stock price.
        K (array): Strike prices.
        sigma (float): Volatility.
        t (float): Current time.
        T (float): Maturity time.
        r (float): Risk-free rate.
    
    Returns:
        np.ndarray: Vega values.
    """
    d1 = (np.log(S_0 / K) + (r + 0.5 * np.power(sigma, 2.0)) * (T - t)) / \
         (sigma * np.sqrt(T - t))
    return S_0 * st.norm.pdf(d1) * np.sqrt(T - t)

def GeneratePathsGBMEuler(NoOfPaths, NoOfSteps, T, r, sigma, S_0):
    """
    Generate GBM paths using Euler discretization.
    
    Args:
        NoOfPaths (int): Number of paths.
        NoOfSteps (int): Number of time steps.
        T (float): Time horizon.
        r (float): Risk-free rate.
        sigma (float): Volatility.
        S_0 (float): Initial stock price.
    
    Returns:
        dict: Dictionary with 'time' and 'S' arrays.
    """
    Z = np.random.normal(0.0, 1.0, [NoOfPaths, NoOfSteps])
    W = np.zeros([NoOfPaths, NoOfSteps+1])
   
    S = np.zeros([NoOfPaths, NoOfSteps+1])
    S[:, 0] = S_0
    
    X = np.zeros([NoOfPaths, NoOfSteps+1])
    X[:, 0] = np.log(S_0)
      
    time = np.zeros([NoOfSteps+1])
        
    dt = T / float(NoOfSteps)
    for i in range(0, NoOfSteps):
        # Making sure that samples from normal have mean 0 and variance 1
        if NoOfPaths > 1:
            Z[:, i] = (Z[:, i] - np.mean(Z[:, i])) / np.std(Z[:, i])
        W[:, i+1] = W[:, i] + np.power(dt, 0.5) * Z[:, i]
        
        X[:, i+1] = X[:, i] + (r - 0.5 * sigma**2.0) * dt + sigma * (W[:, i+1] - W[:, i])
        time[i+1] = time[i] + dt
        
    paths = {"time": time, "S": np.exp(X)}
    return paths


def EUOptionPriceFromMCPathsGeneralized(CP, S, K, T, r):
    """
    Price European options from Monte Carlo paths.
    
    Args:
        CP (OptionType): CALL or PUT.
        S (np.ndarray): Terminal stock prices.
        K (array): Strike prices.
        T (float): Time to maturity.
        r (float): Risk-free rate.
    
    Returns:
        np.ndarray: Option prices.
    """
    result = np.zeros([len(K), 1])
    if CP == OptionType.CALL:
        for (idx, k) in enumerate(K):
            result[idx] = np.exp(-r * T) * np.mean(np.maximum(S - k, 0.0))
    elif CP == OptionType.PUT:
        for (idx, k) in enumerate(K):
            result[idx] = np.exp(-r * T) * np.mean(np.maximum(k - S, 0.0))
    return result


def PathwiseDelta(S0, S, K, r, T):
    """
    Compute delta using the pathwise method.
    
    The pathwise estimator for call delta is:
        Δ = E[exp(-rT) * (S(T)/S0) * 1_{S(T)>K}]
    
    Args:
        S0 (float): Initial stock price.
        S (np.ndarray): Stock price paths.
        K (float): Strike price.
        r (float): Risk-free rate.
        T (float): Time to maturity.
    
    Returns:
        float: Pathwise delta estimate.
    
    Note:
        The pathwise method differentiates the payoff function directly,
        providing an unbiased estimator with typically lower variance than
        finite differences.
    """
    temp1 = S[:, -1] > K
    return np.exp(-r * T) * np.mean(S[:, -1] / S0 * temp1)


def PathwiseVega(S0, S, sigma, K, r, T):
    """
    Compute vega using the pathwise method.
    
    The pathwise estimator for call vega is:
        ν = E[exp(-rT) * S(T) * (1/σ)(log(S(T)/S0) - (r+0.5σ²)T) * 1_{S(T)>K}]
    
    Args:
        S0 (float): Initial stock price.
        S (np.ndarray): Stock price paths.
        sigma (float): Volatility.
        K (float): Strike price.
        r (float): Risk-free rate.
        T (float): Time to maturity.
    
    Returns:
        float: Pathwise vega estimate.
    
    Note:
        Vega measures sensitivity to volatility changes. The pathwise method
        exploits the smooth dependence of GBM paths on σ.
    """
    temp1 = S[:, -1] > K
    temp2 = 1.0 / sigma * S[:, -1] * (np.log(S[:, -1] / S0) - (r + 0.5 * sigma**2.0) * T)
    return np.exp(-r * T) * np.mean(temp1 * temp2)


def mainCalculation():
    """
    Demonstrate pathwise method convergence for delta and vega estimation.
    
    Compares pathwise Monte Carlo estimates against analytical Black-Scholes
    Greeks, showing convergence as the number of paths increases.
    """
    CP = OptionType.CALL
    S0 = 1
    r = 0.06
    sigma = 0.3
    T = 1
    K = np.array([S0])
    t = 0.0

    NoOfSteps = 1000
    delta_Exact = BS_Delta(CP, S0, K, sigma, t, T, r)
    vega_Exact = BS_Vega(S0, K, sigma, t, T, r)
    
    print(f"=== Pathwise Method Convergence Analysis ===")
    print(f"Option: {CP.name}, S0={S0}, K={K[0]}, σ={sigma}, T={T}, r={r}")
    print(f"Analytical Delta: {delta_Exact[0][0]:.6f}")
    print(f"Analytical Vega:  {vega_Exact:.6f}\n")
    
    NoOfPathsV = np.round(np.linspace(5, 1000, 50))
    deltaPathWiseV = np.zeros(len(NoOfPathsV))
    vegaPathWiseV = np.zeros(len(NoOfPathsV))
    
    for (idx, nPaths) in enumerate(NoOfPathsV):
        if idx % 10 == 0:
            print(f'Running simulation with {int(nPaths)} paths')
        np.random.seed(3)
        paths1 = GeneratePathsGBMEuler(int(nPaths), NoOfSteps, T, r, sigma, S0)
        S = paths1["S"]
        delta_pathwise = PathwiseDelta(S0, S, K, r, T)
        deltaPathWiseV[idx] = delta_pathwise
        
        vega_pathwise = PathwiseVega(S0, S, sigma, K, r, T)
        vegaPathWiseV[idx] = vega_pathwise
    
    # Delta convergence plot
    plt.figure(1, figsize=(12, 5))
    plt.grid(True, alpha=0.3)
    plt.plot(NoOfPathsV, deltaPathWiseV, '.-r', linewidth=1.5, markersize=6, label='Pathwise estimate')
    plt.plot(NoOfPathsV, delta_Exact * np.ones([len(NoOfPathsV), 1]), 'b--', 
             linewidth=2, label=f'Analytical ({delta_Exact[0][0]:.4f})')
    plt.xlabel('Number of Paths')
    plt.ylabel('Delta (∂V/∂S)')
    plt.title('Pathwise Delta Convergence')
    plt.legend()
    
    # Vega convergence plot
    plt.figure(2, figsize=(12, 5))
    plt.grid(True, alpha=0.3)
    plt.plot(NoOfPathsV, vegaPathWiseV, '.-r', linewidth=1.5, markersize=6, label='Pathwise estimate')
    plt.plot(NoOfPathsV, vega_Exact * np.ones([len(NoOfPathsV), 1]), 'b--', 
             linewidth=2, label=f'Analytical ({vega_Exact:.4f})')
    plt.xlabel('Number of Paths')
    plt.ylabel('Vega (∂V/∂σ)')
    plt.title('Pathwise Vega Convergence')
    plt.legend()
    
    # Final error report
    final_delta_error = abs(deltaPathWiseV[-1] - delta_Exact[0][0])
    final_vega_error = abs(vegaPathWiseV[-1] - vega_Exact)
    print(f"\n=== Final Results (N={int(NoOfPathsV[-1])} paths) ===")
    print(f"Delta error: {final_delta_error:.6f}")
    print(f"Vega error:  {final_vega_error:.6f}")
    print("\n✓ Pathwise method provides unbiased Greek estimates")


if __name__ == "__main__":
    mainCalculation()