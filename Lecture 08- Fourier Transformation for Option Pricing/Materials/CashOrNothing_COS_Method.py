#%%
"""
Cash-or-Nothing Option Pricing using the COS Method.

Cash-or-nothing (digital) options pay a fixed cash amount if the option ends
in-the-money, otherwise nothing. This module prices them using the COS method.

Created on Thu Jan 16 2019
@author: Lech A. Grzelak
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as st
import time


def CashOrNothingPriceCOSMthd(cf, CP, S0, r, tau, K, N, L):
    """
    Price cash-or-nothing digital options using COS method.
    
    Cash-or-nothing payoffs:
    - Call: pays K if S(T) > K, else 0
    - Put: pays K if S(T) < K, else 0
    
    Args:
        cf (callable): Characteristic function.
        CP (str): 'c' for call, 'p' for put.
        S0 (float): Initial stock price.
        r (float): Risk-free rate.
        tau (float): Time to maturity.
        K (list/array): Strike prices (also cash payouts).
        N (int): Number of expansion terms.
        L (float): Truncation parameter.
    
    Returns:
        np.ndarray: Cash-or-nothing option prices.
    
    Note:
        Digital options are highly sensitive to parameters near strike.
        They are used for speculation and structured products.
    """
    K = np.array(K).reshape([len(K), 1])
    i = np.complex(0.0, 1.0)
    
    x0 = np.log(S0 / K)
    
    # Truncation domain
    a = 0.0 - L * np.sqrt(tau)
    b = 0.0 + L * np.sqrt(tau)
    
    k = np.linspace(0, N-1, N).reshape([N, 1])
    u = k * np.pi / (b - a)

    # Coefficients for cash-or-nothing payoff
    H_k = CashOrNothingCoefficients(CP, a, b, k)
       
    mat = np.exp(i * np.outer((x0 - a), u))
    temp = cf(u) * H_k 
    temp[0] = 0.5 * temp[0]
    
    value = np.exp(-r * tau) * K * np.real(mat.dot(temp))
         
    return value


def CashOrNothingCoefficients(CP, a, b, k):
    """
    Compute Fourier coefficients for cash-or-nothing payoff.
    
    Args:
        CP (str): 'c' for call, 'p' for put.
        a, b (float): Truncation bounds.
        k (np.ndarray): Summation indices.
    
    Returns:
        np.ndarray: Fourier coefficients H_k.
    """
    if str(CP).lower() == "c" or str(CP).lower() == "1":
        c = 0.0
        d = b
        coef = Chi_Psi(a, b, c, d, k)
        Psi_k = coef["psi"]
        if a < b and b < 0.0:
            H_k = np.zeros([len(k), 1])
        else:
            H_k = 2.0 / (b - a) * Psi_k
        
    elif str(CP).lower() == "p" or str(CP).lower() == "-1":
        c = a
        d = 0.0
        coef = Chi_Psi(a, b, c, d, k)
        Psi_k = coef["psi"]
        H_k = 2.0 / (b - a) * Psi_k
    
    return H_k


def Chi_Psi(a, b, c, d, k):
    """
    Compute auxiliary functions for digital payoff integration.
    
    Args:
        a, b (float): Truncation bounds.
        c, d (float): Integration bounds.
        k (np.ndarray): Summation indices.
    
    Returns:
        dict: Dictionary with 'chi' and 'psi' arrays.
    """
    psi = np.sin(k * np.pi * (d - a) / (b - a)) - np.sin(k * np.pi * (c - a) / (b - a))
    psi[1:] = psi[1:] * (b - a) / (k[1:] * np.pi)
    psi[0] = d - c
    
    chi = 1.0 / (1.0 + np.power((k * np.pi / (b - a)), 2.0))
    expr1 = np.cos(k * np.pi * (d - a) / (b - a)) * np.exp(d) - \
            np.cos(k * np.pi * (c - a) / (b - a)) * np.exp(c)
    expr2 = k * np.pi / (b - a) * np.sin(k * np.pi * (d - a) / (b - a)) - \
            k * np.pi / (b - a) * np.sin(k * np.pi * (c - a) / (b - a)) * np.exp(c)
    chi = chi * (expr1 + expr2)
    
    value = {"chi": chi, "psi": psi}
    return value


def BS_Cash_Or_Nothing_Price(CP, S_0, K, sigma, tau, r):
    """
    Analytical Black-Scholes price for cash-or-nothing options.
    
    Args:
        CP (str): 'c' for call, 'p' for put.
        S_0 (float): Initial stock price.
        K (list/array): Strike/payout amounts.
        sigma (float): Volatility.
        tau (float): Time to maturity.
        r (float): Risk-free rate.
    
    Returns:
        np.ndarray: Cash-or-nothing prices.
    """
    cp = str(CP).lower()
    K = np.array(K).reshape([len(K), 1])
    d1 = (np.log(S_0 / K) + (r + 0.5 * np.power(sigma, 2.0)) * tau) / \
         float(sigma * np.sqrt(tau))
    d2 = d1 - sigma * np.sqrt(tau)
    if cp == "c" or cp == "1":
        value = K * np.exp(-r * tau) * st.norm.cdf(d2)
    elif cp == "p" or cp == "-1":
        value = K * np.exp(-r * tau) * (1.0 - st.norm.cdf(d2))
    return value


def mainCalculation():
    """
    Demonstrate cash-or-nothing option pricing with COS method.
    
    Tests convergence and performance for digital options.
    """
    i = np.complex(0.0, 1.0)
    
    CP = "p"
    S0 = 100.0
    r = 0.05
    tau = 0.1
    sigma = 0.2
    K = [120]
    N = [40, 60, 80, 100, 120, 140]
    L = 6
    
    # Characteristic function for GBM
    cf = lambda u: np.exp((r - 0.5 * np.power(sigma, 2.0)) * i * u * tau - 0.5 *
                          np.power(sigma, 2.0) * np.power(u, 2.0) * tau)
    
    # High-accuracy reference value
    val_COS_Exact = CashOrNothingPriceCOSMthd(cf, CP, S0, r, tau, K, np.power(2, 14), L)
    print(f"\n=== Cash-or-Nothing {CP.upper()} Option ===")
    print(f"S0={S0}, K={K[0]}, σ={sigma}, τ={tau}, r={r}")
    print(f"Reference value (N=2^14): {val_COS_Exact[0][0]:.8f}\n")
    
    # Convergence analysis
    NoOfIterations = 1000
    print("=== Convergence Analysis ===")
    errors = []
    times = []
    
    for n in N:
        time_start = time.time()
        for k in range(0, NoOfIterations):
            val_COS = CashOrNothingPriceCOSMthd(cf, CP, S0, r, tau, K, n, L)[0]
        time_stop = time.time()
        
        error = val_COS[0] - val_COS_Exact[0]
        avg_time = (time_stop - time_start) / float(NoOfIterations)
        errors.append(error[0])
        times.append(avg_time * 1000)
        
        print(f"N={n:3d}: Error={error[0]:10.2E}, Time={avg_time*1000:.4f} ms")
    
    # Visualization
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
    
    ax1.semilogy(N, np.abs(errors), 'bo-', linewidth=2, markersize=8)
    ax1.grid(True, alpha=0.3)
    ax1.set_xlabel('Number of Terms N')
    ax1.set_ylabel('Absolute Error')
    ax1.set_title('COS Method Convergence for Digital Options')
    
    ax2.plot(N, times, 'ro-', linewidth=2, markersize=8)
    ax2.grid(True, alpha=0.3)
    ax2.set_xlabel('Number of Terms N')
    ax2.set_ylabel('Computation Time (ms)')
    ax2.set_title('Computational Efficiency')
    
    plt.tight_layout()
    
    print("\n✓ COS method provides fast and accurate pricing for digital options")


if __name__ == "__main__":
    mainCalculation()
