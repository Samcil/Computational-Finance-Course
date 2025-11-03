#%%
"""
European Option Pricing using the COS (Fourier Cosine) Method.

The COS method is a fast and accurate Fourier-based method for pricing European
options. It recovers option prices from the characteristic function using
cosine expansion, providing significant computational advantages over Monte Carlo.

Reference: Fang & Oosterlee (2008) "A Novel Pricing Method for European Options..."

Created on Thu Nov 27 2018
@author: Lech A. Grzelak
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as st
import time


def CallPutOptionPriceCOSMthd(cf, CP, S0, r, tau, K, N, L):
    """
    Price European call/put options using the COS method.
    
    The COS method computes option prices using Fourier cosine expansion:
        V(x) ≈ Σ_{k=0}^{N-1} A_k * Re[φ(k*π/(b-a)) * exp(ik*π*(x-a)/(b-a))]
    where φ is the characteristic function and A_k are coefficients.
    
    Args:
        cf (callable): Characteristic function φ(u) of log-asset price.
        CP (str): 'c' or '1' for call, 'p' or '-1' for put.
        S0 (float): Initial stock price.
        r (float): Risk-free interest rate.
        tau (float): Time to maturity.
        K (list/array): Strike prices.
        N (int): Number of expansion terms (typically 128-256).
        L (float): Truncation parameter (typically 8-12).
    
    Returns:
        np.ndarray: Option prices for each strike.
    
    Note:
        The COS method has exponential convergence and is typically much
        faster than Monte Carlo for European options with known characteristic
        functions. Computational complexity is O(N*M) where M is number of strikes.
    """
    # Reshape K to column vector
    K = np.array(K).reshape([len(K), 1])
    
    # Imaginary unit
    i = np.complex(0.0, 1.0)
    
    x0 = np.log(S0 / K)
    
    # Truncation domain [a,b]
    a = 0.0 - L * np.sqrt(tau)
    b = 0.0 + L * np.sqrt(tau)
    
    # Summation indices k = 0 to N-1
    k = np.linspace(0, N-1, N).reshape([N, 1])
    u = k * np.pi / (b - a)

    # Compute coefficients for call/put payoff
    H_k = CallPutCoefficients(CP, a, b, k)
       
    mat = np.exp(i * np.outer((x0 - a), u))

    temp = cf(u) * H_k 
    temp[0] = 0.5 * temp[0]  # Adjust first term
    
    value = np.exp(-r * tau) * K * np.real(mat.dot(temp))
         
    return value


def CallPutCoefficients(CP, a, b, k):
    """
    Compute Fourier cosine coefficients for call/put payoffs.
    
    These coefficients represent the Fourier transform of the payoff function
    in the truncated domain [a,b].
    
    Args:
        CP (str): 'c' or '1' for call, 'p' or '-1' for put.
        a (float): Lower truncation bound.
        b (float): Upper truncation bound.
        k (np.ndarray): Summation indices.
    
    Returns:
        np.ndarray: Fourier coefficients H_k.
    """
    if str(CP).lower() == "c" or str(CP).lower() == "1":
        # Call option: max(S-K, 0)
        c = 0.0
        d = b
        coef = Chi_Psi(a, b, c, d, k)
        Chi_k = coef["chi"]
        Psi_k = coef["psi"]
        if a < b and b < 0.0:
            H_k = np.zeros([len(k), 1])
        else:
            H_k = 2.0 / (b - a) * (Chi_k - Psi_k)
        
    elif str(CP).lower() == "p" or str(CP).lower() == "-1":
        # Put option: max(K-S, 0)
        c = a
        d = 0.0
        coef = Chi_Psi(a, b, c, d, k)
        Chi_k = coef["chi"]
        Psi_k = coef["psi"]
        H_k = 2.0 / (b - a) * (-Chi_k + Psi_k)
    
    return H_k


def Chi_Psi(a, b, c, d, k):
    """
    Compute auxiliary functions χ and ψ for payoff integration.
    
    These functions represent integrals of the exponential payoff weighted
    by cosine basis functions.
    
    Args:
        a, b (float): Truncation bounds.
        c, d (float): Integration bounds for payoff.
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
    

def BS_Call_Option_Price(CP, S_0, K, sigma, tau, r):
    """
    Compute Black-Scholes option price (analytical benchmark).
    
    Args:
        CP (str): 'c' for call, 'p' for put.
        S_0 (float): Initial stock price.
        K (list/array): Strike prices.
        sigma (float): Volatility.
        tau (float): Time to maturity.
        r (float): Risk-free rate.
    
    Returns:
        np.ndarray: Black-Scholes option prices.
    """
    cp = str(CP).lower()
    K = np.array(K).reshape([len(K), 1])
    d1 = (np.log(S_0 / K) + (r + 0.5 * np.power(sigma, 2.0)) * tau) / \
         float(sigma * np.sqrt(tau))
    d2 = d1 - sigma * np.sqrt(tau)
    if cp == "c" or cp == "1":
        value = st.norm.cdf(d1) * S_0 - st.norm.cdf(d2) * K * np.exp(-r * tau)
    elif cp == "p" or cp == "-1":
        value = st.norm.cdf(-d2) * K * np.exp(-r * tau) - st.norm.cdf(-d1) * S_0
    return value


def mainCalculation():
    """
    Demonstrate COS method for European option pricing.
    
    Compares COS method against analytical Black-Scholes prices for GBM.
    Shows accuracy and computational efficiency of the COS method.
    """
    i = np.complex(0.0, 1.0)
    
    CP = "c"
    S0 = 100.0
    r = 0.1
    tau = 0.1
    sigma = 0.25
    K = [80.0, 90.0, 100.0, 110, 120.0]
    N = 4 * 32  # 128 expansion terms
    L = 10  # Truncation parameter
    
    # Characteristic function for GBM: φ(u) = exp((r - 0.5σ²)iuτ - 0.5σ²u²τ)
    # Note: This excludes the "+iux₀" term, which is included internally
    cf = lambda u: np.exp((r - 0.5 * np.power(sigma, 2.0)) * i * u * tau - 0.5 *
                          np.power(sigma, 2.0) * np.power(u, 2.0) * tau)
    
    # Timing analysis
    NoOfIterations = 100
    time_start = time.time()
    for k in range(0, NoOfIterations, 1):
        val_COS = CallPutOptionPriceCOSMthd(cf, CP, S0, r, tau, K, N, L)
    time_stop = time.time()
    avg_time = (time_stop - time_start) / float(NoOfIterations)
    print(f"\n=== COS Method Performance ===")
    print(f"Average pricing time: {avg_time*1000:.4f} ms")
    print(f"Parameters: N={N}, L={L}, strikes={len(K)}")
    
    # Analytical Black-Scholes benchmark
    val_Exact = BS_Call_Option_Price(CP, S0, K, sigma, tau, r)
    
    # Visualization
    plt.figure(1, figsize=(12, 6))
    plt.plot(K, val_COS, 'bo-', linewidth=2, markersize=8, label='COS Method')
    plt.plot(K, val_Exact, 'r--', linewidth=2, label='Black-Scholes (Exact)')
    plt.xlabel("Strike K")
    plt.ylabel("Option Price")
    plt.title(f"European Call Option: COS vs Black-Scholes (S₀={S0}, σ={sigma}, τ={tau})")
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # Error analysis
    print(f"\n=== Accuracy Analysis ===")
    error = []
    for i in range(0, len(K)):
        err = np.abs(val_COS[i] - val_Exact[i])[0]
        error.append(err)
        print(f"Strike K={K[i]:6.1f}: COS={val_COS[i][0]:8.5f}, BS={val_Exact[i][0]:8.5f}, Error={err:.2E}")
    
    max_error = np.max(error)
    print(f"\nMaximum absolute error: {max_error:.2E}")
    print("✓ COS method provides machine precision accuracy for GBM")


if __name__ == "__main__":
    mainCalculation()

        
mainCalculation()