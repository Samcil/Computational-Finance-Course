"""
Forward Start Options Pricing under Heston Model Using COS Method.

Forward start options are options where the strike is determined at a future
date T₁, typically set at-the-money: K = S(T₁). The payoff at maturity T₂ is:
    Payoff = max(S(T₂) - S(T₁)(1 + k), 0)
where k is the strike adjustment (often k = 0 for ATM forward start).

Under Heston model, the challenge is that V(T₁) is stochastic, requiring:
    1. Integration over V(T₁) distribution (noncentral chi-squared)
    2. Conditional characteristic function given V(T₁)

The Heston dynamics:
    dS(t)/S(t) = rdt + sqrt(V(t))dW₂(t)
    dV(t) = κ(v̄ - V(t))dt + γsqrt(V(t))dW₁(t)
    Corr(dW₁, dW₂) = ρdt

Key insight: V(T₁)|v₀ follows scaled noncentral χ² distribution, enabling
semi-analytical pricing via COS method with conditional CF.

This module:
    - Prices forward start options via COS method
    - Analyzes implied volatility surface
    - Studies impact of time-to-start (T₁) and time-to-maturity (T₂)

Author: Lech A. Grzelak
Created: Jan 2, 2019
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as st
import enum
import scipy.optimize as optimize

# Set i = imaginary number
i = np.complex(0.0, 1.0)


class OptionType(enum.Enum):
    """Enumeration for option types."""
    CALL = 1.0
    PUT = -1.0

def CallPutOptionPriceCOSMthd_FrwdStart(cf, CP, r, T1, T2, K, N, L):
    """
    Price forward start options using COS method.
    
    Forward start option: Strike K_eff = S(T₁)(1 + k) determined at T₁.
    Payoff at T₂: max(S(T₂) - K_eff, 0) = max(S(T₂)/S(T₁) - (1+k), 0)S(T₁)
    
    By normalizing with S(T₁), we price option on S(T₂)/S(T₁) with strike (1+k).
    
    Args:
        cf (function): Characteristic function for log(S(T₂)/S(T₁)) | V(T₁).
        CP (OptionType): CALL or PUT.
        r (float): Risk-free rate.
        T1 (float): Time when strike is determined.
        T2 (float): Option maturity.
        K (np.ndarray): Strike adjustments k (actual strike = S(T₁)(1+k)).
        N (int): COS expansion terms.
        L (float): Truncation domain size.
    
    Returns:
        np.ndarray: Forward start option prices.
    
    Note:
        Discounting is from T₂ (option maturity), not T₁.
        Characteristic function must be conditional on V(T₁).
    """
    tau = T2 - T1
    if K is not np.array:
        K = np.array(K).reshape([len(K), 1])
    
    # Adjust strike (forward start specific)
    K = K + 1.0
    
    i = np.complex(0.0, 1.0)
    x0 = np.log(1.0 / K)
    
    # Truncation domain
    a = 0.0 - L * np.sqrt(tau)
    b = 0.0 + L * np.sqrt(tau)
    
    # Fourier series
    k = np.linspace(0, N - 1, N).reshape([N, 1])
    u = k * np.pi / (b - a)
    
    # Payoff coefficients
    H_k = CallPutCoefficients(CP, a, b, k)
    mat = np.exp(i * np.outer((x0 - a), u))
    temp = cf(u) * H_k
    temp[0] = 0.5 * temp[0]
    value = np.exp(-r * T2) * K * np.real(mat.dot(temp))
    return value

def CallPutCoefficients(CP, a, b, k):
    """
    Compute Fourier coefficients for call/put payoffs in COS method.
    
    Args:
        CP (OptionType): CALL or PUT.
        a, b (float): Truncation domain bounds.
        k (np.ndarray): Fourier indices.
    
    Returns:
        np.ndarray: Payoff coefficients H_k.
    """
    if CP == OptionType.CALL:
        c = 0.0
        d = b
        coef = Chi_Psi(a, b, c, d, k)
        Chi_k = coef["chi"]
        Psi_k = coef["psi"]
        if a < b and b < 0.0:
            H_k = np.zeros([len(k), 1])
        else:
            H_k = 2.0 / (b - a) * (Chi_k - Psi_k)
    elif CP == OptionType.PUT:
        c = a
        d = 0.0
        coef = Chi_Psi(a, b, c, d, k)
        Chi_k = coef["chi"]
        Psi_k = coef["psi"]
        H_k = 2.0 / (b - a) * (-Chi_k + Psi_k)
    
    return H_k


def Chi_Psi(a, b, c, d, k):
    """
    Compute auxiliary functions χ and ψ for COS method payoff integration.
    
    Returns:
        dict: {'chi': χ_k, 'psi': ψ_k}
    """
    psi = (np.sin(k * np.pi * (d - a) / (b - a)) -
           np.sin(k * np.pi * (c - a) / (b - a)))
    psi[1:] = psi[1:] * (b - a) / (k[1:] * np.pi)
    psi[0] = d - c
    
    chi = 1.0 / (1.0 + np.power((k * np.pi / (b - a)), 2.0))
    expr1 = (np.cos(k * np.pi * (d - a) / (b - a)) * np.exp(d) -
             np.cos(k * np.pi * (c - a) / (b - a)) * np.exp(c))
    expr2 = (k * np.pi / (b - a) * np.sin(k * np.pi * (d - a) / (b - a)) -
             k * np.pi / (b - a) * np.sin(k * np.pi * (c - a) / (b - a)) * np.exp(c))
    chi = chi * (expr1 + expr2)
    
    value = {"chi": chi, "psi": psi}
    return value
    
def BS_Call_Option_Price_FrwdStart(K, sigma, T1, T2, r):
    """
    Black-Scholes price for forward start call option.
    
    Forward start BS formula (ATM forward start):
        C = e^(-rT₁)[N(d₁) - (1+k)e^(-r(T₂-T₁))N(d₂)]
    
    Args:
        K (float or np.ndarray): Strike adjustment k.
        sigma (float): Volatility.
        T1 (float): Forward start time.
        T2 (float): Maturity.
        r (float): Risk-free rate.
    
    Returns:
        np.ndarray: Forward start option prices.
    """
    if K is list:
        K = np.array(K).reshape([len(K), 1])
    K = K + 1.0
    tau = T2 - T1
    d1 = (np.log(1.0 / K) + (r + 0.5 * np.power(sigma, 2.0)) * tau) / (sigma * np.sqrt(tau))
    d2 = d1 - sigma * np.sqrt(tau)
    value = np.exp(-r * T1) * st.norm.cdf(d1) - st.norm.cdf(d2) * K * np.exp(-r * T2)
    return value


def ImpliedVolatility_FrwdStart(marketPrice, K, T1, T2, r):
    """
    Calculate implied volatility for forward start option.
    
    Uses two-step process:
        1. Grid search for initial guess
        2. Newton's method for refinement
    
    Args:
        marketPrice (float): Observed option price.
        K (float): Strike adjustment k.
        T1 (float): Forward start time.
        T2 (float): Maturity.
        r (float): Risk-free rate.
    
    Returns:
        float: Implied volatility.
    """
    # Grid search
    sigmaGrid = np.linspace(0, 2, 200)
    optPriceGrid = BS_Call_Option_Price_FrwdStart(K, sigmaGrid, T1, T2, r)
    sigmaInitial = np.interp(marketPrice, optPriceGrid, sigmaGrid)
    print(f"Initial volatility = {sigmaInitial:.6f}")
    
    # Newton's method refinement
    func = lambda sigma: np.power(BS_Call_Option_Price_FrwdStart(K, sigma, T1, T2, r) - marketPrice, 1.0)
    impliedVol = optimize.newton(func, sigmaInitial, tol=1e-15)
    print(f"Final volatility = {impliedVol:.6f}")
    return impliedVol

def ChFHestonModelForwardStart(r, T1, T2, kappa, gamma, vbar, v0, rho):
    """
    Characteristic function for forward start option under Heston model.
    
    The CF integrates over V(T₁) distribution (noncentral χ²) and uses
    conditional CF for log(S(T₂)/S(T₁)) given V(T₁).
    
    Key components:
        - A(u), C(u): Standard Heston CF terms for period [T₁, T₂]
        - c̄(T₁) = γ²(1-e^(-κT₁))/(4κ): Scaling of V(T₁) distribution
        - δ = 4κv̄/γ²: Degrees of freedom
        - κ̄(T₁) = 4κv₀e^(-κT₁)/(γ²(1-e^(-κT₁))): Non-centrality parameter
    
    The CF has form:
        φ(u) = E[e^(iuX)] = exp(term1(u)) · (1/(1-2C(u)c̄))^(δ/2)
    where:
        term1 = A(u) + C(u)c̄κ̄/(1-2C(u)c̄)
        term2 = (1/(1-2C(u)c̄))^(δ/2)
    
    Args:
        r (float): Risk-free rate.
        T1 (float): Forward start time.
        T2 (float): Maturity.
        kappa (float): Mean reversion speed.
        gamma (float): Vol-of-vol.
        vbar (float): Long-term variance.
        v0 (float): Initial variance.
        rho (float): Correlation.
    
    Returns:
        function: Characteristic function φ(u).
    
    Note:
        This CF accounts for stochastic V(T₁) via integration over
        its noncentral chi-squared distribution.
    """
    i = np.complex(0.0, 1.0)
    tau = T2 - T1
    
    # Standard Heston CF components for [T₁, T₂]
    D1 = lambda u: np.sqrt(np.power(kappa - gamma * rho * i * u, 2) + 
                           (u * u + i * u) * gamma * gamma)
    g = lambda u: ((kappa - gamma * rho * i * u - D1(u)) / 
                   (kappa - gamma * rho * i * u + D1(u)))
    C = lambda u: ((1.0 - np.exp(-D1(u) * tau)) / 
                   (gamma * gamma * (1.0 - g(u) * np.exp(-D1(u) * tau))) *
                   (kappa - gamma * rho * i * u - D1(u)))
    
    # Note: Discounting -r*tau performed in COS method
    A = lambda u: (r * i * u * tau + 
                   kappa * vbar * tau / (gamma * gamma) * (kappa - gamma * rho * i * u - D1(u)) -
                   2 * kappa * vbar / (gamma * gamma) * 
                   np.log((1.0 - g(u) * np.exp(-D1(u) * tau)) / (1.0 - g(u))))
    
    # V(T₁) distribution parameters
    c_bar = lambda t1, t2: gamma * gamma / (4.0 * kappa) * (1.0 - np.exp(-kappa * (t2 - t1)))
    delta = 4.0 * kappa * vbar / (gamma * gamma)  # Degrees of freedom
    kappa_bar = lambda t1, t2: (4.0 * kappa * v0 * np.exp(-kappa * (t2 - t1)) / 
                                (gamma * gamma * (1.0 - np.exp(-kappa * (t2 - t1)))))
    
    # Forward start CF: integrate over V(T₁)
    term1 = lambda u: A(u) + C(u) * c_bar(0.0, T1) * kappa_bar(0.0, T1) / (1.0 - 2.0 * C(u) * c_bar(0.0, T1))
    term2 = lambda u: np.power(1.0 / (1.0 - 2.0 * C(u) * c_bar(0.0, T1)), 0.5 * delta)
    cf = lambda u: np.exp(term1(u)) * term2(u)
    return cf

def mainCalculation():
    """
    Analyze Heston forward start option implied volatility.
    
    Two scenarios explored:
        1. Fixed T₁, varying T₂-T₁ (fixed start, varying maturity)
        2. Fixed T₂-T₁, varying T₁ (varying start, fixed tenor)
    
    Expected patterns:
        - Volatility smile/skew due to stochastic volatility
        - Smile shape depends on correlation ρ (leverage effect)
        - Longer periods generally show flatter smiles
    """
    CP = OptionType.CALL
    r = 0.00
    
    # Scenario 1: Fixed T₁, varying tenor (T₂ - T₁)
    TMat1 = [[1.0, 3.0], [2.0, 4.0], [3.0, 5.0], [4.0, 6.0]]
    # Scenario 2: Fixed tenor, varying T₁
    TMat2 = [[1.0, 2.0], [1.0, 3.0], [1.0, 4.0], [1.0, 5.0]]
    
    K = np.linspace(-0.4, 4.0, 50)
    K = np.array(K).reshape([len(K), 1])
    
    N = 500  # COS expansion terms
    L = 10   # Truncation domain
    
    # Heston model parameters
    kappa = 0.6    # Mean reversion speed
    gamma = 0.2    # Vol-of-vol
    vbar = 0.1     # Long-term variance
    rho = -0.5     # Correlation (negative for leverage effect)
    v0 = 0.05      # Initial variance
    
    # Figure 1: Fixed start time, varying maturity
    plt.figure(1, figsize=(10, 6))
    plt.grid(True, alpha=0.3)
    plt.xlabel('Strike Adjustment k', fontsize=11)
    plt.ylabel('Implied Volatility (%)', fontsize=11)
    plt.title('Heston Forward Start: Fixed T₁, Varying Tenor (T₂-T₁)', fontsize=12)
    legend = []
    
    for T_pair in TMat1:
        T1 = T_pair[0]
        T2 = T_pair[1]
        cf = ChFHestonModelForwardStart(r, T1, T2, kappa, gamma, vbar, v0, rho)
        valCOS = CallPutOptionPriceCOSMthd_FrwdStart(cf, CP, r, T1, T2, K, N, L)
        
        IV = np.zeros([len(K), 1])
        for idx in range(0, len(K)):
            IV[idx] = ImpliedVolatility_FrwdStart(valCOS[idx], K[idx], T1, T2, r)
        
        plt.plot(K, IV * 100.0, linewidth=2, marker='o', markersize=4)
        legend.append(f'T₁={T1}, T₂={T2} (tenor={T2-T1})')
    
    plt.legend(legend, fontsize=9)
    plt.tight_layout()
    
    # Figure 2: Fixed tenor, varying start time
    plt.figure(2, figsize=(10, 6))
    plt.grid(True, alpha=0.3)
    plt.xlabel('Strike Adjustment k', fontsize=11)
    plt.ylabel('Implied Volatility (%)', fontsize=11)
    plt.title('Heston Forward Start: Fixed Tenor, Varying T₁', fontsize=12)
    legend = []
    
    for T_pair in TMat2:
        T1 = T_pair[0]
        T2 = T_pair[1]
        cf = ChFHestonModelForwardStart(r, T1, T2, kappa, gamma, vbar, v0, rho)
        valCOS = CallPutOptionPriceCOSMthd_FrwdStart(cf, CP, r, T1, T2, K, N, L)
        
        IV = np.zeros([len(K), 1])
        for idx in range(0, len(K)):
            IV[idx] = ImpliedVolatility_FrwdStart(valCOS[idx], K[idx], T1, T2, r)
        
        plt.plot(K, IV * 100.0, linewidth=2, marker='s', markersize=4)
        legend.append(f'T₁={T1}, T₂={T2} (tenor={T2-T1})')
    
    plt.legend(legend, fontsize=9)
    plt.tight_layout()
    
    print("\n" + "=" * 70)
    print("HESTON FORWARD START OPTIONS IMPLIED VOLATILITY ANALYSIS")
    print("=" * 70)
    print("Model Parameters:")
    print(f"  Heston: κ={kappa}, γ={gamma}, v̄={vbar}, ρ={rho}, v₀={v0}")
    print(f"  Market: r={r}")
    print("-" * 70)
    print("Strike Convention:")
    print("  K = strike adjustment, Actual strike = S(T₁)(1 + k)")
    print("  K=0 → ATM forward start option")
    print("-" * 70)
    print("Analysis:")
    print("  Figure 1: Impact of varying maturity T₂ (fixed start T₁)")
    print("    - Longer tenor → potentially flatter smile")
    print("  Figure 2: Impact of varying start time T₁ (fixed tenor)")
    print("    - Different T₁ → different V(T₁) distributions")
    print("-" * 70)
    print("Key Insight:")
    print("  Forward start options embed uncertainty in both:")
    print("    1. Future strike level S(T₁)")
    print("    2. Future variance V(T₁)")
    print("  Both sources contribute to implied volatility smile/skew")
    print("=" * 70 + "\n")


if __name__ == "__main__":
    mainCalculation()