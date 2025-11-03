"""
Bates Model Implied Volatility Analysis Using COS Method.

The Bates model extends the Heston stochastic volatility model by adding
log-normal jumps in the asset price, combining features of both models:

Stock price dynamics:
    dS(t)/S(t) = (r - λμ)dt + sqrt(V(t))dW₂(t) + (e^J - 1)dN(t)

Variance dynamics:
    dV(t) = κ(v̄ - V(t))dt + γsqrt(V(t))dW₁(t)

where:
    - W₁, W₂ are correlated Brownian motions (correlation ρ)
    - N(t) is a Poisson process with intensity λ (xiP)
    - J ~ N(μ_J, σ_J²) are i.i.d. log-jump sizes
    - λμ = λ(e^(μ_J + 0.5σ_J²) - 1) is the jump compensator

This module:
    - Prices options using COS method with Bates characteristic function
    - Computes implied volatilities from option prices
    - Analyzes impact of jump parameters (λ, μ_J, σ_J) on volatility smile

Author: Lech A. Grzelak
Created: Nov 30, 2018
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
    
def CallPutOptionPriceCOSMthd(cf, CP, S0, r, tau, K, N, L):
    """
    Price European options using COS (Fourier-cosine) method.
    
    The COS method recovers option price from characteristic function via:
        V(x,t) ≈ e^(-rτ) Σ Re[φ(k)H_k exp(ik(x-a))]
    
    Args:
        cf (function): Characteristic function φ(u).
        CP (OptionType): CALL or PUT.
        S0 (float): Initial stock price.
        r (float): Risk-free rate.
        tau (float): Time to maturity.
        K (np.ndarray): Strike prices.
        N (int): Number of expansion terms.
        L (float): Truncation domain size ([−L√τ, L√τ]).
    
    Returns:
        np.ndarray: Option prices.
    """
    if K is not np.array:
        K = np.array(K).reshape([len(K), 1])
    
    i = np.complex(0.0, 1.0)
    x0 = np.log(S0 / K)
    
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
    value = np.exp(-r * tau) * K * np.real(mat.dot(temp))
    return value


def CallPutCoefficients(CP, a, b, k):
    """
    Compute Fourier coefficients for call/put payoffs.
    
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
    Compute auxiliary functions χ and ψ for COS method.
    
    These functions integrate payoff components over truncation domain.
    
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


def BS_Call_Option_Price(CP, S_0, K, sigma, tau, r):
    """
    Black-Scholes option price (for implied volatility calculation).
    
    Args:
        CP (OptionType): CALL or PUT.
        S_0 (float): Current stock price.
        K (float or np.ndarray): Strike price(s).
        sigma (float): Volatility.
        tau (float): Time to maturity.
        r (float): Risk-free rate.
    
    Returns:
        np.ndarray: Option price(s).
    """
    if K is list:
        K = np.array(K).reshape([len(K), 1])
    d1 = (np.log(S_0 / K) + (r + 0.5 * np.power(sigma, 2.0)) * tau) / (sigma * np.sqrt(tau))
    d2 = d1 - sigma * np.sqrt(tau)
    
    if CP == OptionType.CALL:
        value = st.norm.cdf(d1) * S_0 - st.norm.cdf(d2) * K * np.exp(-r * tau)
    elif CP == OptionType.PUT:
        value = st.norm.cdf(-d2) * K * np.exp(-r * tau) - st.norm.cdf(-d1) * S_0
    return value

def ImpliedVolatility(CP, marketPrice, K, T, S_0, r):
    """
    Calculate implied volatility from market price using Newton's method.
    
    Two-step approach:
        1. Grid search for initial guess (interpolation)
        2. Newton's method for precise solution
    
    Args:
        CP (OptionType): CALL or PUT option type.
        marketPrice (float): Observed market price.
        K (float): Strike price.
        T (float): Time to maturity.
        S_0 (float): Initial stock price.
        r (float): Risk-free rate.
    
    Returns:
        float: Implied volatility σ_imp such that BS(σ_imp) = market price.
    
    Note:
        For Bates model, implied volatility reflects both stochastic
        volatility and jump effects, creating characteristic smile/skew.
    """
    # Step 1: Grid search for initial guess
    sigmaGrid = np.linspace(0, 2, 200)
    optPriceGrid = BS_Call_Option_Price(CP, S_0, K, sigmaGrid, T, r)
    sigmaInitial = np.interp(marketPrice, optPriceGrid, sigmaGrid)
    print(f"Initial volatility = {sigmaInitial:.6f}")
    
    # Step 2: Newton's method for refinement
    func = lambda sigma: np.power(BS_Call_Option_Price(CP, S_0, K, sigma, T, r) - marketPrice, 1.0)
    impliedVol = optimize.newton(func, sigmaInitial, tol=1e-10)
    print(f"Final volatility = {impliedVol:.6f}")
    return impliedVol

def ChFBatesModel(r, tau, kappa, gamma, vbar, v0, rho, xiP, muJ, sigmaJ):
    """
    Characteristic function for the Bates model.
    
    The Bates model combines Heston stochastic volatility with Merton jumps.
    
    Characteristic function: φ(u) = exp(A(u) + C(u)v₀)
    where:
        - A(u) = A_Heston(u) + A_jump(u)
        - C(u) = coefficient for initial variance v₀
        - A_jump accounts for jump compensator and MGF of jumps
    
    Args:
        r (float): Risk-free rate.
        tau (float): Time to maturity.
        kappa (float): Mean reversion speed of variance.
        gamma (float): Volatility of variance (vol-of-vol).
        vbar (float): Long-term mean variance.
        v0 (float): Initial variance.
        rho (float): Correlation between price and variance.
        xiP (float): Jump intensity λ.
        muJ (float): Mean log-jump size μ_J.
        sigmaJ (float): Jump volatility σ_J.
    
    Returns:
        function: Characteristic function φ(u).
    
    Note:
        Jump terms:
            - Compensator: -λiuτ(e^(μ_J + 0.5σ_J²) - 1)
            - MGF: λτ(e^(iuμ_J - 0.5σ_J²u²) - 1)
    """
    i = np.complex(0.0, 1.0)
    D1 = lambda u: np.sqrt(np.power(kappa - gamma * rho * i * u, 2) + 
                           (u * u + i * u) * gamma * gamma)
    g = lambda u: ((kappa - gamma * rho * i * u - D1(u)) / 
                   (kappa - gamma * rho * i * u + D1(u)))
    C = lambda u: ((1.0 - np.exp(-D1(u) * tau)) / 
                   (gamma * gamma * (1.0 - g(u) * np.exp(-D1(u) * tau))) * 
                   (kappa - gamma * rho * i * u - D1(u)))
    
    # Heston component (excluding discounting -r*tau, done in COS method)
    AHes = lambda u: (r * i * u * tau + 
                      kappa * vbar * tau / (gamma * gamma) * (kappa - gamma * rho * i * u - D1(u)) -
                      2 * kappa * vbar / (gamma * gamma) * 
                      np.log((1.0 - g(u) * np.exp(-D1(u) * tau)) / (1.0 - g(u))))
    
    # Jump component
    A = lambda u: (AHes(u) - xiP * i * u * tau * (np.exp(muJ + 0.5 * sigmaJ * sigmaJ) - 1.0) +
                   xiP * tau * (np.exp(i * u * muJ - 0.5 * sigmaJ * sigmaJ * u * u) - 1.0))
    
    # Full Bates characteristic function
    cf = lambda u: np.exp(A(u) + C(u) * v0)
    return cf

def mainCalculation():
    """
    Analyze Bates model implied volatility surface.
    
    Explores impact of jump parameters on implied volatility smile/skew:
        1. Jump intensity (λ/xiP): Controls frequency of jumps
        2. Mean jump size (μ_J): Controls direction (negative → skew)
        3. Jump volatility (σ_J): Controls jump size dispersion
    
    Expected patterns:
        - Higher λ → steeper smile
        - Negative μ_J → pronounced left skew (crash premium)
        - Higher σ_J → wider smile
    """
    CP = OptionType.CALL
    S0 = 100.0
    r = 0.0
    tau = 1.0
    
    K = np.linspace(40, 180, 10)
    K = np.array(K).reshape([len(K), 1])
    
    N = 1000  # COS method expansion terms
    L = 6     # Truncation domain size
    
    # Bates model base parameters
    kappa = 1.2      # Mean reversion speed
    gamma = 0.05     # Vol-of-vol
    vbar = 0.05      # Long-term variance
    rho = -0.75      # Correlation (negative for leverage effect)
    v0 = vbar        # Initial variance
    muJ = 0.0        # Mean log-jump size
    sigmaJ = 0.2     # Jump volatility
    xiP = 0.1        # Jump intensity
    
    # Effect of Jump Intensity (xiP/λ)
    plt.figure(1, figsize=(10, 6))
    plt.grid(True, alpha=0.3)
    plt.xlabel('Strike K', fontsize=11)
    plt.ylabel('Implied Volatility (%)', fontsize=11)
    plt.title('Bates Model: Impact of Jump Intensity λ on Implied Volatility', fontsize=12)
    xiPV = [0.01, 0.1, 0.2, 0.3]
    legend = []
    
    for xiPTemp in xiPV:
        # Bates characteristic function
        cf = ChFBatesModel(r, tau, kappa, gamma, vbar, v0, rho, xiPTemp, muJ, sigmaJ)
        
        # Option prices via COS method
        valCOS = CallPutOptionPriceCOSMthd(cf, CP, S0, r, tau, K, N, L)
        
        # Implied volatilities
        IV = np.zeros([len(K), 1])
        for idx in range(0, len(K)):
            IV[idx] = ImpliedVolatility(CP, valCOS[idx], K[idx], tau, S0, r)
        
        plt.plot(K, IV * 100.0, linewidth=2, marker='o', markersize=6)
        legend.append(f'λ = {xiPTemp}')
    
    plt.legend(legend, fontsize=10)
    plt.tight_layout()
    
    # Effect of Mean Jump Size (μ_J)
    plt.figure(2, figsize=(10, 6))
    plt.grid(True, alpha=0.3)
    plt.xlabel('Strike K', fontsize=11)
    plt.ylabel('Implied Volatility (%)', fontsize=11)
    plt.title('Bates Model: Impact of Mean Jump Size μ_J on Implied Volatility', fontsize=12)
    muJPV = [-0.5, -0.25, 0, 0.25]
    legend = []
    
    for muJTemp in muJPV:
        # Bates characteristic function
        cf = ChFBatesModel(r, tau, kappa, gamma, vbar, v0, rho, xiP, muJTemp, sigmaJ)
        
        # Option prices via COS method
        valCOS = CallPutOptionPriceCOSMthd(cf, CP, S0, r, tau, K, N, L)
        
        # Implied volatilities
        IV = np.zeros([len(K), 1])
        for idx in range(0, len(K)):
            IV[idx] = ImpliedVolatility(CP, valCOS[idx], K[idx], tau, S0, r)
        
        plt.plot(K, IV * 100.0, linewidth=2, marker='s', markersize=6)
        legend.append(f'μ_J = {muJTemp}')
    
    plt.legend(legend, fontsize=10)
    plt.tight_layout()
    
    # Effect of Jump Volatility (σ_J)
    plt.figure(3, figsize=(10, 6))
    plt.grid(True, alpha=0.3)
    plt.xlabel('Strike K', fontsize=11)
    plt.ylabel('Implied Volatility (%)', fontsize=11)
    plt.title('Bates Model: Impact of Jump Volatility σ_J on Implied Volatility', fontsize=12)
    sigmaJV = [0.01, 0.15, 0.2, 0.25]
    legend = []
    
    for sigmaJTemp in sigmaJV:
        # Bates characteristic function
        cf = ChFBatesModel(r, tau, kappa, gamma, vbar, v0, rho, xiP, muJ, sigmaJTemp)
        
        # Option prices via COS method
        valCOS = CallPutOptionPriceCOSMthd(cf, CP, S0, r, tau, K, N, L)
        
        # Implied volatilities
        IV = np.zeros([len(K), 1])
        for idx in range(0, len(K)):
            IV[idx] = ImpliedVolatility(CP, valCOS[idx], K[idx], tau, S0, r)
        
        plt.plot(K, IV * 100.0, linewidth=2, marker='^', markersize=6)
        legend.append(f'σ_J = {sigmaJTemp}')
    
    plt.legend(legend, fontsize=10)
    plt.tight_layout()
    
    print("\n" + "=" * 70)
    print("BATES MODEL IMPLIED VOLATILITY ANALYSIS")
    print("=" * 70)
    print("Model Parameters:")
    print(f"  Heston: κ={kappa}, γ={gamma}, v̄={vbar}, ρ={rho}, v₀={v0}")
    print(f"  Base Jumps: λ={xiP}, μ_J={muJ}, σ_J={sigmaJ}")
    print(f"  Market: S₀={S0}, r={r}, T={tau}")
    print("-" * 70)
    print("Analysis:")
    print("  Figure 1: Jump intensity λ effect")
    print("    - Higher λ → steeper smile (more frequent jumps)")
    print("  Figure 2: Mean jump size μ_J effect")
    print("    - Negative μ_J → left skew (downside jump risk)")
    print("  Figure 3: Jump volatility σ_J effect")
    print("    - Higher σ_J → wider smile (larger jump uncertainty)")
    print("=" * 70 + "\n")


if __name__ == "__main__":
    mainCalculation()