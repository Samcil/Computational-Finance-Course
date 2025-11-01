"""
Delta Hedging Jump-Diffusion Models with Black-Scholes Delta.

This module demonstrates the challenges of delta hedging when the underlying
follows a jump-diffusion (Merton) model but the hedge is based on Black-Scholes
delta, which assumes continuous paths.

The Merton jump-diffusion model:
    dS(t)/S(t) = (r - λ(e^(μ_J + 0.5σ_J²) - 1))dt + σdW(t) + (e^J - 1)dN(t)

where:
    - N(t) is a Poisson process with intensity λ (xiP)
    - J ~ N(μ_J, σ_J²) are i.i.d. log-jump sizes
    - W(t) is standard Brownian motion

Hedging Strategy:
    - Rebalance portfolio continuously using Black-Scholes delta
    - P&L = V(t-dt)e^(r·dt) - (Δ_new - Δ_old)S(t) - Payoff + Δ_final·S(T)
    
The P&L distribution shows:
    - Mean close to zero (unbiased hedging)
    - Fat tails due to jump risk
    - Variance higher than in BS model without jumps

Author: Lech A. Grzelak
Created: Dec 12, 2018
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as st
import enum
from mpl_toolkits import mplot3d
from scipy.interpolate import RegularGridInterpolator


class OptionType(enum.Enum):
    """Enumeration for option types."""
    CALL = 1.0
    PUT = -1.0

def GeneratePathsMerton(NoOfPaths, NoOfSteps, S0, T, xiP, muJ, sigmaJ, r, sigma):
    """
    Generate stock price paths under Merton jump-diffusion model.
    
    The Merton model combines geometric Brownian motion with compound Poisson jumps:
        dX(t) = (r - λμ - σ²/2)dt + σdW(t) + JdN(t)
    where:
        X(t) = log(S(t))
        λμ = xiP·(e^(μ_J + 0.5σ_J²) - 1) is the jump compensator
        J ~ N(μ_J, σ_J²) are log-jump sizes
        N(t) ~ Poisson(λt) counts number of jumps
    
    Args:
        NoOfPaths (int): Number of Monte Carlo paths.
        NoOfSteps (int): Number of time steps.
        S0 (float): Initial stock price.
        T (float): Time to maturity (years).
        xiP (float): Jump intensity (λ, expected jumps per year).
        muJ (float): Mean log-jump size (μ_J).
        sigmaJ (float): Std dev of log-jump size (σ_J).
        r (float): Risk-free interest rate.
        sigma (float): Diffusion volatility.
    
    Returns:
        dict: Dictionary with keys:
            - 'time': np.ndarray of time points
            - 'X': np.ndarray of log-price paths
            - 'S': np.ndarray of stock price paths
    
    Note:
        Jump compensator ensures martingale property under risk-neutral measure.
    """
    X = np.zeros([NoOfPaths, NoOfSteps + 1])
    S = np.zeros([NoOfPaths, NoOfSteps + 1])
    time = np.zeros([NoOfSteps + 1])
    
    dt = T / float(NoOfSteps)
    X[:, 0] = np.log(S0)
    S[:, 0] = S0
    
    # Expectation E(e^J) for J~N(muJ, sigmaJ^2)
    EeJ = np.exp(muJ + 0.5 * sigmaJ * sigmaJ)
    ZPois = np.random.poisson(xiP * dt, [NoOfPaths, NoOfSteps])
    Z = np.random.normal(0.0, 1.0, [NoOfPaths, NoOfSteps])
    J = np.random.normal(muJ, sigmaJ, [NoOfPaths, NoOfSteps])
    
    for i in range(0, NoOfSteps):
        # Ensure samples have mean 0 and variance 1
        if NoOfPaths > 1:
            Z[:, i] = (Z[:, i] - np.mean(Z[:, i])) / np.std(Z[:, i])
        
        # Merton model dynamics with jump compensator
        X[:, i + 1] = (X[:, i] + (r - xiP * (EeJ - 1) - 0.5 * sigma * sigma) * dt +
                       sigma * np.sqrt(dt) * Z[:, i] + J[:, i] * ZPois[:, i])
        time[i + 1] = time[i] + dt
    
    S = np.exp(X)
    paths = {"time": time, "X": X, "S": S}
    return paths

def GeneratePathsGBM(NoOfPaths, NoOfSteps, T, r, sigma, S_0):
    """
    Generate GBM paths (for comparison with Merton model).
    
    GBM dynamics: dX(t) = (r - σ²/2)dt + σdW(t), where X(t) = log(S(t))
    
    Args:
        NoOfPaths (int): Number of Monte Carlo paths.
        NoOfSteps (int): Number of time steps.
        T (float): Time to maturity (years).
        r (float): Risk-free interest rate.
        sigma (float): Volatility.
        S_0 (float): Initial stock price.
    
    Returns:
        dict: Dictionary with 'time' and 'S' arrays.
    """
    Z = np.random.normal(0.0, 1.0, [NoOfPaths, NoOfSteps])
    X = np.zeros([NoOfPaths, NoOfSteps + 1])
    W = np.zeros([NoOfPaths, NoOfSteps + 1])
    time = np.zeros([NoOfSteps + 1])
    
    X[:, 0] = np.log(S_0)
    dt = T / float(NoOfSteps)
    
    for i in range(0, NoOfSteps):
        # Ensure samples have mean 0 and variance 1
        if NoOfPaths > 1:
            Z[:, i] = (Z[:, i] - np.mean(Z[:, i])) / np.std(Z[:, i])
        W[:, i + 1] = W[:, i] + np.power(dt, 0.5) * Z[:, i]
        X[:, i + 1] = X[:, i] + (r - 0.5 * sigma * sigma) * dt + sigma * (W[:, i + 1] - W[:, i])
        time[i + 1] = time[i] + dt
    
    # Compute stock price from log-price
    S = np.exp(X)
    paths = {"time": time, "S": S}
    return paths

def BS_Call_Put_Option_Price(CP, S_0, K, sigma, t, T, r):
    """
    Calculate Black-Scholes option price.
    
    Args:
        CP (OptionType): CALL or PUT.
        S_0 (float): Current stock price.
        K (list): Strike price(s).
        sigma (float): Volatility.
        t (float): Current time.
        T (float): Maturity time.
        r (float): Risk-free rate.
    
    Returns:
        np.ndarray: Option price(s).
    """
    K = np.array(K).reshape([len(K), 1])
    d1 = ((np.log(S_0 / K) + (r + 0.5 * np.power(sigma, 2.0)) * (T - t)) / 
          (sigma * np.sqrt(T - t)))
    d2 = d1 - sigma * np.sqrt(T - t)
    
    if CP == OptionType.CALL:
        value = st.norm.cdf(d1) * S_0 - st.norm.cdf(d2) * K * np.exp(-r * (T - t))
    elif CP == OptionType.PUT:
        value = st.norm.cdf(-d2) * K * np.exp(-r * (T - t)) - st.norm.cdf(-d1) * S_0
    return value


def BS_Delta(CP, S_0, K, sigma, t, T, r):
    """
    Calculate Black-Scholes delta (∂V/∂S).
    
    Delta measures sensitivity of option price to stock price changes.
    For hedging, we hold Δ shares of stock for each option sold.
    
    Args:
        CP (OptionType): CALL or PUT.
        S_0 (float): Current stock price.
        K (list): Strike price(s).
        sigma (float): Volatility.
        t (float): Current time.
        T (float): Maturity time.
        r (float): Risk-free rate.
    
    Returns:
        np.ndarray: Delta value(s).
            - Call: Δ ∈ [0, 1]
            - Put: Δ ∈ [-1, 0]
    
    Note:
        When t ≈ T, numerical issues may arise. Code handles t ≥ T edge case.
    """
    # Handle edge case where time grid slightly exceeds maturity
    if t - T > 10e-20 and T - t < 10e-7:
        t = T
    
    K = np.array(K).reshape([len(K), 1])
    d1 = ((np.log(S_0 / K) + (r + 0.5 * np.power(sigma, 2.0)) * (T - t)) / 
          (sigma * np.sqrt(T - t)))
    
    if CP == OptionType.CALL:
        value = st.norm.cdf(d1)
    elif CP == OptionType.PUT:
        value = st.norm.cdf(d1) - 1.0
    return value


def mainCalculation():
    """
    Main calculation demonstrating delta hedging under jump-diffusion dynamics.
    
    Strategy:
        1. Sell call option and delta-hedge using BS delta
        2. Rebalance continuously (5000 steps per year)
        3. Track P&L evolution for all paths
        4. Analyze final P&L distribution
    
    Expected Result:
        - Mean P&L ≈ 0 (unbiased hedging)
        - Non-zero variance due to jump risk (model misspecification)
        - Fatter tails than in pure BS world
    """
    NoOfPaths = 1000
    NoOfSteps = 5000
    T = 1.0
    r = 0.1
    sigma = 0.2      # Diffusion volatility
    xiP = 1.0        # Jump intensity (1 jump per year on average)
    muJ = 0.0        # Mean log-jump size
    sigmaJ = 0.25    # Jump volatility
    s0 = 1.0
    K = [0.95]
    CP = OptionType.CALL
    
    np.random.seed(7)
    Paths = GeneratePathsMerton(NoOfPaths, NoOfSteps, s0, T, xiP, muJ, sigmaJ, r, sigma)
    time = Paths["time"]
    S = Paths["S"]
    
    # Lambda functions for option pricing and delta
    C = lambda t, K, S0: BS_Call_Put_Option_Price(CP, S0, K, sigma, t, T, r)
    Delta = lambda t, K, S0: BS_Delta(CP, S0, K, sigma, t, T, r)
    
    # Initialize portfolio: short option + delta shares
    PnL = np.zeros([NoOfPaths, NoOfSteps + 1])
    delta_init = Delta(0.0, K, s0)
    PnL[:, 0] = C(0.0, K, s0) - delta_init * s0
    
    CallM = np.zeros([NoOfPaths, NoOfSteps + 1])
    CallM[:, 0] = C(0.0, K, s0)
    DeltaM = np.zeros([NoOfPaths, NoOfSteps + 1])
    DeltaM[:, 0] = Delta(0, K, s0)
    
    # Dynamic hedging loop
    for i in range(1, NoOfSteps + 1):
        dt = time[i] - time[i - 1]
        delta_old = Delta(time[i - 1], K, S[:, i - 1])
        delta_curr = Delta(time[i], K, S[:, i])
        
        # P&L evolution: accrue interest and adjust hedge
        PnL[:, i] = PnL[:, i - 1] * np.exp(r * dt) - (delta_curr - delta_old) * S[:, i]
        CallM[:, i] = C(time[i], K, S[:, i])
        DeltaM[:, i] = Delta(time[i], K, S[:, i])
    
    # Final settlement: pay option payoff and liquidate hedge
    PnL[:, -1] = PnL[:, -1] - np.maximum(S[:, -1] - K, 0) + DeltaM[:, -1] * S[:, -1]
    
    # Visualization: Single path evolution
    path_id = 10
    plt.figure(1, figsize=(12, 6))
    plt.plot(time, S[path_id, :], 'b-', linewidth=2, label='Stock Price', alpha=0.8)
    plt.plot(time, CallM[path_id, :], 'r--', linewidth=2, label='Call Price')
    plt.plot(time, DeltaM[path_id, :], 'g-.', linewidth=1.5, label='Delta')
    plt.plot(time, PnL[path_id, :], 'm-', linewidth=2, label='P&L', alpha=0.7)
    plt.legend(fontsize=11, loc='best')
    plt.grid(True, alpha=0.3)
    plt.xlabel('Time (years)', fontsize=11)
    plt.ylabel('Value', fontsize=11)
    plt.title(f'Delta Hedging Under Merton Model (Path {path_id})', fontsize=12)
    plt.tight_layout()
    
    # P&L distribution histogram
    plt.figure(2, figsize=(10, 6))
    plt.hist(PnL[:, -1], bins=100, alpha=0.7, edgecolor='black')
    plt.grid(True, alpha=0.3)
    plt.xlim([-0.1, 0.1])
    plt.xlabel('Final P&L', fontsize=11)
    plt.ylabel('Frequency', fontsize=11)
    plt.title('Distribution of Final P&L\n(Delta Hedging with Jumps)', fontsize=12)
    plt.axvline(x=0, color='r', linestyle='--', linewidth=2, label='Zero P&L')
    plt.axvline(x=np.mean(PnL[:, -1]), color='g', linestyle='--', linewidth=2, 
                label=f'Mean P&L = {np.mean(PnL[:, -1]):.4f}')
    plt.legend(fontsize=10)
    plt.tight_layout()
    
    # Summary statistics
    print("\n" + "=" * 70)
    print("DELTA HEDGING WITH JUMP-DIFFUSION MODEL")
    print("=" * 70)
    print(f"Model: Merton Jump-Diffusion")
    print(f"  σ (diffusion) = {sigma}, λ (jump intensity) = {xiP}, "
          f"μ_J = {muJ}, σ_J = {sigmaJ}")
    print(f"Hedging: Black-Scholes Delta (continuous rebalancing)")
    print(f"Option: Call, K = {K[0]}, T = {T}, S₀ = {s0}")
    print("-" * 70)
    print(f"Number of paths: {NoOfPaths}")
    print(f"Rebalancing steps: {NoOfSteps}")
    print("-" * 70)
    print(f"Final P&L Statistics:")
    print(f"  Mean P&L:     {np.mean(PnL[:, -1]):>10.6f}")
    print(f"  Std Dev P&L:  {np.std(PnL[:, -1]):>10.6f}")
    print(f"  Min P&L:      {np.min(PnL[:, -1]):>10.6f}")
    print(f"  Max P&L:      {np.max(PnL[:, -1]):>10.6f}")
    print("-" * 70)
    print(f"Sample Path {path_id}:")
    print(f"  S₀:              {s0:.6f}")
    print(f"  S(T):            {S[path_id, -1]:.6f}")
    print(f"  Option Payoff:   {np.maximum(S[path_id, -1] - K[0], 0.0):.6f}")
    print(f"  P&L(T-dt):       {PnL[path_id, -2]:.6f}")
    print(f"  Final P&L:       {PnL[path_id, -1]:.6f}")
    print("=" * 70)
    print("\nNote: Jump risk causes P&L variance even with continuous hedging.")
    print("BS delta cannot fully hedge jump-diffusion model (model risk).")
    print("=" * 70 + "\n")


if __name__ == "__main__":
    mainCalculation()