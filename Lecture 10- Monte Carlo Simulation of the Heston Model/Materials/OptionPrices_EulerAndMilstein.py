"""
Monte Carlo Option Pricing: Euler vs Milstein Scheme Convergence Analysis.

This module compares the convergence properties of Euler and Milstein discretization
schemes for pricing European and Cash-or-Nothing options under GBM.

The GBM dynamics: dS(t) = rS(t)dt + σS(t)dW(t)

Discretization Schemes:
    - Euler: S(t+dt) = S(t) + rS(t)dt + σS(t)ΔW
    - Milstein: S(t+dt) = S(t) + rS(t)dt + σS(t)ΔW + 0.5σ²S(t)((ΔW)² - dt)

For European options:
    - Both schemes show good convergence as paths increase
    - Errors decrease with Monte Carlo convergence rate O(1/√N)

For Digital (Cash-or-Nothing) options:
    - Payoff discontinuity at strike causes slower convergence
    - Milstein typically outperforms Euler due to higher order approximation

Author: Lech A. Grzelak
Created: Jan 20, 2019
"""
import numpy as np
import scipy.stats as st
import enum

# Set i = imaginary number
i = np.complex(0.0, 1.0)


# This class defines puts and calls
class OptionType(enum.Enum):
    CALL = 1.0
    PUT = -1.0
    
def BS_Call_Option_Price(CP, S_0, K, sigma, tau, r):
    """
    Calculate European option price using Black-Scholes formula.
    
    Black-Scholes formula:
        Call: C = S₀N(d₁) - Ke^(-rτ)N(d₂)
        Put: P = Ke^(-rτ)N(-d₂) - S₀N(-d₁)
    where:
        d₁ = [ln(S₀/K) + (r + σ²/2)τ] / (σ√τ)
        d₂ = d₁ - σ√τ
    
    Args:
        CP (OptionType): CALL or PUT option type.
        S_0 (float): Initial stock price.
        K (list or np.ndarray): Strike price(s).
        sigma (float): Volatility (annualized).
        tau (float): Time to maturity (years).
        r (float): Risk-free interest rate.
    
    Returns:
        np.ndarray: Option price(s).
    """
    K = np.array(K).reshape([len(K), 1])
    d1 = (np.log(S_0 / K) + (r + 0.5 * np.power(sigma, 2.0)) * tau) / (sigma * np.sqrt(tau))
    d2 = d1 - sigma * np.sqrt(tau)
    
    if CP == OptionType.CALL:
        value = st.norm.cdf(d1) * S_0 - st.norm.cdf(d2) * K * np.exp(-r * tau)
    elif CP == OptionType.PUT:
        value = st.norm.cdf(-d2) * K * np.exp(-r * tau) - st.norm.cdf(-d1) * S_0
    return value

def GeneratePathsGBMEuler(NoOfPaths, NoOfSteps, T, r, sigma, S_0):
    """
    Generate GBM paths using Euler-Maruyama discretization.
    
    Euler scheme for GBM: dS(t) = rS(t)dt + σS(t)dW(t)
        S(t+dt) = S(t) + rS(t)dt + σS(t)ΔW
    
    The Euler scheme has strong convergence order 0.5.
    
    Args:
        NoOfPaths (int): Number of Monte Carlo paths.
        NoOfSteps (int): Number of time steps.
        T (float): Time to maturity (years).
        r (float): Risk-free interest rate.
        sigma (float): Volatility (annualized).
        S_0 (float): Initial stock price.
    
    Returns:
        dict: Dictionary with keys:
            - 'time': np.ndarray of time points
            - 'S': np.ndarray of stock price paths
    """
    Z = np.random.normal(0.0, 1.0, [NoOfPaths, NoOfSteps])
    W = np.zeros([NoOfPaths, NoOfSteps + 1])
    
    # Euler Approximation
    S1 = np.zeros([NoOfPaths, NoOfSteps + 1])
    S1[:, 0] = S_0
    
    time = np.zeros([NoOfSteps + 1])
    dt = T / float(NoOfSteps)
    
    for i in range(0, NoOfSteps):
        # Ensure samples have mean 0 and variance 1
        if NoOfPaths > 1:
            Z[:, i] = (Z[:, i] - np.mean(Z[:, i])) / np.std(Z[:, i])
        W[:, i + 1] = W[:, i] + np.power(dt, 0.5) * Z[:, i]
        
        S1[:, i + 1] = S1[:, i] + r * S1[:, i] * dt + sigma * S1[:, i] * (W[:, i + 1] - W[:, i])
        time[i + 1] = time[i] + dt
    
    paths = {"time": time, "S": S1}
    return paths

def BS_Cash_Or_Nothing_Price(CP, S_0, K, sigma, tau, r):
    """
    Calculate Cash-or-Nothing (Digital) option price using Black-Scholes.
    
    Cash-or-Nothing payoff:
        Call: K if S(T) > K, else 0
        Put: K if S(T) ≤ K, else 0
    
    Black-Scholes formula:
        Call: C = Ke^(-rτ)N(d₂)
        Put: P = Ke^(-rτ)(1 - N(d₂))
    where d₂ = [ln(S₀/K) + (r - σ²/2)τ] / (σ√τ)
    
    Args:
        CP (OptionType): CALL or PUT option type.
        S_0 (float): Initial stock price.
        K (list or np.ndarray): Strike price(s).
        sigma (float): Volatility (annualized).
        tau (float): Time to maturity (years).
        r (float): Risk-free interest rate.
    
    Returns:
        np.ndarray: Digital option price(s).
    
    Note:
        Digital options have discontinuous payoff, making Monte Carlo
        convergence slower than for standard European options.
    """
    K = np.array(K).reshape([len(K), 1])
    d1 = (np.log(S_0 / K) + (r + 0.5 * np.power(sigma, 2.0)) * tau) / (sigma * np.sqrt(tau))
    d2 = d1 - sigma * np.sqrt(tau)
    
    if CP == OptionType.CALL:
        value = K * np.exp(-r * tau) * st.norm.cdf(d2)
    elif CP == OptionType.PUT:
        value = K * np.exp(-r * tau) * (1.0 - st.norm.cdf(d2))
    return value

def GeneratePathsGBMMilstein(NoOfPaths, NoOfSteps, T, r, sigma, S_0):
    """
    Generate GBM paths using Milstein discretization scheme.
    
    Milstein scheme for GBM: dS(t) = rS(t)dt + σS(t)dW(t)
        S(t+dt) = S(t) + rS(t)dt + σS(t)ΔW + 0.5σ²S(t)((ΔW)² - dt)
    
    The Milstein scheme includes the Itô correction term (ΔW)² - dt,
    achieving strong convergence order 1.0 (better than Euler's 0.5).
    
    Args:
        NoOfPaths (int): Number of Monte Carlo paths.
        NoOfSteps (int): Number of time steps.
        T (float): Time to maturity (years).
        r (float): Risk-free interest rate.
        sigma (float): Volatility (annualized).
        S_0 (float): Initial stock price.
    
    Returns:
        dict: Dictionary with keys:
            - 'time': np.ndarray of time points
            - 'S': np.ndarray of stock price paths
    
    Note:
        Milstein scheme uses derivative of diffusion coefficient:
            ∂(σS)/∂S = σ, leading to correction term 0.5σ²S((ΔW)² - dt)
    """
    Z = np.random.normal(0.0, 1.0, [NoOfPaths, NoOfSteps])
    W = np.zeros([NoOfPaths, NoOfSteps + 1])
    
    # Milstein Approximation
    S1 = np.zeros([NoOfPaths, NoOfSteps + 1])
    S1[:, 0] = S_0
    
    time = np.zeros([NoOfSteps + 1])
    dt = T / float(NoOfSteps)
    
    for i in range(0, NoOfSteps):
        # Ensure samples have mean 0 and variance 1
        if NoOfPaths > 1:
            Z[:, i] = (Z[:, i] - np.mean(Z[:, i])) / np.std(Z[:, i])
        W[:, i + 1] = W[:, i] + np.power(dt, 0.5) * Z[:, i]
        
        # Milstein scheme with Itô correction
        S1[:, i + 1] = (S1[:, i] + r * S1[:, i] * dt + 
                        sigma * S1[:, i] * (W[:, i + 1] - W[:, i]) +
                        0.5 * sigma**2.0 * S1[:, i] * 
                        (np.power((W[:, i + 1] - W[:, i]), 2) - dt))
        time[i + 1] = time[i] + dt
    
    paths = {"time": time, "S": S1}
    return paths

def EUOptionPriceFromMCPaths(CP, S, K, T, r):
    """
    Price European option from Monte Carlo terminal stock prices.
    
    Args:
        CP (OptionType): CALL or PUT option type.
        S (np.ndarray): Terminal stock prices from MC simulation.
        K (float): Strike price.
        T (float): Time to maturity.
        r (float): Risk-free rate.
    
    Returns:
        float: Discounted expected payoff.
    """
    if CP == OptionType.CALL:
        return np.exp(-r * T) * np.mean(np.maximum(S - K, 0.0))
    elif CP == OptionType.PUT:
        return np.exp(-r * T) * np.mean(np.maximum(K - S, 0.0))


def CashofNothingPriceFromMCPaths(CP, S, K, T, r):
    """
    Price Cash-or-Nothing option from Monte Carlo terminal stock prices.
    
    Args:
        CP (OptionType): CALL or PUT option type.
        S (np.ndarray): Terminal stock prices from MC simulation.
        K (float): Strike price (also the cash payoff amount).
        T (float): Time to maturity.
        r (float): Risk-free rate.
    
    Returns:
        float: Discounted expected cash payoff.
    """
    if CP == OptionType.CALL:
        return np.exp(-r * T) * K * np.mean((S > K))
    elif CP == OptionType.PUT:
        return np.exp(-r * T) * K * np.mean((S <= K))

def mainCalculation():
    """
    Main calculation comparing Euler and Milstein schemes for option pricing.
    
    Tests convergence for:
        1. European Call Options (smooth payoff)
        2. Cash-or-Nothing Options (discontinuous payoff)
    
    Demonstrates that Milstein generally provides better accuracy due to
    higher-order approximation of the stochastic differential equation.
    """
    CP = OptionType.CALL
    T = 1
    r = 0.06
    sigma = 0.3
    S_0 = 5
    K = [S_0]
    NoOfSteps = 1000
    
    # Number of paths to test
    NoOfPathsV = [100, 1000, 5000, 10000]
    
    # European Call Option Pricing
    exactPrice = BS_Call_Option_Price(CP, S_0, K, sigma, T, r)[0]
    print("\n" + "=" * 80)
    print("EUROPEAN CALL OPTION PRICING")
    print("=" * 80)
    print(f"Parameters: S₀={S_0}, K={K[0]}, r={r}, σ={sigma}, T={T}")
    print(f"Analytical (Black-Scholes) Price: {exactPrice:.6f}")
    print("-" * 80)
    print(f"{'Paths':<10} {'Euler Price':<15} {'Milstein Price':<17} "
          f"{'Euler Error':<15} {'Milstein Error':<15}")
    print("-" * 80)
    
    for NoOfPathsTemp in NoOfPathsV:
        np.random.seed(1)
        PathsEuler = GeneratePathsGBMEuler(NoOfPathsTemp, NoOfSteps, T, r, sigma, S_0)
        np.random.seed(1)
        PathsMilstein = GeneratePathsGBMMilstein(NoOfPathsTemp, NoOfSteps, T, r, sigma, S_0)
        S_Euler = PathsEuler["S"]
        S_Milstein = PathsMilstein["S"]
        priceEuler = EUOptionPriceFromMCPaths(CP, S_Euler[:, -1], K, T, r)
        priceMilstein = EUOptionPriceFromMCPaths(CP, S_Milstein[:, -1], K, T, r)
        errorEuler = priceEuler - exactPrice
        errorMilstein = priceMilstein - exactPrice
        
        print(f"{NoOfPathsTemp:<10} {priceEuler:<15.6f} {priceMilstein:<17.6f} "
              f"{errorEuler:<15.6e} {errorMilstein:<15.6e}")
    
    # Cash-or-Nothing Option Pricing
    exactPrice = BS_Cash_Or_Nothing_Price(CP, S_0, K, sigma, T, r)
    print("\n" + "=" * 80)
    print("CASH-OR-NOTHING (DIGITAL) OPTION PRICING")
    print("=" * 80)
    print(f"Parameters: S₀={S_0}, K={K[0]}, r={r}, σ={sigma}, T={T}")
    print(f"Analytical (Black-Scholes) Price: {exactPrice[0]:.6f}")
    print("-" * 80)
    print(f"{'Paths':<10} {'Euler Price':<15} {'Milstein Price':<17} "
          f"{'Euler Error':<15} {'Milstein Error':<15}")
    print("-" * 80)
    
    for NoOfPathsTemp in NoOfPathsV:
        np.random.seed(1)
        PathsEuler = GeneratePathsGBMEuler(NoOfPathsTemp, NoOfSteps, T, r, sigma, S_0)
        np.random.seed(1)
        PathsMilstein = GeneratePathsGBMMilstein(NoOfPathsTemp, NoOfSteps, T, r, sigma, S_0)
        S_Euler = PathsEuler["S"]
        S_Milstein = PathsMilstein["S"]
        priceEuler = CashofNothingPriceFromMCPaths(CP, S_Euler[:, -1], K[0], T, r)
        priceMilstein = CashofNothingPriceFromMCPaths(CP, S_Milstein[:, -1], K[0], T, r)
        errorEuler = priceEuler - exactPrice
        errorMilstein = priceMilstein - exactPrice
        
        print(f"{NoOfPathsTemp:<10} {priceEuler:<15.6f} {priceMilstein:<17.6f} "
              f"{errorEuler:<15.6e} {errorMilstein:<15.6e}")
    
    print("\n" + "=" * 80)
    print("OBSERVATIONS:")
    print("- European options: Both schemes converge well (smooth payoff)")
    print("- Digital options: Slower convergence due to discontinuous payoff")
    print("- Milstein typically more accurate due to higher-order approximation")
    print("- Convergence rate: O(1/√N) for Monte Carlo (path-dependent)")
    print("=" * 80 + "\n")


if __name__ == "__main__":
    mainCalculation()