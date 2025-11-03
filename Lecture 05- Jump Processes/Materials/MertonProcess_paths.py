#%%
"""
Merton Jump-Diffusion Process Path Generation.

This module implements the Merton (1976) jump-diffusion model, which extends
the Black-Scholes model by adding jump discontinuities to the stock price
dynamics. The model combines a geometric Brownian motion with a compound
Poisson process.

Created on Thu Jan 11 2019
@author: Lech A. Grzelak
"""
import numpy as np
import matplotlib.pyplot as plt


def GeneratePathsMerton(NoOfPaths, NoOfSteps, S0, T, xiP, muJ, sigmaJ, r, sigma):
    """
    Generate paths for the Merton jump-diffusion process.
    
    The Merton model describes the stock price dynamics as:
        dS(t)/S(t) = r*dt + σ*dW(t) + (e^J - 1)*dN(t)
    where:
        - W(t) is a Brownian motion
        - N(t) is a Poisson process with intensity λ (xiP)
        - J ~ N(μ_J, σ_J²) are i.i.d. jump sizes
    
    The drift is adjusted to maintain the martingale property under Q-measure:
        drift = r - λ*(E[e^J] - 1) - 0.5*σ²
    
    Args:
        NoOfPaths (int): Number of Monte Carlo paths to simulate.
        NoOfSteps (int): Number of time steps for discretization.
        S0 (float): Initial stock price.
        T (float): Time horizon.
        xiP (float): Jump intensity λ (expected number of jumps per unit time).
        muJ (float): Mean of log-jump size J.
        sigmaJ (float): Standard deviation of log-jump size J.
        r (float): Risk-free interest rate.
        sigma (float): Diffusion volatility parameter.
    
    Returns:
        dict: A dictionary containing:
            - 'time' (np.ndarray): Time grid of shape (NoOfSteps+1,).
            - 'X' (np.ndarray): Log-price paths of shape (NoOfPaths, NoOfSteps+1).
            - 'S' (np.ndarray): Stock price paths of shape (NoOfPaths, NoOfSteps+1).
    
    Note:
        The martingale correction E[e^J] = exp(μ_J + 0.5*σ_J²) ensures
        the discounted stock price is a martingale under the risk-neutral measure.
    """
    # Create empty matrices
    X = np.zeros([NoOfPaths, NoOfSteps+1])
    S = np.zeros([NoOfPaths, NoOfSteps+1])
    time = np.zeros([NoOfSteps+1])
                
    dt = T / float(NoOfSteps)
    X[:, 0] = np.log(S0)
    S[:, 0] = S0
    
    # Expectation E(e^J) for J~N(muJ, sigmaJ^2) - martingale correction
    EeJ = np.exp(muJ + 0.5 * sigmaJ * sigmaJ)
    
    # Generate random variables
    ZPois = np.random.poisson(xiP * dt, [NoOfPaths, NoOfSteps])
    Z = np.random.normal(0.0, 1.0, [NoOfPaths, NoOfSteps])
    J = np.random.normal(muJ, sigmaJ, [NoOfPaths, NoOfSteps])
    
    for i in range(0, NoOfSteps):
        # Making sure that samples from normal have mean 0 and variance 1
        if NoOfPaths > 1:
            Z[:, i] = (Z[:, i] - np.mean(Z[:, i])) / np.std(Z[:, i])
        
        # Merton SDE discretization with martingale correction
        X[:, i+1] = X[:, i] + (r - xiP * (EeJ - 1) - 0.5 * sigma * sigma) * dt + \
                    sigma * np.sqrt(dt) * Z[:, i] + J[:, i] * ZPois[:, i]
        time[i+1] = time[i] + dt
        
    S = np.exp(X)
    paths = {"time": time, "X": X, "S": S}
    return paths

def mainCalculation():
    """
    Main calculation routine for Merton jump-diffusion process visualization.
    
    Generates and plots sample paths showing both the log-price process X(t)
    and the stock price process S(t). The jumps are clearly visible as
    discontinuous changes in the paths.
    """
    NoOfPaths = 25
    NoOfSteps = 500
    T = 5
    xiP = 1        # Jump intensity: 1 jump per year on average
    muJ = 0        # Mean log-jump size
    sigmaJ = 0.7   # Std dev of log-jump size
    sigma = 0.2    # Diffusion volatility

    S0 = 100
    r = 0.05
    Paths = GeneratePathsMerton(NoOfPaths, NoOfSteps, S0, T, xiP, muJ, sigmaJ, r, sigma)
    timeGrid = Paths["time"]
    X = Paths["X"]
    S = Paths["S"]
           
    plt.figure(1)
    plt.plot(timeGrid, np.transpose(X), alpha=0.7)   
    plt.grid()
    plt.xlabel("Time")
    plt.ylabel("X(t) = log(S(t))")
    plt.title("Merton Jump-Diffusion: Log-Price Process")
    
    plt.figure(2)
    plt.plot(timeGrid, np.transpose(S), alpha=0.7)   
    plt.grid()
    plt.xlabel("Time")
    plt.ylabel("S(t)")
    plt.title(f"Merton Jump-Diffusion: Stock Price (λ={xiP}, σ_J={sigmaJ})")


if __name__ == "__main__":
    mainCalculation()