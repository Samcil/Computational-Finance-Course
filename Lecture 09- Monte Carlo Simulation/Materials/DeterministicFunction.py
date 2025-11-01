#%%
"""
Monte Carlo Integration of Deterministic Functions.

This module demonstrates two Monte Carlo methods for computing definite integrals
of deterministic functions:
1. Hit-or-Miss Monte Carlo: Uses acceptance-rejection sampling
2. Sample Mean Monte Carlo: Uses the law of large numbers

Created on Thu Nov 27 2018
@author: Lech A. Grzelak
"""
import numpy as np
import matplotlib.pyplot as plt


def ComputeIntegral1(NoOfSamples, a, b, c, d, g):
    """
    Compute integral using hit-or-miss Monte Carlo method.
    
    This method estimates the integral ∫_a^b g(x) dx by:
    1. Generating uniform random points (x_i, y_i) in rectangle [a,b] × [c,d]
    2. Counting points where y_i < g(x_i) (hits)
    3. Estimating integral as: (# hits / # total) × area of rectangle
    
    Args:
        NoOfSamples (int): Number of random samples to generate.
        a (float): Lower bound of integration domain.
        b (float): Upper bound of integration domain.
        c (float): Lower bound of y-range for sampling rectangle.
        d (float): Upper bound of y-range for sampling rectangle.
        g (callable): Function to integrate.
    
    Returns:
        float: Monte Carlo estimate of the integral.
    
    Note:
        Requires that c ≤ g(x) ≤ d for all x in [a,b] for accuracy.
        Convergence rate: O(1/√N)
    """
    x_i = np.random.uniform(a, b, NoOfSamples)
    y_i = np.random.uniform(c, d, NoOfSamples)
    
    p_i = g(x_i) > y_i
    p = np.sum(p_i) / NoOfSamples
    
    integral = p * (b - a) * (d - c)
    
    plt.figure(1)
    plt.plot(x_i, y_i, '.r', alpha=0.3, label='Sample points')
    plt.plot(x_i, g(x_i), '.b', alpha=0.5, label='Function g(x)')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Hit-or-Miss Monte Carlo Integration')
    plt.legend()
    plt.grid(True)
    
    return integral


def ComputeIntegral2(NoOfSamples, a, b, g):
    """
    Compute integral using sample mean Monte Carlo method.
    
    This method estimates the integral ∫_a^b g(x) dx by:
    1. Generating uniform random samples x_i in [a,b]
    2. Computing sample mean: mean(g(x_i))
    3. Estimating integral as: (b-a) × mean(g(x_i))
    
    Args:
        NoOfSamples (int): Number of random samples to generate.
        a (float): Lower bound of integration.
        b (float): Upper bound of integration.
        g (callable): Function to integrate.
    
    Returns:
        float: Monte Carlo estimate of the integral.
    
    Note:
        More efficient than hit-or-miss method for smooth functions.
        Based on: ∫_a^b g(x) dx = (b-a) × E[g(X)] where X ~ Uniform(a,b)
        Convergence rate: O(1/√N)
    """
    x_i = np.random.uniform(a, b, NoOfSamples)
       
    integral = (b - a) * np.mean(g(x_i))
    
    return integral


def main():
    """
    Main calculation routine demonstrating Monte Carlo integration.
    
    Computes ∫_0^1 e^x dx using both Monte Carlo methods and compares
    with the analytical solution e^1 - e^0 = e - 1 ≈ 1.71828.
    """
    NoOfSamples = 100000
    
    a = 0.0
    b = 1.0
    c = 0.0
    d = 3.0  # Upper bound should be > max(g(x)) on [a,b]
    
    g = lambda x: np.exp(x)
    
    # Method 1: Hit-or-Miss Monte Carlo
    output1 = ComputeIntegral1(NoOfSamples, a, b, c, d, g)
    print(f'Integral from Hit-or-Miss Monte Carlo: {output1:.6f}')
    
    # Method 2: Sample Mean Monte Carlo
    output2 = ComputeIntegral2(NoOfSamples, a, b, g)
    print(f'Integral from Sample Mean Monte Carlo: {output2:.6f}')
    
    # Analytical solution
    analytical = np.exp(b) - np.exp(a)
    print(f'Integral computed analytically: {analytical:.6f}')
    
    # Error analysis
    print(f'\nError in Method 1: {abs(output1 - analytical):.6f}')
    print(f'Error in Method 2: {abs(output2 - analytical):.6f}')


if __name__ == "__main__":
    main()
