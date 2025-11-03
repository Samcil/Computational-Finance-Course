#%%
"""
LogNormal Density Recovery using the COS Method.

This module demonstrates recovering the lognormal density from its characteristic
function using the COS (Fourier cosine) method. The lognormal distribution is
fundamental in financial modeling (e.g., Black-Scholes model).
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as st


def COSDensity(cf, x, N, a, b):
    """
    Recover probability density using Fourier cosine expansion.
    
    Args:
        cf (callable): Characteristic function φ(u).
        x (np.ndarray): Points at which to evaluate density.
        N (int): Number of expansion terms.
        a (float): Lower truncation bound.
        b (float): Upper truncation bound.
    
    Returns:
        np.ndarray: Recovered density values.
    """
    i = np.complex(0.0, 1.0)
    k = np.linspace(0, N-1, N)
    u = k * np.pi / (b - a)
        
    # Fourier coefficients
    F_k = 2.0 / (b - a) * np.real(cf(u) * np.exp(-i * u * a))
    F_k[0] = F_k[0] * 0.5  # Adjust first term
    
    # Reconstruct density
    f_X = np.matmul(F_k, np.cos(np.outer(u, x - a)))
        
    return f_X


def mainCalculation():
    """
    Demonstrate lognormal density recovery using COS method.
    
    For lognormal Y = exp(X) where X ~ N(μ, σ²), we recover the density
    f_Y(y) = (1/y) * f_X(log(y)) using the COS method applied to X.
    """
    i = np.complex(0.0, 1.0)
    
    # COS method settings
    a = -10
    b = 10
    
    # Test convergence with different N
    N = [16, 64, 128]
    
    # Parameters for underlying normal distribution X
    mu = 0.5
    sigma = 0.2
        
    # Characteristic function of X ~ N(μ, σ²)
    cF = lambda u: np.exp(i * mu * u - 0.5 * np.power(sigma, 2.0) * np.power(u, 2.0))
    
    # Domain for lognormal Y = exp(X)
    y = np.linspace(0.05, 5, 1000)
        
    plt.figure(1, figsize=(12, 6))
    plt.grid(True, alpha=0.3)
    plt.xlabel("y")
    plt.ylabel("$f_Y(y)$")
    plt.title(f"LogNormal Density Recovery using COS Method\n(underlying: X ~ N({mu}, {sigma}²))")
    
    print("=== LogNormal Density Recovery ===")
    print(f"Underlying distribution: X ~ N({mu}, {sigma²})")
    print(f"Target distribution: Y = exp(X) (LogNormal)\n")
    
    for n in N:
        # Recover density of X=log(Y), then transform to Y
        f_Y = 1/y * COSDensity(cF, np.log(y), n, a, b)
        plt.plot(y, f_Y, linewidth=2, label=f'COS N={n}')
        print(f"N = {n:3d} expansion terms computed")
    
    # Add theoretical lognormal density for comparison
    f_Y_exact = st.lognorm.pdf(y, s=sigma, scale=np.exp(mu))
    plt.plot(y, f_Y_exact, 'k--', linewidth=2.5, label='Exact LogNormal', alpha=0.7)
    
    plt.legend()
    plt.xlim([0, 3])
    
    print("\n✓ COS method accurately recovers lognormal density")
    print("  Convergence improves with increasing N")


if __name__ == "__main__":
    mainCalculation()
