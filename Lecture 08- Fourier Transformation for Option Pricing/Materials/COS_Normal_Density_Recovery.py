#%%
"""
Density Recovery using the COS (Fourier Cosine) Method for Normal Distribution.

The COS method recovers probability density functions from their characteristic
functions using Fourier cosine expansion. This provides a fast and accurate
alternative to numerical integration methods.

Reference: Fang & Oosterlee (2008) "A Novel Pricing Method for European Options..."
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as st


def COSDensity(cf, x, N, a, b):
    """
    Recover probability density using Fourier cosine expansion.
    
    The COS method approximates the density f_X(x) using:
        f_X(x) ≈ Σ_{k=0}^{N-1} F_k * cos(k*π*(x-a)/(b-a))
    where F_k are Fourier coefficients computed from the characteristic function.
    
    Args:
        cf (callable): Characteristic function φ(u) = E[exp(iuX)].
        x (np.ndarray): Points at which to evaluate the density.
        N (int): Number of expansion terms (higher N = better accuracy).
        a (float): Lower truncation bound.
        b (float): Upper truncation bound.
    
    Returns:
        np.ndarray: Recovered density values at points x.
    
    Note:
        The truncation domain [a,b] should contain most of the probability mass.
        Typically a = μ - L*σ, b = μ + L*σ with L ∈ [8, 12].
    """
    i = np.complex(0.0, 1.0)  # Imaginary unit
    k = np.linspace(0, N-1, N)
    u = k * np.pi / (b - a)
        
    # Compute Fourier coefficients F_k
    F_k = 2.0 / (b - a) * np.real(cf(u) * np.exp(-i * u * a))
    F_k[0] = F_k[0] * 0.5  # Adjust first term (cosine series convention)
    
    # Reconstruct density using cosine expansion
    f_X = np.matmul(F_k, np.cos(np.outer(u, x - a)))
        
    return f_X


def mainCalculation():
    """
    Demonstrate density recovery for standard normal distribution.
    
    Tests COS method convergence with different numbers of expansion terms,
    comparing recovered density against exact normal PDF.
    """
    i = np.complex(0.0, 1.0)
    
    # COS method settings
    a = -10.0
    b = 10.0
    
    # Test with increasing number of expansion terms
    N = [2**x for x in range(2, 7, 1)]  # [4, 8, 16, 32, 64]
    
    # Normal distribution parameters
    mu = 0.0
    sigma = 1.0 
    
    # Characteristic function of N(μ, σ²): φ(u) = exp(iμu - 0.5σ²u²)
    cF = lambda u: np.exp(i * mu * u - 0.5 * np.power(sigma, 2.0) * np.power(u, 2.0))
    
    # Domain for density evaluation
    x = np.linspace(-10.0, 10, 1000)
    f_XExact = st.norm.pdf(x, mu, sigma)
    
    plt.figure(1, figsize=(12, 6))
    plt.grid(True, alpha=0.3)
    plt.xlabel("x")
    plt.ylabel("$f_X(x)$")
    plt.title("Density Recovery using COS Method: Standard Normal Distribution")
    
    # Plot exact density
    plt.plot(x, f_XExact, 'k-', linewidth=3, label='Exact PDF', alpha=0.7)
    
    # Recover density with different N and compute errors
    colors = ['r', 'g', 'b', 'm', 'c']
    print("=== COS Method Convergence Analysis ===")
    print(f"Distribution: N({mu}, {sigma**2})")
    print(f"Truncation domain: [{a}, {b}]\n")
    
    for idx, n in enumerate(N):
        f_X = COSDensity(cF, x, n, a, b)
        error = np.max(np.abs(f_X - f_XExact))
        print(f"N = {n:3d}: Max error = {error:.6e}")
        
        plt.plot(x, f_X, linestyle='--', color=colors[idx % len(colors)], 
                linewidth=1.5, label=f'COS N={n}')
    
    plt.legend()
    plt.xlim([-4, 4])
    
    print("\n✓ COS method provides exponential convergence")
    print("  Error decreases exponentially with increasing N")


if __name__ == "__main__":
    mainCalculation()
