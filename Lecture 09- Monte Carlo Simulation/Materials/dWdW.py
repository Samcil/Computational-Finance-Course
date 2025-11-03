#%%
"""
Convergence Analysis of Brownian Motion Increment Properties.

This module analyzes the convergence of the mean and variance of squared
Brownian motion increments (ΔW)² with respect to the number of discretization
intervals. This is important for understanding the quadratic variation of
Brownian motion.

Created on Jan 16 2019
@author: Lech A. Grzelak
"""
import numpy as np
import matplotlib.pyplot as plt


def mainCalculation():
    """
    Analyze convergence of E[(ΔW)²] and Var[(ΔW)²] as discretization refines.
    
    This function demonstrates that:
    1. E[(W(t_{i+1}) - W(t_i))²] → Δt as the number of intervals increases
    2. Var[(W(t_{i+1}) - W(t_i))²] → 2(Δt)² as the number of intervals increases
    
    These properties are fundamental to the theory of stochastic calculus and
    illustrate the quadratic variation of Brownian motion.
    
    Theoretical values:
        - E[(ΔW)²] = Δt
        - Var[(ΔW)²] = 2(Δt)²
    """
    NoOfPaths = 50000
    
    mV = []  # Mean values
    vV = []  # Variance values
    dtV = []  # Delta t values
    T = 1.0
    
    for m in range(2, 60, 1):
        t1 = 1.0 * T / m
        t2 = 2.0 * T / m
        dt = t2 - t1
        dtV.append(dt)
        
        # Generate Brownian motion at two consecutive time points
        W_t1 = np.sqrt(t1) * np.random.normal(0.0, 1.0, [NoOfPaths, 1])
        W_t2 = W_t1 + np.sqrt(dt) * np.random.normal(0.0, 1.0, [NoOfPaths, 1])
        
        # Compute squared increment
        X = np.power(W_t2 - W_t1, 2.0)
        mV.append(np.mean(X))
        vV.append(np.var(X))
    
    plt.figure(1, figsize=(10, 6))
    plt.plot(range(2, 60, 1), mV, 'b-', linewidth=2, label='E[(ΔW)²]')
    plt.plot(range(2, 60, 1), vV, 'r--', linewidth=2, label='Var[(ΔW)²]')
    plt.plot(range(2, 60, 1), dtV, 'g:', linewidth=2, label='Δt (theoretical mean)')
    plt.plot(range(2, 60, 1), [2*dt**2 for dt in dtV], 'm:', linewidth=2, 
             label='2(Δt)² (theoretical variance)')
    plt.grid()
    plt.xlabel('Number of discretization intervals (m)')
    plt.ylabel('Value')
    plt.title('Convergence of Brownian Motion Increment Properties')
    plt.legend()
    
    print(f'Final mean E[(ΔW)²] = {mV[-1]:.6f}, expected Δt = {dtV[-1]:.6f}')
    print(f'Final variance Var[(ΔW)²] = {vV[-1]:.8f}, expected 2(Δt)² = {2*dtV[-1]**2:.8f}')


if __name__ == "__main__":
    mainCalculation()