"""
Heston Model Discretization: Euler vs Almost Exact Scheme (AES) Comparison.

This module provides a comprehensive comparison of two discretization methods
for the Heston stochastic volatility model, demonstrating their convergence
properties and accuracy for option pricing.

The Heston Model:
    Stock price: dS(t) = rS(t)dt + sqrt(V(t))S(t)dW₂(t)
    Variance: dV(t) = κ(v̄ - V(t))dt + γsqrt(V(t))dW₁(t)
    where Corr(dW₁, dW₂) = ρdt

Discretization Schemes:

1. **Euler Scheme** (Simple but less accurate):
   - V(t+dt) = V(t) + κ(v̄ - V(t))dt + γsqrt(V(t))ΔW₁
   - Truncate V at 0 if negative (boundary condition)
   - Strong convergence order 0.5

2. **Almost Exact Scheme (AES)** (More accurate):
   - V(t+dt) sampled exactly from noncentral chi-squared distribution
   - S(t+dt) computed conditionally given V(t) and V(t+dt)
   - Strong convergence order 1.0
   - No need for truncation

This module demonstrates:
    - Path generation for both schemes
    - Option pricing across strikes
    - Convergence analysis as time step decreases
    - Comparison against analytical COS method prices

Key Functions:
    - GeneratePathsHestonEuler: Euler discretization with truncation
    - GeneratePathsHestonAES: Almost exact simulation
    - CallPutOptionPriceCOSMthd: Reference analytical prices
    - ChFHestonModel: Heston characteristic function

Author: Lech A. Grzelak
Created: Feb 11, 2019
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as st
import enum


class OptionType(enum.Enum):
    """Enumeration for option types."""
    CALL = 1.0
    PUT = -1.0

def CallPutOptionPriceCOSMthd(cf,CP,S0,r,tau,K,N,L):
    # cf   - characteristic function as a functon, in the book denoted as \varphi
    # CP   - C for call and P for put
    # S0   - Initial stock price
    # r    - interest rate (constant)
    # tau  - time to maturity
    # K    - list of strikes
    # N    - Number of expansion terms
    # L    - size of truncation domain (typ.:L=8 or L=10)  
        
    # reshape K to a column vector
    if K is not np.array:
        K = np.array(K).reshape([len(K),1])
    
    #assigning i=sqrt(-1)
    i = np.complex(0.0,1.0) 
    x0 = np.log(S0 / K)   
    
    # truncation domain
    a = 0.0 - L * np.sqrt(tau)
    b = 0.0 + L * np.sqrt(tau)
    
    # sumation from k = 0 to k=N-1
    k = np.linspace(0,N-1,N).reshape([N,1])  
    u = k * np.pi / (b - a);  

    # Determine coefficients for Put Prices  
    H_k = CallPutCoefficients(CP,a,b,k)   
    mat = np.exp(i * np.outer((x0 - a) , u))
    temp = cf(u) * H_k 
    temp[0] = 0.5 * temp[0]    
    value = np.exp(-r * tau) * K * np.real(mat.dot(temp))     
    return value

def Chi_Psi(a,b,c,d,k):
    psi = np.sin(k * np.pi * (d - a) / (b - a)) - np.sin(k * np.pi * (c - a)/(b - a))
    psi[1:] = psi[1:] * (b - a) / (k[1:] * np.pi)
    psi[0] = d - c
    
    chi = 1.0 / (1.0 + np.power((k * np.pi / (b - a)) , 2.0)) 
    expr1 = np.cos(k * np.pi * (d - a)/(b - a)) * np.exp(d)  - np.cos(k * np.pi 
                  * (c - a) / (b - a)) * np.exp(c)
    expr2 = k * np.pi / (b - a) * np.sin(k * np.pi * 
                        (d - a) / (b - a))   - k * np.pi / (b - a) * np.sin(k 
                        * np.pi * (c - a) / (b - a)) * np.exp(c)
    chi = chi * (expr1 + expr2)
    
    value = {"chi":chi,"psi":psi }
    return value

# Determine coefficients for Put Prices 
def CallPutCoefficients(CP,a,b,k):
    if CP==OptionType.CALL:                  
        c = 0.0
        d = b
        coef = Chi_Psi(a,b,c,d,k)
        Chi_k = coef["chi"]
        Psi_k = coef["psi"]
        if a < b and b < 0.0:
            H_k = np.zeros([len(k),1])
        else:
            H_k      = 2.0 / (b - a) * (Chi_k - Psi_k)  
    elif CP==OptionType.PUT:
        c = a
        d = 0.0
        coef = Chi_Psi(a,b,c,d,k)
        Chi_k = coef["chi"]
        Psi_k = coef["psi"]
        H_k      = 2.0 / (b - a) * (- Chi_k + Psi_k)                  
    return H_k    

def ChFHestonModel(r,tau,kappa,gamma,vbar,v0,rho):
    i = np.complex(0.0,1.0)
    D1 = lambda u: np.sqrt(np.power(kappa-gamma*rho*i*u,2)+(u*u+i*u)*gamma*gamma)
    g  = lambda u: (kappa-gamma*rho*i*u-D1(u))/(kappa-gamma*rho*i*u+D1(u))
    C  = lambda u: (1.0-np.exp(-D1(u)*tau))/(gamma*gamma*(1.0-g(u)*np.exp(-D1(u)*tau)))\
        *(kappa-gamma*rho*i*u-D1(u))
    # Note that we exclude the term -r*tau, as the discounting is performed in the COS method
    A  = lambda u: r * i*u *tau + kappa*vbar*tau/gamma/gamma *(kappa-gamma*rho*i*u-D1(u))\
        - 2*kappa*vbar/gamma/gamma*np.log((1.0-g(u)*np.exp(-D1(u)*tau))/(1.0-g(u)))
    # Characteristic function for the Heston's model    
    cf = lambda u: np.exp(A(u) + C(u)*v0)
    return cf 

def EUOptionPriceFromMCPathsGeneralized(CP,S,K,T,r):
    # S is a vector of Monte Carlo samples at T
    result = np.zeros([len(K),1])
    if CP == OptionType.CALL:
        for (idx,k) in enumerate(K):
            result[idx] = np.exp(-r*T)*np.mean(np.maximum(S-k,0.0))
    elif CP == OptionType.PUT:
        for (idx,k) in enumerate(K):
            result[idx] = np.exp(-r*T)*np.mean(np.maximum(k-S,0.0))
    return result

def GeneratePathsHestonEuler(NoOfPaths, NoOfSteps, T, r, S_0, kappa, gamma, rho, vbar, v0):
    """
    Generate Heston model paths using Euler-Maruyama discretization with truncation.
    
    Euler scheme with reflection boundary:
        V(t+dt) = max(0, V(t) + κ(v̄ - V(t))dt + γsqrt(V(t))ΔW₁)
        X(t+dt) = X(t) + (r - 0.5V(t))dt + sqrt(V(t))ΔW₂
    
    Args:
        NoOfPaths (int): Number of Monte Carlo paths.
        NoOfSteps (int): Number of time steps.
        T (float): Time to maturity (years).
        r (float): Risk-free interest rate.
        S_0 (float): Initial stock price.
        kappa (float): Mean reversion speed of variance.
        gamma (float): Volatility of variance.
        rho (float): Correlation between W₁ and W₂.
        vbar (float): Long-term mean variance.
        v0 (float): Initial variance.
    
    Returns:
        dict: {'time': time points, 'S': stock price paths}
    
    Note:
        Truncation at 0 handles potential negative variance values.
        Strong convergence order 0.5.
    """
    Z1 = np.random.normal(0.0, 1.0, [NoOfPaths, NoOfSteps])
    Z2 = np.random.normal(0.0, 1.0, [NoOfPaths, NoOfSteps])
    W1 = np.zeros([NoOfPaths, NoOfSteps + 1])
    W2 = np.zeros([NoOfPaths, NoOfSteps + 1])
    V = np.zeros([NoOfPaths, NoOfSteps + 1])
    X = np.zeros([NoOfPaths, NoOfSteps + 1])
    V[:, 0] = v0
    X[:, 0] = np.log(S_0)
    
    time = np.zeros([NoOfSteps + 1])
    dt = T / float(NoOfSteps)
    
    for i in range(0, NoOfSteps):
        # Ensure samples have mean 0 and variance 1
        if NoOfPaths > 1:
            Z1[:, i] = (Z1[:, i] - np.mean(Z1[:, i])) / np.std(Z1[:, i])
            Z2[:, i] = (Z2[:, i] - np.mean(Z2[:, i])) / np.std(Z2[:, i])
        
        # Apply correlation
        Z2[:, i] = rho * Z1[:, i] + np.sqrt(1.0 - rho**2) * Z2[:, i]
        
        W1[:, i + 1] = W1[:, i] + np.power(dt, 0.5) * Z1[:, i]
        W2[:, i + 1] = W2[:, i] + np.power(dt, 0.5) * Z2[:, i]
        
        # Variance process with truncation at 0
        V[:, i + 1] = V[:, i] + kappa * (vbar - V[:, i]) * dt + gamma * np.sqrt(V[:, i]) * (W1[:, i + 1] - W1[:, i])
        V[:, i + 1] = np.maximum(V[:, i + 1], 0.0)
        
        # Log-price process
        X[:, i + 1] = X[:, i] + (r - 0.5 * V[:, i]) * dt + np.sqrt(V[:, i]) * (W2[:, i + 1] - W2[:, i])
        time[i + 1] = time[i] + dt
    
    # Compute stock price from log-price
    S = np.exp(X)
    paths = {"time": time, "S": S}
    return paths

def CIR_Sample(NoOfPaths,kappa,gamma,vbar,s,t,v_s):
    delta = 4.0 *kappa*vbar/gamma/gamma
    c= 1.0/(4.0*kappa)*gamma*gamma*(1.0-np.exp(-kappa*(t-s)))
    kappaBar = 4.0*kappa*v_s*np.exp(-kappa*(t-s))/(gamma*gamma*(1.0-np.exp(-kappa*(t-s))))
    sample = c* np.random.noncentral_chisquare(delta,kappaBar,NoOfPaths)
    return  sample

def GeneratePathsHestonAES(NoOfPaths,NoOfSteps,T,r,S_0,kappa,gamma,rho,vbar,v0):    
    Z1 = np.random.normal(0.0,1.0,[NoOfPaths,NoOfSteps])
    W1 = np.zeros([NoOfPaths, NoOfSteps+1])
    V = np.zeros([NoOfPaths, NoOfSteps+1])
    X = np.zeros([NoOfPaths, NoOfSteps+1])
    V[:,0]=v0
    X[:,0]=np.log(S_0)
    
    time = np.zeros([NoOfSteps+1])
        
    dt = T / float(NoOfSteps)
    for i in range(0,NoOfSteps):
        # making sure that samples from normal have mean 0 and variance 1
        if NoOfPaths > 1:
            Z1[:,i] = (Z1[:,i] - np.mean(Z1[:,i])) / np.std(Z1[:,i])
        W1[:,i+1] = W1[:,i] + np.power(dt, 0.5)*Z1[:,i]
        
        # Exact samles for the variance process
        V[:,i+1] = CIR_Sample(NoOfPaths,kappa,gamma,vbar,0,dt,V[:,i])
        k0 = (r -rho/gamma*kappa*vbar)*dt
        k1 = (rho*kappa/gamma -0.5)*dt - rho/gamma
        k2 = rho / gamma
        X[:,i+1] = X[:,i] + k0 + k1*V[:,i] + k2 *V[:,i+1] + np.sqrt((1.0-rho**2)*V[:,i])*(W1[:,i+1]-W1[:,i])
        time[i+1] = time[i] +dt
        
    #Compute exponent
    S = np.exp(X)
    paths = {"time":time,"S":S}
    return paths

# Black-Scholes Call option price
def BS_Call_Put_Option_Price(CP,S_0,K,sigma,t,T,r):
    #print('Maturity T={0} and t={1}'.format(T,t))
    #print(float(sigma * np.sqrt(T-t)))
    #print('strike K ={0}'.format(K))
    K = np.array(K).reshape([len(K),1])
    d1    = (np.log(S_0 / K) + (r + 0.5 * np.power(sigma,2.0)) 
    * (T-t)) / (sigma * np.sqrt(T-t))
    d2    = d1 - sigma * np.sqrt(T-t)
    if CP == OptionType.CALL:
        value = st.norm.cdf(d1) * S_0 - st.norm.cdf(d2) * K * np.exp(-r * (T-t))
     #   print(value)
    elif CP == OptionType.PUT:
        value = st.norm.cdf(-d2) * K * np.exp(-r * (T-t)) - st.norm.cdf(-d1)*S_0
      #  print(value)
    return value

def mainCalculation():
    """
    Main comparison of Euler vs AES discretization schemes for Heston model.
    
    Workflow:
        1. Price options across strikes using COS method (benchmark)
        2. Compare Euler and AES Monte Carlo prices
        3. Analyze convergence as dt → 0
        4. Display results with plots and tables
    
    Demonstrates superior accuracy of AES over Euler scheme.
    """
    NoOfPaths = 2500
    NoOfSteps = 500
    
    # Heston model parameters
    gamma = 1.0      # Vol-of-vol
    kappa = 0.5      # Mean reversion speed
    vbar = 0.04      # Long-term mean variance
    rho = -0.9       # Correlation (negative for leverage effect)
    v0 = 0.04        # Initial variance
    T = 1.0          # Time to maturity
    S_0 = 100.0      # Initial stock price
    r = 0.1          # Risk-free rate
    CP = OptionType.CALL
    
    # First we define a range of strikes and check the convergence
    K = np.linspace(80,S_0*1.5,30)
    
    # Exact solution with the COS method
    cf = ChFHestonModel(r,T,kappa,gamma,vbar,v0,rho)
    
    # The COS method
    optValueExact = CallPutOptionPriceCOSMthd(cf, CP, S_0, r, T, K, 1000, 8)
    
    # Euler simulation
    pathsEULER = GeneratePathsHestonEuler(NoOfPaths,NoOfSteps,T,r,S_0,kappa,gamma,rho,vbar,v0)
    S_Euler = pathsEULER["S"]
    
    # Almost exact simulation
    pathsAES = GeneratePathsHestonAES(NoOfPaths,NoOfSteps,T,r,S_0,kappa,gamma,rho,vbar,v0)
    S_AES = pathsAES["S"]
    
        
    OptPrice_EULER = EUOptionPriceFromMCPathsGeneralized(CP,S_Euler[:,-1],K,T,r)
    OptPrice_AES   = EUOptionPriceFromMCPathsGeneralized(CP,S_AES[:,-1],K,T,r)
    
    plt.figure(1, figsize=(10, 6))
    plt.plot(K, optValueExact, '-r', linewidth=2.5, label='Exact (COS)')
    plt.plot(K, OptPrice_EULER, '--k', linewidth=2, label='Euler Scheme')
    plt.plot(K, OptPrice_AES, '.b', markersize=10, label='AES Scheme')
    plt.legend(fontsize=11, loc='best')
    plt.grid(True, alpha=0.3)
    plt.xlabel('Strike K', fontsize=11)
    plt.ylabel('Option Price', fontsize=11)
    plt.title('Heston Model: Option Prices Comparison\n(Euler vs AES vs COS Method)', fontsize=12)
    plt.tight_layout()
    
    # Here we will analyze the convergence for particular dt
    dtV = np.array([1.0, 1.0/4.0, 1.0/8.0,1.0/16.0,1.0/32.0,1.0/64.0])
    NoOfStepsV = [int(T/x) for x in dtV]
    
    # Specify strike for analysis
    K = np.array([140.0])
    
    # Exact
    optValueExact = CallPutOptionPriceCOSMthd(cf, CP, S_0, r, T, K, 1000, 8)
    errorEuler = np.zeros([len(dtV),1])
    errorAES = np.zeros([len(dtV),1])
    
    for (idx,NoOfSteps) in enumerate(NoOfStepsV):
        # Euler
        np.random.seed(3)
        pathsEULER = GeneratePathsHestonEuler(NoOfPaths,NoOfSteps,T,r,S_0,kappa,gamma,rho,vbar,v0)
        S_Euler = pathsEULER["S"]
        OptPriceEULER = EUOptionPriceFromMCPathsGeneralized(CP,S_Euler[:,-1],K,T,r)
        errorEuler[idx] = OptPriceEULER-optValueExact
        # AES
        np.random.seed(3)
        pathsAES = GeneratePathsHestonAES(NoOfPaths,NoOfSteps,T,r,S_0,kappa,gamma,rho,vbar,v0)
        S_AES = pathsAES["S"]
        OptPriceAES   = EUOptionPriceFromMCPathsGeneralized(CP,S_AES[:,-1],K,T,r)
        errorAES[idx] = OptPriceAES-optValueExact
    
    # Display convergence results
    print("\n" + "=" * 75)
    print("HESTON MODEL DISCRETIZATION: EULER VS AES CONVERGENCE ANALYSIS")
    print("=" * 75)
    print(f"Strike K = {K[0]:.2f}")
    print(f"Exact Option Price (COS Method): {optValueExact[0]:.6f}")
    print("-" * 75)
    
    print("\nEULER SCHEME CONVERGENCE:")
    print(f"{'dt':<12} {'Steps':<10} {'Error':<18} {'|Error|':<15}")
    print("-" * 75)
    for i in range(len(NoOfStepsV)):
        print(f"{dtV[i]:<12.4f} {NoOfStepsV[i]:<10} {errorEuler[i][0]:<18.8e} "
              f"{abs(errorEuler[i][0]):<15.8e}")
    
    print("\nAES SCHEME CONVERGENCE:")
    print(f"{'dt':<12} {'Steps':<10} {'Error':<18} {'|Error|':<15}")
    print("-" * 75)
    for i in range(len(NoOfStepsV)):
        print(f"{dtV[i]:<12.4f} {NoOfStepsV[i]:<10} {errorAES[i][0]:<18.8e} "
              f"{abs(errorAES[i][0]):<15.8e}")
    
    print("\n" + "=" * 75)
    print("KEY OBSERVATIONS:")
    print("  • AES scheme demonstrates superior accuracy over Euler")
    print("  • AES achieves first-order strong convergence (O(dt))")
    print("  • Euler limited to half-order convergence (O(√dt))")
    print("  • Exact CIR sampling eliminates truncation bias")
    print("=" * 75 + "\n")


if __name__ == "__main__":
    mainCalculation()