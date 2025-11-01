# R Port Implementation Plan

## Overview
Porting 31 Python computational finance scripts to R following tidyverse and R package development best practices.

## Design Principles

### 1. Tidyverse Philosophy
- **Pipeable functions**: All functions work with `|>` pipe operator
- **Data-first arguments**: Data/paths argument comes first when appropriate
- **Consistent returns**: Return tibbles for tabular data, lists for complex structures
- **Tidy data**: One observation per row, one variable per column

### 2. Naming Conventions
- **Functions**: `snake_case` (e.g., `generate_gbm_paths`, `calculate_call_price`)
- **Variables**: `snake_case` (e.g., `n_paths`, `time_grid`, `stock_price`)
- **Files**: `snake_case.R` (e.g., `gbm_paths.R`, `heston_model.R`)
- **Meaningful names**: `n_paths` instead of `NoOfPaths`, `interest_rate` instead of `r`

### 3. Modular Structure
- **Separation of concerns**: Simulation, pricing, plotting in separate functions
- **Reusable components**: Common utilities in shared modules
- **Clear interfaces**: Well-defined inputs/outputs
- **Composability**: Functions that work together naturally

### 4. Documentation Standards (roxygen2)
```r
#' Generate Geometric Brownian Motion Paths
#'
#' Simulates stock price paths under GBM dynamics using Euler-Maruyama discretization.
#' The process follows dS(t) = r*S(t)*dt + sigma*S(t)*dW(t).
#'
#' @param n_paths Integer. Number of Monte Carlo paths to simulate.
#' @param n_steps Integer. Number of time steps for discretization.
#' @param maturity Numeric. Time to maturity in years.
#' @param interest_rate Numeric. Risk-free interest rate (annualized).
#' @param volatility Numeric. Volatility parameter (annualized).
#' @param initial_price Numeric. Initial stock price at t=0.
#' @param seed Integer. Random seed for reproducibility. Default is 123.
#'
#' @return A tibble with columns:
#'   \describe{
#'     \item{path_id}{Path identifier (1 to n_paths)}
#'     \item{time}{Time point}
#'     \item{stock_price}{Simulated stock price}
#'     \item{log_price}{Log of stock price (ABM process)}
#'   }
#'
#' @examples
#' # Simulate 100 paths over 1 year
#' paths <- generate_gbm_paths(
#'   n_paths = 100,
#'   n_steps = 252,
#'   maturity = 1.0,
#'   interest_rate = 0.05,
#'   volatility = 0.2,
#'   initial_price = 100
#' )
#'
#' @export
generate_gbm_paths <- function(n_paths, n_steps, maturity,
                                interest_rate, volatility, initial_price,
                                seed = 123) {
  # Implementation
}
```

### 5. Visualization Standards
- Use **ggplot2** for all visualizations
- Return plot objects (not print them)
- Consistent theme and styling
- Proper labels, titles, legends

## Implementation Phases

### Phase 1: Core Building Blocks (Priority 1)

#### 1.1 GBM/ABM Path Generation
**File**: `R/gbm_paths.R`
- `generate_gbm_paths()` - Main path generation
- `generate_abm_paths()` - ABM specifically
- `plot_paths()` - Visualization helper

**Python source**: `Lecture 03/GBM_ABM_paths.py`

#### 1.2 Poisson Process
**File**: `R/poisson_process.R`
- `generate_poisson_process()` - Poisson path generation
- `generate_compound_poisson()` - With jump sizes
- `plot_jump_process()` - Visualization

**Python source**: `Lecture 05/PoissonProcess_paths.py`

#### 1.3 Correlated Brownian Motion
**File**: `R/correlated_brownian_motion.R`
- `generate_correlated_bm()` - Multivariate BM
- `plot_correlation_paths()` - 2D/3D visualization

**Python source**: `Lecture 07/CorrelatedBM.py`

#### 1.4 CIR Process
**File**: `R/cir_process.R`
- `generate_cir_paths()` - CIR variance process
- `generate_cir_exact()` - Exact simulation (noncentral chi-squared)
- `plot_cir_paths()` - Visualization

**Python sources**: 
- `Lecture 10/CIR_paths_Exact.py`
- `Lecture 10/CIR_paths_boundary.py`

#### 1.5 Black-Scholes Pricing
**File**: `R/black_scholes.R`
- `bs_call_price()` - Call option price
- `bs_put_price()` - Put option price
- `bs_delta()` - Delta Greek
- `bs_gamma()` - Gamma Greek
- `bs_vega()` - Vega Greek

**Python sources**: Used across multiple files

### Phase 2: Advanced Models (Priority 2)

_Status: Complete (2025-11-04)_

#### 2.1 Heston Model
**Files**: 
- `R/heston_simulation.R` - Path generation
- `R/heston_pricing.R` - Option pricing via COS

**Functions**:
- `generate_heston_paths_euler()`
- `generate_heston_paths_aes()` - Almost Exact Scheme
- `heston_characteristic_function()`
- `price_heston_option_cos()`

**Python sources**:
- `Lecture 10/HestonModelDiscretization.py`
- `Lecture 10/CIR_ExactSimulation.py`
- `Lecture 10/OptionPrices_EulerAndMilstein.py`

#### 2.2 Merton Jump-Diffusion
**File**: `R/merton_model.R`
- `generate_merton_paths()` - Jump-diffusion simulation
- `merton_characteristic_function()`

**Python source**: `Lecture 05/MertonProcess_paths.py`

#### 2.3 Bates Model
**File**: `R/bates_model.R`
- `generate_bates_paths()` - Heston + jumps
- `bates_characteristic_function()`
- `bates_implied_volatility()`

**Python source**: `Lecture 12/BatesImpliedVolatility.py`

#### 2.4 COS Method
**File**: `R/cos_method.R`
- `cos_call_put_price()` - Generic COS pricing
- `cos_density_recovery()` - Density from CF
- `cos_coefficients()` - Fourier coefficients
- `chi_psi_functions()` - Auxiliary functions

**Python sources**:
- `Lecture 08/CallPut_COS_Method.py`
- `Lecture 08/COS_Normal_Density_Recovery.py`
- `Lecture 08/COS_LogNormal_Density_Recovery.py`
- `Lecture 08/CashOrNothing_COS_Method.py`

#### 2.5 Convergence Analysis
**File**: `R/convergence_analysis.R`
- `euler_convergence_study()` - Euler scheme convergence
- `milstein_convergence_study()` - Milstein scheme convergence
- `plot_convergence_rates()` - Visualization

**Python sources**:
- `Lecture 09/EulerConvergence_GBM.py`
- `Lecture 09/MilsteinConvergence_GBM.py`

> Deliverables implemented: Heston/Merton/Bates specs with COS pricing, convergence utilities, and associated tests are in place. Proceeding to Phase 3 feature work.

### Phase 3: Greeks & Hedging (Priority 3)

_Status: Complete (2025-11-04)_

#### 3.1 Implied Volatility
**File**: `R/implied_volatility.R`
- `implied_volatility_call()` - Newton-Raphson solver
- `implied_volatility_put()` - Put IV
- `plot_volatility_smile()` - IV smile visualization

**Python source**: `Lecture 04/ImpliedVolatility.py`

#### 3.2 Pathwise Sensitivities
**File**: `R/pathwise_sensitivities.R`
- `pathwise_delta()` - Pathwise delta estimation
- `pathwise_vega()` - Pathwise vega estimation
- `compare_with_finite_diff()` - Validation

**Python source**: `Lecture 11/PathwiseSens_DeltaVega.py`

#### 3.3 Delta Hedging
**Files**:
- `R/delta_hedging_bs.R` - BS delta hedging
- `R/delta_hedging_jumps.R` - Hedging with jumps

> Deliverables implemented: Implied volatility solvers, volatility smile plotting, pathwise sensitivities with finite-difference validation, and documentation updates. Discrete-time Black-Scholes hedging with transaction costs and jump-diffusion diagnostics now implemented with accompanying regression tests.

**Python sources**:
- `Lecture 11/BS_Hedging.py`
- `Lecture 11/HedgingWithJumps.py`

### Phase 4: Exotic Options (Priority 4)

_Status: Complete (2025-11-05)_

#### 4.1 Asian Options
**File**: `R/asian_options.R`
- `price_asian_call()` - Monte Carlo pricing
- `price_asian_put()` - Put options
- `asian_variance_reduction()` - Control variates

> Progress 2025-11-04: Implemented Monte Carlo pricing with antithetic variance reduction and added regression tests.

**Python source**: `Lecture 13/AsianOption.py`

#### 4.2 Barrier Options
**File**: `R/barrier_options.R`
- `price_barrier_option()` - Up/down, in/out pricing
- `barrier_hit_probability()` - Analysis

> Progress 2025-11-04: Added Monte Carlo pricing for barrier structures with complementarity checks and barrier hit probability estimator plus regression tests.

**Python source**: `Lecture 13/DigitalPayoffs_CostReduction.py` (partial)

#### 4.3 Digital Options
**File**: `R/digital_options.R`
- `price_digital_call()` - Cash-or-nothing
- `price_digital_put()` - Put digital
- `digital_cos_method()` - COS pricing

**Python source**: `Lecture 08/CashOrNothing_COS_Method.py`

> Progress 2025-11-05: Added Black-Scholes closed-form pricing, COS density integration, probability clamping, and regression tests spanning analytic parity and Fourier accuracy.

#### 4.4 Forward Start Options
**File**: `R/forward_start_options.R`
- `price_forward_start_heston()` - Forward start under Heston
- `forward_start_characteristic_function()`

**Python source**: `Lecture 12/HestonForwardStart2.py`

> Progress 2025-11-05: Implemented Heston forward-start characteristic function, COS pricing pipeline, and Monte Carlo regression tests; Phase 4 feature set complete.

### Phase 5: Package Infrastructure (Priority 5)

#### 5.1 Package Setup Files
- `DESCRIPTION` - Package metadata, dependencies
- `NAMESPACE` - Exported functions (auto-generated by roxygen2)
- `.Rbuildignore` - Files to exclude from package
- `LICENSE` - MIT or GPL-3

#### 5.2 Utility Functions
**File**: `R/utils.R`
- `create_time_grid()` - Time grid generation
- `standardize_normal()` - Normal random variables
- `format_paths_tibble()` - Convert to tidy format

**File**: `R/plotting_utils.R`
- `theme_compfinance()` - Custom ggplot2 theme
- `create_path_plot()` - Standard path visualization
- `create_convergence_plot()` - Convergence plot template

#### 5.3 Testing
**Directory**: `tests/testthat/`
- `test-gbm_paths.R` - Test GBM simulation
- `test-black_scholes.R` - Test BS formulas
- `test-cos_method.R` - Test COS pricing
- etc.

#### 5.4 Documentation
**Directory**: `vignettes/`
- `introduction.Rmd` - Getting started
- `basic_simulations.Rmd` - GBM, Poisson, etc.
- `heston_model.Rmd` - Stochastic volatility
- `cos_method.Rmd` - Fourier pricing
- `hedging_strategies.Rmd` - Greeks and hedging

#### 5.5 Data Documentation
**File**: `R/data.R`
- Document any included datasets
- Example market data for testing

## Key Differences from Python

### 1. Data Structures
- **Python**: NumPy arrays (matrices)
- **R**: Tibbles (tidy data frames) for analysis, matrices for computation

### 2. Random Number Generation
- **Python**: `np.random.normal()`, `np.random.seed()`
- **R**: `rnorm()`, `set.seed()`, `withr::with_seed()` for local seeds

### 3. Plotting
- **Python**: matplotlib
- **R**: ggplot2 with consistent theme

### 4. Linear Algebra
- **Python**: NumPy
- **R**: Base R matrices, MASS::mvrnorm for multivariate

### 5. Function Organization
- **Python**: Scripts with `if __name__ == "__main__":`
- **R**: Exported package functions + unexported helpers

## Dependencies

### Required R Packages
```r
Imports:
  tibble (>= 3.0.0),
  dplyr (>= 1.0.0),
  ggplot2 (>= 3.3.0),
  purrr (>= 0.3.0),
  tidyr (>= 1.1.0),
  rlang (>= 0.4.0),
  MASS,
  stats

Suggests:
  testthat (>= 3.0.0),
  knitr,
  rmarkdown,
  withr
```

## Implementation Timeline

### Week 1: Foundation (Phase 1)
- Set up package structure
- Implement core path generation (GBM, ABM, Poisson, CIR)
- Create plotting utilities
- Write initial tests

### Week 2: Advanced Models (Phase 2)
- Implement Heston model
- Implement Merton and Bates models
- Implement COS method
- Convergence analysis tools

### Week 3: Greeks & Hedging (Phase 3)
- Black-Scholes Greeks
- Implied volatility
- Pathwise sensitivities
- Hedging strategies

### Week 4: Exotic Options (Phase 4)
- Asian options
- Barrier options
- Digital options
- Forward start options

### Week 5: Finalization (Phase 5)
- Complete testing suite
- Write vignettes
- Generate documentation
- Package checks (R CMD check)

## Success Criteria

1. ✅ All 31 Python scripts ported to idiomatic R
2. ✅ 100% roxygen2 documentation coverage
3. ✅ All functions follow tidyverse principles
4. ✅ Passes R CMD check with no errors/warnings
5. ✅ Comprehensive test suite (>80% coverage)
6. ✅ Complete vignettes for all major functionality
7. ✅ Consistent naming and code style
8. ✅ DRY principles applied (no code duplication)

## Notes

- Focus on **correctness** over performance initially
- Use **tidyverse** approaches even if base R is faster
- Prioritize **readability** and **maintainability**
- Include **examples** in all documentation
- Write **tests** alongside implementation
- Keep Python files unchanged (separate R implementation)
