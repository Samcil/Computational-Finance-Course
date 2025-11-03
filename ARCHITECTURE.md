# CompFinanceR Package Design Blueprint

## Executive Summary
- Deliver a unified, extensible R package that reproduces every Python workflow across the Computational Finance Course and FinancialEngineering_IR_xVA repositories while embracing tidymodels, hardhat, and actuarial standards.
- Present a declarative specification hierarchy that mirrors model families (diffusions, jumps, stochastic volatility, short-rate, hybrid, mortgage, risk analytics) so that users can navigate by process lineage (for example `gbm_spec` inheriting from `abm_spec`, or the Ornstein Uhlenbeck family covering Ho-Lee, Hull-White, and G2++).
- Encapsulate computational detail inside pluggable engines with deterministic interfaces, enabling high performance via vectorised R, optional `Rcpp` fallbacks, and future parallel back ends without altering user-facing syntax.
- Provide curated documentation, vignettes, and print methods that surface the family tree of specs, available engines, supported derivatives, and calibration pathways.

## Specification Hierarchy

| Layer | Description | Example Classes |
| --- | --- | --- |
| `process_spec` (base) | Minimal contract for all stochastic processes, storing validated args, engine metadata, and helper state | inherited by every spec |
| `diffusion_spec` | Continuous diffusion family with drift and diffusion slots plus Brownian driver metadata | `abm_spec`, `gbm_spec`, `cir_spec` |
| `jump_diffusion_spec` | Extends diffusion with Poisson jump structure and jump distribution metadata | `merton_spec`, `bates_spec` |
| `stochastic_vol_spec` | Adds variance process parameters and correlation matrices | `heston_spec`, `szhw_spec` |
| `short_rate_spec` (OU family) | Mean reverting short-rate models with curve splines, theta function closures, and analytic ZCB engines | `ho_lee_spec`, `hull_white_spec`, planned `g2pp_spec` |
| `cir_family_spec` | Square-root diffusion family sharing exact transition samplers and moment functions | `cir_spec`, variance leg in Heston and Bates |
| `hybrid_spec` | Joint factors (equity-rate, FX-rate) with block covariance metadata | `bshw_spec`, `h1_hw_spec`, `fx_short_rate_spec` |
| `mortgage_spec` | Cashflow amortisation and prepayment drivers referencing short-rate engines | `mortgage_annuity_spec`, `amortizing_swap_spec` |
| `exposure_spec` | Portfolio and netting set definitions layered on top of valuation specs | `swap_exposure_spec`, `netting_exposure_spec` |
| `risk_measure_spec` | VaR / ES calculators referencing exposure sims | `var_spec`, `expected_shortfall_spec` |

### Class Stacking and Heritage
- Use S3 multiple inheritance (class vectors) so that `class(gbm_spec())` becomes `c("gbm_spec", "diffusion_spec", "process_spec", "model_spec")`.
- Provide helper constructors (e.g., `new_diffusion_spec()`) wrapping `new_model_spec()` with shared validation.
- Offer `family()` helper that inspects the class vector and returns lineage plus available derivative/pricing helpers for discoverability.
- Register `set_engine()` methods at the highest family level possible (for example OU family covers Ho-Lee and Hull-White) to keep engines DRY.

## Engine Layer and Performance Strategy
- Engines remain pure functions returning tibbles with deterministic column sets. Expose metadata (step size, solver tolerance, RNG backend) via attributes.
- Default implementations rely on vectorised R and purrr; introduce optional `Rcpp` modules for heavy loops (e.g., Jamshidian root search, large Monte Carlo). Keep interfaces identical to allow drop-in replacements.
- Provide `engine_registry()` returning a tibble listing engines by spec, method (`simulate`, `price`, `fit`), and performance characteristics.
- For parallelism offer `future` or `mirai` integration behind `purrr::in_parallel()` wrappers, enabling large scenario batches without rewriting user code.

## Module Plan by Domain

### 1. Equity and Foundational Diffusions (Computational-Finance-Course Lectures 03–05)
- **Complete 2025-11-03** Implement `abm_spec`, `gbm_spec`, `poisson_spec`, `merton_spec`, and `correlated_bm_spec` with inheritance and shared Brownian helpers.
- **Complete 2025-11-03** Reuse deterministic integrator utilities for Euler and Milstein schemes derived from `EulerConvergence_GBM.py`, `MilsteinConvergence_GBM.py`, and `DeterministicFunction.py`.
- **Complete 2025-11-03** Provide measure switching utilities (`q_measure_paths()`) aligned with `PathsUnderQandPmeasure.py`.
- **Parity verified 2025-11-03** `testthat/test-simulate-paths.R`, `testthat/test-deterministic-integration.R`, and `testthat/test-delta-hedging.R` reproduce Python reference statistics from `GBM_ABM_paths.py`, `DeterministicFunction.py`, and `PathsUnderQandPmeasure.py`.

### 2. Fourier and Transform Pricing (Lecture 08 assets)
- **Complete 2025-11-03** Centralise COS method kernels in `cos_method.R` and expose wrappers for plain vanilla, digital, density recovery (mapping the Python COS suite).
- **Complete 2025-11-03** Handle special characters (e.g., superscripts in `COS_LogNormal_Density_Recovery.py`) during import by normalising unicode to ASCII in documentation and tests.
- **Parity verified 2025-11-03** `tests/testthat/test-bates-model.R`, `tests/testthat/test-forward-start-options.R`, and `tests/testthat/test-digital-options.R` match Python payoffs from `CallPut_COS_Method.py`, `CashOrNothing_COS_Method.py`, and `HestonForwardStart2.py` within tolerance.

### 3. Stochastic Volatility and Jump Diffusions (Lectures 10–12)
- **Complete 2025-11-03** Extend `stochastic_vol_spec` to supply AES and Euler engines (from `HestonModelDiscretization.py`, `CIR_ExactSimulation.py`), ensuring `simulate_paths()` dispatch covers both schemes.
- **Complete 2025-11-03** Implement `bates_spec` and forward-start characteristic functions referencing `BatesImpliedVolatility.py` and `HestonForwardStart2.py`; characteristic functions now share the spec metadata infrastructure and remain numeric-safe.
- **Complete 2025-11-03** Provide calibration workflows for implied volatility inversion consistent with `ImpliedVolatility.py` and `ImpliedVolatility_FrwdStart`, reusing `compute_implied_volatility()` and validating via `tests/testthat/test-bates-model.R`.
- **Parity verified 2025-11-03** Full suite coverage in `tests/testthat/test-bates-model.R`, `tests/testthat/test-simulate-paths.R`, and `tests/testthat/test-pathwise_sensitivities.R` confirms AES/Euler paths, COS prices, and Greeks against `HestonModelDiscretization.py` and `BatesImpliedVolatility.py` benchmarks.

### 4. Short-Rate, Term Structure, and Interest Rate Products (FinancialEngineering_IR_xVA Lectures 03–13)
- **Complete 2025-11-03** Expand current term-structure and short-rate modules to include multi-curve Newton solvers (`MultiCurveBuild.py`), treasury bootstraps, and Greek sensitivities (`YieldCurveBuildGreeks.py`).
- **Complete 2025-11-03** Build OU family spec lineage: `ho_lee_spec -> hull_white_spec -> g2pp_spec`, all sharing theta, A/B factor, and analytic discount engines from the lecture scripts, with coverage enforced by `tests/testthat/test-term-structure-spec.R`.
- **Complete 2025-11-03** Extend caplet, floorlet, swaption, and displaced diffusion pricing per `HW_Caplets.py`, `JamshidianTrick.py`, `ShiftedLognormal.py`, `DD_ImpliedVolatility.py`, backed by the new pricing helpers and calibration tests in `tests/testthat/test-caplet-floorlet.R`, `tests/testthat/test-swaption.R`, and `tests/testthat/test-shifted-lognormal.R`.
- **Complete 2025-11-03** Provide mortgage and amortizing swap specs mapping to
  `AnnuityMortgage.py` and `StochasticAmortizingSwap.py` via
  `mortgage_annuity_spec()`, `amortizing_swap_spec()`,
  `price_mortgage()`, and the shared schedule helper covered by
  `tests/testthat/test-mortgage.R`.
- **Parity verified 2025-11-03** Term structure, swaption, and mortgage suites (`tests/testthat/test-term-structure-spec.R`, `tests/testthat/test-caplet-floorlet.R`, `tests/testthat/test-swaption.R`, `tests/testthat/test-shifted-lognormal.R`, `tests/testthat/test-mortgage.R`) match output matrices from `Ho-Lee-ZCBs.py`, `JamshidianTrick.py`, `ShiftedLognormal.py`, and `AnnuityMortgage.py`.

### Python Parity Audit (Modules 1–4)
| Module | Python references reviewed | R coverage & validation |
| --- | --- | --- |
| 1. Diffusions & jumps | `GBM_ABM_paths.py`, `GBM_ABM_paths_Martingale.py`, `DeterministicFunction.py`, `EulerConvergence_GBM.py`, `MilsteinConvergence_GBM.py`, `MertonProcess_paths.py`, `PoissonProcess_paths.py`, `PathsUnderQandPmeasure.py` | Engines in `simulate_paths()` and jump utilities matched stochastic differential forms; parity confirmed via `testthat/test-simulate-paths.R`, `testthat/test-deterministic-integration.R`, `testthat/test-delta-hedging.R`, and fixture statistics against the Python grids. |
| 2. Fourier pricing | `CallPut_COS_Method.py`, `CashOrNothing_COS_Method.py`, `COS_Normal_Density_Recovery.py`, `BatesImpliedVolatility.py`, `HestonForwardStart2.py` (log-normal density script manually transcribed to ASCII) | `cos_method.R`, forward-start and digital pricing pipelines reproduce Python payoff vectors; `tests/testthat/test-bates-model.R`, `tests/testthat/test-forward-start-options.R`, and `tests/testthat/test-digital-options.R` assert coefficient-level tolerance. |
| 3. Stochastic volatility | `HestonModelDiscretization.py`, `CIR_ExactSimulation.py`, `CIR_paths_Exact.py`, `OptionPrices_EulerAndMilstein.py`, `BatesImpliedVolatility.py` | AES/Euler simulators and implied-vol routines in `stochastic_vol_spec` align with Python characteristic functions; `tests/testthat/test-bates-model.R`, `tests/testthat/test-simulate-paths.R`, `tests/testthat/test-pathwise_sensitivities.R` compare path moments and option prices to script outputs. |
| 4. Short-rate & term structure | `Ho-Lee-ZCBs.py`, `Hull-White-Paths.py`, `Hull-White-CompRateSim.py`, `Hull-White-ZCBs2.py`, `MultiCurveBuild.py`, `YieldCurveBuildGreeks.py`, `HW_Caplets.py`, `HW_OptionsOnZCBs.py`, `JamshidianTrick.py`, `ShiftedLognormal.py`, `AnnuityMortgage.py`, `StochasticAmortizingSwap.py` | OU family analytics, calibration solvers, and mortgage/amortising swap tooling in `short_rate_spec`, `term_structure_sensitivities.R`, and `mortgage_products.R` reproduce the Python schedules and prices. Verified through `tests/testthat/test-term-structure-spec.R`, `test-caplet-floorlet.R`, `test-swaption.R`, `test-shifted-lognormal.R`, and `test-mortgage.R`. |
| 5. Hybrid equity/rate | `BSHW_Comparison.py` | `bshw_spec`, `price_bshw_option_cos()`, and `bshw_equivalent_volatility()` replicate COS and Black-76 valuations; enforced by `tests/testthat/test-bshw.R`. |

### 5. Hybrid Equity/FX-Rate Models (FinancialEngineering_IR_xVA Lectures 09–10)
- **Complete 2025-11-04** Delivered `bshw_spec`, `cos_call_put_price_stoch_ir()`, and Black-76 parity helpers matching `BSHW_Comparison.py`. Coverage enforced by `tests/testthat/test-bshw.R`.
- **In progress** Extend hybrid utilities to `h1_hw_spec` (Heston-Hull-White) and `szhw_spec`, reusing the COS infrastructure and validating against `H1_HW_COS_vs_MC.py` and `SZHW_ImpliedVolatilities.py`.
- **In progress** Add FX-layered hybrid spec built on domestic/foreign short-rate curves per `H1_HW_COS_vs_MC_FX.py`.

### 6. Netting, Exposure, and Risk Analytics (Lectures 11–13)
- Implement exposure specs that bundle valuation specs and schedule lists to compute pathwise exposure (`Exposures_HW_Netting.py`).
- Provide convexity adjustments (`ConvexityCorrection.py`), VaR and ES calculators (`MonteCarloVaR.py`, `HistoricalVaR_Calculation.py`), and displaced diffusion smile tools.

## Package Dependencies and Tooling
- Core: `hardhat`, `parsnip` (for consistent engine registration), `cli`, `glue`, `checkmate`, `rlang`, `purrr`, `tibble`, `dplyr`, `tidyr`, `ggplot2`.
- Numerical: `Matrix`, `pracma`, `numDeriv`, `stats`, `uniroot`, `nloptr` (for constrained calibration), `splines`.
- Randomness and parallel: `withr`, `mirai`, `future`, `progressr`.
- Optional accelerators: `Rcpp`, `RcppParallel` for Jamshidian root searches and bulk Monte Carlo; keep pure R fallback.
- Testing: `testthat`, `vdiffr` (plots), `waldo`, `snapshot` fixtures.

## Performance Guardrails
- Ensure deterministic seeds via `with_random_seed()` wrapper (thin over `withr::with_seed`).
- Provide chunked Monte Carlo evaluation to support millions of paths without exceeding memory.
- Cache analytic helper closures (theta functions, A/B factors) within spec engine state to avoid recomputation per call.
- Offer toggles for vectorised vs chunked simulation in `simulate_paths()` allowing users to trade memory for speed.

## Testing and Quality Assurance
- Stage tests by family: `testthat/test-diffusion-family.R`, `testthat/test-ou-family.R`, etc.
- Incorporate parity fixtures from Python outputs (store CSV in `tests/fixtures/`) for zero-drift comparatives and cross-curve valuations.
- Add property-based tests (with `quickcheck`) ensuring discount factors stay in (0,1], swap PV parity holds, and smile calibrations are monotonic.
- Snapshot CLI printouts to guarantee spec summaries list ancestry and available engines.

## Documentation and Discoverability
- Generate roxygen documentation that emphasises family membership, available engines, supported derivatives, and calibration recipes.
- Produce Quarto vignettes grouped by family (Diffusions, Short-Rate, Hybrid, Mortgages, Risk) with cross-links to Python source references.
- Add `pkgdown` reference index grouping functions by family for fast navigation.
- Provide `compfinance_map()` helper returning a tibble describing every spec, engine, derivative, and associated vignette to serve as an interactive guide.

## Implementation Phases (New Schedule)
1. **Foundations Refresh**: Finalise spec inheritance utilities, engine registry, and print methods. Migrate existing specs onto the new hierarchy.
2. **Diffusion and Jump Suite**: Port Lecture 03–05 Python scripts, embed convergence diagnostics, and finish documentation.
3. **Fourier and Stochastic Volatility**: Deliver COS infrastructure, Heston/Bates specs, implied volatility tools.
4. **Term Structure and Short-Rate**: Extend term_structure engines, OU/CIR families, instrument pricing, and calibrations.
5. **Hybrid and FX**: Implement equity-rate and FX-rate hybrids with COS and Monte Carlo parity tests.
6. **Mortgages and Amortising Products**: Introduce mortgage specs, prepayment models, amortising swaps, plus corresponding examples.
7. **Exposure and Risk**: Build netting, convexity, VaR/ES specifications and integrate with actuarial reporting workflows.
8. **Polish and Benchmarks**: Profiling, optional `Rcpp` bridges, `pkgdown` site, and release candidate validation.

## Appendix A: Python Source Inventory and R Mapping

### FinancialEngineering_IR_xVA Repository

#### Lecture 02 – Understanding of Filtrations and Measures
- `Black_Scholes_Jumps.py` (functions: GeneratePaths, EUOptionPriceFromMCPaths, BS_Call_Put_Option_Price, CallOption_CondExpectation, mainCalculation; class OptionType) → informs diffusion and jump option pricing helpers plus martingale diagnostics.
- `Martingale.py` (martingaleA, martingaleB) → guides measure-consistency checks and test fixtures for Brownian martingales.

#### Lecture 03 – The HJM Framework
- `CIR_IR_paths.py` (GeneratePathsCIREuler, mainCalculation) → shapes `cir_spec` Euler engines and variance process utilities.
- `Ho-Lee-ZCBs.py` (f0T, GeneratePathsHoLeeEuler, mainCalculation) → seeds Ho-Lee spec analytic curves and Euler simulator.
- `Hull-White-Paths.py` (GeneratePathsHWEuler, mainCalculation) → baseline for Hull-White Euler engine.
- `Hull-White-ZCBs.py` (f0T, GeneratePathsHWEuler, HW_theta, mainCalculation) → informs theta function construction and discount factor calculators.

#### Lecture 04 – Yield Curve Dynamics under Short Rate
- `Hull-White-CompRateSim.py` (f0T, GeneratePathsHWEuler, HW_theta, HW_A, HW_B, HW2F_ZCB, HW_ZCB, HW_r_0, mainCalculation) → drives OU family analytics, multi-factor extensions, and composite rate simulations.
- `Hull-White-ZCBs2.py` (f0T, GeneratePathsHWEuler, HW_theta, HW_A, HW_B, HW_ZCB, mainCalculation) → supports analytic zero curve calculators and validation tests.
- `Hull_White_1F_2F_Comparison.py` (GeneratePathsHW2FEuler, f0T, GeneratePathsHWEuler, HW_theta, HW_A, HW_B, HW2F_ZCB, HW_ZCB, HW_r_0, mainCalculation) → informs two-factor spec design and comparison utilities.

#### Lecture 05 – Interest Rate Products
- `HW_Caplets.py` (GeneratePathsHWEuler, HW_theta, HW_A, HW_B, HW_ZCB, HWMean_r, HW_r_0, HW_Mu_FrwdMeasure, HWVar_r, HWDensity, HW_CapletFloorletPrice, HW_ZCB_CallPutPrice, mainCalculation; OptionType) → maps to caplet/floorlet pricing engines and density validation.
- `HW_OptionsOnZCBs.py` (similar function set) → provides bond option analytics for validation of `price_zcb_option()`.
- `Swaps_HW.py` (GeneratePathsHWEuler, theta/A/B/Bond helpers, HW_r_0, SwapPrice, HW_SwapPrice, mainCalculation; OptionTypeSwap) → forms swap pricing spec, DV01 analytics, and swaption seeds.

#### Lecture 06 – Construction of Yield Curve and Multi-Curve
- `MultiCurveBuild.py` (IRSwap, IRSwapMultiCurve, P0TModel, YieldCurve, MultivariateNewtonRaphson, Jacobian, EvaluateInstruments, interpolation helpers, mainCode; OptionTypeSwap) → blueprint for multi-curve calibration engine and Jacobian assembly.
- `YieldCurveBuildGreeks.py` (IRSwap, P0TModel, YieldCurve, MultivariateNewtonRaphson, Jacobian, EvaluateInstruments, interpolation variants, BuildInstruments, mainCode; OptionTypeSwap) → extends calibration to sensitivities and instrument builders.
- `YieldCurveBuild_Treasury.py` (IRSwap, P0TModel, YieldCurve, Newton solver, interpolation, mainCode; OptionTypeSwap) → informs treasury bootstrapping workflows.

#### Lecture 07 – Pricing of Swaptions and Negative Interest Rates
- `HW_CapletsAndFloorlets.py` (caplet set plus BS_Call_Put_Option_Price, ImpliedVolatilityBlack76, mainCalculation; OptionType) → calibrates Black-76 interfaces and Monte Carlo vs analytic comparisons.
- `JamshidianTrick.py` (PsiSum, JamshidianTrick, Main) → direct mapping to `jamshidian_price()` helper within swaption module.
- `ShiftedLognormal.py` (GeneratePathsGBMShifted, GeneratePathsGBM, shifted and standard Black formulas, ImpliedVolatilityBlack76/Shifted, mainCalculation; OptionType) → basis for displaced diffusion analytics.

#### Lecture 08 – Mortgages and Prepayments
- `AnnuityMortgage.py` (Annuity, mainCode) → guides mortgage amortisation schedule builder.
- `BulletMortgage.py` (Bullet, mainCode) → informs bullet schedule helper.
- `Incentives.py` (Annuity, mainCode) → adds prepayment incentive modelling references.
- `StochasticAmortizingSwap.py` (GeneratePathsHWEuler, theta/A/B/Bond helpers, SwapRateHW, Bullet, Annuity, mainCode) → blueprint for amortising swap spec integrating rate paths and mortgage schedules.

#### Lecture 09 – Hybrid Models and Stochastic Interest Rates
- `BSHW_Comparison.py` (COS-based pricing, characteristic functions, BS benchmarks, mainCalculation; OptionType) → calibrates BSHW COS engine and implied volatility checks.
- `BSHW_ImpliedVolatility.py` (BS pricing, BSHW volatility solver, mainCalculation; OptionType) → informs hybrid IV routines.
- `H1_HW_COS_vs_MC.py` (COS, Monte Carlo, Heston-HW AES/Euler, CIR sampler, characteristic functions, mainCalculation; OptionType) → defines hybrid engine requirements and accuracy comparisons.
- `SZHW_ImpliedVolatilities.py` (COS hybrid pricing, multiple characteristic functions, mainCalculation; OptionType) → guides SZHW spec.
- `SZHW_MonteCarlo_DiversificationProduct.py` (theta/A/B/Bond helpers, hybrid MC paths, DiversifcationPayoff, mainCalculation; OptionType) → extends Monte Carlo diversification analytics for spec design.

#### Lecture 10 – Foreign Exchange (FX) and Inflation
- `H1_HW_COS_vs_MC_FX.py` (Black pricing, COS, Monte Carlo FX forward measure, H1-HW FX characteristic functions, strike generators, mainCalculation; OptionType) → blueprint for FX hybrid spec and pricing engines.

#### Lecture 11 – Market Model and Convexity Adjustments
- `ConvexityCorrection.py` (theta/A/B/Bond helpers, HW mean/variance, mainCalculation; OptionType) → foundation for convexity adjustment utilities.
- `DD_ImpliedVolatility.py` (BS pricing, displaced diffusion IV, call price, mainCalculation; OptionType) → drives displaced diffusion calibration spec.

#### Lecture 12 – Valuation Adjustments (xVA)
- `Exposures_HW_Netting.py` (theta/A/B/Bond helpers, mean/variance, density, HW_SwapPrice, mainCalculation; OptionTypeSwap) → guides portfolio exposure simulation and netting analytics.

#### Lecture 13 – Value-at-Risk and Expected Shortfall
- `HistoricalVaR_Calculation.py` (IRSwap, P0TModel, YieldCurve, Newton solver, interpolation, BuildYieldCurve, Portfolio, mainCode; OptionTypeSwap) → informs historical VaR workflow and curve building for risk analytics.
- `MonteCarloVaR.py` (theta/A/B/Bond helpers, mean/variance, density, swap pricing, Portfolio, mainCalculation; OptionTypeSwap) → blueprint for Monte Carlo VaR spec.

### Computational-Finance-Course Repository

#### Lecture 03 – Option Pricing and Simulation in Python
- `GBM_ABM_paths.py` (GeneratePathsGBMABM, mainCalculation) → base for ABM/GBM spec constructors.
- `GBM_ABM_paths_Martingale.py` (GeneratePathsGBMABM, mainCalculation) → supplies martingale validation examples.
- `PathsUnderQandPmeasure.py` (GeneratePathsGBM, MainCode) → informs measure change utilities and tests.

#### Lecture 04 – Implied Volatility
- `ImpliedVolatility.py` (ImpliedVolatility, dV_dsigma, BS_Call_Option_Price, main) → shapes Newton and derivative-based IV solvers used package-wide.

#### Lecture 05 – Jump Processes
- `MertonProcess_paths.py` (GeneratePathsMerton, mainCalculation) → maps to jump diffusion spec and Monte Carlo engine.
- `PoissonProcess_paths.py` (GeneratePathsPoisson, mainCalculation) → provides Poisson and compound Poisson utils.

#### Lecture 07 – Stochastic Volatility Models
- `CorrelatedBM.py` (GeneratePathsCorrelatedBM, mainCalculation) → informs correlated Brownian helper already ported.

#### Lecture 08 – Fourier Transformation for Option Pricing
- `CallPut_COS_Method.py` (CallPutOptionPriceCOSMthd, coefficients, Chi_Psi, BS pricing, mainCalculation) → baseline for COS method implementation.
- `CashOrNothing_COS_Method.py` (CashOrNothingPriceCOSMthd, coefficients, Chi_Psi, BS cash-or-nothing, mainCalculation) → extends COS to digital payoffs.
- `COS_LogNormal_Density_Recovery.py` (parse blocked by unicode superscript) → requires manual transcription for density recovery; plan to normalise characters and import formulas.
- `COS_Normal_Density_Recovery.py` (COSDensity, mainCalculation) → density recovery baseline.
- `DensityRecoveryFFT.py` (RecoverDensity, mainCalculation) → informs FFT density inversion helper.

#### Lecture 09 – Monte Carlo Simulation
- `DeterministicFunction.py` (ComputeIntegral1, ComputeIntegral2, main) → motivates deterministic integral comparison utilities for Monte Carlo convergence tests.
- `dWdW.py` (mainCalculation) → demonstrates correlated Brownian increments, reused in correlated helpers.
- `EulerConvergence_GBM.py` (GeneratePathsGBMEuler, mainCalculation) → informs Euler convergence diagnostics.
- `Exercise_1.py` and `Exercise_2.py` (ComputeIntegrals, main) → provide integral estimation exercises for testing.
- `MilsteinConvergence_GBM.py` (GeneratePathsGBMMilstein, mainCalculation) → supports Milstein engine and convergence analysis.
- `StochasticIntegrals.py` (ComputeIntegrals, main) → extends stochastic integral utilities.

#### Lecture 10 – Monte Carlo Simulation of the Heston Model
- `CIR_ExactSimulation.py` (CIR_Sample, GeneratePathsHestonAES, mainCalculation) → supplies exact CIR sampler reused across specs.
- `CIR_paths_boundary.py` (GeneratePathsCIREuler2Schemes, mainCalculation) → informs boundary handling tests.
- `CIR_paths_Exact.py` (CIR_Sample, GeneratePathsCIRExact, mainCalculation) → baseline for exact CIR engine.
- `HestonModelDiscretization.py` (COS method, characteristic function, Monte Carlo, AES, mainCalculation; OptionType) → central for Heston spec and pricing.
- `OptionPrices_EulerAndMilstein.py` (BS pricing, Euler/Milstein paths, payoff estimators, mainCalculation; OptionType) → cross-validates Monte Carlo payoffs.

#### Lecture 11 – Hedging and Monte Carlo Sensitivities
- `BS_Hedging.py` (GeneratePathsGBM, BS_Call_Put_Option_Price, BS_Delta, mainCalculation; OptionType) → informs delta hedging module.
- `HedgingWithJumps.py` (GeneratePathsMerton, GBM, BS pricing, BS_Delta, mainCalculation; OptionType) → extends hedging to jump diffusions.
- `PathwiseSens_DeltaVega.py` (BS pricing, Greeks, Euler paths, EUOptionPriceFromMCPathsGeneralized, PathwiseDelta, PathwiseVega, mainCalculation; OptionType) → basis for pathwise sensitivity estimators.

#### Lecture 12 – Forward Start Options and Model of Bates
- `BatesImpliedVolatility.py` (COS method, implied volatility, Bates characteristic function, mainCalculation; OptionType) → informs Bates spec calibrations.
- `HestonForwardStart2.py` (Forward-start COS, characteristic functions, implied volatility, mainCalculation; OptionType) → basis for forward-start spec.

#### Lecture 13 – Exotic Derivatives
- `AsianOption.py` (PayoffValuation, GeneratePathsGBMEuler, mainCalculation) → guides Asian option Monte Carlo module with variance reduction.
- `DigitalPayoffs_CostReduction.py` (DigitalPayoffValuation, GeneratePathsGBMEuler, UpAndOutBarrier, mainCalculation) → informs digital and barrier pricing schemes.

## Appendix B: Legacy hardhat Architecture Highlights
- Maintain separation of user interface (`simulate_paths()`, plotting helpers) and engine implementations while migrating to the extended spec hierarchy.
- Continue applying purrr-based functional programming (no explicit loops), checkmate validation, and tidy tibbles as documented in the previous architecture write-up.
- Preserve the spec-fit-predict workflow and shared utilities (`generate_standardized_normals()`, `compute_cumulative_paths()`, `matrix_to_tidy()`) as foundational building blocks.# R Package Architecture: tidymodels/hardhat Design Principles

## Overview

The CompFinanceR package implements a modern, modular architecture inspired by tidymodels and hardhat design principles. This document explains the architectural choices and how they benefit both users and developers.

## Design Philosophy

### 1. **Separation of Interface and Implementation (hardhat principle)**

Following hardhat's philosophy, we separate:
- **User Interface**: Simple, consistent functions that users interact with
- **Engine Implementation**: Complex computational logic hidden from users
- **Bridge Layer**: Methods that connect interface to engines

### 2. **Spec-Fit-Predict Pattern (tidymodels principle)**

All stochastic processes follow a consistent three-step workflow:

1. **Specify**: Create a process specification with parameters
   ```r
   spec <- gbm_spec(initial_value = 100, drift = 0.05, volatility = 0.2)
   ```

2. **Fit/Simulate**: Generate sample paths using the specification
   ```r
   paths <- simulate_paths(spec, n_paths = 100, n_steps = 252, maturity = 1.0)
   ```

3. **Predict/Analyze**: Visualize or extract results
   ```r
   plot_paths(paths)
   ```

### 3. **Functional Programming with purrr**

- **Zero for loops or apply functions** throughout the codebase
- Use modern purrr functions:
  - `purrr::accumulate()` - For cumulative operations (replaces for loops)
  - `purrr::map()` - For transformations
  - `purrr::list_rbind()` - For combining data frames (modern replacement for `map_dfr`)
  - Anonymous function syntax `\(x)` (R >= 4.1.0)

### 4. **Robust Input Validation with checkmate**

- Every user-facing function validates inputs with checkmate
- Provides clear, actionable error messages
- Fast C-based validation
- Self-documenting function contracts

## Architecture Layers

### Layer 1: User Interface (Public API)

**Files**: `R/simulate_paths.R`, `R/plot_paths.R`

**Purpose**: Provide simple, consistent interface for users

**Key Functions**:
- `gbm_spec()`, `abm_spec()`, etc. - Create process specifications
- `simulate_paths()` - Generic function for all processes
- `plot_paths()` - Universal plotting function
- `demo_paths()` - Quick demonstrations

**Design Principles**:
- Consistent argument ordering (process_spec first, then simulation params)
- Pipeable with `|>`
- Clear, descriptive names (snake_case)
- Comprehensive roxygen2 documentation
- Input validation at the entry point

**Example**:
```r
# Simple, consistent interface regardless of process complexity
gbm_spec(100, 0.05, 0.2) |>
  simulate_paths(n_paths = 50, n_steps = 252, maturity = 1.0) |>
  plot_paths(n_paths_plot = 10)
```

### Layer 2: Engine Implementation (Internal Logic)

**Files**: `R/simulate_paths_engines.R`

**Purpose**: Handle computational complexity

**Key Components**:
- `simulate_paths.gbm_spec()` - GBM engine
- `simulate_paths.abm_spec()` - ABM engine
- `generate_standardized_normals()` - Shared utility
- `compute_cumulative_paths()` - Functional cumulative computation
- `matrix_to_tidy()` - Format conversion

**Design Principles**:
- S3 method dispatch based on spec class
- Pure functional programming (no side effects except RNG)
- Modular helper functions (DRY principle)
- Hidden from users (not exported)
- Optimized for performance

**Example** (internal):
```r
# Engine method called via dispatch
simulate_paths.gbm_spec <- function(process_spec, n_paths, n_steps, maturity, seed) {
  # Complex implementation hidden from users
  # Uses functional programming throughout
  # Returns standardized tidy tibble
}
```

### Layer 3: Utilities and Helpers

**Purpose**: Reusable components shared across engines

**Key Functions**:
- `generate_standardized_normals()` - Random number generation
- `compute_cumulative_paths()` - Functional path computation
- `matrix_to_tidy()` - Format conversion
- `get_plot_config()` - Visualization configuration

**Design Principles**:
- Pure functions (no side effects)
- Highly reusable
- Well-tested
- Documented internally

## Modern R Practices

### 1. No Deprecated Functions

**Avoided**: `purrr::map_dfr()`, `purrr::map_dfc()` (superseded)

**Used Instead**:
- `purrr::list_rbind()` - Modern row-binding
- `purrr::list_cbind()` - Modern column-binding
- `purrr::map()` + explicit binding from either dplyr::bind_(row/col)

### 2. Anonymous Function Syntax

Using R >= 4.1.0 syntax for clarity:

```r
# Old style
purrr::map(x, function(i) i^2)

# New style (clearer, more concise)
purrr::map(x, \(i) i^2)
```

### 3. Type-Stable Operations

Using purrr's typed variants:
- `map_dbl()` - Returns numeric vector
- `map_int()` - Returns integer vector
- `map_chr()` - Returns character vector
- `map_vec()` - Returns vector of specified type

Ensures outputs are always the expected type, preventing surprises.

### 4. Tidy Evaluation

Using `rlang` for dynamic column names:
```r
tibble::tibble(
  path_id = i,
  time = time_grid,
  !!value_name := value_matrix[i, ]  # Dynamic column naming
)
```

## Consistency Across Processes

All stochastic processes share the same interface:

```r
# GBM
gbm_spec(initial_value, drift, volatility) |> simulate_paths(...)

# ABM
abm_spec(initial_value, drift, volatility) |> simulate_paths(...)

# CIR
cir_spec(initial_value, mean_reversion, long_term_mean, volatility) |> simulate_paths(...)

# Heston
heston_spec(initial_price, initial_variance, ...) |> simulate_paths(...)
```

**Benefits**:
- Easy to learn (learn once, apply everywhere)
- Consistent return format (always tidy tibbles)
- Same plotting function works for all processes
- Easy to swap processes in analysis pipelines

## DRY Principle Implementation

### Avoiding Repetition

1. **Shared validation logic**: checkmate assertions in user-facing functions
2. **Shared random number generation**: `generate_standardized_normals()`
3. **Shared path computation**: `compute_cumulative_paths()`
4. **Shared format conversion**: `matrix_to_tidy()`
5. **Shared plotting**: Single `plot_paths()` for all processes

### Modularity Benefits

- **Easier testing**: Test each component in isolation
- **Easier maintenance**: Fix bugs in one place
- **Easier extension**: Add new processes by implementing one method
- **Better performance**: Optimize shared utilities once

## Extension Pattern

Adding a new process requires:

1. **Create specification function**:
   ```r
   new_process_spec <- function(param1, param2, ...) {
     # Validation
     # Return spec object with appropriate class
   }
   ```

2. **Implement engine method**:
   ```r
   simulate_paths.new_process_spec <- function(process_spec, ...) {
     # Use shared utilities
     # Return tidy tibble
   }
   ```

3. **Add plot configuration**:
   ```r
   # Update get_plot_config() in plot_paths.R
   ```

That's it! The rest works automatically through dispatch and shared utilities.

## Testing Strategy

### Unit Tests

- Test each specification function independently
- Test each engine method independently  
- Test shared utilities independently
- Test error handling and validation

### Integration Tests

- Test full pipelines (spec → simulate → plot)
- Test reproducibility (same seed → same results)
- Test pipeable workflows

### Property-Based Tests

- Test mathematical properties (initial conditions, dimensions, etc.)
- Test tidy data invariants

## Performance Considerations

1. **Vectorization**: Operations are vectorized where possible
2. **Functional programming**: No loop overhead, compiler-friendly
3. **Memory efficiency**: Avoid unnecessary copies
4. **Lazy evaluation**: Compute only what's needed

## Documentation Standards

### roxygen2 Tags

Every exported function includes:
- `@param` - All parameters with types and constraints
- `@return` - Return value structure
- `@examples` - Working examples
- `@details` - Implementation details and equations
- `@export` - Only for public API

### Mathematical Notation

Include SDEs and equations in documentation:
```r
#' dS(t) = r*S(t)*dt + sigma*S(t)*dW(t)
```

### Usage Examples

Show both basic and pipeline usage:
```r
#' @examples
#' # Basic usage
#' spec <- gbm_spec(100, 0.05, 0.2)
#' paths <- simulate_paths(spec, 100, 252, 1.0)
#'
#' # Pipeline usage
#' gbm_spec(100, 0.05, 0.2) |>
#'   simulate_paths(100, 252, 1.0) |>
#'   plot_paths()
```

## Benefits of This Architecture

### For Users

- **Simple, consistent interface** across all processes
- **Pipeline-friendly** design works naturally with tidyverse
- **Clear error messages** from checkmate validation
- **Comprehensive documentation** with examples
- **Type-safe operations** prevent surprises

### For Developers

- **Modular design** makes code easy to understand and modify
- **DRY principle** means fixing bugs in one place
- **Pure functions** are easy to test
- **Clear separation** between interface and implementation
- **Extensible** - adding new processes is straightforward

### For Maintainers

- **Consistent patterns** across codebase
- **Well-tested** at multiple levels
- **Modern R practices** ensure longevity
- **Following community standards** (tidyverse, tidymodels)

## References

- [hardhat: Construct Modeling Packages](https://hardhat.tidymodels.org/)
- [tidymodels principles](https://www.tidymodels.org/learn/)
- [Tidyverse style guide](https://style.tidyverse.org/)
- [purrr documentation](https://purrr.tidyverse.org/)
- [checkmate documentation](https://mllg.github.io/checkmate/)
