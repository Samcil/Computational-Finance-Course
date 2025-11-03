# R Port Implementation Plan

## Mission
- Deliver a unified, tidymodels-aligned R toolbox that mirrors the functionality of the Computational Finance Python notebooks while emphasising composable specifications, reusable engines, and reproducible calibration workflows.
- Focus near-term effort on the interest-rate stack (curve construction, short-rate dynamics, swap/caplet analytics, exposure profiling) to close gaps with the FinancialEngineering_IR_xVA reference material.

## Architectural Guardrails
- Model specifications expose immutable parameter lists and defer heavy lifting to pluggable engines (`term_structure_spec`, `short_rate_spec`, `*_spec` families).
- `fit()`, `simulate_paths()`, `predict()`, and `augment()` follow hardhat blueprints, returning tibbles with explicit column conventions.
- Pricing and calibration functions remain pure (no hidden state), accept specs first, and return tidy metric/value pairs for downstream summarisation.
- Shared math utilities (discount splines, RNG helpers, Fourier kernels) live in dedicated modules to maintain DRY principles.
- Tests, examples, and vignettes accompany new features before features are marked complete.

## Current Implementation Snapshot

### Core Framework
| File | Role | Status | Alignment Notes |
| --- | --- | --- | --- |
| spec_utils.R | Minimal tidymodels-style spec constructors & engine wiring | Stable | Adopted by all spec types; extend to record version metadata before GA |
| simulate_paths.R | User-facing dispatcher for stochastic processes | Stable | Covers all current specs; add input schema checks for list-column outputs |
| simulate_paths_engines.R | Engine implementations for GBM/ABM/Gaussian short rate | Stable | Needs extension hooks for variance reduction and multi-factor engines |
| random_utils.R | Deterministic RNG helpers (`with_random_seed`, standardised normals) | Stable | Consider migrating to withr for nested seeds |
| plot_paths.R | ggplot2 path visualisation | Stable | Add facets for multi-asset specs and swap exposure outputs |
| augment_generics.R | Augment generic registration | Stable | Ensure every fit object registers an augment method before release |
| fit_generics.R | Fit generic registration | Stable | Expand documentation with examples for interest-rate fits |
| price_generics.R | Price generic registration | Stable | Reconcile naming with swap/caplet helpers (currently direct functions) |

### Process & Model Specs
| File | Role | Status | Alignment Notes |
| --- | --- | --- | --- |
| bm_specs.R | ABM/GBM spec constructors | Stable | Engines support Euler schemes; add variance-reduced engine option |
| correlated_bm_process.R | Multivariate Brownian motion spec and simulator | Stable | Provides tidy outputs; align covariance validation with Python Cholesky demos |
| poisson_process.R | Poisson and compound Poisson specs | Stable | Integrate jump size distribution controls per Python lecture |
| cir_process.R | CIR process simulator & analytics | Stable | Update docs to highlight reuse inside Heston/CIR calibration |
| heston_model.R | Heston spec, simulators, CF | Stable | Ensure COS helpers reuse shared Fourier code |
| bates_model.R | Bates (Heston + jumps) spec & simulator | Stable | Add calibration vignette referencing lecture material |
| merton_model.R | Jump-diffusion spec & simulator | Stable | Document linkage with delta hedging jump example |
| shifted_lognormal.R | Shifted lognormal transform utilities | Stable | Review integration into swaption smile workflows |

### Fourier, Numerical & Analytics
| File | Role | Status | Alignment Notes |
| --- | --- | --- | --- |
| cos_method.R | COS pricing & density recovery | Stable | Shares kernels with digital options; extend unit tests to Ho-Lee bond options |
| fft_density_recovery.R | FFT-based density inversion | Stable | Needs example showcasing calibration residuals |
| monte_carlo_integration.R | Generic MC integration helpers | Stable | Reuse in swap exposure Monte Carlo to avoid duplicate accumulation logic |
| stochastic_calculus.R | Ito calculus utilities & exact moments | Stable | Reference in documentation for convergence studies |
| convergence_analysis.R | Euler/Milstein convergence diagnostics | Stable | Align output schema with new plotting conventions |

### Interest-Rate Stack
| File | Role | Status | Alignment Notes |
| --- | --- | --- | --- |
| term_structure_spec.R | Tidymodels-style deterministic curve spec | Stable | Only direct engine implemented; add bootstrapped/Newton engines for multi-curve builds |
| term_structure_engines.R | Direct term-structure calibration engine | Stable | Extend with OIS/LIBOR bootstrapping and Newton-Raphson solver (from Python `MultiCurveBuild.py`) |
| discount_curve_utils.R | Spline-based discount/forward utilities | Stable | Add monotonic enforcement & boundary extrapolation used in Python scripts |
| short_rate_spec.R | Ho-Lee/Hull-White spec & simulator | Stable | Currently Ho-Lee/Hull-White only; plan OU/CIR families per strategic hierarchy |
| simulate_paths_engines.R | Gaussian short-rate Monte Carlo engine | Stable | Needs modular drift hooks for multi-curve theta functions |
| short_rate_derivatives.R | Bond options & Jamshidian swaption pricing | Stable | Introduce analytics for multi-factor specs once available |
| caplet_pricing.R | Caplet/floorlet analytics (bond-option mapping) | Stable | Extend to support displaced diffusion (python Lecture 05 extension) |
| short_rate_calibration.R | Caplet/swaption volatility calibration wrappers | Stable | Currently 1-D volatility search; expand to joint (sigma, a) estimation |
| swap_utils.R | Schedule validation helpers | Stable | Add support for stubs and amortising notionals |
| swap_pricing.R | Swap PV/DV01 & schedule builder | Stable | Introduce floating leg day-count conventions to mirror Python schedule builder |
| swap_exposure.R | Monte Carlo swap exposure & PFE summary | Stable | Integrate cholesky-driven multi-curve simulation when two curve specs exist |

### Pricing, Greeks & Hedging
| File | Role | Status | Alignment Notes |
| --- | --- | --- | --- |
| black_scholes.R | Closed-form pricing & Greeks | Stable | Already used by delta hedging; document integration points |
| asian_options.R | Monte Carlo Asian pricing & control variates | Stable | Add link to pathwise sensitivities for Greeks |
| barrier_options.R | Barrier option pricing & diagnostics | Stable | Extend to include rebate support per Python examples |
| digital_options.R | Digital pricing & COS bridge | Stable | Provide ties to term-structure discounting |
| forward_start_options.R | Forward-start pricing under Heston | Stable | Align with python `HestonForwardStart2.py` parameter naming |
| implied_volatility.R | Newton IV solver & smile plotting | Stable | Reuse in swaption calibration reporting |
| pathwise_sensitivities.R | Pathwise Greeks utilities | Stable | Integrate with swap exposure outputs |
| delta_hedging_bs.R | Discrete-time BS hedging backtest | Stable | Factor out portfolio accounting for reuse |
| delta_hedging_jumps.R | Hedging under jumps | Stable | Align random seeds with `random_utils` |

### Package Infrastructure
| File | Role | Status | Alignment Notes |
| --- | --- | --- | --- |
| price_generics.R | Generic registration | Stable | Expand coverage to interest-rate instruments |
| fit_generics.R | Generic registration | Stable | Already used by term structure & short-rate fits |
| augment_generics.R | Generic registration | Stable | Ensure consistent messaging |
| zzz.R | Package hooks & namespace init | Stable | Monitor dependency loading |

## Reference Python Baseline
- `FinancialEngineering_IR_xVA/MultiCurveBuild.py`: Newton-based multi-curve calibration with swap pricing residuals.
- `.../affine_diffusion_materials/QUICK_DEMO.R` & `test_fourier_features.R`: Bond/swap valuation pipelines used for parity checks.
- Lectures 05–12 notebooks: Caplet/swaption analytics, Jamshidian implementations, xVA exposure demos.
- Cholesky demos (e.g., `cholesky_decomposition_demo.R`) inform correlated factor handling.

## Gap Analysis
- Term-structure module lacks piecewise bootstrapping, forward extrapolation, and multi-curve linkage (OIS vs LIBOR) showcased in Python.
- Short-rate calibration searches only volatility; Python workflows solve simultaneously for `(sigma, a)` and optionally initial forward curve shifts.
- No abstraction yet for process inheritance (e.g., OU family) or multiple-factor Gaussian models (G2++), both present in lecture notes.
- Swap exposure simulator assumes single-curve discounting; Python examples support discount/forward curve separation and CSA adjustments.
- Documentation and vignettes have not been refreshed post engine refactor; plan requires a new interest-rate vignette aligned with the strategic port plan.

## Architecture Blueprint
- **Term-Structure Engines**: Introduce `bootstrap_piecewise`, `newton_multi_curve`, and `spline_smooth` engines, each conforming to a shared contract returning discount/zero/forward columns with provenance metadata.
- **Short-Rate Model Hierarchy**: Establish `gaussian_short_rate_spec` (base), `ho_lee_spec`, `hull_white_spec`, and future `g2pp_spec`, using inheritance helpers to register model families and share theta calculations.
- **Calibration Orchestration**: Expand `fit.short_rate_spec()` to accept recipe-style lists of instrument books (caps, swaptions, bonds) and run constrained optimisations via `stats::optim()` or `nloptr`, mirroring Python solver flows.
- **Pricing & Exposure Layer**: Normalise outputs (`metric`, `value`, `unit`) across swaps, caps, swaptions, and exposures; ensure `price_*` functions optionally attach Greeks and scenario metadata when fed path-level inputs.
- **Testing & Diagnostics**: Adopt snapshot-based regression tests for calibration outputs, comparing against Python benchmarks stored under `tests/fixtures/`.

## Implementation Roadmap
- **Phase A – Interest-Rate Foundations (in flight)**
  - Add comprehensive tests for `term_structure_spec` and `short_rate_spec` covering Ho-Lee/Hull-White analytics (`tests/testthat/test-term-structure.R`, `test-short-rate.R`).
  - Build bootstrapping and Newton engines in `term_structure_engines.R`, including dual-curve calibration mirroring `MultiCurveBuild.py`.
  - Refactor `short_rate_calibration.R` to support joint parameter estimation and return broom-like tidiers (`augment()`, `tidy()`).
- **Phase B – Multi-Curve Simulation & Pricing**
  - Extend `short_rate_spec` engine state to store multiple discount/forward curves and adapt Monte Carlo to pull curve-specific theta functions.
  - Update `swap_pricing.R` and `swap_exposure.R` to read forward curves separately from discount curves, enabling FRA-style analytics.
  - Implement curve-consistent caplet/swaption pricing tests referencing Python outputs.
- **Phase C – Advanced Models & Risk**
  - Introduce `g2pp_spec` and correlated Gaussian factor engines, using Cholesky utilities from the Python lectures.
  - Add displaced-diffusion caplet pricing and smile fitting modules to align with Lecture 05 content.
  - Expand exposure analytics to compute CVA/DVA style metrics leveraging Monte Carlo results.
- **Phase D – Documentation & Packaging**
  - Publish an interest-rate vignette demonstrating calibration → simulation → pricing flow.
  - Refresh README and architecture docs to reflect engine-based design.
  - Run `devtools::document()`, `devtools::test()`, `devtools::check()` gates and capture results in `REFACTORING_STATUS.md`.

## Testing & Validation Strategy
- Extend `tests/testthat` with parity tests that ingest saved Python benchmark CSVs for curves, caplets, and swaptions.
- Add high-level integration tests covering multi-curve calibration, swap PV consistency, and exposure quantile stability.
- Use property-based tests (via `quickcheck` or manual loops) for monotonic discount curves and Jamshidian root bracketing.

## Documentation & Examples
- Update roxygen examples to showcase spec creation → fit → predict workflows for both term structures and short-rate models.
- Provide notebooks or vignettes comparing R outputs directly against Python results, highlighting any tolerances.
- Maintain `REFACTORING_SUMMARY.md` with progress bullets per roadmap phase.

## Dependencies & Tooling
- Core imports remain `tibble`, `dplyr`, `purrr`, `tidyr`, `rlang`, `cli`, `checkmate`, `hardhat`.
- Evaluate `Matrix`/`pracma` for multi-curve solvers, `numDeriv` for calibration Jacobians, and `future`/`furrr` for parallel Monte Carlo once deterministic baselines are stabilised.

## Success Criteria
- Multi-curve calibration reproduces Python benchmarks within tolerance (<0.5 bp on zero curve, <1 bp on swap PVs).
- Short-rate calibration supports joint `(sigma, a)` estimation with convergence diagnostics and broom tidiers.
- Swap exposure engine generates EE/PFE curves consistent across single-curve and multi-curve modes.
- Package passes `R CMD check` without errors/warnings/notes and all new modules have targeted tests and documentation.
