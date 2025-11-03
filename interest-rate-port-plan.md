---
title: "Interest Rate Models R Port Plan"
author: "GitHub Copilot"
date: "2025-02-18"
output: html_document
---

## Background
- Translate Lech Grzelak's Python interest-rate notebooks into an R package that already exposes `spec` and `simulate_paths()` dispatchers for diffusion models.
- Preserve actuarial focus: short-rate dynamics, term-structure construction, and derivatives on rates (caps, floors, swaps, swaptions).
- Adopt tidymodels and hardhat idioms so short-rate processes behave like model specifications with `fit()`, `augment()`, and `predict()` workflows.

## Architectural Alignment
- **Model specs**: Introduce `rate_curve_spec()` and `short_rate_spec()` constructors that call `hardhat::new_model_spec()` with `set_engine()` for analytic or Monte Carlo paths. Store validated parameter tibbles and metadata rather than base lists.
- **Forge step**: Use `hardhat::mold()` to accept calibration quotes (`curve_points`, `cap_quotes`, `swaption_vols`) and ensure tidy column names before fitting.
- **Engines**: Separate engines for deterministic curves (e.g., Newton solver vs spline) and stochastic simulation (Euler, exact, analytic). Follow existing `simulate_paths_engines.R` naming, but add rate-specific helpers.
- **Outputs**: Return tibbles with explicit identifiers (`path_id`, `tenor`, `discount_factor`, `rate_path`) to facilitate ggplot2 plots and broom-style augmentations.
- **Error handling**: Continue to use `rlang::abort()` with tidy error classes, validating via `checkmate` where numeric bounds are critical (e.g., positive volatilities, mean-reversion > 0).

## Consistent API Surface
- **Spec constructors**: Keep all specs in dedicated files (`short_rate_spec.R`, upcoming `hybrid_specs.R`, `mortgage_specs.R`) inheriting from `process_spec` so that `simulate_paths()` dispatch remains uniform. Register print methods with cli messaging for discoverability.
- **Simulation dispatch**: Route every stochastic routine through `simulate_paths()`; refactor `simulate_paths.short_rate_spec()` to delegate Gaussian innovations to a shared helper in `simulate_paths_engines.R` (reusing `generate_standardized_normals()` and `compute_cumulative_paths()`). Prevent ad-hoc RNG loops when VaR or exposure modules need rate paths.
- **Pricing generics**: Extend the existing `price_*` generic family—`price_zcb()` already lives in `price_generics.R`; add `price_swaption()`/`price_caplet()` generics rather than bespoke functions per module so developers plug new engines without rewriting callers. Continue returning `tibble(metric, value)` payloads for broom compatibility.
- **Calibration contract**: Have every calibration use `fit()` methods that wrap helper optimisers. Reuse `short_rate_calibration.R` patterns (blueprint storage, `short_rate_fit` class) for forthcoming displaced-diffusion and hybrid calibrations so `augment()` and `predict()` work identically across models.
- **Predict/augment**: Guarantee that `predict(<*_fit>)` exposes `type = c("discount", "zero", "forward", "exposure")` where appropriate, leaning on `augment()` for quote diagnostics. Any new model must implement both to pass API compliance checks.

## Shared Utility Layer & DRY Guardrails
- Centralise schedule validation logic (`start`, `end`, `pay_time`, `accrual_fraction`, `notional`) into `validate_swap_schedule()` inside `swap_pricing.R` (or new `swap_utils.R`) and reuse it in `price_swap()`, `price_swaption()`, `simulate_swap_exposure()`, and future mortgage/amortizing helpers.
- Move Jamshidian-specific helpers (`bond_option_price_short_rate()`, `find_short_rate_root()`) into an internal `short_rate_utilities.R` so caplet/floorlet, swaption, and VaR modules call the same math primitives.
- Promote discount/forward spline builders (`make_log_discount_spline()`, `make_forward_function()`) to an internal module shared by both `short_rate_spec()` and upcoming multi-curve specs to avoid recomputing splines in multiple files.
- Add a `with_random_seed()` wrapper (thin proxy over `withr::with_seed`) to `simulate_paths_engines.R` and use it across all engines for consistent RNG handling and easy audit logging.
- Replace repeated `%>%` style column selections with dedicated tibble formatters (e.g., `as_discount_curve_tbl()`, `as_quote_tbl()`) so future instruments inherit the same cleaning rules.
## Module Inventory And Port Targets

### Term-Structure Foundations
- **Python sources**: `MultiCurveBuild.py`, `YieldCurveBuildGreeks.py`, `YieldCurveBuild_Treasury.py`.
- **R deliverables**:
  - `R/term_structure_spec.R`: `rate_curve_spec(curve_type = c("ois", "libor"), ...)`, `set_engine("nelder_mead")`, `set_engine("spline")`.
  - `fit.term_structure_spec()`: calibrates bootstrapped discount factors using `hardhat::mold()` on input quotes.
  - `augment.term_structure_fit()`: returns tidy discount curves; includes Greek sensitivities by automatic perturbation (align with `YieldCurveBuildGreeks.py`).
  - `predict.term_structure_fit(type = c("discount", "zero", "forward"))`.
  - Utilities: `solve_discount_factor_newton()`, `piecewise_log_fwd()` in `R/term_structure_engines.R`.
  - Tests: `tests/testthat/test-term-structure-spec.R` covering bootstrapping convergence, Jacobian sign, and sensitivity deltas.
  - Vignette section: "Term Structure Calibration" using sample US Treasury data (Quarto chunk).

### Short-Rate Dynamics
- **Python sources**: `CIR_IR_paths.py`, `Ho-Lee-ZCBs.py`, `Hull-White-Paths.py`, `Hull-White-ZCBs.py`, `Hull-White-CompRateSim.py`, `Hull-White-ZCBs2.py`, `Hull_White_1F_2F_Comparison.py`.
- **R deliverables**:
  - `R/short_rate_specs.R`: `short_rate_spec(model = "ho_lee", mean_reversion, vol, shift = 0)`; variants for `hull_white_1f`, `hull_white_2f`, and `cir` with parameter validation.
  - `set_engine()` combinations: `"analytic_zcb"`, `"monte_carlo"`, `"fourier"` (for future extensions).
  - `simulate_paths.short_rate_spec()`: dispatches to `simulate_short_rate_paths()` using Brownian increments already implemented in `simulate_paths_engines.R`, but extends to multifactor correlation via `cov_cholesky()` drawn from Python examples.
  - `price_zcb(short_rate_fit, maturities)`: analytic formulas for Ho-Lee and Hull-White (closed form) plus Monte Carlo fallback (averaging `exp(-\int r dt)`).
  - `compare_one_two_factor()` replicating `Hull_White_1F_2F_Comparison.py`, returning tibble with columns `tenor`, `price_1f`, `price_2f`.
  - Tests: scenario snapshots verifying that analytic and simulated ZCB prices agree within tolerance; regression tests using saved seeds.
  - Roxygen examples showing integration with `term_structure_spec()` for initial curve.

### Caplets, Floorlets, And Swaptions
- **Python sources**: `HW_Caplets.py`, `HW_CapletsAndFloorlets.py`, `HW_OptionsOnZCBs.py`, `ShiftedLognormal.py`, `JamshidianTrick.py`.
- **R deliverables**:
  - `R/caplet_pricing.R`: expand existing `price_caplet()`/`price_floorlet()` into formal generics that reuse `bond_option_price_short_rate()` from `short_rate_derivatives.R`. Add Monte Carlo engines that call the shared short-rate path helper instead of duplicating Euler code.
  - `R/swaption_pricing.R`: keep Jamshidian implementation but wire it through a `price_swaption.short_rate_spec()` method so VaR, calibration, and exposure modules call the same entry point.
  - `R/shifted_lognormal.R`: `shifted_lognormal_call()` and `implied_vol_shifted()` replicating Python implied volatility solver via `stats::uniroot` with fallback to `bbmle` for stability; ensure calibration helpers call `fit(<short_rate_spec>, method = "displaced_diffusion")` once introduced.
  - Introduce reusable coupon schedule helpers (`swap_cashflow_schedule()` already in place) and export a miniature `schedule_spec` for amortizing legs so mortgage functions piggyback on identical table structures.
  - Tests: calibrate to deterministic cases (flat term structure) to match closed-form results; property tests for monotonicity in strike shifts and parity between analytic vs MC paths using shared engines.

### Swap Valuation And Exposures
- **Python source**: `Swaps_HW.py`.
- **R deliverables**:
  - `R/swap_pricing.R`: `price_swap(short_rate_fit, swap_schedule, payoff = c("payer", "receiver"))` returning tibble with PV, DV01, forward rates.
  - `simulate_swap_exposure()` to compute expected exposure profiles under Monte Carlo, returning tidy data for CVA integration (aligns with actuarial risk metrics).
  - Provide `plot_exposure_profile()` using ggplot2.
  - Tests verifying PV parity vs direct discounting and regression tests on exposure percentiles.

### Mortgage And Amortizing Instruments
- **Python sources**: `AnnuityMortgage.py`, `BulletMortgage.py`, `Incentives.py`, `StochasticAmortizingSwap.py`.
- **R deliverables**:
  - `R/mortgage_cashflows.R`: `generate_mortgage_schedule()` returning tibble with outstanding balance, principal, and interest by period; internally call `validate_swap_schedule()` for consistency with swap tooling.
  - `R/prepayment_models.R`: `simulate_prepayment_paths()` implementing incentive-based prepayment dynamics tied to Hull-White rates via the shared short-rate path helper; provide deterministic fallbacks for validation.
  - `R/amortizing_swap.R`: `price_amortizing_swap()` leveraging mortgage schedules and `price_swap()` for floating-leg valuation, ensuring DRY exposure of discount factors.
  - Tests covering balance reconciliation, deterministic amortization parity, and sensitivity of prepayment incentives to rate shocks using the unified schedule utilities.

### Hybrid Equity-Rate Models
- **Python sources**: `BSHW_COS_vs_MC.py`, `BSHW_Euler_Milstein.py`, `BSHW_AE.py`, `SZHW_COS_vs_MC.py`, `H1-HW_COS_vs_MC.py`, `H1-HW_COS_vs_MC_FX.py`.
- **R deliverables**:
  - `R/hybrid_specs.R`: `hybrid_spec(model = c("bshw", "h1_hw", "szhw"), ...)` extending existing `spec` framework with joint equity-rate characteristic functions. Store correlation matrices once and reuse the multi-factor Cholesky helper that will also serve FX models.
  - Engines for `"cos"` and `"monte_carlo"`, reusing `cos_method.R` and the Gaussian short-rate path helper to avoid separate Euler loops.
  - `price_hybrid_option()` supporting equity and FX payoffs with rate coupling via a shared `price_hybrid.short_rate_spec()` method that funnels through pricing generics.
  - Tests confirming convergence between Euler/AES schemes, COS accuracy, and correlation sanity checks against Python benchmarks using identical seeds for both engines.

### FX And Cross-Currency Extensions
- **Python source**: `H1-HW_COS_vs_MC_FX.py`.
- **R deliverables**:
  - `R/fx_rate_models.R`: `fx_short_rate_spec()` capturing domestic/foreign short rate dynamics with shared volatility shocks; reuse the same discount/forward spline builders promoted to the utility layer.
  - `price_fx_option()` integrating hybrid characteristic functions to deliver domestic-discounted FX option prices through the generic pricing interface.
  - Calibration helper to align domestic vs foreign curves and enforce no-arbitrage drift adjustments by calling the existing `fit()` infrastructure with multi-curve blueprints.
  - Tests ensuring put-call parity across domestic and foreign discount curves plus cross-checks versus equity-rate hybrid outputs to prevent duplicated COS code.

### Convexity And Smile Adjustments
- **Python sources**: `ConvexityCorrection.py`, `DD_ImpliedVolatility.py`.
- **R deliverables**:
  - `R/convexity_adjustments.R`: `convexity_correction_forward_rate()` returning tibble with corrected forwards for different tenors under Hull-White dynamics, reusing `b_factor()`/`theta_integral()` from the shared utilities.
  - `R/displaced_diffusion.R`: `displaced_diffusion_iv()` and `calibrate_displacement()` providing volatility smile fitting to cap/floor quotes; integrate with the generic `fit()` pathway so displaced diffusion updates produce `short_rate_fit` objects.
  - Tests benchmarking corrections against analytic expectations and verifying displaced diffusion calibration monotonicity while avoiding redundant numerical integration code.

### Netting Set Exposure Analytics
- **Python source**: `Exposures_HW_Netting.py`.
- **R deliverables**:
  - `R/netting_exposure.R`: `simulate_netting_exposure()` generating pathwise exposures for portfolios of swaps/swaptions by composing `simulate_paths()` and `price_swap()`; ensure shared helper functions compute positive exposure and discounting so swap exposure and VaR codebases stay in sync.
  - `summarise_exposure_distribution()` returning tidy percentiles for CVA/DVA workflows, using the same quantile labelling helper defined once in the utility layer.
  - ggplot helper `plot_netting_exposure()` for exposure term structure visualization leveraging the existing plotting theme from `plot_paths.R`.
  - Tests validating positive exposure aggregation, netting benefits, and regression baselines aligned with Python reference outputs using the shared schedule validator.

### Value-At-Risk And Expected Shortfall
- **Python source**: `MonteCarloVaR.py`.
- **R deliverables**:
  - `R/var_es.R`: `simulate_var(short_rate_fit, portfolio_definition, horizon, alpha)` computing VaR/ES via Hull-White Monte Carlo with tidy outputs (`statistic`, `value`). Reuse `simulate_swap_exposure()` under the hood for portfolio revaluation to remain DRY.
  - Portfolio helpers leveraging existing swap pricing to accommodate heterogeneous notional/coupon mixes while sharing schedule validation and pricing utilities with the swaps module.
  - Diagnostic plotting functions comparing analytic ZCB vs simulated proxies and histogram overlays of P&L, mirroring the Python script, and built on top of `plot_paths.R` utilities.
  - Tests covering quantile stability across seeds and ES consistency checks (ES ≤ VaR for loss-defined conventions) in addition to regression tests checking parity between VaR computed from exposures vs direct revaluation.

### Multi-Curve And Hybrid Extensions
- **Python sources**: `affine_diffusion_materials/QUICK_DEMO.R`, `test_fourier_features.R` (for future Fourier enhancements).
- **R deliverables** (Stretch goals): plan to expose multi-curve bootstrapping with OIS-LIBOR spreads, but schedule for later milestone once base plan complete.

## Implementation Roadmap
- **Milestone A (Week 1)**: Ship `term_structure_spec` and `short_rate_spec` with Ho-Lee analytic ZCBs. Include tests and minimal vignette example.
- **Milestone B (Week 2)**: Add Hull-White 1F Monte Carlo engine, path simulation tests, and analytic vs MC parity checks.
- **Milestone C (Week 3)**: Deliver caplet/floorlet pricing plus shifted lognormal implied volatility tools. Provide calibration vignette chunk.
- **Milestone D (Week 4)**: Implement Jamshidian swaption pricing and swap exposure simulations. Expand documentation with actuarial CVA context.
- **Milestone E (Week 5)**: Polish multi-curve utilities, finalize vignettes, and run package diagnostics (`devtools::check()`).
- **Milestone F (Week 6)**: Deliver mortgage cashflow engines, amortizing swap pricing, and netting exposure analytics with visual toolset.
- **Milestone G (Week 7)**: Ship hybrid equity-rate and FX specifications with COS engines plus convexity/displaced diffusion adjustments.
- **Milestone H (Week 8)**: Implement VaR/ES Monte Carlo workflows, integrate outputs into risk vignettes, and broaden regression tests for exposure distributions.
- **Milestone 0 (Pre-Week 1)**: Consolidate shared utilities (`swap` schedule validation, Gaussian path helper, discount spline builders) and wire `simulate_paths.short_rate_spec()` plus calibration methods to use the unified layer before new features land.

## Testing Strategy
- Use `testthat` contextual blocks keyed to each spec (`test_that("Ho-Lee analytic matches MC", {...})`).
- Store lightweight reference curves in `tests/fixtures/` for deterministic tests; avoid writing plots.
- Leverage `vdiffr` snapshots for ggplot outputs (optional) and `waldo::compare()` for tibble comparisons.
- Integrate property-based tests with `hypothesisr` or `quickcheck` to ensure discount factors remain in (0,1] across parameter draws.

## Documentation And Vignettes
- Extend `vignettes/monte-carlo-diagnostics.qmd` with section on rate models.
- Create new vignette `vignettes/short-rate-models.qmd` summarizing Ho-Lee and Hull-White specs with Quarto code chunks.
- Draft calibration walkthrough referencing actuarial reporting requirements (e.g., Solvency II yield curves).
- Update `R_PORT_PLAN.md` cross-reference table to include interest-rate milestones once initial functions merge.

## Dependencies And Tooling
- Add imports: `hardhat`, `parsnip` (for consistent `set_engine()`), `pracma` (matrix exponentials), `Matrix` (sparse handling for two-factor covariance), `stats` (PCHIP interpolation fallback), `rootSolve` (Newton solvers).
- Ensure random seeds use `withr::with_seed()` inside tests for determinism.
- Evaluate `fansi` or `cli` for progress messaging consistent with existing `simulate_paths()` output.

## Risks And Mitigations
- **Analytic vs Monte Carlo drift**: Mitigate by unit tests on known closed-form values and by exposing integration step size as argument with default per Python prototypes.
- **Performance**: R loops may be slower than NumPy; plan to vectorize integrals with matrix algebra and consider `Rcpp` fallback if profiling justifies.
- **Regulatory sensitivity**: Provide audit trails by logging calibration residuals and saving them in returned objects (tidy tibble with columns `instrument`, `quote`, `model_price`).

## Next Actions
1. Extract shared helpers (`validate_swap_schedule()`, Gaussian path generator, discount spline utilities) and refactor existing short-rate modules to call them.
2. Prototype Ho-Lee analytic ZCB pricing in R and validate against Python script using shared seed data, ensuring the new helper layer is exercised end-to-end.
3. Draft vignette outline to maintain alignment with actuarial reporting use cases, highlighting the spec–fit–predict pattern for interest-rate models.

## Progress Update (2025-11-02)
- `term_structure_spec()` and `short_rate_spec()` now conform to the hardhat
  spec pattern, with purrr-based simulations, analytic Ho-Lee and Hull-White
  discount factors, and deterministic calibration helpers.
- Companion tests in `tests/testthat/test-term-structure-spec.R` and
  `tests/testthat/test-short-rate-spec.R` validate curve preparation, analytic
  versus Monte Carlo parity, and Hull-White mean-reversion constraints.
- Package metadata includes `withr` for deterministic RNG handling; raw
  vector-valued closures are preserved for efficiency in analytic routines.
- Added `swap_cashflow_schedule()` and `price_swap()` with regression tests to
  deliver payer/receiver PV and DV01 analytics against calibrated curves.
- New derivative layer delivers `price_caplet()`, `price_floorlet()`, and
  `price_swaption()` via Jamshidian decomposition plus shifted lognormal helpers
  with standalone tests for deterministic limits and implied volatility checks.
- Implemented `calibrate_short_rate_volatility()` with deterministic caplet
  fixtures, regression coverage, and a calibration print helper to keep the API
  aligned with actuarial usage.
- Added `calibrate_short_rate_swaption_volatility()` using Jamshidian pricing,
  lifted weighted regression tests, and expanded the calibration suite to cover
  swaption quotes alongside caplets.
- Introduced `fit(short_rate_spec)` with hardhat blueprints so short-rate
  calibrations return tidy `short_rate_fit` objects compatible with
  `augment()`, enabling tabular comparisons across multiple models.

## Upcoming Focus
1. Draft vignette sections covering caplet/floorlet analytics, Jamshidian
  swaption calibration, and shifted lognormal implied-volatility tooling.
2. Integrate mortgage, exposure, and VaR outputs into the actuarial CVA
  vignette, highlighting deterministic regression fixtures fed by the shared
  helper layer.
3. Build comparison workflows where multiple short-rate fits live in a single
  tibble row-per-model table, facilitating residual diagnostics, upcoming
  multi-curve extensions, and hybrid equity-rate calibration studies while
  exercising the consolidated generics (`simulate_paths()`, `price_*()`, `fit()`).
