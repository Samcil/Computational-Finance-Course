test_that("H1-HW characteristic function properties hold", {
  curve <- tibble::tibble(
    tenor = seq(0, 10),
    discount_factor = exp(-0.02 * tenor)
  )

  short_rate <- short_rate_spec(
    model = "hull_white",
    volatility = 0.01,
    mean_reversion = 0.2,
    curve = curve
  )

  spec <- h1_hw_spec(
    spot = 100,
    short_rate = short_rate,
    initial_variance = 0.04,
    variance_mean_reversion = 1.5,
    variance_long_term = 0.04,
    variance_volatility = 0.3,
    correlation_eq_var = -0.5,
    correlation_eq_rate = 0.25
  )

  state <- spec$method$engine_state
  u <- seq(-2, 2, length.out = 9)

  cf_values <- h1_hw_characteristic_function(
    u = u,
    maturity = 1,
    theta_fun = state$theta_fun,
    initial_rate = state$initial_rate,
    mean_reversion = state$mean_reversion,
    short_rate_vol = state$short_rate_vol,
    initial_variance = state$heston$initial_variance,
    variance_mean_reversion = state$heston$mean_reversion,
    variance_long_term = state$heston$long_term,
    variance_volatility = state$heston$volatility,
    correlation_eq_var = state$correlation_eq_var,
    correlation_eq_rate = state$correlation_eq_rate
  )

  expect_equal(cf_values[u == 0], 1 + 0i, tolerance = 1e-10)
  expect_equal(cf_values[length(u)], Conj(cf_values[1]), tolerance = 1e-10)
  expect_true(all(Mod(cf_values) <= 1 + 1e-8))
})


test_that("H1-HW COS pricing produces sensible call/put ordering", {
  curve <- tibble::tibble(
    tenor = seq(0, 10),
    discount_factor = exp(-0.025 * tenor)
  )

  short_rate <- short_rate_spec(
    model = "hull_white",
    volatility = 0.015,
    mean_reversion = 0.3,
    curve = curve
  )

  spec <- h1_hw_spec(
    spot = 95,
    short_rate = short_rate,
    initial_variance = 0.05,
    variance_mean_reversion = 2.0,
    variance_long_term = 0.04,
    variance_volatility = 0.35,
    correlation_eq_var = -0.4,
    correlation_eq_rate = 0.15
  )

  strikes <- c(80, 95, 110)
  prices <- price_h1_hw_option_cos(
    process_spec = spec,
    strikes = strikes,
    maturity = 1.5,
    n_terms = 256,
    truncation = 10,
    integration_points = 1500
  )

  expect_true(all(diff(prices$call_price) <= 0))
  expect_true(all(diff(prices$put_price) >= 0))
  parity <- prices$call_price - prices$put_price
  discount_factor <- spec$method$engine_state$discount_fun(1.5)
  expect_equal(parity, spec$spot - discount_factor * strikes, tolerance = 1e-4)
})


test_that("Schoebel-Zhu Hull-White characteristic function is well-behaved", {
  curve <- tibble::tibble(
    tenor = seq(0, 10),
    discount_factor = exp(-0.015 * tenor)
  )

  short_rate <- short_rate_spec(
    model = "hull_white",
    volatility = 0.012,
    mean_reversion = 0.25,
    curve = curve
  )

  spec <- szhw_spec(
    spot = 120,
    short_rate = short_rate,
    initial_volatility = 0.2,
    vol_mean_reversion = 2.2,
    vol_long_term = 0.18,
    vol_volatility = 0.4,
    correlation_eq_vol = -0.3,
    correlation_rate_vol = -0.2,
    correlation_eq_rate = 0.1
  )

  cf_values <- szhw_characteristic_function(
    u = c(-1.5, 0, 1.5),
    maturity = 2.0,
    discount_fun = spec$method$engine_state$discount_fun,
    initial_rate = spec$method$engine_state$initial_rate,
    mean_reversion = spec$method$engine_state$mean_reversion,
    short_rate_vol = spec$method$engine_state$short_rate_vol,
    initial_volatility = spec$method$engine_state$ou$initial,
    vol_mean_reversion = spec$method$engine_state$ou$mean_reversion,
    vol_long_term = spec$method$engine_state$ou$long_term,
    vol_volatility = spec$method$engine_state$ou$volatility,
    correlation_eq_vol = spec$method$engine_state$correlation_eq_vol,
    correlation_rate_vol = spec$method$engine_state$correlation_rate_vol,
    correlation_eq_rate = spec$method$engine_state$correlation_eq_rate,
    integration_points = 1200
  )

  expect_equal(cf_values[2], 1 + 0i, tolerance = 1e-10)
  expect_equal(cf_values[1], Conj(cf_values[3]), tolerance = 1e-10)

  strikes <- c(100, 120, 140)
  prices <- price_szhw_option_cos(
    process_spec = spec,
    strikes = strikes,
    maturity = 2.0,
    n_terms = 256,
    truncation = 10,
    integration_points = 1200
  )

  expect_true(all(prices$call_price > 0))
  expect_true(all(prices$put_price > 0))
})


test_that("FX H1-HW pricing respects arbitrage bounds", {
  curve_d <- tibble::tibble(
    tenor = seq(0, 10),
    discount_factor = exp(-0.02 * tenor)
  )
  curve_f <- tibble::tibble(
    tenor = seq(0, 10),
    discount_factor = exp(-0.01 * tenor)
  )

  domestic <- short_rate_spec(
    model = "hull_white",
    volatility = 0.01,
    mean_reversion = 0.25,
    curve = curve_d
  )

  foreign <- short_rate_spec(
    model = "hull_white",
    volatility = 0.008,
    mean_reversion = 0.2,
    curve = curve_f
  )

  spec <- fx_h1_hw_spec(
    spot_fx = 1.2,
    domestic_short_rate = domestic,
    foreign_short_rate = foreign,
    initial_variance = 0.04,
    variance_mean_reversion = 1.4,
    variance_long_term = 0.05,
    variance_volatility = 0.35,
    correlation_eq_var = -0.45,
    correlation_eq_domestic = 0.25,
    correlation_eq_foreign = -0.1,
    correlation_var_domestic = -0.3,
    correlation_var_foreign = 0.2,
    correlation_domestic_foreign = 0.5
  )

  strikes <- c(1.0, 1.2, 1.4)
  maturity <- 1.25
  prices <- price_fx_h1_hw_option_cos(
    process_spec = spec,
    strikes = strikes,
    maturity = maturity,
    n_terms = 256,
    truncation = 8,
    integration_points = 1200
  )

  discount_d <- spec$method$engine_state$domestic$discount_fun(maturity)
  discount_f <- spec$method$engine_state$foreign$discount_fun(maturity)
  forward <- spec$spot_fx * discount_f / discount_d

  expect_true(all(prices$call_price >= 0))
  expect_true(all(prices$call_price <= discount_d * forward + 1e-6))
  expect_true(all(diff(prices$call_price) <= 1e-8))
})


test_that("SZHW Monte Carlo simulation is reproducible and numeraire-consistent", {
  curve <- tibble::tibble(
    tenor = seq(0, 5),
    discount_factor = exp(-0.02 * tenor)
  )

  short_rate <- short_rate_spec(
    model = "hull_white",
    volatility = 0.01,
    mean_reversion = 0.3,
    curve = curve
  )

  spec <- szhw_spec(
    spot = 90,
    short_rate = short_rate,
    initial_volatility = 0.25,
    vol_mean_reversion = 1.5,
    vol_long_term = 0.2,
    vol_volatility = 0.35,
    correlation_eq_vol = -0.25,
    correlation_rate_vol = 0.3,
    correlation_eq_rate = -0.1
  )

  n_paths <- 7
  n_steps <- 24
  maturity <- 2

  paths_a <- simulate_paths(
    process_spec = spec,
    n_paths = n_paths,
    n_steps = n_steps,
    maturity = maturity,
    seed = 321
  )

  paths_b <- simulate_paths(
    process_spec = spec,
    n_paths = n_paths,
    n_steps = n_steps,
    maturity = maturity,
    seed = 321
  )

  paths_c <- simulate_paths(
    process_spec = spec,
    n_paths = n_paths,
    n_steps = n_steps,
    maturity = maturity,
    seed = 123
  )

  expect_equal(paths_a, paths_b)
  expect_false(isTRUE(all.equal(paths_a$spot, paths_c$spot)))
  expect_equal(nrow(paths_a), n_paths * (n_steps + 1))
  expect_true(all(paths_a$discount_factor > 0))
  expect_true(all(paths_a$money_market > 0))
  expect_true(max(abs(paths_a$money_market * paths_a$discount_factor - 1)) < 1e-8)
})


test_that("SZHW diversification value matches manual payoff replication", {
  curve <- tibble::tibble(
    tenor = seq(0, 6),
    discount_factor = exp(-0.03 * tenor)
  )

  short_rate <- short_rate_spec(
    model = "hull_white",
    volatility = 0.012,
    mean_reversion = 0.25,
    curve = curve
  )

  spec <- szhw_spec(
    spot = 100,
    short_rate = short_rate,
    initial_volatility = 0.22,
    vol_mean_reversion = 1.1,
    vol_long_term = 0.18,
    vol_volatility = 0.28,
    correlation_eq_vol = -0.35,
    correlation_rate_vol = 0.2,
    correlation_eq_rate = -0.15
  )

  omega <- c(0, 0.25, 0.6, 1)
  maturity <- 1.5
  settlement <- 2
  seed <- 987
  n_paths <- 200
  n_steps <- 36

  function_values <- szhw_diversification_value(
    process_spec = spec,
    omega = omega,
    maturity = maturity,
    settlement = settlement,
    n_paths = n_paths,
    n_steps = n_steps,
    seed = seed
  )

  mc_paths <- simulate_paths(
    process_spec = spec,
    n_paths = n_paths,
    n_steps = n_steps,
    maturity = maturity,
    seed = seed
  )

  terminal_slice <- mc_paths |>
    dplyr::filter(abs(time - maturity) < 1e-10)

  sr_state <- short_rate_state(spec$args$short_rate)
  discount_settlement <- sr_state$discount_fun(settlement)

  zero_coupon_T_T1 <- price_zcb(
    object = spec$args$short_rate,
    maturities = rep(settlement, nrow(terminal_slice)),
    valuation_time = maturity,
    short_rate = terminal_slice$short_rate
  )

  bond_leg <- zero_coupon_T_T1 / discount_settlement
  equity_leg <- terminal_slice$spot / spec$spot
  manual <- vapply(
    omega,
    function(w) {
      payoff <- w * equity_leg + (1 - w) * bond_leg
      mean(pmax(payoff, 0) / terminal_slice$money_market)
    },
    numeric(1)
  )

  expect_equal(function_values$omega, omega)
  expect_equal(function_values$present_value, manual, tolerance = 1e-10)
})
