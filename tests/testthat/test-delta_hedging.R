test_that("delta hedging without intermediate rebalancing matches static replication", {
  option_spec <- black_scholes_spec(
    option_type = "call",
    strike = 100,
    maturity = 1,
    risk_free_rate = 0.02
  )

  spot_paths <- tibble::tibble(
    path_id = 1L,
    time = c(0, 1),
    stock_price = c(100, 110)
  )

  hedged <- simulate_delta_hedge_bs(
    spot_paths = spot_paths,
    option_spec = option_spec,
    volatility = 0.2,
    hedge_steps = 1,
    transaction_cost = 0
  )

  bs_metrics <- price_options(option_spec, spot = 100, volatility = 0.2)
  delta0 <- bs_metrics$delta
  price0 <- bs_metrics$price
  expected_cash0 <- price0 - delta0 * 100
  expected_cash_final <- expected_cash0 * exp(0.02)
  expected_portfolio <- expected_cash_final + delta0 * 110

  expect_equal(nrow(hedged), 2L)
  expect_equal(hedged$portfolio_value[2], as.numeric(expected_portfolio), tolerance = 1e-8)
  payoff <- max(110 - 100, 0)
  expect_equal(hedged$hedging_error[2], hedged$portfolio_value[2] - payoff, tolerance = 1e-10)
  expect_true(all(is.na(hedged$hedging_error[1])))
})


test_that("transaction costs reduce cash by proportional notional", {
  option_spec <- black_scholes_spec(
    option_type = "put",
    strike = 105,
    maturity = 1,
    risk_free_rate = 0,
    dividend_yield = 0
  )

  spot_paths <- tibble::tibble(
    path_id = 1L,
    time = c(0, 0.5, 1),
    stock_price = c(100, 105, 102)
  )

  hedged_no_cost <- simulate_delta_hedge_bs(
    spot_paths = spot_paths,
    option_spec = option_spec,
    volatility = 0.25,
    hedge_steps = 2,
    transaction_cost = 0
  )

  hedged_cost <- simulate_delta_hedge_bs(
    spot_paths = spot_paths,
    option_spec = option_spec,
    volatility = 0.25,
    hedge_steps = 2,
    transaction_cost = 0.01
  )

  delta_series <- hedged_no_cost$delta
  delta_changes <- c(delta_series[1], diff(delta_series))
  trade_mask <- abs(delta_changes) > 1e-10
  expected_cost <- sum(abs(delta_changes[trade_mask]) * spot_paths$stock_price[trade_mask]) * 0.01

  bond_diff <- hedged_no_cost$bond_position - hedged_cost$bond_position
  expect_equal(bond_diff[3], expected_cost, tolerance = 1e-8)
})


test_that("Black-Scholes delta hedging keeps errors small under GBM", {
  gbm_paths <- simulate_paths(
    gbm_spec(initial_value = 100, drift = 0.02, volatility = 0.2),
    n_paths = 256,
    n_steps = 128,
    maturity = 1,
    seed = 42
  )

  option_spec <- black_scholes_spec(
    option_type = "call",
    strike = 100,
    maturity = 1,
    risk_free_rate = 0.02
  )

  hedged <- simulate_delta_hedge_bs(
    spot_paths = gbm_paths,
    option_spec = option_spec,
    volatility = 0.2,
    hedge_steps = 128,
    transaction_cost = 0
  )

  terminal_errors <- hedged |>
    dplyr::filter(.data$time == 1) |>
    dplyr::pull(.data$hedging_error)

  expect_lt(abs(mean(terminal_errors)), 0.25)
  expect_lt(stats::sd(terminal_errors), 1)
})


test_that("jump risk inflates hedging error dispersion", {
  gbm_paths <- simulate_paths(
    gbm_spec(initial_value = 100, drift = 0.02, volatility = 0.2),
    n_paths = 256,
    n_steps = 128,
    maturity = 1,
    seed = 99
  )

  merton_paths <- simulate_paths(
    merton_spec(
      initial_price = 100,
      risk_free_rate = 0.02,
      volatility = 0.2,
      jump_intensity = 0.5,
      jump_mean = -0.1,
      jump_sd = 0.2
    ),
    n_paths = 256,
    n_steps = 128,
    maturity = 1,
    seed = 99
  )

  option_spec <- black_scholes_spec(
    option_type = "call",
    strike = 100,
    maturity = 1,
    risk_free_rate = 0.02
  )

  hedged_gbm <- simulate_delta_hedge_bs(
    spot_paths = gbm_paths,
    option_spec = option_spec,
    volatility = 0.2,
    hedge_steps = 128,
    transaction_cost = 0
  )

  hedged_jump <- simulate_delta_hedge_jumps(
    spot_paths = merton_paths,
    option_spec = option_spec,
    diffusion_volatility = 0.2,
    hedge_steps = 128,
    transaction_cost = 0
  )

  err_gbm <- hedged_gbm |>
    dplyr::filter(.data$time == 1) |>
    dplyr::pull(.data$hedging_error)
  err_jump <- hedged_jump |>
    dplyr::filter(.data$time == 1) |>
    dplyr::pull(.data$hedging_error)

  expect_gt(mean(abs(err_jump)), mean(abs(err_gbm)))
  expect_gt(stats::sd(err_jump), stats::sd(err_gbm))
})
