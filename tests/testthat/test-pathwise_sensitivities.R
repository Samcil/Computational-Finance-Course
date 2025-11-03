test_that("pathwise_delta matches analytic delta for call options", {
  maturity <- 1
  risk_free_rate <- 0.03
  volatility <- 0.2
  initial_price <- 100
  strikes <- c(90, 100, 110)

  spec <- gbm_spec(
    initial_value = initial_price,
    drift = risk_free_rate,
    volatility = volatility
  )

  paths <- simulate_paths(
    process_spec = spec,
    n_paths = 5000,
    n_steps = 32,
    maturity = maturity,
    seed = 101
  )

  delta_tbl <- pathwise_delta(
    paths,
    strikes = strikes,
    risk_free_rate = risk_free_rate,
    option_type = "call",
    initial_price = initial_price
  )

  analytic_delta <- purrr::map_dbl(
    strikes,
    ~ price_options(
      black_scholes_spec(
        option_type = "call",
        strike = .x,
        maturity = maturity,
        risk_free_rate = risk_free_rate
      ),
      spot = initial_price,
      volatility = volatility
    )$delta
  )

  expect_equal(delta_tbl$delta_pathwise, analytic_delta, tolerance = 0.08)
})


test_that("pathwise_vega matches analytic vega for put options", {
  maturity <- 0.75
  risk_free_rate <- 0.025
  volatility <- 0.18
  initial_price <- 95
  strikes <- c(90, 95, 100)

  spec <- gbm_spec(
    initial_value = initial_price,
    drift = risk_free_rate,
    volatility = volatility
  )

  paths <- simulate_paths(
    process_spec = spec,
    n_paths = 5000,
    n_steps = 32,
    maturity = maturity,
    seed = 202
  )

  vega_tbl <- pathwise_vega(
    paths,
    strikes = strikes,
    risk_free_rate = risk_free_rate,
    option_type = "put",
    volatility = volatility,
    initial_price = initial_price
  )

  analytic_vega <- purrr::map_dbl(
    strikes,
    ~ price_options(
      black_scholes_spec(
        option_type = "put",
        strike = .x,
        maturity = maturity,
        risk_free_rate = risk_free_rate
      ),
      spot = initial_price,
      volatility = volatility
    )$vega
  )

  expect_equal(vega_tbl$vega_pathwise, analytic_vega, tolerance = 0.12)
})


test_that("compare_with_finite_diff agrees with analytic Greeks", {
  maturity <- 1
  risk_free_rate <- 0.04
  volatility <- 0.22
  initial_price <- 120
  strikes <- c(100, 110, 120)

  spec <- gbm_spec(
    initial_value = initial_price,
    drift = risk_free_rate,
    volatility = volatility
  )

  comparison <- compare_with_finite_diff(
    process_spec = spec,
    strikes = strikes,
    maturity = maturity,
    n_paths = 6000,
    n_steps = 32,
    risk_free_rate = risk_free_rate,
    volatility = volatility,
    option_type = "call",
    bump_spot = 0.5,
    bump_vol = 5e-4,
    seed = 303
  )

  expect_setequal(
    names(comparison),
    c(
      "strike",
      "delta_pathwise",
      "vega_pathwise",
      "delta_fd",
      "vega_fd",
      "delta_analytic",
      "vega_analytic"
    )
  )

  expect_equal(comparison$delta_pathwise, comparison$delta_analytic, tolerance = 0.1)
  expect_equal(comparison$delta_fd, comparison$delta_analytic, tolerance = 0.12)
  expect_equal(comparison$vega_pathwise, comparison$vega_analytic, tolerance = 0.15)
  expect_equal(comparison$vega_fd, comparison$vega_analytic, tolerance = 0.18)
})
