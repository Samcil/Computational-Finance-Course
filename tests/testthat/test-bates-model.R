test_that("bates_spec creates valid specification", {
  spec <- bates_spec(
    initial_price = 100,
    initial_variance = 0.04,
    risk_free_rate = 0.02,
    dividend_yield = 0.01,
    mean_reversion = 1.5,
    long_term_variance = 0.04,
    vol_of_vol = 0.5,
    correlation = -0.6,
    jump_intensity = 0.7,
    jump_mean = -0.1,
    jump_sd = 0.25,
    scheme = "aes"
  )

  expect_s3_class(spec, "bates_spec")
  expect_true(inherits(spec, "process_spec"))
  expect_equal(spec$jump_intensity, 0.7)
  expect_equal(spec$scheme, "aes")
})

test_that("simulate_paths.bates_spec generates paths for both schemes", {
  spec <- bates_spec(
    initial_price = 95,
    initial_variance = 0.05,
    risk_free_rate = 0.015,
    dividend_yield = 0,
    mean_reversion = 1.2,
    long_term_variance = 0.05,
    vol_of_vol = 0.35,
    correlation = -0.4,
    jump_intensity = 0.6,
    jump_mean = -0.08,
    jump_sd = 0.2,
    scheme = "aes"
  )

  aes_paths <- simulate_paths(spec, n_paths = 6, n_steps = 50, maturity = 1, seed = 123)
  expect_s3_class(aes_paths, "tbl_df")
  expect_equal(length(unique(aes_paths$path_id)), 6)
  expect_equal(length(unique(aes_paths$time)), 51)
  expect_equal(attr(aes_paths, "process_type"), "bates")
  expect_equal(attr(aes_paths, "scheme"), "aes")

  euler_paths <- simulate_paths(spec, n_paths = 6, n_steps = 50, maturity = 1, seed = 123, scheme = "euler")
  expect_equal(attr(euler_paths, "scheme"), "euler")
  expect_equal(length(unique(euler_paths$path_id)), 6)
})

test_that("bates_characteristic_function collapses to Heston when jumps are absent", {
  maturity <- 1
  strikes <- c(90, 100, 110)

  base_heston <- heston_spec(
    initial_price = 100,
    initial_variance = 0.04,
    risk_free_rate = 0.02,
    dividend_yield = 0,
    mean_reversion = 1.4,
    long_term_variance = 0.04,
    vol_of_vol = 0.3,
    correlation = -0.5,
    scheme = "aes"
  )

  bates_no_jumps <- bates_spec(
    initial_price = base_heston$initial_price,
    initial_variance = base_heston$initial_variance,
    risk_free_rate = base_heston$risk_free_rate,
    dividend_yield = base_heston$dividend_yield,
    mean_reversion = base_heston$mean_reversion,
    long_term_variance = base_heston$long_term_variance,
    vol_of_vol = base_heston$vol_of_vol,
    correlation = base_heston$correlation,
    jump_intensity = 0,
    jump_mean = 0,
    jump_sd = 0,
    scheme = "aes"
  )

  cf_bates <- bates_characteristic_function(bates_no_jumps, maturity)
  cf_heston <- function(u) {
    heston_characteristic_function(
      u = u,
      maturity = maturity,
      initial_price = base_heston$initial_price,
      initial_variance = base_heston$initial_variance,
      risk_free_rate = base_heston$risk_free_rate,
      dividend_yield = base_heston$dividend_yield,
      mean_reversion = base_heston$mean_reversion,
      long_term_variance = base_heston$long_term_variance,
      vol_of_vol = base_heston$vol_of_vol,
      correlation = base_heston$correlation
    )
  }

  u_grid <- seq(0, 5, length.out = 11)
  expect_equal(cf_bates(u_grid), cf_heston(u_grid), tolerance = 1e-10)

  bates_prices <- cos_call_put_price(
    cf = cf_bates,
    option_type = "call",
    spot = base_heston$initial_price,
    risk_free_rate = base_heston$risk_free_rate,
    maturity = maturity,
    strikes = strikes,
    n_terms = 256,
    truncation = 10
  )$price

  heston_prices <- price_heston_option_cos(
    base_heston,
    strikes = strikes,
    maturity = maturity,
    n_terms = 256,
    truncation = 10
  )$call_price

  expect_equal(bates_prices, heston_prices, tolerance = 1e-6)
})

test_that("bates_implied_volatility recovers option prices", {
  spec <- bates_spec(
    initial_price = 105,
    initial_variance = 0.05,
    risk_free_rate = 0.01,
    dividend_yield = 0,
    mean_reversion = 1.3,
    long_term_variance = 0.05,
    vol_of_vol = 0.4,
    correlation = -0.35,
    jump_intensity = 0.5,
    jump_mean = -0.05,
    jump_sd = 0.25,
    scheme = "aes"
  )

  strikes <- c(90, 100, 110)
  iv_tbl <- bates_implied_volatility(
    process_spec = spec,
    maturity = 1,
    strikes = strikes,
    option_type = "call",
    n_terms = 256,
    truncation = 10,
    vol_interval = c(1e-4, 1.5)
  )

  expect_false(any(is.na(iv_tbl$implied_volatility)))
  expect_true(all(iv_tbl$implied_volatility > 0))

  recovered_prices <- purrr::map2_dbl(
    strikes,
    iv_tbl$implied_volatility,
    ~ price_options(
      black_scholes_spec(
        option_type = "call",
        strike = .x,
        maturity = 1,
        risk_free_rate = spec$risk_free_rate,
        dividend_yield = spec$dividend_yield
      ),
      spot = spec$initial_price,
      volatility = .y
    )$price
  )

  expect_equal(recovered_prices, iv_tbl$option_price, tolerance = 1e-6)
})
