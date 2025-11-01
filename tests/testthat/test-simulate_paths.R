capture_cli_output <- function(expr) {
  messages <- testthat::capture_messages(force(expr))
  if (length(messages) == 0) {
    return("")
  }
  paste0(messages, collapse = "")
}

test_that("gbm_spec creates valid specification", {
  spec <- gbm_spec(
    initial_value = 100,
    drift = 0.05,
    volatility = 0.2
  )

  expect_s3_class(spec, "gbm_spec")
  expect_s3_class(spec, "process_spec")
  expect_equal(spec$initial_value, 100)
  expect_equal(spec$drift, 0.05)
  expect_equal(spec$volatility, 0.2)
  expect_equal(spec$process_type, "gbm")
})


test_that("abm_spec creates valid specification", {
  spec <- abm_spec(
    initial_value = 0,
    drift = 0.03,
    volatility = 0.15
  )

  expect_s3_class(spec, "abm_spec")
  expect_s3_class(spec, "process_spec")
  expect_equal(spec$initial_value, 0)
  expect_equal(spec$drift, 0.03)
  expect_equal(spec$volatility, 0.15)
  expect_equal(spec$process_type, "abm")
})


test_that("gbm_spec validates inputs", {
  expect_error(
    gbm_spec(initial_value = -10, drift = 0.05, volatility = 0.2),
    "initial_value"
  )

  expect_error(
    gbm_spec(initial_value = 100, drift = 0.05, volatility = -0.2),
    "volatility"
  )

  expect_error(
    gbm_spec(initial_value = 100, drift = Inf, volatility = 0.2),
    "drift"
  )
})


test_that("simulate_paths works with GBM spec", {
  spec <- gbm_spec(100, 0.05, 0.2)

  paths <- simulate_paths(
    process_spec = spec,
    n_paths = 10,
    n_steps = 50,
    maturity = 1.0,
    seed = 123
  )

  # Check structure
  expect_s3_class(paths, "tbl_df")
  expect_true(all(c("path_id", "time", "stock_price") %in% names(paths)))

  # Check dimensions
  expect_equal(nrow(paths), 10 * 51) # 10 paths * 51 time points
  expect_equal(length(unique(paths$path_id)), 10)
  expect_equal(length(unique(paths$time)), 51)

  # Check initial condition
  initial_prices <- paths |>
    dplyr::filter(time == 0) |>
    dplyr::pull(stock_price)
  expect_equal(initial_prices, rep(100, 10))

  # Check attributes
  expect_equal(attr(paths, "process_type"), "gbm")
  expect_s3_class(attr(paths, "spec"), "gbm_spec")
})


test_that("simulate_paths works with ABM spec", {
  spec <- abm_spec(0, 0.03, 0.15)

  paths <- simulate_paths(
    process_spec = spec,
    n_paths = 5,
    n_steps = 100,
    maturity = 1.0,
    seed = 456
  )

  # Check structure
  expect_s3_class(paths, "tbl_df")
  expect_true(all(c("path_id", "time", "value") %in% names(paths)))

  # Check dimensions
  expect_equal(nrow(paths), 5 * 101)
  expect_equal(length(unique(paths$path_id)), 5)

  # Check initial condition
  initial_values <- paths |>
    dplyr::filter(time == 0) |>
    dplyr::pull(value)
  expect_equal(initial_values, rep(0, 5))

  # Check attributes
  expect_equal(attr(paths, "process_type"), "abm")
})


test_that("simulate_paths is reproducible with seed", {
  spec <- gbm_spec(100, 0.05, 0.2)

  paths1 <- simulate_paths(spec, 10, 50, 1.0, seed = 789)
  paths2 <- simulate_paths(spec, 10, 50, 1.0, seed = 789)

  expect_equal(paths1, paths2)
})


test_that("simulate_paths validates inputs", {
  spec <- gbm_spec(100, 0.05, 0.2)

  expect_error(
    simulate_paths(spec, n_paths = 0, n_steps = 50, maturity = 1.0),
    "n_paths"
  )

  expect_error(
    simulate_paths(spec, n_paths = 10, n_steps = -5, maturity = 1.0),
    "n_steps"
  )

  expect_error(
    simulate_paths(spec, n_paths = 10, n_steps = 50, maturity = -1.0),
    "maturity"
  )

  expect_error(
    simulate_paths("not a spec", n_paths = 10, n_steps = 50, maturity = 1.0),
    "process_spec"
  )
})


test_that("pipeline workflow works end-to-end", {
  # This tests the entire tidyverse pipeline
  result <- gbm_spec(100, 0.05, 0.2) |>
    simulate_paths(n_paths = 20, n_steps = 100, maturity = 1.0) |>
    dplyr::filter(path_id <= 5) |>
    dplyr::select(path_id, time, stock_price)

  expect_s3_class(result, "tbl_df")
  expect_equal(length(unique(result$path_id)), 5)
  expect_true(all(result$path_id <= 5))
})


test_that("spec print methods work", {
  gbm <- gbm_spec(100, 0.05, 0.2)
  abm <- abm_spec(0, 0.03, 0.15)

  old_opts <- options(
    cli.unicode = FALSE,
    cli.num_colors = 1,
    cli.width = 80
  )
  on.exit(options(old_opts), add = TRUE)

  gbm_output <- capture_cli_output(print(gbm))
  abm_output <- capture_cli_output(print(abm))

  expect_match(gbm_output, "Geometric", fixed = TRUE)
  expect_match(abm_output, "Arithmetic", fixed = TRUE)
})


# Poisson Process Tests -----------------------------------------------------

test_that("poisson_spec creates valid specification", {
  spec <- poisson_spec(intensity = 1.0, initial_value = 0)

  expect_s3_class(spec, "poisson_spec")
  expect_s3_class(spec, "process_spec")
  expect_equal(spec$intensity, 1.0)
  expect_equal(spec$initial_value, 0)
})


test_that("poisson_spec validates parameters", {
  expect_error(
    poisson_spec(intensity = -1.0),
    "intensity"
  )

  expect_error(
    poisson_spec(intensity = NA),
    "intensity"
  )
})


test_that("poisson_spec print method works", {
  spec <- poisson_spec(intensity = 2.0)

  old_opts <- options(
    cli.unicode = FALSE,
    cli.num_colors = 1,
    cli.width = 80
  )
  on.exit(options(old_opts), add = TRUE)

  poisson_output <- capture_cli_output(print(spec))

  expect_match(poisson_output, "Poisson", fixed = TRUE)
  expect_match(poisson_output, "martingale", ignore.case = TRUE)
})


test_that("simulate_paths.poisson_spec returns correct structure", {
  spec <- poisson_spec(intensity = 1.0)
  paths <- simulate_paths(spec, n_paths = 10, n_steps = 100, maturity = 10)

  expect_s3_class(paths, "tbl_df")
  expect_named(paths, c("path_id", "time", "count", "compensated_count"))
  expect_equal(length(unique(paths$path_id)), 10)
  expect_equal(nrow(paths), 10 * 101) # 10 paths * 101 time points
})


test_that("poisson process starts at initial value", {
  spec <- poisson_spec(intensity = 1.0, initial_value = 5)
  paths <- simulate_paths(spec, n_paths = 10, n_steps = 100, maturity = 10)

  initial_counts <- paths[paths$time == 0, "count", drop = TRUE]
  expect_true(all(initial_counts == 5))
})


test_that("compensated poisson process has zero mean", {
  spec <- poisson_spec(intensity = 1.0)
  paths <- simulate_paths(spec, n_paths = 1000, n_steps = 100, maturity = 10, seed = 42)

  # At final time, compensated process should have mean ≈ 0
  final_compensated <- paths[paths$time == 10, "compensated_count", drop = TRUE]
  mean_final <- mean(final_compensated)

  # Should be close to 0 (within 0.2 for 1000 paths)
  expect_lt(abs(mean_final), 0.2)
})


test_that("poisson intensity affects jump frequency", {
  # Low intensity
  spec_low <- poisson_spec(intensity = 0.5)
  paths_low <- simulate_paths(spec_low, n_paths = 100, n_steps = 100, maturity = 10, seed = 123)
  final_low <- paths_low[paths_low$time == 10, "count", drop = TRUE]
  mean_low <- mean(final_low)

  # High intensity
  spec_high <- poisson_spec(intensity = 2.0)
  paths_high <- simulate_paths(spec_high, n_paths = 100, n_steps = 100, maturity = 10, seed = 123)
  final_high <- paths_high[paths_high$time == 10, "count", drop = TRUE]
  mean_high <- mean(final_high)

  # Higher intensity should give more jumps
  expect_gt(mean_high, mean_low)

  # Should be close to theoretical E[N(T)] = λT
  expect_lt(abs(mean_low - 0.5 * 10), 2) # ≈ 5
  expect_lt(abs(mean_high - 2.0 * 10), 3) # ≈ 20
})


test_that("poisson process pipeline works", {
  result <- poisson_spec(1.5) |>
    simulate_paths(n_paths = 20, n_steps = 200, maturity = 15) |>
    dplyr::filter(path_id <= 5) |>
    dplyr::select(path_id, time, count, compensated_count)

  expect_s3_class(result, "tbl_df")
  expect_equal(length(unique(result$path_id)), 5)
  expect_true(all(result$path_id <= 5))
  expect_true("count" %in% names(result))
  expect_true("compensated_count" %in% names(result))
})


# CIR Process Tests ---------------------------------------------------------

test_that("cir_spec creates valid specification", {
  spec <- cir_spec(0.05, mean_reversion = 3, long_term_mean = 0.04, volatility = 0.25)
  expect_s3_class(spec, "cir_spec")
  expect_s3_class(spec, "process_spec")
  expect_equal(spec$initial_value, 0.05)
  expect_equal(spec$mean_reversion, 3)
})

test_that("simulate_paths.cir_spec returns tidy output", {
  spec <- cir_spec(0.05, 2.5, 0.04, 0.15)
  paths <- simulate_paths(spec, n_paths = 20, n_steps = 200, maturity = 1, seed = 321)
  expect_s3_class(paths, "tbl_df")
  expect_named(paths, c("path_id", "time", "variance"))
  expect_equal(length(unique(paths$path_id)), 20)
  expect_equal(length(unique(paths$time)), 201)
  expect_equal(attr(paths, "process_type"), "cir")
})


# Correlated Brownian Motion Tests -----------------------------------------

test_that("correlated_bm_spec validates structure", {
  cov_matrix <- matrix(c(0.04, 0.02, 0.02, 0.09), nrow = 2)
  spec <- correlated_bm_spec(
    initial_values = c(100, 95),
    drift = c(0.05, 0.04),
    covariance_matrix = cov_matrix,
    component_names = c("Asset_A", "Asset_B")
  )
  expect_s3_class(spec, "correlated_bm_spec")
  expect_equal(spec$component_names, c("Asset_A", "Asset_B"))
})

test_that("simulate_paths.correlated_bm_spec returns tidy multivariate output", {
  cov_matrix <- matrix(c(0.04, 0.02, 0.02, 0.09), nrow = 2)
  spec <- correlated_bm_spec(
    initial_values = c(100, 95),
    drift = c(0.05, 0.04),
    covariance_matrix = cov_matrix,
    component_names = c("Asset_A", "Asset_B")
  )
  paths <- simulate_paths(spec, n_paths = 5, n_steps = 50, maturity = 1, seed = 99)
  expect_s3_class(paths, "tbl_df")
  expect_true(all(c("path_id", "time", "component", "value") %in% names(paths)))
  expect_equal(sort(unique(paths$component)), c("Asset_A", "Asset_B"))
  expect_equal(attr(paths, "process_type"), "correlated_bm")
})


# Heston Model Tests -------------------------------------------------------

test_that("heston_spec creates valid specification", {
  spec <- heston_spec(
    initial_price = 100,
    initial_variance = 0.04,
    risk_free_rate = 0.02,
    mean_reversion = 1.5,
    long_term_variance = 0.04,
    vol_of_vol = 0.6,
    correlation = -0.7
  )
  expect_s3_class(spec, "heston_spec")
  expect_s3_class(spec, "process_spec")
  expect_equal(spec$scheme, "euler")
  expect_equal(spec$initial_price, 100)
})

test_that("simulate_paths.heston_spec returns joint price and variance", {
  spec <- heston_spec(
    initial_price = 100,
    initial_variance = 0.04,
    risk_free_rate = 0.02,
    mean_reversion = 1.3,
    long_term_variance = 0.04,
    vol_of_vol = 0.5,
    correlation = -0.6
  )
  paths <- simulate_paths(spec, n_paths = 8, n_steps = 120, maturity = 1, seed = 77)
  expect_s3_class(paths, "tbl_df")
  expect_true(all(c("path_id", "time", "stock_price", "variance") %in% names(paths)))
  expect_equal(length(unique(paths$path_id)), 8)
  expect_equal(length(unique(paths$time)), 121)
  expect_equal(attr(paths, "process_type"), "heston")
})

test_that("simulate_paths.heston_spec supports AES scheme", {
  spec <- heston_spec(
    initial_price = 90,
    initial_variance = 0.09,
    risk_free_rate = 0.015,
    mean_reversion = 1.2,
    long_term_variance = 0.08,
    vol_of_vol = 0.7,
    correlation = -0.5,
    scheme = "aes"
  )

  paths <- simulate_paths(spec, n_paths = 6, n_steps = 60, maturity = 0.5, seed = 101)
  expect_s3_class(paths, "tbl_df")
  expect_true(all(paths$variance >= 0))
  expect_equal(attr(paths, "process_type"), "heston")
  expect_equal(attr(paths, "scheme"), "aes")
})

test_that("simulate_paths can override Heston scheme via arguments", {
  spec <- heston_spec(
    initial_price = 100,
    initial_variance = 0.05,
    risk_free_rate = 0.02,
    mean_reversion = 1.4,
    long_term_variance = 0.04,
    vol_of_vol = 0.6,
    correlation = -0.5,
    scheme = "euler"
  )

  override_paths <- simulate_paths(
    spec,
    n_paths = 3,
    n_steps = 30,
    maturity = 0.5,
    seed = 2024,
    scheme = "aes"
  )

  expect_s3_class(override_paths, "tbl_df")
  expect_equal(attr(override_paths, "scheme"), "aes")
  expect_equal(spec$scheme, "euler")
  expect_error(
    simulate_paths(spec, n_paths = 2, n_steps = 10, maturity = 0.1, scheme = "invalid"),
    "scheme"
  )
})

test_that("generate_heston_paths convenience wrappers match simulate_paths", {
  base_spec <- heston_spec(
    initial_price = 105,
    initial_variance = 0.05,
    risk_free_rate = 0.01,
    mean_reversion = 1.8,
    long_term_variance = 0.04,
    vol_of_vol = 0.6,
    correlation = -0.4
  )

  direct_paths <- simulate_paths(base_spec, n_paths = 4, n_steps = 40, maturity = 0.5, seed = 555)
  wrapper_paths <- generate_heston_paths_euler(
    n_paths = 4,
    n_steps = 40,
    maturity = 0.5,
    initial_price = 105,
    initial_variance = 0.05,
    risk_free_rate = 0.01,
    mean_reversion = 1.8,
    long_term_variance = 0.04,
    vol_of_vol = 0.6,
    correlation = -0.4,
    seed = 555
  )

  expect_equal(direct_paths, wrapper_paths)

  aes_paths <- generate_heston_paths_aes(
    n_paths = 3,
    n_steps = 30,
    maturity = 0.5,
    initial_price = 95,
    initial_variance = 0.07,
    risk_free_rate = 0.02,
    mean_reversion = 1.1,
    long_term_variance = 0.06,
    vol_of_vol = 0.5,
    correlation = -0.3,
    seed = 99
  )

  expect_s3_class(aes_paths, "tbl_df")
  expect_equal(attr(aes_paths, "scheme"), "aes")
  expect_true(all(aes_paths$variance >= 0))
})

test_that("heston_characteristic_function returns valid complex values", {
  cf_vals <- heston_characteristic_function(
    u = c(0, 0.5, 1),
    maturity = 1,
    initial_price = 100,
    initial_variance = 0.04,
    risk_free_rate = 0.02,
    dividend_yield = 0.01,
    mean_reversion = 1.5,
    long_term_variance = 0.04,
    vol_of_vol = 0.6,
    correlation = -0.7
  )

  expect_type(cf_vals, "complex")
  expect_equal(cf_vals[1], 1 + 0i)
  expect_true(all(Mod(cf_vals) <= 1.1))
})

test_that("price_heston_option_cos enforces put call parity", {
  spec <- heston_spec(
    initial_price = 100,
    initial_variance = 0.04,
    risk_free_rate = 0.02,
    dividend_yield = 0.01,
    mean_reversion = 1.5,
    long_term_variance = 0.04,
    vol_of_vol = 0.5,
    correlation = -0.6
  )

  strikes <- c(90, 100, 110)
  prices <- price_heston_option_cos(spec, strikes = strikes, maturity = 1, n_terms = 512, truncation = 12)

  expect_s3_class(prices, "tbl_df")
  expect_true(all(prices$call_price > 0))
  expect_true(all(prices$put_price > 0))

  parity_lhs <- prices$call_price - prices$put_price
  parity_rhs <- spec$initial_price * exp(-spec$dividend_yield * 1) - strikes * exp(-spec$risk_free_rate * 1)
  expect_lt(max(abs(parity_lhs - parity_rhs)), 1e-2)
})


# Merton Jump-Diffusion Tests ---------------------------------------------

test_that("merton_spec creates valid specification", {
  spec <- merton_spec(
    initial_price = 120,
    risk_free_rate = 0.03,
    volatility = 0.25,
    jump_intensity = 1.2,
    jump_mean = -0.05,
    jump_sd = 0.2
  )
  expect_s3_class(spec, "merton_spec")
  expect_s3_class(spec, "process_spec")
  expect_equal(spec$jump_intensity, 1.2)
})

test_that("simulate_paths.merton_spec returns jump-diffusion paths", {
  spec <- merton_spec(
    initial_price = 100,
    risk_free_rate = 0.02,
    volatility = 0.2,
    jump_intensity = 0.8,
    jump_mean = -0.1,
    jump_sd = 0.3
  )
  paths <- simulate_paths(spec, n_paths = 12, n_steps = 150, maturity = 1, seed = 88)
  expect_s3_class(paths, "tbl_df")
  expect_true(all(c("path_id", "time", "stock_price") %in% names(paths)))
  expect_equal(length(unique(paths$path_id)), 12)
  expect_equal(length(unique(paths$time)), 151)
  expect_equal(attr(paths, "process_type"), "merton")
})


# Black-Scholes Pricing Tests ----------------------------------------------

test_that("price_options.black_scholes_spec returns prices and Greeks", {
  spec <- black_scholes_spec("call", strike = 100, maturity = 1, risk_free_rate = 0.02)
  pricing <- price_options(spec, spot = c(95, 100), volatility = c(0.15, 0.2))
  expect_s3_class(pricing, "tbl_df")
  expect_true(all(c("spot", "volatility", "price", "delta", "gamma", "vega", "theta", "rho") %in% names(pricing)))
  expect_equal(nrow(pricing), 4)
  expect_true(all(pricing$price > 0))
})
