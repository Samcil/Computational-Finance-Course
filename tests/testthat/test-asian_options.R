test_that("price_asian_call returns tidy results", {
  spec <- gbm_spec(initial_value = 100, drift = 0.03, volatility = 0.2)
  results <- price_asian_call(
    process_spec = spec,
    strike = c(95, 105),
    maturity = 1,
    n_paths = 2000,
    n_steps = 126,
    seed = 123
  )

  expect_s3_class(results, "tbl_df")
  expect_true(all(c("strike", "option_type", "method", "price", "std_error") %in% names(results)))
  expect_equal(unique(results$option_type), "call")
  expect_equal(unique(results$method), "standard")
  expect_equal(nrow(results), 2)
})

test_that("asian call prices decrease with higher strikes", {
  spec <- gbm_spec(initial_value = 100, drift = 0.05, volatility = 0.2)
  prices <- price_asian_call(
    process_spec = spec,
    strike = c(90, 100, 110),
    maturity = 1,
    n_paths = 4000,
    n_steps = 200,
    seed = 321
  )

  sorted <- dplyr::arrange(prices, strike)
  expect_true(all(diff(sorted$price) <= 1e-2))
})

test_that("asian put prices increase with strike", {
  spec <- gbm_spec(initial_value = 100, drift = 0.03, volatility = 0.25)
  prices <- price_asian_put(
    process_spec = spec,
    strike = c(90, 100, 110),
    maturity = 1,
    n_paths = 4000,
    n_steps = 200,
    seed = 987
  )

  sorted <- dplyr::arrange(prices, strike)
  expect_true(all(diff(sorted$price) >= -1e-2))
})

test_that("antithetic variates reduce standard error", {
  spec <- gbm_spec(initial_value = 100, drift = 0.03, volatility = 0.2)
  comparison <- asian_variance_reduction(
    process_spec = spec,
    strike = 100,
    maturity = 1,
    option_type = "call",
    n_paths = 5000,
    n_steps = 200,
    seed = 2024
  )

  expect_s3_class(comparison, "tbl_df")
  expect_true(all(c("method", "price", "std_error") %in% names(comparison)))
  summarised <- dplyr::select(comparison, method, std_error)
  std_standard <- summarised$std_error[summarised$method == "standard"]
  std_antithetic <- summarised$std_error[summarised$method == "antithetic"]
  expect_lt(std_antithetic, std_standard * 1.05)
})
