test_that("price_barrier_option returns expected columns", {
  spec <- gbm_spec(initial_value = 100, drift = 0.03, volatility = 0.2)
  results <- price_barrier_option(
    process_spec = spec,
    strike = c(95, 105),
    barrier = c(110, 120),
    maturity = 1,
    barrier_type = "up-and-out",
    option_type = "call",
    n_paths = 3000,
    n_steps = 150,
    seed = 101
  )

  expect_s3_class(results, "tbl_df")
  expect_true(all(c(
    "strike",
    "barrier",
    "option_type",
    "barrier_type",
    "price",
    "std_error",
    "hit_probability",
    "n_paths",
    "n_steps"
  ) %in% names(results)))
  expect_equal(unique(results$option_type), "call")
  expect_equal(unique(results$barrier_type), "up-and-out")
  expect_equal(nrow(results), 4)
})

test_that("up-and-out call with barrier below spot has zero value", {
  spec <- gbm_spec(initial_value = 100, drift = 0.02, volatility = 0.2)
  results <- price_barrier_option(
    process_spec = spec,
    strike = 100,
    barrier = 90,
    maturity = 1,
    barrier_type = "up-and-out",
    option_type = "call",
    n_paths = 2000,
    n_steps = 100,
    seed = 202
  )

  expect_equal(results$price, 0)
  expect_equal(results$hit_probability, 1)
})

test_that("up-and-in plus up-and-out approximates vanilla call", {
  spec <- gbm_spec(initial_value = 100, drift = 0.03, volatility = 0.25)
  strike <- 100
  barrier <- 130
  maturity <- 1

  bs_spec <- black_scholes_spec(
    option_type = "call",
    strike = strike,
    maturity = maturity,
    risk_free_rate = spec$drift
  )
  vanilla_price <- price_options(
    bs_spec,
    spot = spec$initial_value,
    volatility = spec$volatility
  )$price

  out_price <- price_barrier_option(
    process_spec = spec,
    strike = strike,
    barrier = barrier,
    maturity = maturity,
    barrier_type = "up-and-out",
    option_type = "call",
    n_paths = 8000,
    n_steps = 200,
    seed = 303
  )$price

  in_price <- price_barrier_option(
    process_spec = spec,
    strike = strike,
    barrier = barrier,
    maturity = maturity,
    barrier_type = "up-and-in",
    option_type = "call",
    n_paths = 8000,
    n_steps = 200,
    seed = 303
  )$price

  expect_lt(abs((in_price + out_price) - vanilla_price), 0.75)
})

test_that("barrier_hit_probability returns sensible values", {
  spec <- gbm_spec(initial_value = 100, drift = 0.02, volatility = 0.2)
  probs <- barrier_hit_probability(
    process_spec = spec,
    barrier = c(80, 130),
    maturity = 1,
    barrier_type = "down",
    n_paths = 3000,
    n_steps = 150,
    seed = 404
  )

  expect_s3_class(probs, "tbl_df")
  expect_equal(nrow(probs), 2)
  expect_true(all(probs$probability >= 0 & probs$probability <= 1))
})
