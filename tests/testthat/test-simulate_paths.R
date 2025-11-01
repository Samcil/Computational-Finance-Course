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
  expect_equal(nrow(paths), 10 * 51)  # 10 paths * 51 time points
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
  
  # Should not error
  expect_output(print(gbm), "Geometric")
  expect_output(print(abm), "Arithmetic")
})
