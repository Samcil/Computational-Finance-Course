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
  
  expect_output(print(spec), "Poisson")
  expect_output(print(spec), "martingale")
})


test_that("simulate_paths.poisson_spec returns correct structure", {
  spec <- poisson_spec(intensity = 1.0)
  paths <- simulate_paths(spec, n_paths = 10, n_steps = 100, maturity = 10)
  
  expect_s3_class(paths, "tbl_df")
  expect_named(paths, c("path_id", "time", "count", "compensated_count"))
  expect_equal(length(unique(paths$path_id)), 10)
  expect_equal(nrow(paths), 10 * 101)  # 10 paths * 101 time points
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
  expect_lt(abs(mean_low - 0.5 * 10), 2)  # ≈ 5
  expect_lt(abs(mean_high - 2.0 * 10), 3)  # ≈ 20
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
