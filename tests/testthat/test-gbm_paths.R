test_that("simulate_paths.gbm_spec returns tidy output", {
  spec <- gbm_spec(initial_value = 100, drift = 0.05, volatility = 0.2)
  paths <- simulate_paths(spec, n_paths = 10, n_steps = 100, maturity = 1, seed = 123)

  expect_s3_class(paths, "tbl_df")
  expect_true(all(c("path_id", "time", "stock_price") %in% names(paths)))
  expect_equal(nrow(paths), 10 * (100 + 1))
  expect_equal(length(unique(paths$path_id)), 10)
  expect_equal(length(unique(paths$time)), 101)
  expect_equal(attr(paths, "process_type"), "gbm")
})

test_that("gbm paths start at initial value", {
  spec <- gbm_spec(initial_value = 120, drift = 0.02, volatility = 0.18)
  paths <- simulate_paths(spec, n_paths = 4, n_steps = 60, maturity = 1, seed = 321)

  start_values <- dplyr::filter(paths, time == 0)
  expect_equal(start_values$stock_price, rep(120, 4))
})

test_that("gbm time grid is uniform and ends at maturity", {
  spec <- gbm_spec(initial_value = 90, drift = 0.03, volatility = 0.25)
  maturity <- 1.5
  n_steps <- 150
  paths <- simulate_paths(spec, n_paths = 3, n_steps = n_steps, maturity = maturity, seed = 999)

  time_grid <- sort(unique(paths$time))
  expect_equal(time_grid[[1]], 0)
  expect_equal(time_grid[[length(time_grid)]], maturity)
  expect_equal(diff(time_grid), rep(maturity / n_steps, n_steps), tolerance = 1e-10)
})

test_that("gbm simulation is reproducible with seed", {
  spec <- gbm_spec(initial_value = 100, drift = 0.04, volatility = 0.15)
  params <- list(process_spec = spec, n_paths = 6, n_steps = 80, maturity = 1.2)

  paths1 <- do.call(simulate_paths, c(params, seed = 42))
  paths2 <- do.call(simulate_paths, c(params, seed = 42))

  expect_equal(paths1$stock_price, paths2$stock_price)
})

test_that("gbm_spec validates inputs", {
  expect_error(gbm_spec(initial_value = -100, drift = 0.05, volatility = 0.2))
  expect_error(gbm_spec(initial_value = 100, drift = 0.05, volatility = -0.2))
})
