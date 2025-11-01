test_that("generate_gbm_abm_paths produces correct dimensions", {
  paths <- generate_gbm_abm_paths(
    n_paths = 10,
    n_steps = 100,
    maturity = 1.0,
    interest_rate = 0.05,
    volatility = 0.2,
    initial_price = 100,
    return_format = "tidy"
  )
  
  # Check it's a tibble
  expect_s3_class(paths, "tbl_df")
  
  # Check dimensions: 10 paths * (100 + 1) time points = 1010 rows
  expect_equal(nrow(paths), 10 * 101)
  
  # Check column names
  expect_named(paths, c("path_id", "time", "log_price", "stock_price"))
})

test_that("generate_gbm_abm_paths matrix format works", {
  paths <- generate_gbm_abm_paths(
    n_paths = 10,
    n_steps = 100,
    maturity = 1.0,
    interest_rate = 0.05,
    volatility = 0.2,
    initial_price = 100,
    return_format = "matrix"
  )
  
  # Check it's a list
  expect_type(paths, "list")
  
  # Check components
  expect_named(paths, c("time", "log_price", "stock_price"))
  
  # Check dimensions
  expect_length(paths$time, 101)
  expect_equal(dim(paths$log_price), c(10, 101))
  expect_equal(dim(paths$stock_price), c(10, 101))
})

test_that("GBM paths start at initial price", {
  initial_price <- 100
  paths <- generate_gbm_abm_paths(
    n_paths = 5,
    n_steps = 50,
    maturity = 1.0,
    interest_rate = 0.05,
    volatility = 0.2,
    initial_price = initial_price,
    return_format = "matrix"
  )
  
  # All paths should start at initial_price
  expect_equal(paths$stock_price[, 1], rep(initial_price, 5))
  
  # Log prices should start at log(initial_price)
  expect_equal(paths$log_price[, 1], rep(log(initial_price), 5))
})

test_that("Time grid is correct", {
  maturity <- 1.5
  n_steps <- 150
  paths <- generate_gbm_abm_paths(
    n_paths = 5,
    n_steps = n_steps,
    maturity = maturity,
    interest_rate = 0.05,
    volatility = 0.2,
    initial_price = 100,
    return_format = "matrix"
  )
  
  # Check time grid starts at 0
  expect_equal(paths$time[1], 0)
  
  # Check time grid ends at maturity
  expect_equal(paths$time[length(paths$time)], maturity)
  
  # Check time step size
  dt <- maturity / n_steps
  expect_equal(diff(paths$time), rep(dt, n_steps), tolerance = 1e-10)
})

test_that("Reproducibility with seed", {
  params <- list(
    n_paths = 10,
    n_steps = 100,
    maturity = 1.0,
    interest_rate = 0.05,
    volatility = 0.2,
    initial_price = 100,
    return_format = "matrix"
  )
  
  paths1 <- do.call(generate_gbm_abm_paths, c(params, seed = 42))
  paths2 <- do.call(generate_gbm_abm_paths, c(params, seed = 42))
  
  # Should be identical with same seed
  expect_equal(paths1$stock_price, paths2$stock_price)
  expect_equal(paths1$log_price, paths2$log_price)
})

test_that("Input validation works", {
  # Negative paths
  expect_error(
    generate_gbm_abm_paths(
      n_paths = -5,
      n_steps = 100,
      maturity = 1.0,
      interest_rate = 0.05,
      volatility = 0.2,
      initial_price = 100
    )
  )
  
  # Non-positive volatility
  expect_error(
    generate_gbm_abm_paths(
      n_paths = 10,
      n_steps = 100,
      maturity = 1.0,
      interest_rate = 0.05,
      volatility = -0.2,
      initial_price = 100
    )
  )
  
  # Non-positive initial price
  expect_error(
    generate_gbm_abm_paths(
      n_paths = 10,
      n_steps = 100,
      maturity = 1.0,
      interest_rate = 0.05,
      volatility = 0.2,
      initial_price = 0
    )
  )
})
