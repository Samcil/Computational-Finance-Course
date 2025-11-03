test_that("euler_convergence_study returns expected RMSE profile", {
  result <- euler_convergence_study(
    initial_price = 100,
    drift = 0.05,
    volatility = 0.2,
    maturity = 1,
    n_paths = 4000,
    step_grid = c(8L, 16L, 32L, 64L),
    seed = 123
  )

  expect_s3_class(result, "tbl_df")
  expect_equal(result$scheme, rep("euler", 4))
  expect_equal(result$n_steps, c(8L, 16L, 32L, 64L))
  expect_equal(
    result$rmse,
    c(1.07076113483952, 0.747591725588231, 0.530532343079937, 0.389592292775292),
    tolerance = 1e-9
  )
  expect_true(all(diff(result$rmse) < 0))
})


test_that("milstein_convergence_study improves upon Euler scheme", {
  euler_result <- euler_convergence_study(
    initial_price = 100,
    drift = 0.05,
    volatility = 0.2,
    maturity = 1,
    n_paths = 4000,
    step_grid = c(8L, 16L, 32L, 64L),
    seed = 321
  )

  milstein_result <- milstein_convergence_study(
    initial_price = 100,
    drift = 0.05,
    volatility = 0.2,
    maturity = 1,
    n_paths = 4000,
    step_grid = c(8L, 16L, 32L, 64L),
    seed = 321
  )

  expect_equal(
    milstein_result$rmse,
    c(0.149279818238376, 0.0735864645682173, 0.0364757906602909, 0.0188459915600748),
    tolerance = 1e-9
  )
  expect_true(all(milstein_result$rmse < euler_result$rmse))
})


test_that("plot_convergence_rates validates inputs and returns ggplot", {
  convergence_data <- dplyr::bind_rows(
    euler_convergence_study(100, 0.05, 0.2, 1, 1000, c(8L, 16L), seed = 42),
    milstein_convergence_study(100, 0.05, 0.2, 1, 1000, c(8L, 16L), seed = 42)
  )

  plot_obj <- plot_convergence_rates(convergence_data)
  expect_s3_class(plot_obj, c("gg", "ggplot"))

  bad_data <- tibble::tibble(scheme = "euler", rmse = 0.1)
  expect_error(plot_convergence_rates(bad_data), "convergence_results must contain columns")
})
