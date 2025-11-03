test_that("trapezoidal and Simpson integrals converge for smooth functions", {
  step_grid <- c(20, 40, 80, 160)
  sine_results <- benchmark_deterministic_integrators(
    g = sin,
    lower = 0,
    upper = pi,
    step_grid = step_grid,
    exact = 2
  )

  sine_errors <- split(sine_results, sine_results$method)
  expect_true(all(vapply(sine_errors, function(df) all(diff(df$abs_error) < 0), logical(1))))
  expect_lt(min(sine_results$abs_error), 1e-7)

  exp_results <- benchmark_deterministic_integrators(
    g = exp,
    lower = 0,
    upper = 1,
    step_grid = step_grid,
    exact = exp(1) - 1
  )

  exp_errors <- split(exp_results, exp_results$method)
  expect_true(all(vapply(exp_errors, function(df) all(diff(df$abs_error) < 0), logical(1))))
  final_simpson <- exp_results$abs_error[exp_results$n_steps == max(step_grid) & exp_results$method == "simpson"]
  expect_lt(final_simpson, 1e-8)
})

test_that("Simpson requires even grid", {
  expect_error(simpson_integral(sin, 0, pi, n_steps = 5))
})

test_that("benchmark filters odd grids for Simpson", {
  res <- benchmark_deterministic_integrators(
    g = sin,
    lower = 0,
    upper = pi,
    step_grid = c(21, 42),
    exact = 2
  )

  expect_setequal(unique(res$method), c("trapezoidal", "simpson"))
  expect_true(all(res$n_steps %in% c(21L, 42L)))
  simpson_rows <- res[res$method == "simpson", , drop = FALSE]
  expect_true(all(simpson_rows$n_steps %% 2 == 0))
})
