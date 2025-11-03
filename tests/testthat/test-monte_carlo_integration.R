test_that("sample-mean Monte Carlo integration is accurate", {
  target <- exp(1) - 1
  res <- mc_integrate_sample_mean(
    g = exp,
    lower = 0,
    upper = 1,
    n_samples = 50000,
    seed = 2024
  )

  expect_equal(res$estimate, target, tolerance = 5e-3)
  expect_lt(res$std_error, 2e-2)
})

test_that("hit-or-miss Monte Carlo integration delivers reasonable accuracy", {
  g <- function(x) exp(x)
  res <- mc_integrate_hit_or_miss(
    g = g,
    lower = 0,
    upper = 1,
    lower_bound = 0,
    upper_bound = 3,
    n_samples = 100000,
    seed = 2025
  )

  target <- exp(1) - 1
  expect_equal(res$estimate, target, tolerance = 1.5e-2)
  expect_lt(res$std_error, 5e-2)
})
