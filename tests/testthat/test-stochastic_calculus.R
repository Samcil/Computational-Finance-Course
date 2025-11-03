test_that("stochastic integrals align with theoretical benchmarks", {
  res <- simulate_stochastic_integrals(
    maturity = 2,
    n_paths = 4000,
    n_steps = 400,
    integrand = base::identity,
    seed = 123
  )

  terminal <- dplyr::filter(res, .data$time == 2)
  brownian_final <- terminal$brownian
  ito_final <- terminal$ito_integral
  riemann_final <- terminal$riemann_integral
  quad_final <- terminal$quadratic_variation

  expect_equal(mean(ito_final), 0, tolerance = 0.02)
  ito_reference <- 0.5 * (brownian_final^2 - quad_final)
  expect_equal(ito_final, ito_reference, tolerance = 1e-10)

  expect_equal(mean(quad_final), 2, tolerance = 0.05)

  expect_equal(mean(riemann_final), 0, tolerance = 0.05)
  expect_equal(stats::var(riemann_final), (2^3) / 3, tolerance = 0.05)
})

test_that("quadratic variation diagnostics recover theoretical moments", {
  diag <- quadratic_variation_diagnostics(
    maturity = 1,
    n_paths = 6000,
    step_grid = c(50, 200, 500),
    seed = 456
  )

  expect_equal(diag$mean_sq_increment, diag$expected_mean, tolerance = 0.01)
  expect_equal(diag$var_sq_increment, diag$expected_variance, tolerance = 0.02)
})
