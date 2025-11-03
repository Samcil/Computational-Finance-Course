test_that("convexity correction collapses when volatility is zero", {
  curve <- tibble::tibble(tenor = 0:5, discount_factor = exp(-0.02 * (0:5)))
  spec_zero <- short_rate_spec(model = "ho_lee", volatility = 0, curve = curve)

  result <- convexity_correction_forward_rate(spec_zero, start = 1, end = 2)

  expect_equal(result$forward_rate, result$adjusted_forward)
  expect_true(all(abs(result$convexity_adjustment) < 1e-12))
  expect_equal(result$present_value, result$discount_factor * result$adjusted_forward)
})

test_that("convexity adjustment is positive for Hull-White inputs", {
  curve <- tibble::tibble(tenor = 0:8, discount_factor = exp(-0.015 * (0:8)))
  spec_hw <- short_rate_spec(
    model = "hull_white",
    volatility = 0.02,
    mean_reversion = 0.1,
    curve = curve
  )

  starts <- c(0.5, 1)
  ends <- c(1, 2)
  result <- convexity_correction_forward_rate(spec_hw, start = starts, end = ends)

  expect_true(all(result$convexity_adjustment > 0))
  expect_equal(result$present_value, result$discount_factor * result$adjusted_forward)
})

test_that("convexity correction validates start and end times", {
  curve <- tibble::tibble(tenor = 0:3, discount_factor = exp(-0.03 * (0:3)))
  spec <- short_rate_spec(model = "ho_lee", volatility = 0.01, curve = curve)

  expect_error(convexity_correction_forward_rate(spec, start = 1, end = 0.5))
})
