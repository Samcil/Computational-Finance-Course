test_that("BSHW COS prices align with Black-76 equivalents", {
  curve <- tibble::tibble(
    tenor = seq(0, 20, by = 1),
    discount_factor = exp(-0.05 * tenor)
  )
  short_rate <- short_rate_spec(
    model = "hull_white",
    volatility = 0.05,
    mean_reversion = 0.1,
    curve = curve
  )
  spec <- bshw_spec(
    spot = 100,
    short_rate = short_rate,
    equity_vol = 0.2,
    correlation = 0.3
  )

  strikes <- seq(60, 140, length.out = 5)
  maturity <- 5

  cos_prices <- price_bshw_option_cos(
    spec,
    strikes = strikes,
    maturity = maturity,
    n_terms = 512,
    truncation = 10,
    integration_points = 1500
  )

  analytic_prices <- price_bshw_option_black76(
    spec,
    strikes = strikes,
    maturity = maturity,
    integration_points = 1500
  )

  expect_equal(cos_prices$call_price, analytic_prices$call_price, tolerance = 1e-4)
  expect_equal(cos_prices$put_price, analytic_prices$put_price, tolerance = 1e-4)
})

test_that("BSHW equivalent volatility remains positive", {
  vol <- bshw_equivalent_volatility(
    maturity = 7,
    equity_vol = 0.25,
    short_rate_vol = 0.03,
    correlation = -0.4,
    mean_reversion = 0.2,
    integration_points = 1200
  )

  expect_gt(vol, 0)
})
