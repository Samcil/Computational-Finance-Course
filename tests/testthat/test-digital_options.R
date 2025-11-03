test_that("digital call closed form matches Black-Scholes formula", {
  spec <- black_scholes_spec(
    option_type = "call",
    strike = 100,
    maturity = 1,
    risk_free_rate = 0.03,
    dividend_yield = 0.01
  )
  prices <- price_digital_call(
    model_spec = spec,
    spot = c(95, 100, 105),
    volatility = c(0.18, 0.22)
  )

  expect_s3_class(prices, "tbl_df")
  expect_true(all(prices$price >= 0 & prices$price <= exp(-0.03)))
  expect_equal(nrow(prices), 6)

  expected <- exp(-0.03) * stats::pnorm(
    (
      log(prices$spot / spec$strike) +
        (spec$risk_free_rate - spec$dividend_yield - 0.5 * prices$volatility^2) *
          spec$maturity
    ) / (prices$volatility * sqrt(spec$maturity))
  )
  expect_equal(prices$price, expected, tolerance = 1e-8)
})

test_that("digital put plus call equals discounted payout", {
  spec <- black_scholes_spec(
    option_type = "call",
    strike = 100,
    maturity = 2,
    risk_free_rate = 0.025,
    dividend_yield = 0
  )
  call_price <- price_digital_call(spec, spot = 100, volatility = 0.2, payout = 2)$price
  put_price <- price_digital_put(spec, spot = 100, volatility = 0.2, payout = 2)$price
  expect_equal(call_price + put_price, 2 * exp(-0.025 * 2), tolerance = 1e-8)
})

test_that("digital cos method aligns with Black-Scholes analytic", {
  spot <- 100
  sigma <- 0.2
  r <- 0.03
  q <- 0
  maturity <- 1
  strikes <- c(90, 100, 110)

  cf_bs <- function(u) {
    drift <- log(spot) + (r - q - 0.5 * sigma^2) * maturity
    exp(1i * u * drift - 0.5 * sigma^2 * maturity * u^2)
  }

  cos_prices <- digital_cos_method(
    cf = cf_bs,
    option_type = "call",
    spot = spot,
    risk_free_rate = r,
    maturity = maturity,
    strikes = strikes,
    n_terms = 256,
    truncation = 8,
    grid_size = 4001
  )

  expected_prices <- exp(-r * maturity) * stats::pnorm(
    (
      log(spot / strikes) + (r - q - 0.5 * sigma^2) * maturity
    ) / (sigma * sqrt(maturity))
  )

  expect_equal(cos_prices$price, expected_prices, tolerance = 5e-3)
})
