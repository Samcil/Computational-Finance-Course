test_that("shifted lognormal price collapses to intrinsic when volatility vanishes", {
  forward <- -0.01
  strike <- 0.0
  shift <- 0.05
  price <- CompFinanceR::shifted_lognormal_price(
    option = "call",
    forward = forward,
    strike = strike,
    volatility = 0,
    maturity = 1,
    shift = shift,
    discount = 0.97
  )
  expect_equal(price, 0.97 * max(forward - strike, 0))
})


test_that("shifted implied volatility reproduces observed price", {
  forward <- -0.005
  strike <- 0
  shift <- 0.02
  vol <- 0.35
  maturity <- 0.5
  discount <- 0.99
  option_price <- CompFinanceR::shifted_lognormal_price(
    option = "put",
    forward = forward,
    strike = strike,
    volatility = vol,
    maturity = maturity,
    shift = shift,
    discount = discount
  )

  implied <- CompFinanceR::implied_volatility_shifted(
    option = "put",
    price = option_price,
    forward = forward,
    strike = strike,
    maturity = maturity,
    shift = shift,
    discount = discount
  )

  expect_equal(implied, vol, tolerance = 1e-6)
})
