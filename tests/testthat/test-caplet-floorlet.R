test_that("caplet degenerates to intrinsic value under zero volatility", {
  curve <- tibble::tibble(tenor = 0:5, discount_factor = exp(-0.03 * tenor))
  spec <- CompFinanceR::short_rate_spec(model = "ho_lee", volatility = 0.0, curve = curve)

  state <- spec$method$engine_state
  reset <- 1
  payment <- 1.5
  accrual <- payment - reset
  forward <- (state$discount_fun(reset) / state$discount_fun(payment) - 1) / accrual
  strike <- forward - 0.01

  caplet_price <- CompFinanceR::price_caplet(
    spec = spec,
    reset = reset,
    payment = payment,
    strike = strike,
    accrual = accrual,
    notional = 1
  ) |>
    dplyr::filter(metric == "price") |>
    dplyr::pull(value)

  expected <- state$discount_fun(payment) * accrual * max(forward - strike, 0)
  expect_equal(caplet_price, expected, tolerance = 1e-10)
})


test_that("floorlet intrinsic value recovered when volatility is zero", {
  curve <- tibble::tibble(tenor = 0:5, discount_factor = exp(-0.02 * tenor))
  spec <- CompFinanceR::short_rate_spec(model = "ho_lee", volatility = 0.0, curve = curve)

  state <- spec$method$engine_state
  reset <- 2
  payment <- 2.5
  accrual <- payment - reset
  forward <- (state$discount_fun(reset) / state$discount_fun(payment) - 1) / accrual
  strike <- forward + 0.005

  floorlet_price <- CompFinanceR::price_floorlet(
    spec = spec,
    reset = reset,
    payment = payment,
    strike = strike,
    accrual = accrual,
    notional = 5
  ) |>
    dplyr::filter(metric == "price") |>
    dplyr::pull(value)

  expected <- 5 * state$discount_fun(payment) * accrual * max(strike - forward, 0)
  expect_equal(floorlet_price, expected, tolerance = 1e-10)
})
