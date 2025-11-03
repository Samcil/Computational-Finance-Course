test_that("implied_volatility_call recovers true volatility", {
  spot <- 100
  maturity <- 1
  risk_free_rate <- 0.02
  dividend_yield <- 0.01
  strikes <- c(90, 100, 110)
  true_vols <- c(0.18, 0.22, 0.25)

  prices <- purrr::map2_dbl(
    strikes,
    true_vols,
    ~ price_options(
      black_scholes_spec(
        option_type = "call",
        strike = .x,
        maturity = maturity,
        risk_free_rate = risk_free_rate,
        dividend_yield = dividend_yield
      ),
      spot = spot,
      volatility = .y
    )$price
  )

  option_tbl <- tibble::tibble(strike = strikes, price = prices)
  iv_result <- implied_volatility_call(
    option_tbl,
    spot = spot,
    maturity = maturity,
    risk_free_rate = risk_free_rate,
    dividend_yield = dividend_yield
  )

  expect_equal(iv_result$implied_volatility, true_vols, tolerance = 1e-4)
})

test_that("implied_volatility_put recovers true volatility", {
  spot <- 100
  maturity <- 0.5
  risk_free_rate <- 0.015
  dividend_yield <- 0.005
  strikes <- c(90, 100, 110)
  true_vols <- c(0.17, 0.21, 0.24)

  prices <- purrr::map2_dbl(
    strikes,
    true_vols,
    ~ price_options(
      black_scholes_spec(
        option_type = "put",
        strike = .x,
        maturity = maturity,
        risk_free_rate = risk_free_rate,
        dividend_yield = dividend_yield
      ),
      spot = spot,
      volatility = .y
    )$price
  )

  option_tbl <- tibble::tibble(strike = strikes, price = prices)
  iv_result <- implied_volatility_put(
    option_tbl,
    spot = spot,
    maturity = maturity,
    risk_free_rate = risk_free_rate,
    dividend_yield = dividend_yield
  )

  expect_equal(iv_result$implied_volatility, true_vols, tolerance = 1e-4)
})

test_that("implied volatility returns NA outside arbitrage bounds", {
  spot <- 100
  maturity <- 1
  risk_free_rate <- 0.02
  option_tbl <- tibble::tibble(
    strike = c(100, 100),
    price = c(0, 150) # below intrinsic and above upper bound
  )

  iv_result <- implied_volatility_call(
    option_tbl,
    spot = spot,
    maturity = maturity,
    risk_free_rate = risk_free_rate
  )

  expect_true(all(is.na(iv_result$implied_volatility)))
})

test_that("plot_volatility_smile produces ggplot object", {
  sample_data <- tibble::tibble(
    strike = c(90, 100, 110),
    implied_volatility = c(0.18, 0.2, 0.19),
    series = c("Jan", "Jan", "Jan")
  )

  smile_plot <- plot_volatility_smile(sample_data, colour_var = "series")
  expect_s3_class(smile_plot, "ggplot")
})
