test_that("forward start characteristic function equals 1 at u = 0", {
  cf_val <- forward_start_characteristic_function(
    u = 0,
    start = 0.5,
    maturity = 1,
    risk_free_rate = 0.01,
    dividend_yield = 0,
    mean_reversion = 1.5,
    long_term_variance = 0.04,
    vol_of_vol = 0.2,
    initial_variance = 0.04,
    correlation = -0.5
  )
  expect_equal(cf_val, 1 + 0i)
})

test_that("forward start COS agrees with Monte Carlo pricing", {
  spec <- heston_spec(
    initial_price = 90,
    initial_variance = 0.09,
    risk_free_rate = 0.015,
    dividend_yield = 0,
    mean_reversion = 1.2,
    long_term_variance = 0.09,
    vol_of_vol = 0.35,
    correlation = -0.4,
    scheme = "aes"
  )
  start <- 0.75
  maturity <- 1.5
  strike_adj <- 0.05

  cos_price <- price_forward_start_heston(
    process_spec = spec,
    start = start,
    maturity = maturity,
    strike_adjustment = strike_adj,
    option_type = "call",
    n_terms = 256,
    truncation = 10
  )$price

  set.seed(123)
  paths <- generate_heston_paths_aes(
    n_paths = 8000,
    n_steps = 360,
    maturity = maturity,
    initial_price = spec$initial_price,
    initial_variance = spec$initial_variance,
    risk_free_rate = spec$risk_free_rate,
    dividend_yield = spec$dividend_yield,
    mean_reversion = spec$mean_reversion,
    long_term_variance = spec$long_term_variance,
    vol_of_vol = spec$vol_of_vol,
    correlation = spec$correlation,
    seed = 42
  )

  payoff_tbl <- paths |>
    dplyr::group_by(.data$path_id) |>
    dplyr::summarise(
      s_start = stock_price[which.min(abs(time - start))],
      s_maturity = stock_price[which.min(abs(time - maturity))],
      .groups = "drop"
    )

  discounted_payoff <- exp(-spec$risk_free_rate * maturity) *
    pmax(payoff_tbl$s_maturity - (1 + strike_adj) * payoff_tbl$s_start, 0)

  mc_price <- mean(discounted_payoff)
  expect_equal(cos_price, mc_price, tolerance = 0.6)
})

test_that("forward start call prices decline with higher strike adjustments", {
  spec <- heston_spec(
    initial_price = 100,
    initial_variance = 0.05,
    risk_free_rate = 0.015,
    dividend_yield = 0.01,
    mean_reversion = 1.4,
    long_term_variance = 0.05,
    vol_of_vol = 0.3,
    correlation = -0.5
  )

  start <- 0.6
  maturity <- 1.4
  strike_adj <- c(0, 0.05, 0.1)

  prices <- price_forward_start_heston(
    process_spec = spec,
    start = start,
    maturity = maturity,
    strike_adjustment = strike_adj,
    option_type = "call",
    n_terms = 256,
    truncation = 10
  )$price

  expect_true(all(diff(prices) < 0))
})

test_that("forward start COS reproduces Python reference up to spot scaling", {
  python_style <- function(spec, start, maturity, strike_adj,
                           option_type = "call", n_terms = 256L,
                           truncation = 10) {
    tau <- maturity - start
    a <- -truncation * sqrt(tau)
    b <- truncation * sqrt(tau)
    k_indices <- seq_len(n_terms) - 1
    u <- k_indices * pi / (b - a)
    payoff_coef <- cos_coefficients(option_type, a, b, k_indices)
    cf_ratio <- function(u_eval) {
      forward_start_characteristic_function(
        u = u_eval,
        start = start,
        maturity = maturity,
        risk_free_rate = spec$risk_free_rate,
        dividend_yield = spec$dividend_yield,
        mean_reversion = spec$mean_reversion,
        long_term_variance = spec$long_term_variance,
        vol_of_vol = spec$vol_of_vol,
        initial_variance = spec$initial_variance,
        correlation = spec$correlation
      )
    }

    strike_multiplier <- 1 + strike_adj
    exponential_matrix <- exp(1i * outer(log(1 / strike_multiplier) - a, u))
    cf_values <- cf_ratio(u)
    cf_values[1] <- cf_values[1] * 0.5
    weights <- payoff_coef * cf_values
    discount_factor <- exp(-spec$risk_free_rate * maturity)
    discount_factor * strike_multiplier * as.numeric(Re(exponential_matrix %*% weights))
  }

  spec <- heston_spec(
    initial_price = 125,
    initial_variance = 0.08,
    risk_free_rate = 0.01,
    dividend_yield = 0.02,
    mean_reversion = 1.1,
    long_term_variance = 0.08,
    vol_of_vol = 0.4,
    correlation = -0.3
  )
  start <- 0.6
  maturity <- 1.4
  strike_adj <- c(-0.05, 0, 0.05)

  python_values <- python_style(
    spec = spec,
    start = start,
    maturity = maturity,
    strike_adj = strike_adj,
    option_type = "call",
    n_terms = 512,
    truncation = 12
  )
  r_prices <- price_forward_start_heston(
    process_spec = spec,
    start = start,
    maturity = maturity,
    strike_adjustment = strike_adj,
    option_type = "call",
    n_terms = 512,
    truncation = 12
  )$price

  scaling <- spec$initial_price * exp(-spec$dividend_yield * start)
  expect_equal(r_prices / scaling, python_values, tolerance = 1e-8)
})
