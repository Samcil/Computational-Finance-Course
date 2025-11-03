portfolio_trade <- function(fixed_rate, type, netting_set) {
  schedule <- swap_cashflow_schedule(start = 0, end = 2, frequency = 1, notional = 1)
  list(
    schedule = schedule,
    fixed_rate = fixed_rate,
    type = type,
    trade_id = paste(type, fixed_rate, sep = "_"),
    netting_set = netting_set
  )
}

test_that("simulate_exposure aggregates trades into netting sets", {
  curve <- tibble::tibble(tenor = 0:4, discount_factor = exp(-0.015 * (0:4)))
  spec <- short_rate_spec(model = "ho_lee", volatility = 0, curve = curve)

  trades <- list(
    portfolio_trade(0.03, "receiver", "netting"),
    portfolio_trade(0.01, "payer", "netting")
  )

  exposure <- exposure_spec(trades = trades, quantiles = c(0.5, 0.95), discount = TRUE)
  result <- simulate_exposure(
    spec = exposure,
    process_spec = spec,
    n_paths = 1,
    n_steps = 10,
    maturity = 2,
    seed = 123,
    keep_paths = TRUE
  )

  expect_s3_class(result, "portfolio_exposure_result")

  trade_paths <- dplyr::arrange(result$paths$trades, .data$time, .data$trade_id)
  netting_paths <- dplyr::arrange(result$paths$netting, .data$time)

  manual_netting <- trade_paths |>
    dplyr::group_by(.data$path_id, .data$time, .data$netting_set) |>
    dplyr::summarise(
      discount_factor = dplyr::first(.data$discount_factor),
      value = sum(.data$value),
      positive_exposure = max(value, 0),
      discounted_positive_exposure = positive_exposure * discount_factor,
      .groups = "drop"
    ) |>
    dplyr::arrange(.data$time)

  expect_equal(netting_paths, manual_netting)

  portfolio_paths <- dplyr::arrange(result$paths$portfolio, .data$time)
  expect_equal(
    portfolio_paths$positive_exposure,
    pmax(manual_netting$value, 0)
  )
  expect_equal(portfolio_paths$discounted_positive_exposure, portfolio_paths$positive_exposure * portfolio_paths$discount_factor)
})

test_that("exposure summaries honour discount flag", {
  curve <- tibble::tibble(tenor = 0:3, discount_factor = exp(-0.02 * (0:3)))
  spec <- short_rate_spec(model = "ho_lee", volatility = 0, curve = curve)
  trades <- list(portfolio_trade(0.025, "receiver", "solo"))

  exposure <- exposure_spec(trades = trades, quantiles = 0.9, discount = FALSE)
  result <- simulate_exposure(
    spec = exposure,
    process_spec = spec,
    n_paths = 1,
    n_steps = 8,
    maturity = 2,
    seed = 42,
    keep_paths = FALSE
  )

  expect_false("expected_discounted_exposure" %in% names(result$portfolio))
  expect_false("expected_discounted_exposure" %in% names(result$trades))
})

test_that("Monte Carlo VaR collapses with deterministic exposure", {
  curve <- tibble::tibble(tenor = 0:4, discount_factor = exp(-0.015 * (0:4)))
  spec <- short_rate_spec(model = "ho_lee", volatility = 0, curve = curve)
  trades <- list(portfolio_trade(0.02, "receiver", "set"))

  risk <- risk_measure_spec(trades = trades, horizon = 2, alpha = 0.05, method = "monte_carlo")
  measures <- simulate_risk_measures(
    spec = risk,
    process_spec = spec,
    n_paths = 2,
    n_steps = 10,
    maturity = 2,
    seed = 900
  )

  expect_true(all(abs(measures$value) < 1e-10))
})

test_that("historical VaR matches shocked curve valuations", {
  curve_rate <- function(r) tibble::tibble(tenor = 0:4, discount_factor = exp(-r * (0:4)))
  base_curve <- curve_rate(0.01)
  up_curve <- curve_rate(0.015)
  down_curve <- curve_rate(0.005)

  trades <- list(
    portfolio_trade(0.03, "receiver", "set"),
    portfolio_trade(0.01, "payer", "set")
  )

  portfolio_value <- function(rate) {
    discount <- function(t) exp(-rate * t)

    trade_value <- function(trade) {
      schedule <- trade$schedule
      df_pay <- discount(schedule$pay_time)
      df_start <- discount(schedule$start)
      df_end <- discount(schedule$end)
      float_leg <- sum(schedule$notional * (df_start - df_end))
      fixed_leg <- sum(schedule$notional * trade$fixed_rate * schedule$accrual_fraction * df_pay)
      if (trade$type == "payer") {
        float_leg - fixed_leg
      } else {
        fixed_leg - float_leg
      }
    }

    sum(vapply(trades, trade_value, numeric(1)))
  }

  rates <- c(0.01, 0.015, 0.005)
  scenario_values <- vapply(rates, portfolio_value, numeric(1))
  baseline <- scenario_values[1]
  losses <- baseline - scenario_values

  expected_var <- stats::quantile(losses, probs = 0.05, names = FALSE)
  expected_es <- mean(losses[losses >= expected_var])

  risk <- risk_measure_spec(trades = trades, horizon = 2, alpha = 0.05, method = "historical")
  measures <- simulate_risk_measures(
    spec = risk,
    curve_scenarios = list(base_curve, up_curve, down_curve),
    base_curve = base_curve
  )

  var_value <- measures |>
    dplyr::filter(.data$metric == "VaR") |>
    dplyr::pull(.data$value)
  es_value <- measures |>
    dplyr::filter(.data$metric == "ES") |>
    dplyr::pull(.data$value)

  expect_equal(var_value, expected_var)
  expect_equal(es_value, expected_es)
})
