test_that("caplet volatility calibration recovers known parameter", {
  curve <- tibble::tibble(tenor = 0:10, discount_factor = exp(-0.03 * tenor))
  true_spec <- CompFinanceR::short_rate_spec(
    model = "ho_lee",
    volatility = 0.015,
    curve = curve
  )

  caplet_grid <- tibble::tibble(
    reset = c(0.5, 1.0, 1.5, 2.0),
    payment = reset + 0.5,
    strike = 0.02
  )

  caplet_quotes <- caplet_grid |>
    dplyr::mutate(
      quote = purrr::pmap_dbl(
        list(reset, payment, strike),
        ~ CompFinanceR::price_caplet(
          spec = true_spec,
          reset = ..1,
          payment = ..2,
          strike = ..3
        )$value
      )
    )

  start_spec <- CompFinanceR::short_rate_spec(
    model = "ho_lee",
    volatility = 0.005,
    curve = curve
  )

  result <- CompFinanceR::calibrate_short_rate_volatility(
    spec = start_spec,
    caplet_quotes = caplet_quotes,
    lower = 1e-4,
    upper = 0.05
  )

  expect_equal(result$spec$args$volatility, true_spec$args$volatility, tolerance = 1e-4)
  expect_lt(max(abs(result$fitted$error)), 5e-4)
})


test_that("caplet calibration respects weights column", {
  curve <- tibble::tibble(tenor = 0:8, discount_factor = exp(-0.025 * tenor))
  true_spec <- CompFinanceR::short_rate_spec(
    model = "ho_lee",
    volatility = 0.012,
    curve = curve
  )

  caplet_quotes <- tibble::tibble(
    reset = c(0.5, 1, 1.5),
    payment = reset + 0.5,
    strike = 0.018,
    weight = c(1, 5, 1)
  ) |>
    dplyr::mutate(
      quote = purrr::pmap_dbl(
        list(reset, payment, strike),
        ~ CompFinanceR::price_caplet(
          spec = true_spec,
          reset = ..1,
          payment = ..2,
          strike = ..3
        )$value
      )
    )

  start_spec <- CompFinanceR::short_rate_spec(
    model = "ho_lee",
    volatility = 0.02,
    curve = curve
  )

  result <- CompFinanceR::calibrate_short_rate_volatility(
    spec = start_spec,
    caplet_quotes = caplet_quotes,
    lower = 1e-4,
    upper = 0.05
  )

  expect_equal(result$fitted$weight, caplet_quotes$weight)
  expect_lt(max(abs(result$fitted$error)), 1e-6)
})


test_that("swaption volatility calibration recovers known parameter", {
  curve <- tibble::tibble(tenor = 0:15, discount_factor = exp(-0.028 * tenor))
  true_spec <- CompFinanceR::short_rate_spec(
    model = "ho_lee",
    volatility = 0.018,
    curve = curve
  )

  schedule <- CompFinanceR::swap_cashflow_schedule(start = 1, end = 6, frequency = 2, notional = 1e6)
  fixed_rate <- 0.03

  swaption_quote <- CompFinanceR::price_swaption(
    spec = true_spec,
    schedule = schedule,
    fixed_rate = fixed_rate,
    type = "payer"
  ) |>
    dplyr::pull(value)

  swaption_quotes <- tibble::tibble(
    schedule = list(schedule),
    fixed_rate = fixed_rate,
    type = "payer",
    quote = swaption_quote
  )

  start_spec <- CompFinanceR::short_rate_spec(
    model = "ho_lee",
    volatility = 0.006,
    curve = curve
  )

  result <- CompFinanceR::calibrate_short_rate_swaption_volatility(
    spec = start_spec,
    swaption_quotes = swaption_quotes,
    lower = 1e-4,
    upper = 0.05
  )

  expect_equal(result$spec$args$volatility, true_spec$args$volatility, tolerance = 1e-4)
  expect_lt(max(abs(result$fitted$error)), 5e-4)
})


test_that("swaption calibration propagates weights", {
  curve <- tibble::tibble(tenor = 0:20, discount_factor = exp(-0.03 * tenor))
  true_spec <- CompFinanceR::short_rate_spec(
    model = "ho_lee",
    volatility = 0.02,
    curve = curve
  )

  schedule_short <- CompFinanceR::swap_cashflow_schedule(start = 1, end = 4, frequency = 2, notional = 5e5)
  schedule_long <- CompFinanceR::swap_cashflow_schedule(start = 2, end = 8, frequency = 1, notional = 1e6)

  fixed_rate_short <- 0.028
  fixed_rate_long <- 0.032

  quotes <- tibble::tibble(
    schedule = list(schedule_short, schedule_long),
    fixed_rate = c(fixed_rate_short, fixed_rate_long),
    type = c("payer", "receiver"),
    weight = c(3, 1)
  ) |>
    dplyr::mutate(
      quote = purrr::pmap_dbl(
        list(schedule, fixed_rate, type),
        ~ CompFinanceR::price_swaption(
          spec = true_spec,
          schedule = ..1,
          fixed_rate = ..2,
          type = ..3
        ) |>
          dplyr::pull(value)
      )
    )

  start_spec <- CompFinanceR::short_rate_spec(
    model = "ho_lee",
    volatility = 0.01,
    curve = curve
  )

  result <- CompFinanceR::calibrate_short_rate_swaption_volatility(
    spec = start_spec,
    swaption_quotes = quotes,
    lower = 1e-4,
    upper = 0.05
  )

  expect_equal(result$fitted$weight, quotes$weight)
  expect_length(result$fitted$schedule, nrow(quotes))
  expect_lt(max(abs(result$fitted$error)), 2e-3)
})


test_that("fit.short_rate_spec calibrates using caplet quotes", {
  curve <- tibble::tibble(tenor = 0:10, discount_factor = exp(-0.03 * tenor))
  true_spec <- CompFinanceR::short_rate_spec(
    model = "ho_lee",
    volatility = 0.014,
    curve = curve
  )

  caplet_quotes <- tibble::tibble(
    reset = c(0.5, 1, 1.5),
    payment = reset + 0.5,
    strike = 0.02
  ) |>
    dplyr::mutate(
      quote = purrr::pmap_dbl(
        list(reset, payment, strike),
        ~ CompFinanceR::price_caplet(
          spec = true_spec,
          reset = ..1,
          payment = ..2,
          strike = ..3
        )$value
      )
    )

  start_spec <- CompFinanceR::short_rate_spec(
    model = "ho_lee",
    volatility = 0.006,
    curve = curve
  )

  fitted <- CompFinanceR::fit(
    object = start_spec,
    data = caplet_quotes,
    method = "caplet_volatility"
  )

  expect_s3_class(fitted, "short_rate_fit")
  expect_equal(fitted$spec$args$volatility, true_spec$args$volatility, tolerance = 1e-4)

  augmented <- CompFinanceR::augment(fitted)
  expect_s3_class(augmented, "tbl_df")
  expect_equal(nrow(augmented), nrow(caplet_quotes))
  expect_lt(max(abs(augmented$error)), 5e-4)
})


test_that("fit.short_rate_spec calibrates using swaption quotes", {
  curve <- tibble::tibble(tenor = 0:12, discount_factor = exp(-0.0275 * tenor))
  true_spec <- CompFinanceR::short_rate_spec(
    model = "ho_lee",
    volatility = 0.017,
    curve = curve
  )

  schedule <- CompFinanceR::swap_cashflow_schedule(start = 1, end = 6, frequency = 2, notional = 1e6)
  swaption_quote <- CompFinanceR::price_swaption(
    spec = true_spec,
    schedule = schedule,
    fixed_rate = 0.031,
    type = "payer"
  ) |>
    dplyr::pull(value)

  quotes <- tibble::tibble(
    schedule = list(schedule),
    fixed_rate = 0.031,
    type = "payer",
    quote = swaption_quote
  )

  start_spec <- CompFinanceR::short_rate_spec(
    model = "ho_lee",
    volatility = 0.008,
    curve = curve
  )

  fitted <- CompFinanceR::fit(
    object = start_spec,
    data = quotes,
    method = "swaption_volatility"
  )

  expect_s3_class(fitted, "short_rate_fit")
  expect_equal(fitted$spec$args$volatility, true_spec$args$volatility, tolerance = 1e-4)

  augmented <- CompFinanceR::augment(fitted)
  expect_equal(nrow(augmented), nrow(quotes))
  expect_lt(max(abs(augmented$error)), 5e-4)
})
