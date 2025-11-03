test_that("par swap exposure collapses under deterministic curve", {
  curve <- tibble::tibble(tenor = 0:5, discount_factor = exp(-0.03 * tenor))
  spec <- CompFinanceR::short_rate_spec(model = "ho_lee", volatility = 0.0, curve = curve)
  schedule <- CompFinanceR::swap_cashflow_schedule(start = 0, end = 3, frequency = 1, notional = 1e6)

  discount_fun <- spec$method$engine_state$discount_fun_vec
  df_start <- discount_fun(schedule[["start"]])
  df_end <- discount_fun(schedule[["end"]])
  float_leg <- sum(schedule[["notional"]] * (df_start - df_end))
  annuity <- sum(schedule[["notional"]] * schedule[["accrual_fraction"]] * discount_fun(schedule[["pay_time"]]))
  par_rate <- float_leg / annuity

  result <- CompFinanceR::simulate_swap_exposure(
    spec = spec,
    schedule = schedule,
    fixed_rate = par_rate,
    type = "payer",
    n_paths = 1,
    n_steps = 30,
    keep_paths = TRUE
  )

  time_zero_value <- result$paths |>
    dplyr::filter(time == 0) |>
    dplyr::pull(swap_value)
  expect_true(abs(time_zero_value) < 1e-8)
  expect_true(all(result$paths$positive_exposure < 1e-8))
  expect_true(all(result$summary$expected_exposure < 1e-8))
  expect_true(all(result$summary$expected_discounted_exposure < 1e-8))
})


test_that("deterministic exposure matches analytic swap PV at time zero", {
  curve <- tibble::tibble(tenor = 0:5, discount_factor = exp(-0.025 * tenor))
  spec <- CompFinanceR::short_rate_spec(model = "ho_lee", volatility = 0.0, curve = curve)
  schedule <- CompFinanceR::swap_cashflow_schedule(start = 0, end = 4, frequency = 1, notional = 2e6)

  discount_fun <- spec$method$engine_state$discount_fun_vec
  df_start <- discount_fun(schedule[["start"]])
  df_end <- discount_fun(schedule[["end"]])
  float_leg <- sum(schedule[["notional"]] * (df_start - df_end))
  annuity <- sum(schedule[["notional"]] * schedule[["accrual_fraction"]] * discount_fun(schedule[["pay_time"]]))
  par_rate <- float_leg / annuity
  rich_rate <- par_rate - 0.002

  analytic_pv <- CompFinanceR::price_swap(spec, schedule, fixed_rate = rich_rate, type = "payer") |>
    dplyr::filter(metric == "pv") |>
    dplyr::pull(value)

  result <- CompFinanceR::simulate_swap_exposure(
    spec = spec,
    schedule = schedule,
    fixed_rate = rich_rate,
    type = "payer",
    n_paths = 1,
    n_steps = 40,
    keep_paths = TRUE
  )

  time_zero <- dplyr::filter(result$paths, time == 0)
  expect_equal(time_zero$swap_value, analytic_pv, tolerance = 1e-8)
  expect_equal(time_zero$positive_exposure, max(analytic_pv, 0), tolerance = 1e-8)
  expect_equal(time_zero$discounted_positive_exposure, time_zero$positive_exposure, tolerance = 1e-8)

  pfe_col <- dplyr::select(result$summary, dplyr::starts_with("pfe_"))
  expect_equal(pfe_col[1, 1, drop = TRUE], max(analytic_pv, 0), tolerance = 1e-8)
})
