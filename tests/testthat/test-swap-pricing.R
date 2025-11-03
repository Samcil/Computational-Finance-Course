test_that("swap_cashflow_schedule constructs expected grid", {
  schedule <- CompFinanceR::swap_cashflow_schedule(start = 0, end = 2, frequency = 4, notional = 5)

  expect_s3_class(schedule, "tbl_df")
  expect_equal(nrow(schedule), 8)
  expect_equal(schedule$start[1], 0)
  expect_equal(schedule$end[8], 2)
  expect_true(all(schedule$accrual_fraction == 0.25))
  expect_true(all(schedule$notional == 5))
})

test_that("par swap has near-zero PV and correct DV01", {
  curve <- tibble::tibble(tenor = seq(0, 5), discount_factor = exp(-0.03 * tenor))
  spec <- CompFinanceR::short_rate_spec(volatility = 0.0, curve = curve)
  schedule <- CompFinanceR::swap_cashflow_schedule(start = 0, end = 5, frequency = 1, notional = 1e6)

  discount_fun <- spec$method$engine_state$discount_fun_vec
  df_pay <- discount_fun(schedule$pay_time)
  df_start <- discount_fun(schedule$start)
  df_end <- discount_fun(schedule$end)

  pv_float <- sum(schedule$notional * (df_start - df_end))
  annuity <- sum(schedule$notional * schedule$accrual_fraction * df_pay)
  par_rate <- pv_float / annuity

  result <- CompFinanceR::price_swap(spec, schedule, fixed_rate = par_rate, type = "payer")

  pv_value <- result |>
    dplyr::filter(metric == "pv") |>
    dplyr::pull(value)
  dv01_value <- result |>
    dplyr::filter(metric == "dv01") |>
    dplyr::pull(value)

  expect_equal(pv_value, 0, tolerance = 1e-6, scale = 1e6)
  expect_equal(dv01_value, annuity * 1e-4, tolerance = 1e-6)
})
