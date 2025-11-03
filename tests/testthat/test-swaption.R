test_that("payer swaption collapses to positive swap PV at exercise when volatility is zero", {
  curve <- tibble::tibble(tenor = 0:5, discount_factor = exp(-0.03 * tenor))
  spec <- CompFinanceR::short_rate_spec(model = "ho_lee", volatility = 0.0, curve = curve)

  schedule <- CompFinanceR::swap_cashflow_schedule(start = 1, end = 4, frequency = 1, notional = 1e6)
  par_rate_components <- CompFinanceR::price_swap(spec, schedule, fixed_rate = 0, type = "payer")
  annuity <- CompFinanceR::price_swap(spec, schedule, fixed_rate = 0, type = "receiver") |>
    dplyr::filter(metric == "dv01") |>
    dplyr::pull(value) * 1e4

  par_rate <- par_rate_components |>
    dplyr::filter(metric == "pv") |>
    dplyr::pull(value) / annuity

  rich_rate <- par_rate - 0.002

  swap_value <- CompFinanceR::price_swap(spec, schedule, fixed_rate = rich_rate, type = "payer") |>
    dplyr::filter(metric == "pv") |>
    dplyr::pull(value)

  swaption_value <- CompFinanceR::price_swaption(spec, schedule, fixed_rate = rich_rate, type = "payer") |>
    dplyr::filter(metric == "price") |>
    dplyr::pull(value)

  expect_equal(swaption_value, max(swap_value, 0), tolerance = 1e-8)
})


test_that("receiver swaption intrinsic value recovered under zero volatility", {
  curve <- tibble::tibble(tenor = 0:6, discount_factor = exp(-0.025 * tenor))
  spec <- CompFinanceR::short_rate_spec(model = "ho_lee", volatility = 0.0, curve = curve)
  schedule <- CompFinanceR::swap_cashflow_schedule(start = 0.5, end = 3.5, frequency = 2, notional = 2e6)

  par_swap <- CompFinanceR::price_swap(spec, schedule, fixed_rate = 0.03, type = "receiver")
  swap_pv <- dplyr::filter(par_swap, metric == "pv") |> dplyr::pull(value)

  swaption_value <- CompFinanceR::price_swaption(spec, schedule, fixed_rate = 0.03, type = "receiver") |>
    dplyr::filter(metric == "price") |>
    dplyr::pull(value)

  expect_equal(swaption_value, max(swap_pv, 0), tolerance = 1e-8)
})
