test_that("short_rate_spec stores curve data", {
  curve <- tibble::tibble(
    tenor = c(0, 1, 2),
    discount_factor = exp(-0.02 * tenor)
  )

  spec <- CompFinanceR::short_rate_spec(
    model = "ho_lee",
    volatility = 0.01,
    curve = curve
  )

  expect_s3_class(spec, "short_rate_spec")
  state <- spec$method$engine_state
  expect_true(all(c("discount_fun", "theta_fun") %in% names(state)))
  expect_equal(state$discount_fun(1), curve$discount_factor[2], tolerance = 1e-10)
})


test_that("simulate_paths for Ho-Lee yields tidy output", {
  curve <- tibble::tibble(
    tenor = c(0, 1, 2),
    discount_factor = exp(-0.025 * tenor)
  )
  spec <- CompFinanceR::short_rate_spec(volatility = 0.01, curve = curve)

  paths <- CompFinanceR::simulate_paths(spec, n_paths = 5, n_steps = 10, maturity = 1, seed = 42)

  expect_true(all(c("short_rate", "discount_factor") %in% names(paths)))
  expect_equal(dplyr::n_distinct(paths$path_id), 5)
  expect_true(all(paths$time >= 0))
  at_time_zero <- paths |>
    dplyr::filter(abs(time) < 1e-12) |>
    dplyr::pull(discount_factor)

  expect_equal(unname(at_time_zero), rep(1, 5))
})


test_that("Monte Carlo discount matches initial curve on average", {
  curve <- tibble::tibble(
    tenor = seq(0, 5, by = 0.5),
    discount_factor = exp(-0.03 * seq(0, 5, by = 0.5))
  )

  spec <- CompFinanceR::short_rate_spec(volatility = 0.01, curve = curve)
  paths <- CompFinanceR::simulate_paths(spec, n_paths = 2000, n_steps = 50, maturity = 2, seed = 123)

  mc_estimate <- paths |>
    dplyr::filter(abs(time - 2) < 1e-8) |>
    dplyr::summarise(estimate = mean(discount_factor)) |>
    dplyr::pull(estimate)

  expected <- CompFinanceR::price_zcb(spec, 2)

  expect_equal(mc_estimate, expected, tolerance = 0.01)
})


test_that("Ho-Lee analytic pricing matches conditional Monte Carlo", {
  curve <- tibble::tibble(
    tenor = seq(0, 5, by = 0.5),
    discount_factor = exp(-0.025 * seq(0, 5, by = 0.5))
  )

  spec <- CompFinanceR::short_rate_spec(volatility = 0.01, curve = curve)
  maturity <- 2
  valuetime <- 1
  target <- 1.5

  paths <- CompFinanceR::simulate_paths(spec, n_paths = 2000, n_steps = 200, maturity = maturity, seed = 321)

  path_t <- paths |>
    dplyr::filter(abs(time - valuetime) < 1e-10) |>
    dplyr::arrange(path_id)

  path_T <- paths |>
    dplyr::filter(abs(time - target) < 1e-10) |>
    dplyr::arrange(path_id)

  ratio <- path_T$discount_factor / path_t$discount_factor

  analytic <- CompFinanceR::price_zcb(
    spec,
    maturities = target,
    valuation_time = valuetime,
    short_rate = path_t$short_rate
  )

  expect_equal(length(analytic), length(ratio))
  expect_equal(mean(analytic), mean(ratio), tolerance = 0.01)
  expect_true(all(analytic > 0))
})


test_that("Hull-White specification requires mean reversion", {
  curve <- tibble::tibble(
    tenor = c(0, 1, 2),
    discount_factor = exp(-0.02 * tenor)
  )

  expect_error(
    CompFinanceR::short_rate_spec(
      model = "hull_white",
      volatility = 0.01,
      curve = curve
    ),
    "mean_reversion"
  )
})


test_that("Hull-White analytic pricing matches conditional Monte Carlo", {
  curve <- tibble::tibble(
    tenor = seq(0, 5, by = 0.5),
    discount_factor = exp(-0.03 * seq(0, 5, by = 0.5))
  )

  spec <- CompFinanceR::short_rate_spec(
    model = "hull_white",
    volatility = 0.01,
    mean_reversion = 0.15,
    curve = curve
  )

  maturity <- 2
  valuetime <- 1
  target <- 1.75

  paths <- CompFinanceR::simulate_paths(
    spec,
    n_paths = 2500,
    n_steps = 200,
    maturity = maturity,
    seed = 456
  )

  path_t <- paths |>
    dplyr::filter(abs(time - valuetime) < 1e-10) |>
    dplyr::arrange(path_id)

  path_T <- paths |>
    dplyr::filter(abs(time - target) < 1e-10) |>
    dplyr::arrange(path_id)

  ratio <- path_T$discount_factor / path_t$discount_factor

  analytic <- CompFinanceR::price_zcb(
    spec,
    maturities = target,
    valuation_time = valuetime,
    short_rate = path_t$short_rate
  )

  expect_equal(length(analytic), length(ratio))
  expect_equal(mean(analytic), mean(ratio), tolerance = 0.015)
  expect_true(all(analytic > 0))
})
