test_that("mortgage annuity schedule fully amortises principal", {
  principal <- 200000
  schedule <- CompFinanceR::mortgage_annuity_schedule(
    principal = principal,
    coupon_rate = 0.025,
    maturity = 25,
    frequency = 12
  )

  expect_s3_class(schedule, "tbl_df")
  expect_equal(
    sum(schedule$principal_component + schedule$prepayment_component),
    principal,
    tolerance = 1e-6
  )
  expect_equal(schedule$outstanding_end[nrow(schedule)], 0, tolerance = 1e-8)
  expect_true(min(schedule$principal_component) >= -1e-8)
  expect_true(all(schedule$prepayment_component == 0))
})


test_that("mortgage valuation is flat at coupon curve", {
  principal <- 150000
  coupon_rate <- 0.03
  maturity <- 20
  frequency <- 12

  mortgage <- CompFinanceR::mortgage_annuity_spec(
    principal = principal,
    coupon_rate = coupon_rate,
    maturity = maturity,
    frequency = frequency
  )

  tenors <- seq(0, maturity, by = 1 / frequency)
  discount_factor <- (1 + coupon_rate / frequency)^(-frequency * tenors)
  curve <- tibble::tibble(
    tenor = tenors,
    discount_factor = discount_factor
  )
  discount_spec <- CompFinanceR::short_rate_spec(
    model = "ho_lee",
    volatility = 0,
    curve = curve
  )

  lender_result <- CompFinanceR::price_mortgage(mortgage, discount_spec, type = "lender")
  borrower_result <- CompFinanceR::price_mortgage(mortgage, discount_spec, type = "borrower")

  pv_lender <- lender_result |>
    dplyr::filter(metric == "pv") |>
    dplyr::pull(value)
  pv_borrower <- borrower_result |>
    dplyr::filter(metric == "pv") |>
    dplyr::pull(value)

  expect_equal(pv_lender, 0, tolerance = 1e-4, scale = principal)
  expect_equal(pv_borrower, -pv_lender, tolerance = 1e-8)
})


test_that("amortizing swap schedule aligns with mortgage and prices flat", {
  principal <- 500000
  coupon_rate <- 0.028
  maturity <- 15
  frequency <- 4

  mortgage <- CompFinanceR::mortgage_annuity_spec(
    principal = principal,
    coupon_rate = coupon_rate,
    maturity = maturity,
    frequency = frequency,
    prepayment = 0.05
  )
  swap_schedule <- CompFinanceR::amortizing_swap_schedule(mortgage)

  expect_equal(nrow(swap_schedule), frequency * maturity)
  expect_equal(
    sum(swap_schedule$principal_payment + swap_schedule$prepayment_payment),
    principal,
    tolerance = 1e-6
  )
  expect_equal(swap_schedule$notional[[1]], principal, tolerance = 1e-6)

  tenors <- seq(0, maturity, by = 1 / frequency)
  curve <- tibble::tibble(
    tenor = tenors,
    discount_factor = exp(-coupon_rate * tenors)
  )
  discount_spec <- CompFinanceR::short_rate_spec(
    model = "ho_lee",
    volatility = 0,
    curve = curve
  )

  discount_fun <- discount_spec$method$engine_state$discount_fun_vec
  df_pay <- discount_fun(swap_schedule$pay_time)
  df_start <- discount_fun(swap_schedule$start)
  df_end <- discount_fun(swap_schedule$end)
  pv_float <- sum(swap_schedule$notional * (df_start - df_end))
  annuity <- sum(swap_schedule$notional * swap_schedule$accrual_fraction * df_pay)
  par_rate <- pv_float / annuity

  swap_spec <- CompFinanceR::amortizing_swap_spec(mortgage)
  swap_price <- CompFinanceR::price_amortizing_swap(
    amortizing_swap = swap_spec,
    discount_spec = discount_spec,
    fixed_rate = par_rate
  )

  pv_value <- swap_price |>
    dplyr::filter(metric == "pv") |>
    dplyr::pull(value)

  expect_equal(pv_value, 0, tolerance = 1e-4, scale = principal)
})


test_that("mortgage schedule matches Python annuity output", {
  principal <- 1000000
  rate_per_period <- 0.05
  periods <- 30
  cpr <- 0.1

  schedule <- CompFinanceR::mortgage_annuity_schedule(
    principal = principal,
    coupon_rate = rate_per_period,
    maturity = periods,
    frequency = 1,
    prepayment = cpr
  )

  python_annuity <- function(rate, notional, periods, cpr_scalar) {
    M <- matrix(0, nrow = periods + 1, ncol = 6)
    M[, 1] <- 0:periods
    M[1, 2] <- notional
    for (t in seq_len(periods)) {
      remaining_periods <- periods - (t - 1)
      installment <- rate * M[t, 2] / (1 - 1 / (1 + rate)^remaining_periods)
      interest_payment <- rate * M[t, 2]
      principal_payment <- installment - interest_payment
      prepayment_payment <- cpr_scalar * (M[t, 2] - principal_payment)

      M[t + 1, 6] <- installment
      M[t + 1, 5] <- interest_payment
      M[t + 1, 4] <- principal_payment
      M[t + 1, 3] <- prepayment_payment
      M[t + 1, 2] <- M[t, 2] - principal_payment - prepayment_payment
    }
    M
  }

  ref <- python_annuity(rate_per_period, principal, periods, cpr)
  expect_equal(schedule$outstanding_end, ref[-1, 2], tolerance = 1e-6)
  expect_equal(schedule$interest_component, ref[-1, 5], tolerance = 1e-6)
  expect_equal(schedule$principal_component, ref[-1, 4], tolerance = 1e-6)
  expect_equal(schedule$prepayment_component, ref[-1, 3], tolerance = 1e-6)
})


test_that("bullet mortgage schedule matches lecture script output", {
  principal <- 750000
  rate_per_period <- 0.04
  periods <- 20
  cpr <- 0.02

  schedule <- CompFinanceR::mortgage_bullet_schedule(
    principal = principal,
    coupon_rate = rate_per_period,
    maturity = periods,
    frequency = 1,
    prepayment = cpr
  )

  bullet_reference <- function(rate, notional, periods, cpr_scalar) {
    M <- matrix(0, nrow = periods + 1, ncol = 6)
    M[, 1] <- 0:periods
    M[1, 2] <- notional
    if (periods >= 2) {
      for (t in 2:periods) {
        M[t, 5] <- rate * M[t - 1, 2]
        M[t, 4] <- 0
        scheduled_outstanding <- M[t - 1, 2] - M[t, 4]
        M[t, 3] <- scheduled_outstanding * cpr_scalar
        M[t, 2] <- scheduled_outstanding - M[t, 3]
        M[t, 6] <- M[t, 4] + M[t, 3] + M[t, 5]
      }
    }
    M[periods + 1, 5] <- rate * M[periods, 2]
    M[periods + 1, 4] <- M[periods, 2]
    M[periods + 1, 6] <- M[periods + 1, 4] + M[periods + 1, 3] + M[periods + 1, 5]
    M
  }

  ref <- bullet_reference(rate_per_period, principal, periods, cpr)

  expect_s3_class(schedule, "tbl_df")
  expect_equal(schedule$outstanding_start, ref[seq_len(periods), 2], tolerance = 1e-6)
  expect_equal(schedule$outstanding_end, ref[-1, 2], tolerance = 1e-6)
  expect_equal(schedule$prepayment_component, ref[-1, 3], tolerance = 1e-6)
  expect_equal(schedule$principal_component, ref[-1, 4], tolerance = 1e-6)
  expect_equal(schedule$interest_component, ref[-1, 5], tolerance = 1e-6)
  expect_equal(schedule$payment, ref[-1, 6], tolerance = 1e-6)
  expect_equal(schedule$outstanding_end[periods], 0, tolerance = 1e-10)
})


test_that("prepayment incentive reproduces logistic lecture curve", {
  spreads <- seq(-0.05, 0.08, length.out = 25)
  base <- 0.04
  amplitude <- 0.1
  slope <- 115
  threshold <- 0.02

  incentive <- CompFinanceR::mortgage_prepayment_incentive(
    rate_spread = spreads,
    base = base,
    amplitude = amplitude,
    slope = slope,
    threshold = threshold
  )

  reference <- base + amplitude / (1 + exp(slope * (threshold - spreads)))

  expect_equal(incentive, reference, tolerance = 1e-12, scale = 1)
  expect_true(all(diff(incentive) >= -1e-12))
  expect_equal(
    CompFinanceR::mortgage_prepayment_incentive(threshold, base, amplitude, slope, threshold),
    base + amplitude / 2,
    tolerance = 1e-12
  )
})
