test_that("term_structure_spec stores configuration", {
  spec <- CompFinanceR::term_structure_spec(curve_type = "ois", quote_type = "zero_rate")
  expect_s3_class(spec, "term_structure_spec")
  expect_equal(spec$args$curve_type, "ois")
  expect_equal(spec$args$quote_type, "zero_rate")
  expect_equal(spec$method$engine, "direct")
})


test_that("fit calibrates discount curve from zero rates", {
  spec <- CompFinanceR::term_structure_spec()
  quotes <- tibble::tibble(
    tenor = c(0.5, 1, 2),
    quote = c(0.02, 0.022, 0.025)
  )

  fit_obj <- CompFinanceR::fit(spec, quotes)
  calibration <- CompFinanceR::augment(fit_obj)

  expected_df <- as.numeric(exp(-quotes$quote * quotes$tenor))
  expect_equal(unname(calibration$discount_factor), expected_df, tolerance = 1e-10)
  expect_equal(calibration$zero_rate, quotes$quote, tolerance = 1e-10)
})


test_that("predict interpolates discount factors", {
  spec <- CompFinanceR::term_structure_spec()
  quotes <- tibble::tibble(
    tenor = c(1, 2),
    quote = c(0.02, 0.03)
  )
  fit_obj <- CompFinanceR::fit(spec, quotes)

  new_tenors <- tibble::tibble(tenor = c(1.5))
  preds <- stats::predict(fit_obj, new_data = new_tenors, type = "discount")

  expected_zero <- stats::approx(
    x = quotes$tenor,
    y = quotes$quote,
    xout = 1.5,
    rule = 2
  )$y
  expected_discount <- exp(-expected_zero * 1.5)

  expect_equal(preds$.pred, expected_discount, tolerance = 1e-8)
})


test_that("bootstrap engine solves par swap discount", {
  spec <- CompFinanceR::term_structure_spec(
    engine = "bootstrap",
    engine_options = list(fixed_leg_frequency = 2)
  )

  quotes <- tibble::tibble(
    tenor = c(0.5, 1.0, 1.5),
    quote = c(0.02, 0.021, 0.025),
    instrument = c("zero_rate", "zero_rate", "par_swap")
  )

  fit_obj <- CompFinanceR::fit(spec, quotes)
  calibration <- CompFinanceR::augment(fit_obj)

  expected_d_0_5 <- exp(-0.02 * 0.5)
  expected_d_1_0 <- exp(-0.021 * 1.0)
  step <- 0.5
  sum_prior <- step * (expected_d_0_5 + expected_d_1_0)
  expected_d_1_5 <- (1 - 0.025 * sum_prior) / (1 + 0.025 * step)

  expect_equal(calibration$discount_factor[calibration$tenor == 0.5], expected_d_0_5, tolerance = 1e-10)
  expect_equal(calibration$discount_factor[calibration$tenor == 1.0], expected_d_1_0, tolerance = 1e-10)
  expect_equal(calibration$discount_factor[calibration$tenor == 1.5], expected_d_1_5, tolerance = 1e-10)
})


test_that("newton engine matches bootstrap solution", {
  quotes <- tibble::tibble(
    tenor = c(0.5, 1.0, 1.5),
    quote = c(0.02, 0.021, 0.025),
    instrument = c("zero_rate", "zero_rate", "par_swap")
  )

  bootstrap_spec <- CompFinanceR::term_structure_spec(
    engine = "bootstrap",
    engine_options = list(fixed_leg_frequency = 2)
  )
  newton_spec <- CompFinanceR::term_structure_spec(
    engine = "newton",
    engine_options = list(fixed_leg_frequency = 2, max_iter = 300)
  )

  boot_fit <- CompFinanceR::fit(bootstrap_spec, quotes)
  newton_fit <- CompFinanceR::fit(newton_spec, quotes)

  boot_curve <- CompFinanceR::augment(boot_fit)
  newton_curve <- CompFinanceR::augment(newton_fit)

  expect_equal(newton_curve$discount_factor, boot_curve$discount_factor, tolerance = 1e-6)
})


test_that("multi-curve newton calibrates discount and projection curves", {
  ois_tenors <- c(0.5, 1.0, 1.5, 2.0)
  ois_zero <- c(0.02, 0.021, 0.022, 0.023)
  ois_discount <- setNames(exp(-ois_zero * ois_tenors), sprintf("%.1f", ois_tenors))

  libor_tenors <- c(0.5, 1.0, 1.5, 2.0)
  libor_zero <- c(0.025, 0.027, 0.028, 0.029)
  libor_discount <- setNames(exp(-libor_zero * libor_tenors), sprintf("%.1f", libor_tenors))

  discount_df_lookup <- function(time) {
    key <- sprintf("%.1f", time)
    ois_discount[[key]]
  }

  projection_df_lookup <- function(time) {
    if (abs(time) < 1e-12) {
      return(1)
    }
    key <- sprintf("%.1f", time)
    libor_discount[[key]]
  }

  compute_par_rate <- function(maturity) {
    payment_times <- seq(0.5, maturity, by = 0.5)
    accrual <- diff(c(0, payment_times))

    discount_leg <- 0
    floating_leg <- 0
    prev_time <- 0

    for (j in seq_along(payment_times)) {
      t_j <- payment_times[j]
      alpha_j <- accrual[j]
      discount_df <- discount_df_lookup(t_j)
      discount_leg <- discount_leg + alpha_j * discount_df

      prev_proj <- if (abs(prev_time) < 1e-12) 1 else projection_df_lookup(prev_time)
      curr_proj <- projection_df_lookup(t_j)
      forward_rate <- (prev_proj / curr_proj - 1) / alpha_j
      floating_leg <- floating_leg + alpha_j * discount_df * forward_rate

      prev_time <- t_j
    }

    floating_leg / discount_leg
  }

  par_swaps <- c(
    compute_par_rate(1.0),
    compute_par_rate(1.5),
    compute_par_rate(2.0)
  )

  market_quotes <- tibble::tribble(
    ~curve, ~tenor, ~instrument, ~quote,
    "ois", 0.5, "zero_rate", 0.02,
    "ois", 1.0, "zero_rate", 0.021,
    "ois", 1.5, "zero_rate", 0.022,
    "ois", 2.0, "zero_rate", 0.023,
    "libor", 0.5, "zero_rate", 0.025,
    "libor", 1.0, "par_swap", par_swaps[1],
    "libor", 1.5, "par_swap", par_swaps[2],
    "libor", 2.0, "par_swap", par_swaps[3]
  )

  spec <- CompFinanceR::term_structure_spec(
    curve_type = "ois",
    engine = "multi_curve_newton",
    engine_options = list(
      instrument_col = "instrument",
      curve_col = "curve",
      discount_curve_id = "ois",
      fixed_leg_frequency = 2,
      max_iter = 400
    )
  )

  fit_obj <- CompFinanceR::fit(spec, market_quotes)
  calibration <- CompFinanceR::augment(fit_obj)

  expect_setequal(unique(calibration$curve), c("ois", "libor"))

  ois_curve <- calibration[calibration$curve == "ois", , drop = FALSE]
  libor_curve <- calibration[calibration$curve == "libor", , drop = FALSE]

  expect_equal(
    ois_curve$discount_factor,
    unname(ois_discount[sprintf("%.1f", ois_curve$tenor)]),
    tolerance = 5e-9
  )

  expect_equal(
    libor_curve$discount_factor,
    unname(libor_discount[sprintf("%.1f", libor_curve$tenor)]),
    tolerance = 1e-6
  )

  zero_pred <- stats::predict(
    fit_obj,
    new_data = tibble::tibble(curve = "libor", tenor = 1.25),
    type = "zero"
  )

  expect_equal(zero_pred$curve, "libor")
  expect_length(zero_pred$.pred, 1)
  expect_true(is.numeric(zero_pred$.pred))
})


test_that("price_swap_term_structure supports multi-curve valuation", {
  ois_tenors <- c(0.5, 1.0, 1.5, 2.0)
  ois_zero <- c(0.02, 0.021, 0.022, 0.023)
  ois_discount <- setNames(exp(-ois_zero * ois_tenors), sprintf("%.1f", ois_tenors))

  libor_tenors <- c(0.5, 1.0, 1.5, 2.0)
  libor_zero <- c(0.025, 0.027, 0.028, 0.029)
  libor_discount <- setNames(exp(-libor_zero * libor_tenors), sprintf("%.1f", libor_tenors))

  discount_df_lookup <- function(time) {
    key <- sprintf("%.1f", time)
    ois_discount[[key]]
  }

  projection_df_lookup <- function(time) {
    if (abs(time) < 1e-12) {
      return(1)
    }
    key <- sprintf("%.1f", time)
    libor_discount[[key]]
  }

  compute_par_rate <- function(maturity) {
    payment_times <- seq(0.5, maturity, by = 0.5)
    accrual <- diff(c(0, payment_times))

    discount_leg <- 0
    floating_leg <- 0
    prev_time <- 0

    for (j in seq_along(payment_times)) {
      t_j <- payment_times[j]
      alpha_j <- accrual[j]
      discount_df <- discount_df_lookup(t_j)
      discount_leg <- discount_leg + alpha_j * discount_df

      prev_proj <- if (abs(prev_time) < 1e-12) 1 else projection_df_lookup(prev_time)
      curr_proj <- projection_df_lookup(t_j)
      forward_rate <- (prev_proj / curr_proj - 1) / alpha_j
      floating_leg <- floating_leg + alpha_j * discount_df * forward_rate

      prev_time <- t_j
    }

    floating_leg / discount_leg
  }

  par_swaps <- c(
    compute_par_rate(1.0),
    compute_par_rate(1.5),
    compute_par_rate(2.0)
  )

  market_quotes <- tibble::tribble(
    ~curve, ~tenor, ~instrument, ~quote,
    "ois", 0.5, "zero_rate", 0.02,
    "ois", 1.0, "zero_rate", 0.021,
    "ois", 1.5, "zero_rate", 0.022,
    "ois", 2.0, "zero_rate", 0.023,
    "libor", 0.5, "zero_rate", 0.025,
    "libor", 1.0, "par_swap", par_swaps[1],
    "libor", 1.5, "par_swap", par_swaps[2],
    "libor", 2.0, "par_swap", par_swaps[3]
  )

  spec <- CompFinanceR::term_structure_spec(
    curve_type = "ois",
    engine = "multi_curve_newton",
    engine_options = list(
      instrument_col = "instrument",
      curve_col = "curve",
      discount_curve_id = "ois",
      fixed_leg_frequency = 2,
      max_iter = 400
    )
  )

  fit_obj <- CompFinanceR::fit(spec, market_quotes)
  schedule <- CompFinanceR::swap_cashflow_schedule(0, 2, frequency = 2, notional = 1)

  par_rate <- par_swaps[3]

  multi_curve_price <- CompFinanceR::price_swap_term_structure(
    fit = fit_obj,
    schedule = schedule,
    fixed_rate = par_rate,
    type = "payer",
    discount_curve = "ois",
    projection_curve = "libor"
  )

  expect_lt(
    abs(multi_curve_price$value[multi_curve_price$metric == "pv"]),
    1e-6
  )

  single_curve_price <- CompFinanceR::price_swap_term_structure(
    fit = fit_obj,
    schedule = schedule,
    fixed_rate = par_rate,
    type = "payer",
    discount_curve = "ois",
    projection_curve = "ois"
  )

  expect_true(
    abs(single_curve_price$value[single_curve_price$metric == "pv"]) > 1e-6
  )
})


test_that("treasury bootstrap recovers coupon bond discounts", {
  true_tenors <- c(0.5, 1.0, 1.5, 2.0)
  true_zero <- c(0.02, 0.021, 0.022, 0.023)
  true_discount <- exp(-true_zero * true_tenors)
  coupon_rate <- 0.03
  frequency <- 2
  coupon_amount <- coupon_rate / frequency

  price_per_unit <- coupon_amount * sum(true_discount[1:3]) + (coupon_amount + 1) * true_discount[4]

  quotes <- tibble::tribble(
    ~tenor, ~quote, ~instrument, ~coupon_rate,
    0.5, true_zero[1], "zero_rate", 0,
    1.0, true_zero[2], "zero_rate", 0,
    1.5, true_zero[3], "zero_rate", 0,
    2.0, price_per_unit, "treasury_price", coupon_rate
  )

  spec <- CompFinanceR::term_structure_spec(
    engine = "treasury_bootstrap",
    engine_options = list(
      instrument_col = "instrument",
      coupon_frequency = frequency,
      face_value = 1
    )
  )

  fit_obj <- CompFinanceR::fit(spec, quotes)
  calibration <- CompFinanceR::augment(fit_obj)

  calibration <- calibration[order(calibration$tenor), , drop = FALSE]

  expect_equal(calibration$discount_factor, true_discount, tolerance = 1e-9)
})


test_that("term_structure_discount_function reproduces calibrated discounts", {
  spec <- CompFinanceR::term_structure_spec()
  quotes <- tibble::tibble(
    tenor = c(0.5, 1, 2, 3),
    quote = c(0.02, 0.021, 0.022, 0.023)
  )

  fit_obj <- CompFinanceR::fit(spec, quotes)
  discount_fun <- CompFinanceR::term_structure_discount_function(fit_obj)

  calibration <- CompFinanceR::augment(fit_obj)
  expect_equal(
    unname(discount_fun(calibration$tenor)),
    unname(calibration$discount_factor),
    tolerance = 1e-10
  )
})


test_that("price_swap_term_structure matches short-rate pricing", {
  spec <- CompFinanceR::term_structure_spec()
  quotes <- tibble::tibble(
    tenor = c(0.5, 1, 1.5, 2),
    quote = c(0.02, 0.022, 0.023, 0.024)
  )

  fit_obj <- CompFinanceR::fit(spec, quotes)
  schedule <- CompFinanceR::swap_cashflow_schedule(0, 2, frequency = 2, notional = 1e6)

  term_structure_price <- CompFinanceR::price_swap_term_structure(
    fit = fit_obj,
    schedule = schedule,
    fixed_rate = 0.023,
    type = "payer"
  )

  sr_spec <- CompFinanceR::short_rate_spec(
    model = "ho_lee",
    volatility = 0,
    curve = fit_obj
  )

  short_rate_price <- CompFinanceR::price_swap(
    spec = sr_spec,
    schedule = schedule,
    fixed_rate = 0.023,
    type = "payer"
  )

  expect_equal(
    term_structure_price$value[term_structure_price$metric == "pv"],
    short_rate_price$value[short_rate_price$metric == "pv"],
    tolerance = 1e-10
  )

  expect_equal(
    term_structure_price$value[term_structure_price$metric == "dv01"],
    short_rate_price$value[short_rate_price$metric == "dv01"],
    tolerance = 1e-10
  )
})


test_that("term_structure_quote_sensitivities computes bump-and-revalue deltas", {
  spec <- CompFinanceR::term_structure_spec(
    engine = "bootstrap",
    engine_options = list(fixed_leg_frequency = 2)
  )

  quotes <- tibble::tribble(
    ~tenor, ~quote, ~instrument,
    0.5, 0.02, "zero_rate",
    1.0, 0.021, "zero_rate",
    1.5, 0.024, "par_swap",
    2.0, 0.025, "par_swap"
  )

  valuation_fun <- function(fit) {
    schedule <- CompFinanceR::swap_cashflow_schedule(0, 2, frequency = 2, notional = 1)
    result <- CompFinanceR::price_swap_term_structure(
      fit = fit,
      schedule = schedule,
      fixed_rate = 0.0235,
      type = "payer"
    )
    result$value[result$metric == "pv"]
  }

  bump <- 1e-5
  sensitivity_tbl <- CompFinanceR::term_structure_quote_sensitivities(
    spec = spec,
    market_data = quotes,
    valuation_fun = valuation_fun,
    bump_size = bump,
    bump_type = "additive"
  )

  expect_s3_class(sensitivity_tbl, "tbl_df")
  expect_equal(nrow(sensitivity_tbl), nrow(quotes))

  base_fit <- CompFinanceR::fit(spec, quotes)
  base_value <- valuation_fun(base_fit)
  expect_true(all(abs(sensitivity_tbl$base_value - base_value) < 1e-12))

  expected <- purrr::map_dbl(seq_len(nrow(quotes)), function(idx) {
    bumped_quotes <- quotes
    bumped_quotes$quote[idx] <- bumped_quotes$quote[idx] + bump
    bumped_fit <- CompFinanceR::fit(spec, bumped_quotes)
    bumped_value <- valuation_fun(bumped_fit)
    (bumped_value - base_value) / bump
  })

  expect_equal(sensitivity_tbl$sensitivity, expected, tolerance = 1e-8)
})
