#' Generate a Swap Cashflow Schedule
#'
#' Construct a tidy tibble describing the payment schedule for a fixed-for-floating
#' interest rate swap. The helper assumes equal spacing between coupon dates and a
#' constant notional unless supplied otherwise.
#'
#' @param start Numeric. Start time of the swap in years.
#' @param end Numeric. Maturity time of the swap in years. Must be greater than
#'   `start`.
#' @param frequency Integer. Number of coupon payments per year (e.g., `2` for
#'   semi-annual, `4` for quarterly).
#' @param notional Numeric vector or scalar giving the notional outstanding over
#'   the life of the swap. A scalar is recycled across periods.
#' @param daycount Numeric vector or scalar representing the accrual fraction per
#'   period. Defaults to `1 / frequency`.
#'
#' @return A tibble with columns `period`, `start`, `end`, `pay_time`,
#'   `accrual_fraction`, and `notional`.
#'
#' @examples
#' schedule <- swap_cashflow_schedule(
#'   start = 0, end = 5, frequency = 2,
#'   notional = 1e6
#' )
#'
#' @export
swap_cashflow_schedule <- function(start,
                                   end,
                                   frequency = 2,
                                   notional = 1,
                                   daycount = NULL) {
  checkmate::assert_number(start, lower = 0)
  checkmate::assert_number(end, lower = start)
  checkmate::assert_int(frequency, lower = 1)
  step <- 1 / frequency
  periods_raw <- (end - start) / step
  if (!checkmate::test_number(periods_raw, lower = 1, finite = TRUE) ||
    abs(periods_raw - round(periods_raw)) > 1e-8) {
    rlang::abort("`end - start` must be a positive multiple of `1 / frequency`")
  }
  periods <- as.integer(round(periods_raw))

  period_starts <- start + step * (seq_len(periods) - 1)
  period_ends <- period_starts + step

  accrual_fraction <- if (is.null(daycount)) {
    rep(step, periods)
  } else {
    checkmate::assert_numeric(daycount, any.missing = FALSE, finite = TRUE)
    if (length(daycount) == 1L) {
      rep(daycount, periods)
    } else if (length(daycount) == periods) {
      as.numeric(daycount)
    } else {
      rlang::abort("`daycount` must have length 1 or match the number of periods")
    }
  }

  checkmate::assert_numeric(notional, any.missing = FALSE, finite = TRUE)
  notional_vec <- if (length(notional) == 1L) {
    rep(notional, periods)
  } else if (length(notional) == periods) {
    as.numeric(notional)
  } else {
    rlang::abort("`notional` must have length 1 or match the number of periods")
  }

  tibble::tibble(
    period = seq_len(periods),
    start = period_starts,
    end = period_ends,
    pay_time = period_ends,
    accrual_fraction = accrual_fraction,
    notional = notional_vec
  )
}


#' Price a Fixed-for-Floating Interest Rate Swap
#'
#' Compute the present value and DV01 of a plain vanilla interest rate swap
#' under a supplied short-rate specification using discount factors from the
#' calibrated curve. Positive PV indicates value to the receiver of the fixed
#' leg.
#'
#' @param spec A `short_rate_spec` created by [short_rate_spec()].
#' @param schedule Tibble produced by [swap_cashflow_schedule()] containing the
#'   swap payment dates.
#' @param fixed_rate Numeric. Fixed coupon rate expressed in decimal form.
#' @param type Character. Either "payer" (pay fixed, receive float) or
#'   "receiver" (receive fixed, pay float).
#'
#' @return A tibble with metrics `pv` and `dv01`.
#'
#' @examples
#' curve <- tibble::tibble(tenor = 0:5, discount_factor = exp(-0.03 * tenor))
#' spec <- short_rate_spec(volatility = 0.0, curve = curve)
#' schedule <- swap_cashflow_schedule(0, 5, frequency = 1, notional = 1e6)
#' price_swap(spec, schedule, fixed_rate = 0.03)
#'
#' @export
price_swap <- function(spec,
                       schedule,
                       fixed_rate,
                       type = c("payer", "receiver")) {
  type <- rlang::arg_match(type)
  checkmate::assert_class(spec, "short_rate_spec")
  schedule <- validate_swap_schedule(schedule)

  checkmate::assert_number(fixed_rate, finite = TRUE)

  state <- short_rate_state(spec)
  discount_fun <- state$discount_fun_vec
  df_pay <- discount_fun(schedule$pay_time)
  df_start <- discount_fun(schedule$start)
  df_end <- discount_fun(schedule$end)

  pv_fixed <- sum(schedule$notional * fixed_rate * schedule$accrual_fraction * df_pay)
  pv_float <- sum(schedule$notional * (df_start - df_end))

  signed_pv <- if (type == "payer") {
    pv_float - pv_fixed
  } else {
    pv_fixed - pv_float
  }

  dv01 <- sum(schedule$notional * schedule$accrual_fraction * df_pay) * 1e-4

  tibble::tibble(
    metric = c("pv", "dv01"),
    value = c(signed_pv, dv01)
  )
}
