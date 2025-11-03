#' Price Caplets Under Short-Rate Models
#'
#' Evaluate a single-period caplet using the analytic transformation to
#' zero-coupon bond options implied by Gaussian short-rate models. The function
#' assumes the supplied `short_rate_spec` exposes the calibrated curve and
#' volatility needed for analytic Hull-White or Ho-Lee pricing.
#'
#' @param spec A [`short_rate_spec`] object.
#' @param reset Numeric scalar. Reset time $T_1$ in years.
#' @param payment Numeric scalar. Payment time $T_2$ in years (must exceed
#'   `reset`).
#' @param strike Numeric scalar. Strike rate expressed in decimal form.
#' @param notional Numeric scalar. Notional of the caplet. Defaults to `1`.
#' @param accrual Numeric scalar. Accrual fraction for the period. Defaults to
#'   `payment - reset` when omitted.
#' @param method Character. Currently only `"analytic"` is supported.
#'
#' @return A tibble with columns `metric = "price"` and the corresponding
#'   numeric `value`.
#' @export
price_caplet <- function(spec,
                         reset,
                         payment,
                         strike,
                         notional = 1,
                         accrual = payment - reset,
                         method = c("analytic")) {
  method <- rlang::arg_match(method)
  value <- caplet_floorlet_price(spec, reset, payment, strike, notional, accrual, payoff = "caplet", method = method)

  tibble::tibble(metric = "price", value = value)
}


#' Price Floorlets Under Short-Rate Models
#'
#' Evaluate a single-period floorlet via its analytic bond-option
#' representation under Gaussian short-rate dynamics.
#'
#' @inheritParams price_caplet
#'
#' @return A tibble with columns `metric = "price"` and `value`.
#' @export
price_floorlet <- function(spec,
                           reset,
                           payment,
                           strike,
                           notional = 1,
                           accrual = payment - reset,
                           method = c("analytic")) {
  method <- rlang::arg_match(method)
  value <- caplet_floorlet_price(spec, reset, payment, strike, notional, accrual, payoff = "floorlet", method = method)

  tibble::tibble(metric = "price", value = value)
}


caplet_floorlet_price <- function(spec,
                                  reset,
                                  payment,
                                  strike,
                                  notional,
                                  accrual,
                                  payoff,
                                  method) {
  checkmate::assert_class(spec, "short_rate_spec")
  checkmate::assert_number(reset, lower = 0, finite = TRUE)
  checkmate::assert_number(payment, lower = 0, finite = TRUE)
  if (payment <= reset) {
    rlang::abort("`payment` must be greater than `reset`")
  }
  checkmate::assert_number(strike, finite = TRUE)
  checkmate::assert_number(notional, finite = TRUE)
  checkmate::assert_number(accrual, lower = 0, finite = TRUE)

  if (!identical(method, "analytic")) {
    rlang::abort("Only analytic pricing is currently implemented")
  }

  delta <- accrual
  k_factor <- 1 + delta * strike
  notional_adj <- notional * k_factor
  bond_strike <- 1 / k_factor
  option_type <- if (identical(payoff, "caplet")) "put" else "call"

  option_value <- bond_option_price_short_rate(
    spec = spec,
    option = option_type,
    strike = bond_strike,
    expiry = reset,
    maturity = payment
  )

  notional_adj * option_value
}
