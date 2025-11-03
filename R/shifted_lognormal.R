#' Shifted Lognormal Option Pricing
#'
#' Evaluate a Black-style option on a shifted underlying, useful for rate
#' instruments that can be negative while maintaining a lognormal diffusion on a
#' shifted process.
#'
#' @param option Character. Either `"call"` or `"put"`.
#' @param forward Numeric scalar. Forward level of the underlying.
#' @param strike Numeric scalar. Strike level of the payoff.
#' @param volatility Numeric scalar. Lognormal volatility (annualised).
#' @param maturity Numeric scalar. Time to maturity in years.
#' @param shift Numeric scalar. Additive shift applied to both forward and
#'   strike. Defaults to `0`.
#' @param discount Numeric scalar. Discount factor applied to the payoff.
#'
#' @return Option price as a numeric scalar.
#' @export
shifted_lognormal_price <- function(option,
                                    forward,
                                    strike,
                                    volatility,
                                    maturity,
                                    shift = 0,
                                    discount = 1) {
  option <- rlang::arg_match(option, c("call", "put"))
  checkmate::assert_number(forward, finite = TRUE)
  checkmate::assert_number(strike, finite = TRUE)
  checkmate::assert_number(volatility, lower = 0, finite = TRUE)
  checkmate::assert_number(maturity, lower = 0, finite = TRUE)
  checkmate::assert_number(shift, finite = TRUE)
  checkmate::assert_number(discount, lower = 0, finite = TRUE)

  adjusted_forward <- forward + shift
  adjusted_strike <- strike + shift

  if (adjusted_forward <= 0 || adjusted_strike <= 0) {
    rlang::abort("`forward + shift` and `strike + shift` must be positive")
  }

  if (volatility < 1e-10 || maturity < 1e-10) {
    intrinsic <- if (option == "call") {
      max(forward - strike, 0)
    } else {
      max(strike - forward, 0)
    }
    return(discount * intrinsic)
  }

  sigma_sqrt_t <- volatility * sqrt(maturity)
  d1 <- (log(adjusted_forward / adjusted_strike) + 0.5 * volatility^2 * maturity) / sigma_sqrt_t
  d2 <- d1 - sigma_sqrt_t

  call_value <- discount * (adjusted_forward * stats::pnorm(d1) - adjusted_strike * stats::pnorm(d2))

  if (option == "call") {
    call_value
  } else {
    call_value - discount * (adjusted_forward - adjusted_strike)
  }
}


#' Shifted Black Implied Volatility
#'
#' Recover the shifted lognormal volatility that reproduces a given option
#' premium.
#'
#' @inheritParams shifted_lognormal_price
#' @param price Numeric scalar. Observed option price.
#' @param tol Numeric scalar. Root-finding tolerance.
#'
#' @return Implied volatility as a numeric scalar.
#' @export
implied_volatility_shifted <- function(option,
                                       price,
                                       forward,
                                       strike,
                                       maturity,
                                       shift = 0,
                                       discount = 1,
                                       tol = 1e-8) {
  option <- rlang::arg_match(option, c("call", "put"))
  checkmate::assert_number(price, lower = 0, finite = TRUE)
  checkmate::assert_number(maturity, lower = 0, finite = TRUE)

  target <- function(vol) {
    shifted_lognormal_price(option, forward, strike, vol, maturity, shift, discount) - price
  }

  if (abs(target(0)) < tol) {
    return(0)
  }

  lower <- 1e-8
  upper <- 5
  value_lower <- target(lower)
  value_upper <- target(upper)
  while (value_lower * value_upper > 0 && upper < 20) {
    upper <- upper * 1.5
    value_upper <- target(upper)
  }

  if (value_lower * value_upper > 0) {
    rlang::abort("Failed to bracket implied volatility")
  }

  stats::uniroot(target, interval = c(lower, upper), tol = tol)$root
}
