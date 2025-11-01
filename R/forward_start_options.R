#' Heston Forward Start Characteristic Function
#'
#' Computes the characteristic function of the forward start log-return
#' \eqn{\log(S(T_2) / S(T_1))} under the Heston stochastic volatility model.
#' The formulation follows Grzelak (2019) and integrates out the stochastic
#' variance \eqn{V(T_1)} using its non-central chi-squared distribution.
#'
#' @param u Numeric vector of real arguments where the characteristic function
#'   is evaluated.
#' @param start Numeric. Forward-start time \eqn{T_1} (years) at which the strike
#'   is set. Must satisfy `0 <= start < maturity`.
#' @param maturity Numeric. Option maturity \eqn{T_2} (years).
#' @param risk_free_rate Numeric. Continuously compounded risk-free rate \eqn{r}.
#' @param dividend_yield Numeric. Continuous dividend yield \eqn{q}. Default is 0.
#' @param mean_reversion Numeric. Mean reversion speed \eqn{\kappa > 0}.
#' @param long_term_variance Numeric. Long-run variance level \eqn{\theta > 0}.
#' @param vol_of_vol Numeric. Volatility of variance \eqn{\xi > 0}.
#' @param initial_variance Numeric. Initial variance \eqn{v_0 \ge 0}.
#' @param correlation Numeric. Correlation \eqn{\rho \in [-1, 1]} between price
#'   and variance Brownian motions.
#'
#' @return Complex vector `cf(u)`.
#' @export
forward_start_characteristic_function <- function(u,
                                                  start,
                                                  maturity,
                                                  risk_free_rate,
                                                  dividend_yield = 0,
                                                  mean_reversion,
                                                  long_term_variance,
                                                  vol_of_vol,
                                                  initial_variance,
                                                  correlation) {
  checkmate::assert_numeric(u, any.missing = FALSE)
  checkmate::assert_number(start, lower = 0, finite = TRUE)
  checkmate::assert_number(maturity, lower = 0, finite = TRUE)
  if (start >= maturity) {
    rlang::abort("`start` must be strictly less than `maturity`")
  }
  checkmate::assert_number(risk_free_rate, finite = TRUE)
  checkmate::assert_number(dividend_yield, finite = TRUE)
  checkmate::assert_number(mean_reversion, lower = .Machine$double.eps, finite = TRUE)
  checkmate::assert_number(long_term_variance, lower = 0, finite = TRUE)
  checkmate::assert_number(vol_of_vol, lower = .Machine$double.eps, finite = TRUE)
  checkmate::assert_number(initial_variance, lower = 0, finite = TRUE)
  checkmate::assert_number(correlation, lower = -1, upper = 1, finite = TRUE)

  tau <- maturity - start
  iu <- 1i * u
  vol_sq <- vol_of_vol^2

  d <- sqrt((mean_reversion - correlation * vol_of_vol * iu)^2 +
      vol_sq * (u^2 + iu))
  g <- (mean_reversion - correlation * vol_of_vol * iu - d) /
    (mean_reversion - correlation * vol_of_vol * iu + d)
  exp_neg_d_tau <- exp(-d * tau)

  C <- (mean_reversion - correlation * vol_of_vol * iu - d) *
    ((1 - exp_neg_d_tau) / (vol_sq * (1 - g * exp_neg_d_tau)))
  A <- (risk_free_rate - dividend_yield) * iu * tau +
    (mean_reversion * long_term_variance / vol_sq) * (
      (mean_reversion - correlation * vol_of_vol * iu - d) * tau -
        2 * log((1 - g * exp_neg_d_tau) / (1 - g))
    )

  if (start <= .Machine$double.eps) {
    c_bar <- 0
    cb_kappa_bar <- initial_variance
  } else {
    exp_neg_kappa_start <- exp(-mean_reversion * start)
    c_bar <- vol_sq * (1 - exp_neg_kappa_start) / (4 * mean_reversion)
    kappa_bar <- (4 * mean_reversion * initial_variance * exp_neg_kappa_start) /
      (vol_sq * (1 - exp_neg_kappa_start))
    cb_kappa_bar <- c_bar * kappa_bar
  }

  delta <- 4 * mean_reversion * long_term_variance / vol_sq
  denom <- 1 - 2 * C * c_bar

  exp(A + C * cb_kappa_bar / denom) * denom^(-0.5 * delta)
}

#' Price Forward Start Options under Heston
#'
#' Prices forward start options where the strike is determined at time `start`
#' as \eqn{K = (1 + k) S(T_1)} and the payoff at `maturity` is evaluated under
#' Heston dynamics. Pricing is performed via the COS method applied to the
#' log-return \eqn{\log(S(T_2) / S(T_1))} using the characteristic function
#' produced by [forward_start_characteristic_function()].
#'
#' @param process_spec A `heston_spec` object describing the model parameters.
#' @param start Numeric. Forward-start time \eqn{T_1} in years.
#' @param maturity Numeric. Final maturity \eqn{T_2} in years.
#' @param strike_adjustment Numeric vector of strike adjustments \eqn{k}. The
#'   effective strike multiple is `1 + strike_adjustment`.
#' @param option_type Character. One of "call" or "put".
#' @param n_terms Integer. Number of COS expansion terms. Default is 256.
#' @param truncation Numeric. Truncation domain size controlling the COS
#'   integration interval. Default is 10.
#'
#' @return A tibble with columns `strike_adjustment`, `strike_multiplier`, and
#'   `price` (time-zero option value).
#'
#' @export
#' @seealso [forward_start_characteristic_function()], [cos_call_put_price()]
price_forward_start_heston <- function(process_spec,
                                       start,
                                       maturity,
                                       strike_adjustment = 0,
                                       option_type = c("call", "put"),
                                       n_terms = 256L,
                                       truncation = 10) {
  checkmate::assert_class(process_spec, "heston_spec")
  checkmate::assert_number(start, lower = 0, finite = TRUE)
  checkmate::assert_number(maturity, lower = 0, finite = TRUE)
  if (start >= maturity) {
    rlang::abort("`start` must be strictly less than `maturity`")
  }
  checkmate::assert_numeric(strike_adjustment, any.missing = FALSE, finite = TRUE)
  option_type <- rlang::arg_match(option_type)
  checkmate::assert_integerish(n_terms, lower = 1, len = 1)
  checkmate::assert_number(truncation, lower = 0, finite = TRUE)

  n_terms <- as.integer(n_terms)
  strike_adjustment <- as.numeric(strike_adjustment)
  strike_multiplier <- 1 + strike_adjustment
  if (any(strike_multiplier <= 0)) {
    rlang::abort("Strike multiplier 1 + strike_adjustment must be positive")
  }

  spec <- process_spec
  cf_ratio <- function(u) {
    forward_start_characteristic_function(
      u = u,
      start = start,
      maturity = maturity,
      risk_free_rate = spec$risk_free_rate,
      dividend_yield = spec$dividend_yield,
      mean_reversion = spec$mean_reversion,
      long_term_variance = spec$long_term_variance,
      vol_of_vol = spec$vol_of_vol,
      initial_variance = spec$initial_variance,
      correlation = spec$correlation
    )
  }

  tau <- maturity - start
  a <- -truncation * sqrt(tau)
  b <- truncation * sqrt(tau)
  k_indices <- seq_len(n_terms) - 1
  u <- k_indices * pi / (b - a)

  payoff_coefficients <- cos_coefficients(option_type, a, b, k_indices)
  cf_values <- cf_ratio(u)
  cf_values[1] <- cf_values[1] * 0.5
  weights <- payoff_coefficients * cf_values

  x0 <- log(1 / strike_multiplier)
  exponential_matrix <- exp(1i * outer(x0 - a, u))
  cos_values <- Re(exponential_matrix %*% weights)

  discount <- exp(-spec$risk_free_rate * maturity)
  normalized_prices <- discount * strike_multiplier * as.numeric(cos_values)
  spot_scaling <- spec$initial_price * exp(-spec$dividend_yield * start)

  tibble::tibble(
    strike_adjustment = strike_adjustment,
    strike_multiplier = strike_multiplier,
    price = spot_scaling * normalized_prices
  )
}
