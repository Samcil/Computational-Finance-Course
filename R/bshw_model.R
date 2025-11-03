#' Black-Scholes Hull-White Hybrid Specification
#'
#' Construct a hybrid specification combining log-normal equity dynamics with
#' a one-factor Hull-White short-rate model. The specification reuses the
#' calibrated short-rate functions to form characteristic functions and pricing
#' routines consistent with the accompanying Python lectures.
#'
#' @param spot Numeric. Current equity spot price.
#' @param short_rate A `short_rate_spec` constructed with `model = "hull_white"`.
#' @param equity_vol Numeric. Instantaneous equity volatility.
#' @param correlation Numeric between -1 and 1 capturing equity/short-rate
#'   correlation.
#'
#' @return An object of class `bshw_spec` inheriting from `hybrid_spec`.
#' @export
#' @examples
#' curve <- tibble::tibble(tenor = c(0, 1, 2, 5, 10), discount_factor = exp(-0.03 * tenor))
#' hw_spec <- short_rate_spec(model = "hull_white", volatility = 0.01, mean_reversion = 0.1, curve = curve)
#' spec <- bshw_spec(spot = 100, short_rate = hw_spec, equity_vol = 0.2, correlation = 0.3)
#' price_bshw_option_cos(spec, strikes = c(90, 100, 110), maturity = 2)
#'
bshw_spec <- function(spot,
                      short_rate,
                      equity_vol,
                      correlation) {
  checkmate::assert_number(spot, lower = .Machine$double.eps, finite = TRUE)
  checkmate::assert_class(short_rate, "short_rate_spec")
  if (!identical(short_rate$args$model, "hull_white")) {
    rlang::abort("`short_rate` must be a Hull-White specification")
  }
  checkmate::assert_number(equity_vol, lower = 0, finite = TRUE)
  checkmate::assert_number(correlation, lower = -1, upper = 1, finite = TRUE)

  mean_reversion <- short_rate$args$mean_reversion
  if (mean_reversion <= 0) {
    rlang::abort("Hull-White mean reversion must be strictly positive")
  }

  state <- short_rate_state(short_rate)

  spec <- new_hybrid_spec(
    class = "bshw_spec",
    args = list(
      spot = spot,
      short_rate = short_rate,
      equity_vol = equity_vol,
      correlation = correlation
    ),
    process_type = "bshw"
  )

  spec <- set_process_metadata(
    spec,
    model_variant = "bshw",
    state_variables = c("equity_price", "short_rate"),
    engines = list(
      price = c("cos", "black76")
    )
  )

  spec$method$engine_state <- list(
    discount_fun = state$discount_fun,
    theta_fun = state$theta_fun_vec,
    initial_rate = short_rate$args$initial_rate,
    mean_reversion = mean_reversion,
    short_rate_vol = short_rate$args$volatility
  )

  spec
}

#' @export
print.bshw_spec <- function(x, ...) {
  cli::cli_h2("Black-Scholes Hull-White Specification")
  cli::cli_dl(c(
    "Spot" = cli::col_cyan(format(x$spot, digits = 6)),
    "Equity volatility" = cli::col_cyan(format(x$equity_vol, digits = 6)),
    "Correlation" = cli::col_cyan(format(x$correlation, digits = 6)),
    "Short-rate mean reversion" = cli::col_cyan(format(x$short_rate$args$mean_reversion, digits = 6)),
    "Short-rate volatility" = cli::col_cyan(format(x$short_rate$args$volatility, digits = 6))
  ))
  cli::cli_text("")
  cli::cli_text("{.strong Lineage:} {format_spec_lineage(x)}")
  invisible(x)
}

#' Black-Scholes Hull-White Characteristic Function
#'
#' Evaluate the characteristic function of the Black-Scholes Hull-White hybrid
#' under the terminal measure. The implementation follows Grzelak and
#' Oosterlee (2019) and matches the Python reference scripts bundled with the
#' course.
#'
#' @param u Numeric vector of integration arguments.
#' @param maturity Numeric. Time to maturity.
#' @param theta_fun Function returning the Hull-White drift adjustment \eqn{\theta(t)}.
#' @param initial_rate Numeric. Short rate at time zero.
#' @param mean_reversion Numeric. Hull-White mean reversion parameter.
#' @param short_rate_vol Numeric. Short-rate volatility.
#' @param equity_vol Numeric. Equity spot volatility.
#' @param correlation Numeric. Correlation between equity and short rate.
#' @param integration_points Integer. Number of trapezoidal steps for the drift
#'   integral.
#'
#' @return Complex vector of characteristic function evaluations.
#' @export
bshw_characteristic_function <- function(u,
                                         maturity,
                                         theta_fun,
                                         initial_rate,
                                         mean_reversion,
                                         short_rate_vol,
                                         equity_vol,
                                         correlation,
                                         integration_points = 2000L) {
  checkmate::assert_numeric(u, any.missing = FALSE)
  checkmate::assert_number(maturity, lower = 0, finite = TRUE)
  checkmate::assert_function(theta_fun)
  checkmate::assert_number(initial_rate, finite = TRUE)
  checkmate::assert_number(mean_reversion, lower = .Machine$double.eps, finite = TRUE)
  checkmate::assert_number(short_rate_vol, lower = 0, finite = TRUE)
  checkmate::assert_number(equity_vol, lower = 0, finite = TRUE)
  checkmate::assert_number(correlation, lower = -1, upper = 1, finite = TRUE)
  checkmate::assert_integerish(integration_points, lower = 1, len = 1)

  integration_points <- as.integer(integration_points)

  lambda <- mean_reversion
  eta <- short_rate_vol
  sigma <- equity_vol
  rho <- correlation

  grid <- seq(0, maturity, length.out = integration_points + 1)
  theta_vals <- theta_fun(pmax(maturity - grid, 0))
  exp_neg <- exp(-lambda * grid)
  step_size <- if (integration_points > 0) maturity / integration_points else 0

  integrals <- bshw_theta_integrals_cpp(u, lambda, theta_vals, exp_neg, step_size)
  iu <- 1i * u

  term1 <- 0.5 * sigma^2 * iu * (iu - 1) * maturity
  term2 <- iu * rho * sigma * eta / lambda * (iu - 1) * (maturity + (exp(-lambda * maturity) - 1) / lambda)
  term3 <- (eta^2 / (4 * lambda^3)) * (1i + u)^2 * (3 + exp(-2 * lambda * maturity) - 4 * exp(-lambda * maturity) - 2 * lambda * maturity)
  term4 <- lambda * integrals
  c_T <- (iu - 1) / lambda * (1 - exp(-lambda * maturity))

  exp(term1 + term2 + term3 + term4 + c_T * initial_rate)
}

#' Price BSHW Options via COS Method
#'
#' Price European calls and puts under the Black-Scholes Hull-White model using
#' the COS expansion. Relies on `cos_call_put_price_stoch_ir()` for the payoff
#' coefficients.
#'
#' @param process_spec A `bshw_spec` created by [bshw_spec()].
#' @param strikes Numeric vector of strikes.
#' @param maturity Numeric. Time to maturity.
#' @param n_terms Integer number of COS terms.
#' @param truncation Numeric truncation width.
#' @param integration_points Integer number of steps for the drift integral.
#'
#' @return Tibble with columns `strike`, `call_price`, and `put_price`.
#' @export
price_bshw_option_cos <- function(process_spec,
                                  strikes,
                                  maturity,
                                  n_terms = 256L,
                                  truncation = 8,
                                  integration_points = 2000L) {
  checkmate::assert_class(process_spec, "bshw_spec")
  checkmate::assert_numeric(strikes, lower = .Machine$double.eps, any.missing = FALSE, finite = TRUE)
  checkmate::assert_number(maturity, lower = 0, finite = TRUE)
  checkmate::assert_integerish(n_terms, lower = 1, len = 1)
  checkmate::assert_number(truncation, lower = 0, finite = TRUE)
  checkmate::assert_integerish(integration_points, lower = 1, len = 1)

  state <- process_spec$method$engine_state
  discount_factor <- state$discount_fun(maturity)

  cf <- function(u) {
    bshw_characteristic_function(
      u = u,
      maturity = maturity,
      theta_fun = state$theta_fun,
      initial_rate = state$initial_rate,
      mean_reversion = state$mean_reversion,
      short_rate_vol = state$short_rate_vol,
      equity_vol = process_spec$equity_vol,
      correlation = process_spec$correlation,
      integration_points = integration_points
    )
  }

  call_prices <- cos_call_put_price_stoch_ir(
    cf = cf,
    option_type = "call",
    spot = process_spec$spot,
    maturity = maturity,
    strikes = strikes,
    discount_factor = discount_factor,
    n_terms = n_terms,
    truncation = truncation
  )

  put_prices <- cos_call_put_price_stoch_ir(
    cf = cf,
    option_type = "put",
    spot = process_spec$spot,
    maturity = maturity,
    strikes = strikes,
    discount_factor = discount_factor,
    n_terms = n_terms,
    truncation = truncation
  )

  tibble::tibble(
    strike = strikes,
    call_price = call_prices$price,
    put_price = put_prices$price
  )
}

#' Closed-Form BSHW Option Prices
#'
#' Computes Black-76 style prices for the BSHW hybrid using the equivalent
#' normal volatility derived from integrating the forward variance curve.
#'
#' @inheritParams price_bshw_option_cos
#' @return Tibble with columns `strike`, `call_price`, and `put_price`.
#' @export
price_bshw_option_black76 <- function(process_spec,
                                      strikes,
                                      maturity,
                                      integration_points = 2000L) {
  checkmate::assert_class(process_spec, "bshw_spec")
  checkmate::assert_numeric(strikes, lower = .Machine$double.eps, any.missing = FALSE, finite = TRUE)
  checkmate::assert_number(maturity, lower = 0, finite = TRUE)
  checkmate::assert_integerish(integration_points, lower = 1, len = 1)

  state <- process_spec$method$engine_state
  discount_factor <- state$discount_fun(maturity)
  forward <- process_spec$spot / discount_factor
  equivalent_vol <- bshw_equivalent_volatility(
    maturity = maturity,
    equity_vol = process_spec$equity_vol,
    short_rate_vol = state$short_rate_vol,
    correlation = process_spec$correlation,
    mean_reversion = state$mean_reversion,
    integration_points = integration_points
  )

  calls <- black_forward_price(forward, strikes, equivalent_vol, maturity)
  puts <- calls - forward + strikes
  tibble::tibble(
    strike = strikes,
    call_price = discount_factor * calls,
    put_price = discount_factor * puts
  )
}

#' Equivalent BSHW Volatility
#'
#' Integrate the instantaneous forward volatility curve implied by the hybrid
#' model to obtain the Black-76 equivalent volatility.
#'
#' @param maturity Numeric maturity.
#' @param equity_vol Equity volatility parameter.
#' @param short_rate_vol Hull-White short-rate volatility.
#' @param correlation Equity/short-rate correlation.
#' @param mean_reversion Hull-White mean reversion.
#' @param integration_points Number of trapezoidal steps.
#'
#' @return Numeric scalar volatility.
#' @export
bshw_equivalent_volatility <- function(maturity,
                                       equity_vol,
                                       short_rate_vol,
                                       correlation,
                                       mean_reversion,
                                       integration_points = 2000L) {
  checkmate::assert_number(maturity, lower = .Machine$double.eps, finite = TRUE)
  checkmate::assert_number(equity_vol, lower = 0, finite = TRUE)
  checkmate::assert_number(short_rate_vol, lower = 0, finite = TRUE)
  checkmate::assert_number(correlation, lower = -1, upper = 1, finite = TRUE)
  checkmate::assert_number(mean_reversion, lower = .Machine$double.eps, finite = TRUE)
  checkmate::assert_integerish(integration_points, lower = 1, len = 1)

  lambda <- mean_reversion
  eta <- short_rate_vol
  sigma <- equity_vol
  rho <- correlation

  integrand <- function(t) {
    br <- (exp(-lambda * (maturity - t)) - 1) / lambda
    forward_sigma_sq <- sigma^2 + eta^2 * br^2 - 2 * rho * sigma * eta * br
    forward_sigma_sq
  }

  variance_integral <- trapezoidal_integral(integrand, lower = 0, upper = maturity, n_steps = integration_points)
  sqrt(variance_integral / maturity)
}

black_forward_price <- function(forward, strike, volatility, maturity) {
  checkmate::assert_numeric(forward, lower = .Machine$double.eps, any.missing = FALSE, finite = TRUE)
  checkmate::assert_numeric(strike, lower = .Machine$double.eps, any.missing = FALSE, finite = TRUE)
  checkmate::assert_number(volatility, lower = 0, finite = TRUE)
  checkmate::assert_number(maturity, lower = 0, finite = TRUE)

  sqrt_t <- sqrt(maturity)
  d1 <- (log(forward / strike) + 0.5 * volatility^2 * maturity) / (volatility * sqrt_t)
  d2 <- d1 - volatility * sqrt_t
  forward * stats::pnorm(d1) - strike * stats::pnorm(d2)
}
