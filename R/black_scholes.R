#' Create Black-Scholes Option Specification
#'
#' Defines a European option contract under the Black-Scholes-Merton model.
#' The specification integrates with `price_options()` to compute prices and
#' Greeks in a pipe-friendly workflow.
#'
#' @param option_type Character. One of `"call"` or `"put"`.
#' @param strike Numeric. Strike price \eqn{K > 0}.
#' @param maturity Numeric. Time to maturity (in years) \eqn{T > 0}.
#' @param risk_free_rate Numeric. Continuously compounded risk-free rate \eqn{r}.
#' @param dividend_yield Numeric. Continuous dividend yield \eqn{q}. Default is 0.
#'
#' @details
#' The Black-Scholes framework assumes log-normal dynamics for the underlying
#' asset price \eqn{S(t)} under the risk-neutral measure:
#'
#' \deqn{dS(t) = (r - q)S(t) dt + \sigma S(t) dW(t)}
#'
#' @return An object of class `black_scholes_spec`.
#'
#' @examples
#' spec <- black_scholes_spec(
#'   option_type = "call",
#'   strike = 100,
#'   maturity = 1,
#'   risk_free_rate = 0.02,
#'   dividend_yield = 0.01
#' )
#'
#' spec |>
#'   price_options(spot = 105, volatility = 0.2)
#'
#' @export
black_scholes_spec <- function(option_type = c("call", "put"),
                               strike,
                               maturity,
                               risk_free_rate,
                               dividend_yield = 0) {
  option_type <- rlang::arg_match(option_type)
  checkmate::assert_number(strike, lower = .Machine$double.eps, finite = TRUE)
  checkmate::assert_number(maturity, lower = .Machine$double.eps, finite = TRUE)
  checkmate::assert_number(risk_free_rate, finite = TRUE)
  checkmate::assert_number(dividend_yield, finite = TRUE)

  structure(
    list(
      option_type = option_type,
      strike = strike,
      maturity = maturity,
      risk_free_rate = risk_free_rate,
      dividend_yield = dividend_yield
    ),
    class = c("black_scholes_spec", "process_spec")
  )
}

#' @export
print.black_scholes_spec <- function(x, ...) {
  cli::cli_h2("Black-Scholes Specification")
  cli::cli_dl(c(
    "Option type" = cli::col_cyan(x$option_type),
    "Strike" = cli::col_green(x$strike),
    "Maturity" = cli::col_blue(x$maturity),
    "Risk-free rate" = cli::col_magenta(x$risk_free_rate),
    "Dividend yield" = cli::col_magenta(x$dividend_yield)
  ))
  cli::cli_text("")
  cli::cli_alert_info("Use {.fn price_options} to compute prices and Greeks")
  invisible(x)
}

#' Price Options and Compute Greeks
#'
#' Generic function for pricing option specifications and returning the full set
#' of analytic Greeks in tidy format.
#'
#' @param model_spec Option specification object.
#' @param ... Additional arguments passed to methods.
#'
#' @return A tibble containing model-dependent outputs.
#'
#' @export
price_options <- function(model_spec, ...) {
  UseMethod("price_options", model_spec)
}

#' Black-Scholes Pricing and Greeks
#'
#' Computes European option price together with Delta, Gamma, Vega, Theta and
#' Rho using closed-form Black-Scholes formulas.
#'
#' @inheritParams price_options
#' @param spot Numeric vector. Spot prices \eqn{S_0 > 0}.
#' @param volatility Numeric vector. Volatilities \eqn{\sigma > 0}.
#'
#' @return A tibble with columns:
#'   \describe{
#'     \item{spot}{Input spot price}
#'     \item{volatility}{Input volatility}
#'     \item{price}{Option price}
#'     \item{delta}{Sensitivity to spot}
#'     \item{gamma}{Second derivative w.r.t. spot}
#'     \item{vega}{Sensitivity to volatility}
#'     \item{theta}{Sensitivity to time (per year)}
#'     \item{rho}{Sensitivity to interest rates}
#'   }
#'
#' @examples
#' black_scholes_spec("put", strike = 100, maturity = 0.5, risk_free_rate = 0.01) |>
#'   price_options(spot = c(95, 100, 105), volatility = c(0.15, 0.2))
#'
#' @export
price_options.black_scholes_spec <- function(model_spec,
                                             spot,
                                             volatility,
                                             ...) {
  checkmate::assert_numeric(spot, lower = .Machine$double.eps, any.missing = FALSE, finite = TRUE)
  checkmate::assert_numeric(volatility, lower = .Machine$double.eps, any.missing = FALSE, finite = TRUE)
  maturity <- model_spec$maturity
  strike <- model_spec$strike
  r <- model_spec$risk_free_rate
  q <- model_spec$dividend_yield
  option_type <- model_spec$option_type

  grid <- tidyr::crossing(
    tibble::tibble(spot = spot),
    tibble::tibble(volatility = volatility)
  )

  forward_factor <- exp(-q * maturity)
  discount_factor <- exp(-r * maturity)
  sqrt_maturity <- sqrt(maturity)

  d1 <- (
    log(grid$spot / strike) + (r - q + 0.5 * grid$volatility^2) * maturity
  ) / (grid$volatility * sqrt_maturity)
  d2 <- d1 - grid$volatility * sqrt_maturity

  pdf_d1 <- stats::dnorm(d1)
  cdf_d1 <- stats::pnorm(d1)
  cdf_minus_d1 <- stats::pnorm(-d1)
  cdf_d2 <- stats::pnorm(d2)
  cdf_minus_d2 <- stats::pnorm(-d2)

  if (option_type == "call") {
    price <- grid$spot * forward_factor * cdf_d1 - strike * discount_factor * cdf_d2
    delta <- forward_factor * cdf_d1
    theta <- -(
      grid$spot * forward_factor * pdf_d1 * grid$volatility / (2 * sqrt_maturity)
    ) - r * strike * discount_factor * cdf_d2 + q * grid$spot * forward_factor * cdf_d1
    rho <- maturity * strike * discount_factor * cdf_d2
  } else {
    price <- strike * discount_factor * cdf_minus_d2 - grid$spot * forward_factor * cdf_minus_d1
    delta <- forward_factor * (cdf_d1 - 1)
    theta <- -(
      grid$spot * forward_factor * pdf_d1 * grid$volatility / (2 * sqrt_maturity)
    ) + r * strike * discount_factor * cdf_minus_d2 - q * grid$spot * forward_factor * cdf_minus_d1
    rho <- -maturity * strike * discount_factor * cdf_minus_d2
  }

  gamma <- forward_factor * pdf_d1 / (grid$spot * grid$volatility * sqrt_maturity)
  vega <- grid$spot * forward_factor * pdf_d1 * sqrt_maturity

  tibble::tibble(
    spot = grid$spot,
    volatility = grid$volatility,
    price = price,
    delta = delta,
    gamma = gamma,
    vega = vega,
    theta = theta,
    rho = rho
  )
}
