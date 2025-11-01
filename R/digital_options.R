#' Black-Scholes Cash-or-Nothing Digital Call
#'
#' Prices a cash-or-nothing digital call option under the
#' Black-Scholes-Merton model using the closed-form solution.
#'
#' @param model_spec A `black_scholes_spec` object created with
#'   [black_scholes_spec()].
#' @param spot Numeric vector of current spot prices.
#' @param volatility Numeric vector of volatilities.
#' @param payout Numeric cash payout delivered when the option finishes in the
#'   money. Default is 1.
#'
#' @return A tibble with columns `spot`, `volatility`, and `price`.
#'
#' @examples
#' bs_spec <- black_scholes_spec(
#'   option_type = "call",
#'   strike = 100,
#'   maturity = 1,
#'   risk_free_rate = 0.02,
#'   dividend_yield = 0.01
#' )
#' price_digital_call(bs_spec, spot = c(95, 105), volatility = c(0.15, 0.2))
#'
#' @export
price_digital_call <- function(model_spec,
                               spot,
                               volatility,
                               payout = 1) {
  digital_option_price(
    model_spec = model_spec,
    spot = spot,
    volatility = volatility,
    payout = payout,
    option_type = "call"
  )
}

#' Black-Scholes Cash-or-Nothing Digital Put
#'
#' Prices a cash-or-nothing digital put option under the
#' Black-Scholes-Merton model using the closed-form solution.
#'
#' @inheritParams price_digital_call
#'
#' @export
price_digital_put <- function(model_spec,
                              spot,
                              volatility,
                              payout = 1) {
  digital_option_price(
    model_spec = model_spec,
    spot = spot,
    volatility = volatility,
    payout = payout,
    option_type = "put"
  )
}

digital_option_price <- function(model_spec,
                                 spot,
                                 volatility,
                                 payout,
                                 option_type) {
  checkmate::assert_class(model_spec, "black_scholes_spec")
  checkmate::assert_numeric(spot, lower = .Machine$double.eps, any.missing = FALSE, finite = TRUE)
  checkmate::assert_numeric(volatility, lower = .Machine$double.eps, any.missing = FALSE, finite = TRUE)
  checkmate::assert_number(payout, lower = 0, finite = TRUE)

  maturity <- model_spec$maturity
  strike <- model_spec$strike
  r <- model_spec$risk_free_rate
  q <- model_spec$dividend_yield

  grid <- tidyr::crossing(
    tibble::tibble(spot = spot),
    tibble::tibble(volatility = volatility)
  )

  discount_factor <- exp(-r * maturity)
  sqrt_maturity <- sqrt(maturity)

  d2 <- (
    log(grid$spot / strike) + (r - q - 0.5 * grid$volatility^2) * maturity
  ) / (grid$volatility * sqrt_maturity)

  if (option_type == "call") {
    price <- payout * discount_factor * stats::pnorm(d2)
  } else {
    price <- payout * discount_factor * stats::pnorm(-d2)
  }

  tibble::tibble(
    spot = grid$spot,
    volatility = grid$volatility,
    price = price
  )
}

#' Digital Option Pricing via COS Density Approximation
#'
#' Approximates cash-or-nothing digital option prices using the COS method by
#' reconstructing the terminal log-price density and integrating the relevant
#' tail probabilities.
#'
#' @param cf Characteristic function of the log-price under the risk-neutral
#'   measure.
#' @param option_type Character string, either "call" or "put".
#' @param spot Numeric scalar for the current spot price.
#' @param risk_free_rate Numeric risk-free rate.
#' @param maturity Numeric time to maturity in years.
#' @param strikes Numeric vector of strikes.
#' @param n_terms Integer number of COS terms used to reconstruct the density.
#' @param truncation Numeric truncation parameter controlling the integration
#'   interval.
#' @param grid_size Integer number of grid points used for numerical
#'   integration. Default is 2001.
#' @param payout Numeric cash payout delivered when the option finishes in the
#'   money. Default is 1.
#' @param lower_bound Optional numeric scalar specifying the lower bound of the
#'   truncated log-price domain. Must be supplied together with `upper_bound`.
#' @param upper_bound Optional numeric scalar specifying the upper bound of the
#'   truncated log-price domain. Must be supplied together with `lower_bound`.
#'
#' @return A tibble with columns `strike` and `price`.
#'
#' @export
#' @seealso [cos_density_recovery()], [cos_call_put_price()]
digital_cos_method <- function(cf,
                               option_type = c("call", "put"),
                               spot,
                               risk_free_rate,
                               maturity,
                               strikes,
                               n_terms = 256L,
                               truncation = 8,
                               grid_size = 2001L,
                               payout = 1,
                               lower_bound = NULL,
                               upper_bound = NULL) {
  option_type <- rlang::arg_match(option_type)
  checkmate::assert_function(cf)
  checkmate::assert_number(spot, lower = .Machine$double.eps, finite = TRUE)
  checkmate::assert_number(risk_free_rate, finite = TRUE)
  checkmate::assert_number(maturity, lower = 0, finite = TRUE)
  checkmate::assert_numeric(strikes, lower = .Machine$double.eps, any.missing = FALSE, finite = TRUE)
  checkmate::assert_integerish(n_terms, lower = 1, len = 1)
  checkmate::assert_number(truncation, lower = 0, finite = TRUE)
  checkmate::assert_integerish(grid_size, lower = 101, len = 1)
  checkmate::assert_number(payout, lower = 0, finite = TRUE)
  if (xor(is.null(lower_bound), is.null(upper_bound))) {
    rlang::abort("`lower_bound` and `upper_bound` must be supplied together")
  }
  if (!is.null(lower_bound) && !is.null(upper_bound) && lower_bound >= upper_bound) {
    rlang::abort("`lower_bound` must be strictly less than `upper_bound`")
  }

  n_terms <- as.integer(n_terms)
  grid_size <- as.integer(grid_size)
  strikes <- as.numeric(strikes)

  if (!is.null(lower_bound) && !is.null(upper_bound)) {
    a <- lower_bound
    b <- upper_bound
  } else {
    half_width <- truncation * sqrt(maturity)
    center <- log(spot)
    a <- center - half_width
    b <- center + half_width
  }
  x_grid <- seq(a, b, length.out = grid_size)

  density_tbl <- cos_density_recovery(
    cf = cf,
    x = x_grid,
    maturity = maturity,
    n_terms = n_terms,
    truncation = truncation,
    lower_bound = a,
    upper_bound = b
  )
  density <- pmax(density_tbl$density, 0)
  dx <- density_tbl$x[2] - density_tbl$x[1]

  total_mass <- sum(density) * dx
  if (total_mass <= 0) {
    rlang::abort("Recovered density has non-positive total mass. Adjust truncation or integration bounds.")
  }
  if (!isTRUE(all.equal(total_mass, 1, tolerance = 1e-6))) {
    density <- density / total_mass
    total_mass <- 1
  }

  cdf <- cumsum(density) * dx
  cdf_fun <- stats::approxfun(density_tbl$x, cdf, method = "linear", rule = 2)

  thresholds <- log(strikes)
  discount_factor <- exp(-risk_free_rate * maturity)

  probabilities <- purrr::map_dbl(
    thresholds,
    \(threshold) {
      cdf_at_threshold <- cdf_fun(threshold)
      if (option_type == "call") {
        prob <- total_mass - cdf_at_threshold
      } else {
        prob <- cdf_at_threshold
      }
      prob <- max(prob, 0)
      min(prob, total_mass)
    }
  )

  tibble::tibble(
    strike = strikes,
    price = discount_factor * payout * probabilities
  )
}
