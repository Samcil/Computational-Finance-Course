#' Asian Option Pricing via Monte Carlo
#'
#' Prices arithmetic-average Asian options using Monte Carlo simulation under a
#' Geometric Brownian Motion (GBM) model. Supports standard sampling and
#' antithetic variates for variance reduction.
#'
#' @param process_spec A `gbm_spec` object describing the underlying asset
#'   dynamics under the risk-neutral measure.
#' @param strike Numeric vector of strike prices.
#' @param maturity Numeric. Time to maturity in years.
#' @param n_paths Integer. Number of Monte Carlo paths to simulate. Default is
#'   1000.
#' @param n_steps Integer. Number of monitoring points along each path. Default
#'   is 252.
#' @param seed Optional integer seed for reproducibility. When supplied, the
#'   same seed is used for standard and antithetic simulations.
#' @param method Character. Sampling approach, either "standard" Monte Carlo or
#'   "antithetic" variates. Default is "standard".
#'
#' @return A tibble with columns `strike`, `option_type`, `method`, `price`,
#'   `std_error`, `n_paths`, `n_steps`, and `maturity`.
#'
#' @examples
#' spec <- gbm_spec(initial_value = 100, drift = 0.03, volatility = 0.2)
#' price_asian_call(spec, strike = c(95, 105), maturity = 1, n_paths = 2000,
#'                  n_steps = 126, seed = 42)
#'
#' @export
price_asian_call <- function(process_spec,
                             strike,
                             maturity,
                             n_paths = 1000L,
                             n_steps = 252L,
                             seed = NULL,
                             method = c("standard", "antithetic")) {
  method <- rlang::arg_match(method)
  asian_option_price(
    option_type = "call",
    process_spec = process_spec,
    strike = strike,
    maturity = maturity,
    n_paths = n_paths,
    n_steps = n_steps,
    seed = seed,
    method = method
  )
}

#' @rdname price_asian_call
#' @export
price_asian_put <- function(process_spec,
                            strike,
                            maturity,
                            n_paths = 10000L,
                            n_steps = 252L,
                            seed = NULL,
                            method = c("standard", "antithetic")) {
  method <- rlang::arg_match(method)
  asian_option_price(
    option_type = "put",
    process_spec = process_spec,
    strike = strike,
    maturity = maturity,
    n_paths = n_paths,
    n_steps = n_steps,
    seed = seed,
    method = method
  )
}

#' Compare Variance Reduction Strategies for Asian Options
#'
#' Convenience wrapper that computes Asian option prices using both standard
#' Monte Carlo and antithetic variates, allowing quick comparison of point
#' estimates and standard errors.
#'
#' @inheritParams price_asian_call
#' @param option_type Character. Option payoff, either "call" or "put".
#'
#' @return A tibble combining results for each variance reduction method.
#'
#' @export
asian_variance_reduction <- function(process_spec,
                                     strike,
                                     maturity,
                                     option_type = c("call", "put"),
                                     n_paths = 1000L,
                                     n_steps = 252L,
                                     seed = NULL) {
  option_type <- rlang::arg_match(option_type)
  standard <- asian_option_price(
    option_type = option_type,
    process_spec = process_spec,
    strike = strike,
    maturity = maturity,
    n_paths = n_paths,
    n_steps = n_steps,
    seed = seed,
    method = "standard"
  )
  antithetic <- asian_option_price(
    option_type = option_type,
    process_spec = process_spec,
    strike = strike,
    maturity = maturity,
    n_paths = n_paths,
    n_steps = n_steps,
    seed = seed,
    method = "antithetic"
  )

  dplyr::bind_rows(standard, antithetic)
}

asian_option_price <- function(option_type,
                               process_spec,
                               strike,
                               maturity,
                               n_paths,
                               n_steps,
                               seed,
                               method) {
  checkmate::assert_class(process_spec, "gbm_spec")
  checkmate::assert_numeric(strike, lower = 0, any.missing = FALSE, finite = TRUE)
  checkmate::assert_number(maturity, lower = 0, finite = TRUE)
  checkmate::assert_integerish(n_paths, lower = 1, len = 1)
  checkmate::assert_integerish(n_steps, lower = 1, len = 1)
  if (!is.null(seed)) {
    checkmate::assert_integerish(seed, len = 1)
    set.seed(as.integer(seed))
  }

  n_paths <- as.integer(n_paths)
  n_steps <- as.integer(n_steps)

  S0 <- process_spec$initial_value
  r <- process_spec$drift
  sigma <- process_spec$volatility

  dt <- maturity / n_steps
  drift_term <- (r - 0.5 * sigma^2) * dt
  vol_term <- sigma * sqrt(dt)

  base_normals <- generate_standardized_normals(n_paths, n_steps)

  simulate_averages <- function(z_matrix) {
    log_increments <- drift_term + vol_term * z_matrix
    log_paths <- compute_cumulative_paths(log(S0), log_increments)
    rowMeans(exp(log_paths))
  }

  averages_pos <- simulate_averages(base_normals)
  payoffs_pos <- asian_payoff_matrix(averages_pos, strike, option_type)

  payoffs <- payoffs_pos

  if (identical(method, "antithetic")) {
    averages_neg <- simulate_averages(-base_normals)
    payoffs_neg <- asian_payoff_matrix(averages_neg, strike, option_type)
    payoffs <- (payoffs_pos + payoffs_neg) / 2
  }

  discount_factor <- exp(-r * maturity)
  sample_size <- nrow(payoffs)
  col_mean <- colMeans(payoffs)
  col_sq_mean <- colMeans(payoffs^2)

  if (sample_size > 1) {
    moment_var <- pmax(col_sq_mean - col_mean^2, 0)
    unbiased_var <- moment_var * sample_size / (sample_size - 1)
    std_error <- sqrt(unbiased_var / sample_size)
  } else {
    std_error <- rep(NA_real_, length(col_mean))
  }

  tibble::tibble(
    strike = strike,
    option_type = option_type,
    method = method,
    price = discount_factor * col_mean,
    std_error = discount_factor * std_error,
    n_paths = sample_size,
    n_steps = n_steps,
    maturity = maturity
  )
}

asian_payoff_matrix <- function(average_prices, strike, option_type) {
  n_paths <- length(average_prices)
  n_strikes <- length(strike)

  avg_matrix <- matrix(average_prices, nrow = n_paths, ncol = n_strikes)
  strike_matrix <- matrix(rep(strike, each = n_paths), nrow = n_paths, ncol = n_strikes)

  if (identical(option_type, "call")) {
    pmax(avg_matrix - strike_matrix, 0)
  } else {
    pmax(strike_matrix - avg_matrix, 0)
  }
}
