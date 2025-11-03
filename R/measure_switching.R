#' Measure Switching Utilities for Diffusions
#'
#' Provides helpers to simulate paths under the physical (P) and risk-neutral
#' (Q) measures using common random numbers together with Radon--Nikodym
#' densities for reweighting. The implementation is aligned with the
#' `PathsUnderQandPmeasure.py` lecture script.
#'
#' @param spec A `gbm_spec` representing dynamics under the physical measure.
#' @param risk_free_rate Numeric risk-free rate used for the risk-neutral
#'   measure.
#' @param n_paths Integer number of Monte Carlo paths.
#' @param n_steps Integer number of Euler steps.
#' @param maturity Numeric time horizon in years.
#' @param dividend_yield Optional continuously compounded dividend yield. Must
#'   be scalar. Defaults to `0`.
#' @param market_price_of_risk Optional scalar market price of risk (lambda).
#'   When `NULL` it is implied from the specification drift via
#'   `(mu - (r - q)) / sigma`.
#' @param seed Integer seed controlling the shared random numbers.
#'
#' @return Tibble with columns `measure`, `path_id`, `time`, `stock_price`, and
#'   `radon_nikodym`. Rows tagged with `measure = "P"` include the Radon--Nikodym
#'   derivative `dQ/dP` for the corresponding path and time. The tibble inherits
#'   the S3 class `measure_paths` and stores attributes `market_price_of_risk`,
#'   `risk_free_rate`, and `dividend_yield`.
#'
#' @examples
#' spec <- gbm_spec(initial_value = 100, drift = 0.08, volatility = 0.2)
#' paths <- q_measure_paths(
#'   spec = spec,
#'   risk_free_rate = 0.03,
#'   n_paths = 32,
#'   n_steps = 12,
#'   maturity = 1,
#'   dividend_yield = 0.01,
#'   seed = 123
#' )
#' dplyr::filter(paths, time == 1)
#'
#' @export
q_measure_paths <- function(spec,
                            risk_free_rate,
                            n_paths,
                            n_steps,
                            maturity,
                            dividend_yield = 0,
                            market_price_of_risk = NULL,
                            seed = 123) {
  checkmate::assert_true(inherits(spec, "gbm_spec"))
  checkmate::assert_number(risk_free_rate, finite = TRUE)
  checkmate::assert_integerish(n_paths, lower = 1, any.missing = FALSE)
  checkmate::assert_integerish(n_steps, lower = 1, any.missing = FALSE)
  checkmate::assert_number(maturity, lower = .Machine$double.eps, finite = TRUE)
  checkmate::assert_number(dividend_yield, finite = TRUE)
  checkmate::assert_integerish(seed, len = 1, any.missing = FALSE)
  if (!is.null(market_price_of_risk)) {
    checkmate::assert_number(market_price_of_risk, finite = TRUE)
  }

  S0 <- spec$initial_value
  mu_p <- spec$drift
  sigma <- spec$volatility

  checkmate::assert_number(S0, lower = .Machine$double.eps, finite = TRUE)
  checkmate::assert_number(sigma, lower = .Machine$double.eps, finite = TRUE)

  lambda <- if (is.null(market_price_of_risk)) {
    (mu_p - (risk_free_rate - dividend_yield)) / sigma
  } else {
    market_price_of_risk
  }

  dt <- maturity / n_steps
  sqrt_dt <- sqrt(dt)
  times <- seq(0, maturity, length.out = n_steps + 1)

  Z <- with_random_seed(seed, {
    generate_standardized_normals(n_paths, n_steps)
  })

  dW_p <- sqrt_dt * Z

  log_inc_p <- (mu_p - 0.5 * sigma^2) * dt + sigma * dW_p
  log_paths_p <- compute_cumulative_paths(log(S0), log_inc_p)
  prices_p <- exp(log_paths_p)

  mu_q <- risk_free_rate - dividend_yield
  dW_q <- dW_p + lambda * dt
  log_inc_q <- (mu_q - 0.5 * sigma^2) * dt + sigma * dW_q
  log_paths_q <- compute_cumulative_paths(log(S0), log_inc_q)
  prices_q <- exp(log_paths_q)

  W_p <- compute_cumulative_paths(0, dW_p)
  time_matrix <- matrix(times, nrow = n_paths, ncol = length(times), byrow = TRUE)
  radon_matrix <- exp(-lambda * W_p - 0.5 * lambda^2 * time_matrix)

  radon_tidy <- matrix_to_tidy(radon_matrix, times, value_name = "radon_nikodym")

  tidy_p <- matrix_to_tidy(prices_p, times, value_name = "stock_price") |>
    dplyr::mutate(measure = "P") |>
    dplyr::left_join(radon_tidy, by = c("path_id", "time"))

  tidy_q <- matrix_to_tidy(prices_q, times, value_name = "stock_price") |>
    dplyr::mutate(measure = "Q", radon_nikodym = 1)

  result <- dplyr::bind_rows(tidy_p, tidy_q) |>
    dplyr::relocate(measure, .before = path_id)

  attr(result, "market_price_of_risk") <- lambda
  attr(result, "risk_free_rate") <- risk_free_rate
  attr(result, "dividend_yield") <- dividend_yield
  class(result) <- unique(c("measure_paths", class(result)))
  result
}


resolve_lambda_grid <- function(market_price_of_risk, times, n_steps) {
  if (is.null(market_price_of_risk)) {
    return(rep(0, n_steps))
  }

  if (is.numeric(market_price_of_risk)) {
    if (length(market_price_of_risk) == 1) {
      return(rep(market_price_of_risk, n_steps))
    }
    checkmate::assert_numeric(market_price_of_risk, len = n_steps, any.missing = FALSE, finite = TRUE)
    return(market_price_of_risk)
  }

  if (is.function(market_price_of_risk)) {
    lambda_vals <- purrr::map_dbl(times, market_price_of_risk)
    checkmate::assert_numeric(lambda_vals, len = n_steps, any.missing = FALSE, finite = TRUE)
    return(lambda_vals)
  }

  rlang::abort("market_price_of_risk must be NULL, numeric scalar/vector, or function of time")
}


#' Measure Switching Utilities for Short-Rate Models
#'
#' Generates parallel simulations of Gaussian short-rate models under the
#' physical (P) and risk-neutral (Q) measures using a shared Brownian driver.
#' The routine returns Radon--Nikodym weights so that expectations computed
#' under the physical measure can be reweighted to match risk-neutral pricing.
#'
#' @param spec A `short_rate_spec` created by [short_rate_spec()].
#' @param n_paths Integer number of Monte Carlo paths.
#' @param n_steps Integer number of time steps.
#' @param maturity Numeric time horizon in years.
#' @param market_price_of_risk Market price of risk input controlling the
#'   Girsanov drift adjustment. Accepts a numeric scalar, numeric vector of
#'   length `n_steps`, or a function of time returning numeric values. Defaults
#'   to `0`, implying that P and Q coincide.
#' @param seed Integer seed for reproducibility.
#'
#' @return Tibble with columns `measure`, `path_id`, `time`, `short_rate`,
#'   `discount_factor`, and `radon_nikodym`. Rows labelled `measure = "P"`
#'   contain the Radon--Nikodym derivative `dQ/dP`. The tibble carries the class
#'   `short_rate_measure_paths` in addition to `measure_paths`.
#'
#' @examples
#' curve <- tibble::tibble(tenor = c(0, 1, 2), discount_factor = exp(-0.03 * tenor))
#' spec <- short_rate_spec(model = "ho_lee", volatility = 0.01, curve = curve)
#' paths <- short_rate_measure_paths(
#'   spec = spec,
#'   n_paths = 1000,
#'   n_steps = 120,
#'   maturity = 1,
#'   market_price_of_risk = 0.3,
#'   seed = 42
#' )
#' dplyr::filter(paths, time == 1, measure == "P")
#'
#' @export
short_rate_measure_paths <- function(spec,
                                     n_paths,
                                     n_steps,
                                     maturity,
                                     market_price_of_risk = 0,
                                     seed = 123) {
  checkmate::assert_true(inherits(spec, "short_rate_spec"))
  checkmate::assert_integerish(n_paths, lower = 1, any.missing = FALSE)
  checkmate::assert_integerish(n_steps, lower = 1, any.missing = FALSE)
  checkmate::assert_number(maturity, lower = .Machine$double.eps, finite = TRUE)
  checkmate::assert_integerish(seed, len = 1, any.missing = FALSE)

  args <- spec$args
  state <- short_rate_state(spec)

  sigma <- args$volatility
  a <- args$mean_reversion

  dt <- maturity / n_steps
  times <- seq(0, maturity, length.out = n_steps + 1)
  step_times <- times[-length(times)]

  lambda_step <- resolve_lambda_grid(market_price_of_risk, step_times, n_steps)

  shocks <- with_random_seed(seed, {
    generate_standardized_normals(n_paths, n_steps)
  })

  theta_values <- state$theta_fun_vec(step_times)

  evolved <- gaussian_measure_switch_cpp(
    shocks = shocks,
    theta = theta_values,
    initial_rate = args$initial_rate,
    mean_reversion = a,
    volatility = sigma,
    dt = dt,
    lambda = lambda_step
  )

  tidy_common <- function(rates, discounts, measure, weights) {
    rates_tidy <- matrix_to_tidy(rates, times, value_name = "short_rate")
    discounts_tidy <- matrix_to_tidy(discounts, times, value_name = "discount_factor")
    out <- dplyr::left_join(
      rates_tidy,
      discounts_tidy,
      by = c("path_id", "time")
    ) |>
      dplyr::mutate(measure = measure, .before = path_id)
    if (!is.null(weights)) {
      weights_tidy <- matrix_to_tidy(weights, times, value_name = "radon_nikodym")
      out <- dplyr::left_join(out, weights_tidy, by = c("path_id", "time"))
    } else {
      out$radon_nikodym <- 1
    }
    out
  }

  tidy_p <- tidy_common(evolved$rates_p, evolved$discounts_p, "P", evolved$radon)
  tidy_q <- tidy_common(evolved$rates_q, evolved$discounts_q, "Q", NULL)

  result <- dplyr::bind_rows(tidy_p, tidy_q)
  grid <- tibble::tibble(time = step_times, lambda = lambda_step)

  attr(result, "market_price_of_risk") <- grid
  attr(result, "process_type") <- "short_rate"
  attr(result, "spec") <- spec
  class(result) <- unique(c("short_rate_measure_paths", "measure_paths", class(result)))
  result
}
