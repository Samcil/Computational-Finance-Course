#' Barrier Option Pricing via Monte Carlo
#'
#' Prices single-asset barrier options under a Geometric Brownian Motion model
#' using Monte Carlo simulation. Supports up/down and in/out structures for call
#' and put payoffs.
#'
#' @param process_spec A `gbm_spec` object describing the underlying asset
#'   dynamics under the risk-neutral measure.
#' @param strike Numeric vector of strike prices.
#' @param barrier Numeric vector of barrier levels. When both `strike` and
#'   `barrier` have length greater than one, all combinations are evaluated.
#' @param maturity Numeric. Time to maturity in years.
#' @param barrier_type Character. One of "up-and-out", "up-and-in",
#'   "down-and-out", or "down-and-in".
#' @param option_type Character. Option payoff, either "call" or "put".
#' @param n_paths Integer. Number of Monte Carlo paths. Default is 1000.
#' @param n_steps Integer. Number of monitoring points along each path. Default
#'   is 252.
#' @param seed Optional integer seed for reproducibility.
#'
#' @return A tibble containing `strike`, `barrier`, `option_type`,
#'   `barrier_type`, `price`, `std_error`, `hit_probability`, `n_paths`, and
#'   `n_steps`.
#'
#' @examples
#' spec <- gbm_spec(initial_value = 100, drift = 0.02, volatility = 0.2)
#' price_barrier_option(
#'   process_spec = spec,
#'   strike = 100,
#'   barrier = 120,
#'   maturity = 1,
#'   barrier_type = "up-and-out",
#'   option_type = "call",
#'   n_paths = 5000,
#'   n_steps = 200,
#'   seed = 123
#' )
#'
#' @export
price_barrier_option <- function(process_spec,
                                 strike,
                                 barrier,
                                 maturity,
                                 barrier_type = c("up-and-out", "up-and-in", "down-and-out", "down-and-in"),
                                 option_type = c("call", "put"),
                                 n_paths = 1000L,
                                 n_steps = 252L,
                                 seed = NULL) {
  barrier_type <- rlang::arg_match(barrier_type)
  option_type <- rlang::arg_match(option_type)

  checkmate::assert_class(process_spec, "gbm_spec")
  checkmate::assert_numeric(strike, lower = 0, any.missing = FALSE, finite = TRUE)
  checkmate::assert_numeric(barrier, lower = 0, any.missing = FALSE, finite = TRUE)
  checkmate::assert_number(maturity, lower = 0, finite = TRUE)
  checkmate::assert_integerish(n_paths, lower = 1, len = 1)
  checkmate::assert_integerish(n_steps, lower = 1, len = 1)
  if (!is.null(seed)) {
    checkmate::assert_integerish(seed, len = 1)
    set.seed(as.integer(seed))
  }

  path_summary <- summarise_path_extrema(process_spec, n_paths, n_steps, maturity)

  param_grid <- tidyr::expand_grid(strike = strike, barrier = barrier)
  discount_factor <- exp(-process_spec$drift * maturity)
  sample_size <- nrow(path_summary)

  results <- purrr::pmap_dfr(
    param_grid,
    \(strike, barrier) {
      hit <- barrier_hit_indicator(
        path_max = path_summary$path_max,
        path_min = path_summary$path_min,
        barrier = barrier,
        barrier_type = barrier_type
      )
      active <- if (endsWith(barrier_type, "out")) !hit else hit

      intrinsic <- option_intrinsic_payoff(
        terminal_price = path_summary$terminal_price,
        strike = strike,
        option_type = option_type
      )
      payoff <- ifelse(active, intrinsic, 0)

      price_stats(
        payoff = payoff,
        discount_factor = discount_factor,
        hit_indicator = hit,
        barrier = barrier,
        strike = strike,
        option_type = option_type,
        barrier_type = barrier_type,
        sample_size = sample_size,
        n_steps = n_steps
      )
    }
  )

  tibble::as_tibble(results)
}

#' Estimate Barrier Hit Probabilities
#'
#' Computes the probability that a barrier is breached prior to maturity under
#' the specified GBM dynamics using Monte Carlo simulation.
#'
#' @inheritParams price_barrier_option
#'
#' @return A tibble containing `barrier`, `barrier_type`, `probability`,
#'   `n_paths`, and `n_steps`.
#'
#' @export
barrier_hit_probability <- function(process_spec,
                                    barrier,
                                    maturity,
                                    barrier_type = c("up", "down"),
                                    n_paths = 1000L,
                                    n_steps = 252L,
                                    seed = NULL) {
  barrier_type <- rlang::arg_match(barrier_type)
  checkmate::assert_class(process_spec, "gbm_spec")
  checkmate::assert_numeric(barrier, lower = 0, any.missing = FALSE, finite = TRUE)
  checkmate::assert_number(maturity, lower = 0, finite = TRUE)
  checkmate::assert_integerish(n_paths, lower = 1, len = 1)
  checkmate::assert_integerish(n_steps, lower = 1, len = 1)
  if (!is.null(seed)) {
    checkmate::assert_integerish(seed, len = 1)
    set.seed(as.integer(seed))
  }

  path_summary <- summarise_path_extrema(process_spec, n_paths, n_steps, maturity)
  sample_size <- nrow(path_summary)

  purrr::map_dfr(
    barrier,
    \(barrier_level) {
      hit <- if (identical(barrier_type, "up")) {
        path_summary$path_max >= barrier_level
      } else {
        path_summary$path_min <= barrier_level
      }

      probability <- mean(hit)
      tibble::tibble(
        barrier = barrier_level,
        barrier_type = barrier_type,
        probability = probability,
        n_paths = sample_size,
        n_steps = n_steps
      )
    }
  )
}

summarise_path_extrema <- function(process_spec, n_paths, n_steps, maturity) {
  paths <- simulate_paths(
    process_spec = process_spec,
    n_paths = n_paths,
    n_steps = n_steps,
    maturity = maturity
  )

  paths |>
    dplyr::group_by(path_id) |>
    dplyr::summarise(
      terminal_price = dplyr::last(stock_price),
      path_max = max(stock_price),
      path_min = min(stock_price),
      .groups = "drop"
    )
}

barrier_hit_indicator <- function(path_max, path_min, barrier, barrier_type) {
  if (startsWith(barrier_type, "up")) {
    path_max >= barrier
  } else {
    path_min <= barrier
  }
}

option_intrinsic_payoff <- function(terminal_price, strike, option_type) {
  if (identical(option_type, "call")) {
    pmax(terminal_price - strike, 0)
  } else {
    pmax(strike - terminal_price, 0)
  }
}

price_stats <- function(payoff,
                        discount_factor,
                        hit_indicator,
                        barrier,
                        strike,
                        option_type,
                        barrier_type,
                        sample_size,
                        n_steps) {
  discounted_payoff <- discount_factor * payoff
  mean_payoff <- mean(discounted_payoff)
  variance <- stats::var(discounted_payoff)
  std_error <- if (sample_size > 1) {
    sqrt(variance / sample_size)
  } else {
    NA_real_
  }

  tibble::tibble(
    strike = strike,
    barrier = barrier,
    option_type = option_type,
    barrier_type = barrier_type,
    price = mean_payoff,
    std_error = std_error,
    hit_probability = mean(hit_indicator),
    n_paths = sample_size,
    n_steps = n_steps
  )
}
