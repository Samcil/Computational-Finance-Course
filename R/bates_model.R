#' Create Bates Model Specification
#'
#' Defines the Bates (1996) stochastic volatility model with jumps, combining
#' Heston variance dynamics with Merton-style normally distributed jumps. The
#' specification integrates with `simulate_paths()` and helper utilities to
#' generate paths, evaluate characteristic functions, and recover Black-Scholes
#' implied volatilities.
#'
#' @param initial_price Numeric. Initial asset price \eqn{S_0 > 0}.
#' @param initial_variance Numeric. Initial variance \eqn{v_0 \ge 0}.
#' @param risk_free_rate Numeric. Continuously compounded risk-free rate \eqn{r}.
#' @param dividend_yield Numeric. Continuous dividend yield \eqn{q}. Default is 0.
#' @param mean_reversion Numeric. Variance mean reversion speed \eqn{\kappa > 0}.
#' @param long_term_variance Numeric. Long-run variance level \eqn{\theta > 0}.
#' @param vol_of_vol Numeric. Volatility of variance \eqn{\xi > 0}.
#' @param correlation Numeric. Instantaneous correlation \eqn{\rho \in [-1, 1]}.
#' @param jump_intensity Numeric. Jump intensity \eqn{\lambda \ge 0}.
#' @param jump_mean Numeric. Mean of log jump sizes \eqn{\mu_J}.
#' @param jump_sd Numeric. Standard deviation of log jump sizes \eqn{\sigma_J > 0}.
#' @param scheme Character. Simulation scheme, one of "aes" (default) or "euler".
#'
#' @return An object of class `bates_spec` inheriting from `process_spec`.
#'
#' @examples
#' spec <- bates_spec(
#'   initial_price = 100,
#'   initial_variance = 0.04,
#'   risk_free_rate = 0.02,
#'   dividend_yield = 0.01,
#'   mean_reversion = 1.5,
#'   long_term_variance = 0.04,
#'   vol_of_vol = 0.5,
#'   correlation = -0.6,
#'   jump_intensity = 0.8,
#'   jump_mean = -0.05,
#'   jump_sd = 0.2
#' )
#'
#' spec |>
#'   simulate_paths(n_paths = 50, n_steps = 200, maturity = 1)
#'
#' @export
bates_spec <- function(initial_price,
                       initial_variance,
                       risk_free_rate,
                       dividend_yield = 0,
                       mean_reversion,
                       long_term_variance,
                       vol_of_vol,
                       correlation,
                       jump_intensity,
                       jump_mean,
                       jump_sd,
                       scheme = c("aes", "euler")) {
  scheme <- rlang::arg_match(scheme)
  checkmate::assert_number(initial_price, lower = .Machine$double.eps, finite = TRUE)
  checkmate::assert_number(initial_variance, lower = 0, finite = TRUE)
  checkmate::assert_number(risk_free_rate, finite = TRUE)
  checkmate::assert_number(dividend_yield, finite = TRUE)
  checkmate::assert_number(mean_reversion, lower = 0, finite = TRUE)
  checkmate::assert_number(long_term_variance, lower = 0, finite = TRUE)
  checkmate::assert_number(vol_of_vol, lower = 0, finite = TRUE)
  checkmate::assert_number(correlation, lower = -1, upper = 1, finite = TRUE)
  checkmate::assert_number(jump_intensity, lower = 0, finite = TRUE)
  checkmate::assert_number(jump_mean, finite = TRUE)
  checkmate::assert_number(jump_sd, lower = 0, finite = TRUE)

  feller <- 2 * mean_reversion * long_term_variance
  if (feller < vol_of_vol^2) {
    cli::cli_alert_warning(
      "Feller condition not met: 2\u03ba\u03b8 = {round(feller, 4)} < \u03be\u00b2 = {round(vol_of_vol^2, 4)}"
    )
    cli::cli_alert_info("Variance process may hit zero; Euler simulation reflects at zero, AES sampling remains non-negative")
  }

  structure(
    list(
      initial_price = initial_price,
      initial_variance = initial_variance,
      risk_free_rate = risk_free_rate,
      dividend_yield = dividend_yield,
      mean_reversion = mean_reversion,
      long_term_variance = long_term_variance,
      vol_of_vol = vol_of_vol,
      correlation = correlation,
      jump_intensity = jump_intensity,
      jump_mean = jump_mean,
      jump_sd = jump_sd,
      scheme = scheme
    ),
    class = c("bates_spec", "process_spec")
  )
}

#' @export
print.bates_spec <- function(x, ...) {
  cli::cli_h2("Bates Model Specification")
  cli::cli_dl(c(
    "Initial price (S0)" = cli::col_cyan("{format(x$initial_price, digits = 6)}"),
    "Initial variance (v0)" = cli::col_cyan("{format(x$initial_variance, digits = 6)}"),
    "Risk-free rate (r)" = cli::col_magenta("{format(x$risk_free_rate, digits = 6)}"),
    "Dividend yield (q)" = cli::col_magenta("{format(x$dividend_yield, digits = 6)}"),
    "Mean reversion (kappa)" = cli::col_green("{format(x$mean_reversion, digits = 6)}"),
    "Long-term variance (theta)" = cli::col_green("{format(x$long_term_variance, digits = 6)}"),
    "Vol of vol (xi)" = cli::col_blue("{format(x$vol_of_vol, digits = 6)}"),
    "Correlation (rho)" = cli::col_yellow("{format(x$correlation, digits = 6)}"),
    "Jump intensity (lambda)" = cli::col_green("{format(x$jump_intensity, digits = 6)}"),
    "Jump mean (mu_J)" = cli::col_yellow("{format(x$jump_mean, digits = 6)}"),
    "Jump SD (sigma_J)" = cli::col_yellow("{format(x$jump_sd, digits = 6)}"),
    "Scheme" = cli::col_cyan(x$scheme)
  ))
  cli::cli_text("")
  cli::cli_alert_info("Use {.fn simulate_paths} to generate jump-diffusion volatility paths")
  invisible(x)
}

#' Simulate Bates Model Paths
#'
#' Monte Carlo simulation of the Bates stochastic volatility model using either
#' the Euler-Maruyama discretisation or Andersen's quadratic exponential scheme
#' for the variance process augmented with normally distributed jumps.
#'
#' @param process_spec A `bates_spec` object.
#' @inheritParams simulate_paths
#' @param scheme Optional scheme override ("aes" or "euler"). Defaults to the
#'   scheme stored on the specification.
#'
#' @return A tibble containing `path_id`, `time`, `stock_price`, and `variance`.
#'
#' @export
simulate_paths.bates_spec <- function(process_spec,
                                      n_paths,
                                      n_steps,
                                      maturity,
                                      seed = 123,
                                      scheme = NULL,
                                      ...) {
  checkmate::assert_integerish(n_paths, lower = 1, len = 1)
  checkmate::assert_integerish(n_steps, lower = 1, len = 1)
  checkmate::assert_number(maturity, lower = 0, finite = TRUE)
  checkmate::assert_integerish(seed, len = 1)
  if (!is.null(scheme)) {
    checkmate::assert_choice(scheme, c("aes", "euler"))
  }

  n_paths <- as.integer(n_paths)
  n_steps <- as.integer(n_steps)
  seed <- as.integer(seed)

  selected_scheme <- if (is.null(scheme)) process_spec$scheme else scheme

  if (identical(selected_scheme, "aes")) {
    simulate_bates_aes(process_spec, n_paths, n_steps, maturity, seed)
  } else if (identical(selected_scheme, "euler")) {
    simulate_bates_euler(process_spec, n_paths, n_steps, maturity, seed)
  } else {
    rlang::abort("Unsupported Bates simulation scheme")
  }
}

simulate_bates_euler <- function(process_spec, n_paths, n_steps, maturity, seed) {
  set.seed(seed)

  dt <- maturity / n_steps
  sqrt_dt <- sqrt(dt)
  time_grid <- seq(0, maturity, length.out = n_steps + 1)

  kappa <- process_spec$mean_reversion
  theta <- process_spec$long_term_variance
  xi <- process_spec$vol_of_vol
  rho <- process_spec$correlation
  r <- process_spec$risk_free_rate
  q <- process_spec$dividend_yield
  S0 <- process_spec$initial_price
  v0 <- process_spec$initial_variance
  lambda <- process_spec$jump_intensity
  jump_mean <- process_spec$jump_mean
  jump_sd <- process_spec$jump_sd
  jump_compensation <- exp(jump_mean + 0.5 * jump_sd^2) - 1

  z_var <- generate_standardized_normals(n_paths, n_steps)
  z_independent <- generate_standardized_normals(n_paths, n_steps)
  dW_v <- sqrt_dt * z_var
  dW_s <- rho * dW_v + sqrt(1 - rho^2) * sqrt_dt * z_independent

  jump_counts <- matrix(
    stats::rpois(n_paths * n_steps, lambda * dt),
    nrow = n_paths,
    ncol = n_steps
  )
  jump_totals <- jump_counts * jump_mean + sqrt(jump_counts) * jump_sd * matrix(
    stats::rnorm(n_paths * n_steps),
    nrow = n_paths,
    ncol = n_steps
  )

  variance_history <- purrr::accumulate(
    .x = seq_len(n_steps),
    .init = rep(v0, n_paths),
    .f = \(state, step_idx) {
      variance_prev <- pmax(state, 0)
      drift <- kappa * (theta - variance_prev) * dt
      diffusion <- xi * sqrt(variance_prev) * dW_v[, step_idx]
      pmax(variance_prev + drift + diffusion, 0)
    }
  )
  variance_matrix <- do.call(cbind, variance_history)

  price_history <- purrr::accumulate(
    .x = seq_len(n_steps),
    .init = rep(S0, n_paths),
    .f = \(state, step_idx) {
      variance_prev <- pmax(variance_matrix[, step_idx], 0)
      drift_component <- (r - q - lambda * jump_compensation - 0.5 * variance_prev) * dt
      diffusion_component <- sqrt(variance_prev) * dW_s[, step_idx]
      jump_component <- jump_totals[, step_idx]
      state * exp(drift_component + diffusion_component + jump_component)
    }
  )
  stock_matrix <- do.call(cbind, price_history)

  paths_tidy <- purrr::map(
    seq_len(n_paths),
    purrr::in_parallel(
      \(path_idx) {
        tibble::tibble(
          path_id = path_idx,
          time = time_grid,
          stock_price = stock_matrix[path_idx, ],
          variance = variance_matrix[path_idx, ]
        )
      },
      time_grid = time_grid,
      stock_matrix = stock_matrix,
      variance_matrix = variance_matrix
    )
  ) |> purrr::list_rbind()

  attr(paths_tidy, "process_type") <- "bates"
  attr(paths_tidy, "spec") <- process_spec
  attr(paths_tidy, "scheme") <- "euler"

  paths_tidy
}

simulate_bates_aes <- function(process_spec, n_paths, n_steps, maturity, seed) {
  set.seed(seed)

  dt <- maturity / n_steps
  time_grid <- seq(0, maturity, length.out = n_steps + 1)

  kappa <- process_spec$mean_reversion
  theta <- process_spec$long_term_variance
  xi <- process_spec$vol_of_vol
  rho <- process_spec$correlation
  r <- process_spec$risk_free_rate
  q <- process_spec$dividend_yield
  S0 <- process_spec$initial_price
  v0 <- process_spec$initial_variance
  lambda <- process_spec$jump_intensity
  jump_mean <- process_spec$jump_mean
  jump_sd <- process_spec$jump_sd
  jump_compensation <- exp(jump_mean + 0.5 * jump_sd^2) - 1

  z1 <- matrix(stats::rnorm(n_paths * n_steps), nrow = n_paths, ncol = n_steps)
  z1 <- standardize_columns(z1)
  delta_w1 <- sqrt(dt) * z1

  jump_counts <- matrix(
    stats::rpois(n_paths * n_steps, lambda * dt),
    nrow = n_paths,
    ncol = n_steps
  )
  jump_totals <- jump_counts * jump_mean + sqrt(jump_counts) * jump_sd * matrix(
    stats::rnorm(n_paths * n_steps),
    nrow = n_paths,
    ncol = n_steps
  )

  variance_matrix <- matrix(0, nrow = n_paths, ncol = n_steps + 1)
  variance_matrix[, 1] <- v0
  log_price_matrix <- matrix(log(S0), nrow = n_paths, ncol = n_steps + 1)

  k0 <- (r - q - lambda * jump_compensation - rho * theta * kappa / xi) * dt
  k1 <- (rho * kappa / xi - 0.5) * dt - rho / xi
  k2 <- rho / xi
  correlation_sqrt <- sqrt(pmax(1 - rho^2, 0))

  bates_state <- purrr::reduce(
    .x = seq_len(n_steps),
    .init = list(
      variance_matrix = variance_matrix,
      log_price_matrix = log_price_matrix
    ),
    .f = \(state, step_idx) {
      current_variance <- pmax(state$variance_matrix[, step_idx], 0)
      next_variance <- cir_exact_sample(current_variance, dt, kappa, theta, xi)
      state$variance_matrix[, step_idx + 1] <- next_variance

      diffusion_term <- correlation_sqrt * sqrt(pmax(current_variance, 0)) * delta_w1[, step_idx]
      state$log_price_matrix[, step_idx + 1] <- state$log_price_matrix[, step_idx] +
        k0 + k1 * current_variance + k2 * next_variance + diffusion_term + jump_totals[, step_idx]

      state
    }
  )

  variance_matrix <- bates_state$variance_matrix
  log_price_matrix <- bates_state$log_price_matrix

  stock_matrix <- exp(log_price_matrix)

  paths_tidy <- purrr::map(
    seq_len(n_paths),
    purrr::in_parallel(
      \(path_idx) {
        tibble::tibble(
          path_id = path_idx,
          time = time_grid,
          stock_price = stock_matrix[path_idx, ],
          variance = variance_matrix[path_idx, ]
        )
      },
      time_grid = time_grid,
      stock_matrix = stock_matrix,
      variance_matrix = variance_matrix
    )
  ) |> purrr::list_rbind()

  attr(paths_tidy, "process_type") <- "bates"
  attr(paths_tidy, "spec") <- process_spec
  attr(paths_tidy, "scheme") <- "aes"

  paths_tidy
}

#' Generate Bates Model Paths
#'
#' Convenience wrapper that builds a Bates specification and calls
#' `simulate_paths()` in a single step.
#'
#' @inheritParams bates_spec
#' @inheritParams simulate_paths
#' @param n_paths Integer. Number of Monte Carlo paths.
#' @param n_steps Integer. Number of time steps per path.
#' @param maturity Numeric. Simulation horizon in years.
#' @param seed Integer. Random seed. Default is 123.
#'
#' @return A tibble of simulated Bates paths.
#' @export
generate_bates_paths <- function(n_paths,
                                  n_steps,
                                  maturity,
                                  initial_price,
                                  initial_variance,
                                  risk_free_rate,
                                  dividend_yield = 0,
                                  mean_reversion,
                                  long_term_variance,
                                  vol_of_vol,
                                  correlation,
                                  jump_intensity,
                                  jump_mean,
                                  jump_sd,
                                  scheme = c("aes", "euler"),
                                  seed = 123) {
  scheme <- rlang::arg_match(scheme)
  spec <- bates_spec(
    initial_price = initial_price,
    initial_variance = initial_variance,
    risk_free_rate = risk_free_rate,
    dividend_yield = dividend_yield,
    mean_reversion = mean_reversion,
    long_term_variance = long_term_variance,
    vol_of_vol = vol_of_vol,
    correlation = correlation,
    jump_intensity = jump_intensity,
    jump_mean = jump_mean,
    jump_sd = jump_sd,
    scheme = scheme
  )

  simulate_paths(
    process_spec = spec,
    n_paths = n_paths,
    n_steps = n_steps,
    maturity = maturity,
    seed = seed,
    scheme = scheme
  )
}

#' Characteristic Function of the Bates Model
#'
#' Returns the characteristic function of the log-return under the Bates model
#' with stochastic volatility and jumps.
#'
#' @param process_spec A `bates_spec` object.
#' @param maturity Numeric. Time horizon in years.
#' @param risk_free_rate Optional override for the risk-free rate.
#' @param dividend_yield Optional override for the dividend yield.
#'
#' @return A function of `u` returning complex-valued characteristic
#'   evaluations.
#' @export
bates_characteristic_function <- function(process_spec,
                                          maturity,
                                          risk_free_rate = process_spec$risk_free_rate,
                                          dividend_yield = process_spec$dividend_yield) {
  checkmate::assert_class(process_spec, "bates_spec")
  checkmate::assert_number(maturity, lower = 0, finite = TRUE)
  checkmate::assert_number(risk_free_rate, finite = TRUE)
  checkmate::assert_number(dividend_yield, finite = TRUE)

  lambda <- process_spec$jump_intensity
  jump_mean <- process_spec$jump_mean
  jump_sd <- process_spec$jump_sd
  jump_compensation <- exp(jump_mean + 0.5 * jump_sd^2) - 1

  base_cf <- function(u) {
    heston_characteristic_function(
      u = u,
      maturity = maturity,
      initial_price = process_spec$initial_price,
      initial_variance = process_spec$initial_variance,
      risk_free_rate = risk_free_rate - lambda * jump_compensation,
      dividend_yield = dividend_yield,
      mean_reversion = process_spec$mean_reversion,
      long_term_variance = process_spec$long_term_variance,
      vol_of_vol = process_spec$vol_of_vol,
      correlation = process_spec$correlation
    )
  }

  function(u) {
    u <- as.complex(u)
    base_cf(u) * exp(lambda * maturity * (exp(1i * u * jump_mean - 0.5 * jump_sd^2 * u^2) - 1))
  }
}

#' Bates Implied Volatility Smile
#'
#' Computes Black-Scholes implied volatilities corresponding to Bates model
#' option prices obtained via the COS method.
#'
#' @inheritParams bates_characteristic_function
#' @param strikes Numeric vector of strike prices.
#' @param option_type Character string, either "call" or "put".
#' @param n_terms Integer number of COS terms. Default is 256.
#' @param truncation Numeric truncation width for the COS method. Default is 8.
#' @param vol_interval Numeric length-2 vector giving the search interval for
#'   implied volatility. Default is `c(1e-4, 3)`.
#' @param tol Numeric tolerance passed to `uniroot()`. Default is 1e-6.
#' @param max_iter Integer maximum iterations for `uniroot()`. Default is 100.
#'
#' @return A tibble with columns `strike`, `option_price`, and
#'   `implied_volatility`.
#' @export
bates_implied_volatility <- function(process_spec,
                                     maturity,
                                     strikes,
                                     option_type = c("call", "put"),
                                     n_terms = 256L,
                                     truncation = 8,
                                     vol_interval = c(1e-4, 3),
                                     tol = 1e-6,
                                     max_iter = 100) {
  option_type <- rlang::arg_match(option_type)
  checkmate::assert_class(process_spec, "bates_spec")
  checkmate::assert_number(maturity, lower = 0, finite = TRUE)
  checkmate::assert_numeric(strikes, lower = .Machine$double.eps, any.missing = FALSE, finite = TRUE)
  checkmate::assert_integerish(n_terms, lower = 1, len = 1)
  checkmate::assert_number(truncation, lower = 0, finite = TRUE)
  checkmate::assert_numeric(vol_interval, len = 2, finite = TRUE, any.missing = FALSE)
  checkmate::assert_number(tol, lower = 0, finite = TRUE)
  checkmate::assert_integerish(max_iter, lower = 1, len = 1)

  cf <- bates_characteristic_function(process_spec, maturity, process_spec$risk_free_rate, process_spec$dividend_yield)
  cos_prices <- cos_call_put_price(
    cf = cf,
    option_type = option_type,
    spot = process_spec$initial_price,
    risk_free_rate = process_spec$risk_free_rate,
    maturity = maturity,
    strikes = strikes,
    n_terms = n_terms,
    truncation = truncation
  )

  bs_forward_factor <- exp(-process_spec$dividend_yield * maturity)
  bs_discount_factor <- exp(-process_spec$risk_free_rate * maturity)
  intrinsic_lower <- if (option_type == "call") {
    pmax(process_spec$initial_price * bs_forward_factor - strikes * bs_discount_factor, 0)
  } else {
    pmax(strikes * bs_discount_factor - process_spec$initial_price * bs_forward_factor, 0)
  }
  upper_bound <- if (option_type == "call") {
    process_spec$initial_price * bs_forward_factor
  } else {
    strikes * bs_discount_factor
  }

  implied_vols <- purrr::map2_dbl(
    cos_prices$price,
    seq_along(strikes),
    \(target_price, idx) {
      strike <- strikes[idx]
      lower_price <- intrinsic_lower[idx]
      upper_price <- upper_bound[idx]
      if (target_price < lower_price - 1e-8 || target_price > upper_price + 1e-8) {
        return(NA_real_)
      }

      bs_spec <- black_scholes_spec(
        option_type = option_type,
        strike = strike,
        maturity = maturity,
        risk_free_rate = process_spec$risk_free_rate,
        dividend_yield = process_spec$dividend_yield
      )

      pricing_difference <- function(vol) {
        price_options(bs_spec, spot = process_spec$initial_price, volatility = vol)$price - target_price
      }

      lower <- vol_interval[1]
      upper <- vol_interval[2]
      f_lower <- pricing_difference(lower)
      f_upper <- pricing_difference(upper)

      if (abs(f_lower) < tol) {
        return(lower)
      }
      if (abs(f_upper) < tol) {
        return(upper)
      }

      if (f_lower * f_upper > 0) {
        expansion <- c(2, 5, 10)
        expansion_results <- purrr::map(
          expansion,
          purrr::in_parallel(
            \(mult) {
              candidate <- vol_interval[2] * mult
              f_candidate <- pricing_difference(candidate)
              tibble::tibble(
                candidate = candidate,
                f_candidate = f_candidate,
                should_update = f_lower * f_candidate <= 0
              )
            },
            vol_interval = vol_interval,
            pricing_difference = pricing_difference,
            f_lower = f_lower
          )
        ) |> purrr::list_rbind()

        valid_candidate <- expansion_results |>
          dplyr::filter(should_update) |>
          dplyr::slice_head(n = 1)

        if (nrow(valid_candidate) > 0) {
          upper <- valid_candidate$candidate
          f_upper <- valid_candidate$f_candidate
        }
      }

      if (f_lower * f_upper > 0) {
        return(NA_real_)
      }

      stats::uniroot(pricing_difference, lower = lower, upper = upper, tol = tol, maxiter = max_iter)$root
    }
  )

  tibble::tibble(
    strike = strikes,
    option_price = cos_prices$price,
    implied_volatility = implied_vols
  )
}
