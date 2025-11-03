#' Create Heston Model Specification
#'
#' Defines the stochastic volatility dynamics introduced by Heston (1993).
#' The specification integrates with `simulate_paths()` to generate joint stock
#' price and variance trajectories in a pipe-friendly workflow.
#'
#' @param initial_price Numeric. Initial asset price \eqn{S_0 > 0}.
#' @param initial_variance Numeric. Initial variance \eqn{v_0 \ge 0}.
#' @param risk_free_rate Numeric. Continuously compounded risk-free rate \eqn{r}.
#' @param dividend_yield Numeric. Continuous dividend yield \eqn{q}. Default is 0.
#' @param mean_reversion Numeric. Mean reversion speed \eqn{\kappa > 0}.
#' @param long_term_variance Numeric. Long-term variance level \eqn{\theta > 0}.
#' @param vol_of_vol Numeric. Volatility of variance \eqn{\xi > 0}.
#' @param correlation Numeric. Instantaneous correlation \eqn{\rho \in [-1, 1]} between
#'   the Brownian motions driving price and variance.
#' @param scheme Character. Simulation scheme, either `"euler"` for Euler-Maruyama
#'   or `"aes"` for the Andersen quadratic exponential scheme.
#'
#' @details
#' The Heston model evolves according to:
#' \deqn{\begin{aligned}
#'   dS(t) &= (r - q) S(t) dt + \sqrt{v(t)} S(t) dW_S(t), \\
#'   dv(t) &= \kappa (\theta - v(t)) dt + \xi \sqrt{v(t)} dW_v(t),
#' \end{aligned}}
#' with correlation \eqn{dW_S(t) dW_v(t) = \rho dt}.
#'
#' **Feller Condition**: \eqn{2\kappa\theta \ge \xi^2} prevents the variance from
#' hitting zero. When the condition is violated the Euler scheme still runs but
#' the variance is reflected at zero during simulation.
#'
#' @return An object of class `heston_spec` inheriting from `model_spec`.
#'
#' @examples
#' spec <- heston_spec(
#'   initial_price = 100,
#'   initial_variance = 0.04,
#'   risk_free_rate = 0.02,
#'   dividend_yield = 0,
#'   mean_reversion = 1.5,
#'   long_term_variance = 0.04,
#'   vol_of_vol = 0.5,
#'   correlation = -0.7
#' )
#'
#' spec |>
#'   simulate_paths(n_paths = 100, n_steps = 252, maturity = 1)
#'
#' @export
heston_spec <- function(initial_price,
                        initial_variance,
                        risk_free_rate,
                        dividend_yield = 0,
                        mean_reversion,
                        long_term_variance,
                        vol_of_vol,
                        correlation,
                        scheme = c("euler", "aes")) {
  scheme <- rlang::arg_match(scheme)
  checkmate::assert_number(initial_price, lower = .Machine$double.eps, finite = TRUE)
  checkmate::assert_number(initial_variance, lower = 0, finite = TRUE)
  checkmate::assert_number(risk_free_rate, finite = TRUE)
  checkmate::assert_number(dividend_yield, finite = TRUE)
  checkmate::assert_number(mean_reversion, lower = 0, finite = TRUE)
  checkmate::assert_number(long_term_variance, lower = 0, finite = TRUE)
  checkmate::assert_number(vol_of_vol, lower = 0, finite = TRUE)
  checkmate::assert_number(correlation, lower = -1, upper = 1, finite = TRUE)
  if (scheme == "aes" && vol_of_vol <= 0) {
    rlang::abort("Volatility of variance must be positive for the AES scheme")
  }

  feller <- 2 * mean_reversion * long_term_variance
  if (feller < vol_of_vol^2) {
    cli::cli_alert_warning(
      "Feller condition not met: 2\u03ba\u03b8 = {round(feller, 4)} < \u03be\u00b2 = {round(vol_of_vol^2, 4)}"
    )
    cli::cli_alert_info("Variance process may hit zero; Euler simulation reflects at zero, AES sampling remains non-negative")
  }

  spec <- new_stochastic_vol_spec(
    class = "heston_spec",
    args = list(
      initial_price = initial_price,
      initial_variance = initial_variance,
      risk_free_rate = risk_free_rate,
      dividend_yield = dividend_yield,
      mean_reversion = mean_reversion,
      long_term_variance = long_term_variance,
      vol_of_vol = vol_of_vol,
      correlation = correlation,
      scheme = scheme
    ),
    process_type = "heston"
  )

  spec <- set_process_metadata(
    spec,
    variance_process = "cir",
    default_scheme = scheme,
    engines = list(
      simulate = c("euler", "aes"),
      price = "cos"
    )
  )

  spec
}

#' @export
print.heston_spec <- function(x, ...) {
  cli::cli_h2("Heston Model Specification")
  cli::cli_text("Process: dS(t) = (r - q) S(t) dt + \u221av(t) S(t) dW_S(t)")
  cli::cli_text("         dv(t) = \u03ba(\u03b8 - v(t)) dt + \u03be \u221av(t) dW_v(t)")
  cli::cli_text("")
  cli::cli_dl(c(
    "Initial price (S0)" = cli::col_cyan("{format(x$initial_price, digits = 6)}"),
    "Initial variance (v0)" = cli::col_cyan("{format(x$initial_variance, digits = 6)}"),
    "Risk-free rate (r)" = cli::col_magenta("{format(x$risk_free_rate, digits = 6)}"),
    "Dividend yield (q)" = cli::col_magenta("{format(x$dividend_yield, digits = 6)}"),
    "Mean reversion (\u03ba)" = cli::col_green("{format(x$mean_reversion, digits = 6)}"),
    "Long-term variance (\u03b8)" = cli::col_green("{format(x$long_term_variance, digits = 6)}"),
    "Vol of vol (\u03be)" = cli::col_blue("{format(x$vol_of_vol, digits = 6)}"),
    "Correlation (\u03c1)" = cli::col_yellow("{format(x$correlation, digits = 6)}"),
    "Scheme" = cli::col_cyan(x$scheme)
  ))
  cli::cli_text("")
  cli::cli_text("{.strong Lineage:} {format_spec_lineage(x)}")
  cli::cli_text("")
  cli::cli_alert_info("Use {.fn simulate_paths} to generate joint price/variance paths")
  invisible(x)
}

#' Simulate Heston Model Paths
#'
#' Generates Monte Carlo sample paths for the Heston stochastic volatility model
#' using either the Euler-Maruyama discretisation or the Andersen Quadratic
#' Exponential scheme (AES) for the variance process.
#'
#' @param process_spec A `heston_spec` object.
#' @inheritParams simulate_paths
#'
#' @return A tibble with columns `path_id`, `time`, `stock_price`, and `variance`.
#'
#' @export
simulate_paths.heston_spec <- function(process_spec,
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
    checkmate::assert_choice(scheme, c("euler", "aes"))
  }

  n_paths <- as.integer(n_paths)
  n_steps <- as.integer(n_steps)
  seed <- as.integer(seed)

  selected_scheme <- if (is.null(scheme)) {
    process_spec$scheme
  } else {
    scheme
  }

  if (identical(selected_scheme, "euler")) {
    simulate_heston_euler(process_spec, n_paths, n_steps, maturity, seed)
  } else if (identical(selected_scheme, "aes")) {
    simulate_heston_aes(process_spec, n_paths, n_steps, maturity, seed)
  } else {
    rlang::abort("Unsupported Heston simulation scheme")
  }
}

simulate_heston_euler <- function(process_spec, n_paths, n_steps, maturity, seed) {
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

  z_var <- generate_standardized_normals(n_paths, n_steps)
  z_independent <- generate_standardized_normals(n_paths, n_steps)
  dW_v <- sqrt_dt * z_var
  dW_s <- rho * dW_v + sqrt(1 - rho^2) * sqrt_dt * z_independent

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
      drift_component <- (r - q - 0.5 * variance_prev) * dt
      diffusion_component <- sqrt(variance_prev) * dW_s[, step_idx]
      state * exp(drift_component + diffusion_component)
    }
  )
  stock_matrix <- do.call(cbind, price_history)

  paths_tidy <- purrr::map(
    seq_len(n_paths),
    \(path_idx) {
      tibble::tibble(
        path_id = path_idx,
        time = time_grid,
        stock_price = stock_matrix[path_idx, ],
        variance = variance_matrix[path_idx, ]
      )
    }
  ) |> purrr::list_rbind()

  attr(paths_tidy, "process_type") <- "heston"
  attr(paths_tidy, "spec") <- process_spec
  attr(paths_tidy, "scheme") <- "euler"

  paths_tidy
}

simulate_heston_aes <- function(process_spec, n_paths, n_steps, maturity, seed) {
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

  if (xi <= 0) {
    rlang::abort("AES scheme requires strictly positive volatility of variance")
  }

  # Brownian increments used for the log-price dynamics
  z1 <- matrix(stats::rnorm(n_paths * n_steps), nrow = n_paths, ncol = n_steps)
  z1 <- standardize_columns(z1)
  delta_w1 <- sqrt(dt) * z1

  variance_matrix <- matrix(0, nrow = n_paths, ncol = n_steps + 1)
  variance_matrix[, 1] <- v0
  log_price_matrix <- matrix(log(S0), nrow = n_paths, ncol = n_steps + 1)

  k0 <- (r - q - rho * theta * kappa / xi) * dt
  k1 <- (rho * kappa / xi - 0.5) * dt - rho / xi
  k2 <- rho / xi
  correlation_sqrt <- sqrt(pmax(1 - rho^2, 0))

  heston_state <- purrr::reduce(
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
        k0 + k1 * current_variance + k2 * next_variance + diffusion_term

      state
    }
  )

  variance_matrix <- heston_state$variance_matrix
  log_price_matrix <- heston_state$log_price_matrix

  stock_matrix <- exp(log_price_matrix)

  paths_tidy <- purrr::map(
    seq_len(n_paths),
    \(path_idx) {
      tibble::tibble(
        path_id = path_idx,
        time = time_grid,
        stock_price = stock_matrix[path_idx, ],
        variance = variance_matrix[path_idx, ]
      )
    }
  ) |> purrr::list_rbind()

  attr(paths_tidy, "process_type") <- "heston"
  attr(paths_tidy, "spec") <- process_spec
  attr(paths_tidy, "scheme") <- "aes"

  paths_tidy
}

cir_exact_sample <- function(current_variance, dt, mean_reversion, long_term_variance, vol_of_vol) {
  kappa <- mean_reversion
  theta <- long_term_variance
  xi <- vol_of_vol

  if (dt <= .Machine$double.eps) {
    return(pmax(current_variance, 0))
  }

  if (kappa <= 0) {
    # Fall back to Euler update when mean reversion is zero
    drift <- kappa * (theta - current_variance) * dt
    diffusion <- xi * sqrt(pmax(current_variance, 0)) * sqrt(dt) * stats::rnorm(length(current_variance))
    return(pmax(current_variance + drift + diffusion, 0))
  }

  c <- (xi^2 * (1 - exp(-kappa * dt))) / (4 * kappa)
  d <- 4 * kappa * theta / (xi^2)
  non_centrality <- 4 * kappa * exp(-kappa * dt) * pmax(current_variance, 0) / (xi^2 * (1 - exp(-kappa * dt)))

  stats::rchisq(length(current_variance), df = d, ncp = non_centrality) * c
}

standardize_columns <- function(matrix_input) {
  column_means <- colMeans(matrix_input)
  centered <- sweep(matrix_input, 2, column_means, FUN = "-")
  column_sd <- sqrt(colSums(centered^2) / pmax(nrow(matrix_input) - 1, 1))
  column_sd[column_sd == 0] <- 1
  sweep(centered, 2, column_sd, FUN = "/")
}

#' Generate Heston Paths with Euler Scheme
#'
#' Convenience wrapper to simulate Heston model paths using the Euler-Maruyama
#' discretisation without manually creating a specification object.
#'
#' @inheritParams heston_spec
#' @inheritParams simulate_paths
#' @param n_paths Integer. Number of Monte Carlo paths.
#' @param n_steps Integer. Number of time steps per path.
#' @param maturity Numeric. Time horizon in years.
#' @param seed Integer. Random seed for reproducibility.
#'
#' @return Tidy tibble of simulated paths containing columns `path_id`, `time`,
#'   `stock_price`, and `variance`.
#' @export
generate_heston_paths_euler <- function(n_paths,
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
                                        seed = 123) {
  spec <- heston_spec(
    initial_price = initial_price,
    initial_variance = initial_variance,
    risk_free_rate = risk_free_rate,
    dividend_yield = dividend_yield,
    mean_reversion = mean_reversion,
    long_term_variance = long_term_variance,
    vol_of_vol = vol_of_vol,
    correlation = correlation,
    scheme = "euler"
  )

  simulate_paths(spec, n_paths = n_paths, n_steps = n_steps, maturity = maturity, seed = seed)
}

#' Generate Heston Paths with Andersen Quadratic Exponential Scheme
#'
#' Simulates Heston model paths using the Andersen (2008) quadratic exponential
#' scheme (AES) which samples the variance process exactly via a noncentral
#' chi-square draw. This approach improves stability when the Feller condition
#' is violated and offers faster convergence for option pricing tasks.
#'
#' @inheritParams generate_heston_paths_euler
#'
#' @return Tidy tibble of simulated paths containing columns `path_id`, `time`,
#'   `stock_price`, and `variance`.
#' @export
generate_heston_paths_aes <- function(n_paths,
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
                                      seed = 123) {
  spec <- heston_spec(
    initial_price = initial_price,
    initial_variance = initial_variance,
    risk_free_rate = risk_free_rate,
    dividend_yield = dividend_yield,
    mean_reversion = mean_reversion,
    long_term_variance = long_term_variance,
    vol_of_vol = vol_of_vol,
    correlation = correlation,
    scheme = "aes"
  )

  simulate_paths(spec, n_paths = n_paths, n_steps = n_steps, maturity = maturity, seed = seed)
}

#' Heston Log-Return Characteristic Function
#'
#' Computes the characteristic function of the Heston (1993) log-return under the
#' risk-neutral measure. The implementation follows the complex-valued closed
#' form presented in the original paper and is vectorised over the argument `u`.
#'
#' @param u Numeric vector of evaluation points.
#' @param maturity Numeric. Time horizon in years.
#' @param initial_price Numeric. Current asset price. Included for interface
#'   consistency; the characteristic function is defined for log-returns and does
#'   not depend on `initial_price`.
#' @param initial_variance Numeric. Initial instantaneous variance.
#' @param risk_free_rate Numeric. Continuously compounded risk-free rate.
#' @param dividend_yield Numeric. Continuous dividend yield.
#' @param mean_reversion Numeric. Mean reversion speed of the variance process.
#' @param long_term_variance Numeric. Long-run variance level.
#' @param vol_of_vol Numeric. Volatility of variance.
#' @param correlation Numeric. Correlation between price and variance Brownian
#'   motions.
#'
#' @return Complex vector containing the characteristic function evaluated at
#'   `u`.
#' @export
heston_characteristic_function <- function(u,
                                           maturity,
                                           initial_price,
                                           initial_variance,
                                           risk_free_rate,
                                           dividend_yield,
                                           mean_reversion,
                                           long_term_variance,
                                           vol_of_vol,
                                           correlation) {
  checkmate::assert_numeric(u, any.missing = FALSE)
  checkmate::assert_number(maturity, lower = 0, finite = TRUE)
  checkmate::assert_number(initial_price, lower = .Machine$double.eps, finite = TRUE)
  checkmate::assert_number(initial_variance, lower = 0, finite = TRUE)
  checkmate::assert_number(risk_free_rate, finite = TRUE)
  checkmate::assert_number(dividend_yield, finite = TRUE)
  checkmate::assert_number(mean_reversion, lower = 0, finite = TRUE)
  checkmate::assert_number(long_term_variance, lower = 0, finite = TRUE)
  checkmate::assert_number(vol_of_vol, lower = .Machine$double.eps, finite = TRUE)
  checkmate::assert_number(correlation, lower = -1, upper = 1, finite = TRUE)

  iu <- 1i * u
  d <- sqrt((correlation * vol_of_vol * iu - mean_reversion)^2 + vol_of_vol^2 * (iu + u^2))
  g <- (mean_reversion - correlation * vol_of_vol * iu - d) /
    (mean_reversion - correlation * vol_of_vol * iu + d)

  exp_term <- exp(iu * (risk_free_rate - dividend_yield) * maturity)
  C <- (mean_reversion * long_term_variance / (vol_of_vol^2)) *
    ((mean_reversion - correlation * vol_of_vol * iu - d) * maturity -
      2 * log((1 - g * exp(-d * maturity)) / (1 - g)))
  D <- ((mean_reversion - correlation * vol_of_vol * iu - d) / (vol_of_vol^2)) *
    ((1 - exp(-d * maturity)) / (1 - g * exp(-d * maturity)))

  exp_term * exp(C + D * initial_variance)
}

#' Price European Options on Heston via COS Method
#'
#' Uses the Fourier Cosine (COS) expansion to price European call and put
#' options under the Heston stochastic volatility model. Requires the companion
#' `cos_method` utilities providing the expansion coefficients.
#'
#' @param process_spec A `heston_spec` object containing the model parameters.
#'   Only the scalar fields are used; the simulation `scheme` is ignored.
#' @param strikes Numeric vector of strike prices.
#' @param maturity Numeric. Option maturity in years.
#' @param n_terms Integer. Number of COS expansion terms. Default is 256.
#' @param truncation Numeric. Truncation width parameter controlling the
#'   integration interval. Default is 8.
#'
#' @return A tibble with columns `strike`, `call_price`, and `put_price`.
#' @export
price_heston_option_cos <- function(process_spec,
                                    strikes,
                                    maturity,
                                    n_terms = 256L,
                                    truncation = 8) {
  checkmate::assert_class(process_spec, "heston_spec")
  checkmate::assert_numeric(strikes, lower = .Machine$double.eps, any.missing = FALSE, finite = TRUE)
  checkmate::assert_number(maturity, lower = 0, finite = TRUE)
  checkmate::assert_integerish(n_terms, lower = 1, len = 1)
  checkmate::assert_number(truncation, lower = 0, finite = TRUE)

  spec <- process_spec
  cf <- function(u) {
    heston_characteristic_function(
      u = u,
      maturity = maturity,
      initial_price = spec$initial_price,
      initial_variance = spec$initial_variance,
      risk_free_rate = spec$risk_free_rate,
      dividend_yield = spec$dividend_yield,
      mean_reversion = spec$mean_reversion,
      long_term_variance = spec$long_term_variance,
      vol_of_vol = spec$vol_of_vol,
      correlation = spec$correlation
    )
  }

  call_prices <- cos_call_put_price(
    cf = cf,
    option_type = "call",
    spot = spec$initial_price,
    risk_free_rate = spec$risk_free_rate,
    maturity = maturity,
    strikes = strikes,
    n_terms = n_terms,
    truncation = truncation
  )

  put_prices <- cos_call_put_price(
    cf = cf,
    option_type = "put",
    spot = spec$initial_price,
    risk_free_rate = spec$risk_free_rate,
    maturity = maturity,
    strikes = strikes,
    n_terms = n_terms,
    truncation = truncation
  )

  tibble::tibble(
    strike = strikes,
    call_price = call_prices$price,
    put_price = put_prices$price
  )
}
