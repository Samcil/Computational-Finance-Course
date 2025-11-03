# Hybrid Models (Equity/Rate/FX)
#
# This module implements the hybrid specifications and characteristic functions
# that extend the core diffusion, stochastic volatility, and short-rate
# infrastructure. All implementations mirror the Python lecture scripts to
# maintain parity while reusing the shared COS kernel and deterministic
# quadrature utilities available in the package.

hypergeom1F1_neg_half <- function(b, z) {
  if (!requireNamespace("gsl", quietly = TRUE)) {
    rlang::abort("The `gsl` package is required to evaluate the hypergeometric function.")
  }

  z <- as.vector(z)
  b <- if (length(b) == 1) rep(b, length(z)) else as.vector(b)
  if (length(z) == 0) {
    return(numeric())
  }
  if (length(b) != length(z)) {
    rlang::abort("`b` must be scalar or match the length of `z`")
  }

  result <- gsl::hyperg_1F1(-0.5, b, z)
  if (any(!is.finite(result))) {
    rlang::abort("`gsl::hyperg_1F1` returned non-finite values; check input parameters.")
  }
  result
}

mean_cir_sqrt <- function(t, kappa, v0, v_bar, sigma) {
  checkmate::assert_numeric(t, any.missing = FALSE, finite = TRUE)
  out <- numeric(length(t))
  near_zero <- abs(t) < 1e-10
  out[near_zero] <- sqrt(v0)

  if (any(!near_zero)) {
    ts <- t[!near_zero]
    delta <- 4 * kappa * v_bar / (sigma^2)
    c_val <- (sigma^2 / (4 * kappa)) * (1 - exp(-kappa * ts))
    kappa_bar <- 4 * kappa * v0 * exp(-kappa * ts) / (sigma^2 * (1 - exp(-kappa * ts)))
    factor <- sqrt(2 * c_val) * gamma((1 + delta) / 2) / gamma(delta / 2)
    hyper <- hypergeom1F1_neg_half(delta / 2, -0.5 * kappa_bar)
    out[!near_zero] <- factor * hyper
  }

  out
}

trapz_complex <- function(x, y) {
  checkmate::assert_numeric(x, any.missing = FALSE, finite = TRUE)
  checkmate::assert_true(length(x) == length(y))
  if (length(x) < 2) {
    return(0 + 0i)
  }
  dx <- diff(x)
  avg <- (y[-1] + y[-length(y)]) / 2
  sum(dx * avg)
}

require_hull_white <- function(short_rate) {
  if (!identical(short_rate$args$model, "hull_white")) {
    rlang::abort("`short_rate` must be a Hull-White specification")
  }
  short_rate_state(short_rate)
}

# ---------------------------------------------------------------------------
# Heston-Hull-White (H1-HW) specification
# ---------------------------------------------------------------------------

#' Heston-Hull-White Hybrid Specification
#'
#' @inheritParams bshw_spec
#' @param initial_variance Initial Heston variance.
#' @param variance_mean_reversion Heston mean reversion (kappa).
#' @param variance_long_term Long-term variance level (theta).
#' @param variance_volatility Volatility of variance (gamma).
#' @param correlation_eq_var Correlation between equity and variance drivers.
#' @param correlation_eq_rate Correlation between equity and short-rate drivers.
#'
#' @return A `h1_hw_spec` inheriting from `hybrid_spec`.
#' @export
h1_hw_spec <- function(spot,
                       short_rate,
                       initial_variance,
                       variance_mean_reversion,
                       variance_long_term,
                       variance_volatility,
                       correlation_eq_var,
                       correlation_eq_rate) {
  checkmate::assert_number(spot, lower = .Machine$double.eps, finite = TRUE)
  checkmate::assert_class(short_rate, "short_rate_spec")
  state <- require_hull_white(short_rate)
  checkmate::assert_number(initial_variance, lower = 0, finite = TRUE)
  checkmate::assert_number(variance_mean_reversion, lower = .Machine$double.eps, finite = TRUE)
  checkmate::assert_number(variance_long_term, lower = 0, finite = TRUE)
  checkmate::assert_number(variance_volatility, lower = .Machine$double.eps, finite = TRUE)
  checkmate::assert_number(correlation_eq_var, lower = -1, upper = 1, finite = TRUE)
  checkmate::assert_number(correlation_eq_rate, lower = -1, upper = 1, finite = TRUE)

  spec <- new_hybrid_spec(
    class = "h1_hw_spec",
    args = list(
      spot = spot,
      short_rate = short_rate,
      initial_variance = initial_variance,
      variance_mean_reversion = variance_mean_reversion,
      variance_long_term = variance_long_term,
      variance_volatility = variance_volatility,
      correlation_eq_var = correlation_eq_var,
      correlation_eq_rate = correlation_eq_rate
    ),
    process_type = "h1_hw"
  )

  spec <- set_process_metadata(
    spec,
    model_variant = "heston_hull_white",
    state_variables = c("equity_price", "short_rate", "variance"),
    engines = list(price = "cos")
  )

  spec$method$engine_state <- list(
    discount_fun = state$discount_fun,
    theta_fun = state$theta_fun_vec,
    initial_rate = short_rate$args$initial_rate,
    mean_reversion = short_rate$args$mean_reversion,
    short_rate_vol = short_rate$args$volatility,
    heston = list(
      initial_variance = initial_variance,
      mean_reversion = variance_mean_reversion,
      long_term = variance_long_term,
      volatility = variance_volatility
    ),
    correlation_eq_var = correlation_eq_var,
    correlation_eq_rate = correlation_eq_rate
  )

  spec
}

#' @export
print.h1_hw_spec <- function(x, ...) {
  cli::cli_h2("Heston-Hull-White Specification")
  cli::cli_dl(c(
    "Spot" = cli::col_cyan(format(x$spot, digits = 6)),
    "Variance (v0)" = cli::col_cyan(format(x$initial_variance, digits = 6)),
    "Variance mean reversion" = cli::col_cyan(format(x$variance_mean_reversion, digits = 6)),
    "Variance long-term level" = cli::col_cyan(format(x$variance_long_term, digits = 6)),
    "Variance volatility" = cli::col_cyan(format(x$variance_volatility, digits = 6)),
    "rho_{xv}" = cli::col_cyan(format(x$correlation_eq_var, digits = 6)),
    "rho_{xr}" = cli::col_cyan(format(x$correlation_eq_rate, digits = 6)),
    "Short-rate mean reversion" = cli::col_cyan(format(x$short_rate$args$mean_reversion, digits = 6)),
    "Short-rate volatility" = cli::col_cyan(format(x$short_rate$args$volatility, digits = 6))
  ))
  cli::cli_text("")
  cli::cli_text("{.strong Lineage:} {format_spec_lineage(x)}")
  invisible(x)
}

#' Heston-Hull-White Characteristic Function
#'
#' Mirrors `ChFH1HWModel` from the Python lectures.
#'
#' @param u Numeric integration vector.
#' @param maturity Time to maturity.
#' @param theta_fun Hull-White theta function.
#' @param initial_rate Initial short rate.
#' @param mean_reversion Hull-White mean reversion.
#' @param short_rate_vol Hull-White volatility (eta).
#' @param initial_variance Heston initial variance.
#' @param variance_mean_reversion Heston kappa.
#' @param variance_long_term Heston theta.
#' @param variance_volatility Heston gamma.
#' @param correlation_eq_var rho_{xv}.
#' @param correlation_eq_rate rho_{xr}.
#' @param integration_points Trapezoid steps for deterministic integrals.
#'
#' @return Complex vector of characteristic function values.
#' @export
h1_hw_characteristic_function <- function(u,
                                          maturity,
                                          theta_fun,
                                          initial_rate,
                                          mean_reversion,
                                          short_rate_vol,
                                          initial_variance,
                                          variance_mean_reversion,
                                          variance_long_term,
                                          variance_volatility,
                                          correlation_eq_var,
                                          correlation_eq_rate,
                                          integration_points = 2000L) {
  checkmate::assert_numeric(u, any.missing = FALSE)
  checkmate::assert_number(maturity, lower = 0, finite = TRUE)
  if (maturity <= 0) {
    return(rep(1 + 0i, length(u)))
  }

  checkmate::assert_function(theta_fun)
  checkmate::assert_number(initial_rate, finite = TRUE)
  checkmate::assert_number(mean_reversion, lower = .Machine$double.eps, finite = TRUE)
  checkmate::assert_number(short_rate_vol, lower = 0, finite = TRUE)
  checkmate::assert_number(initial_variance, lower = 0, finite = TRUE)
  checkmate::assert_number(variance_mean_reversion, lower = .Machine$double.eps, finite = TRUE)
  checkmate::assert_number(variance_long_term, lower = 0, finite = TRUE)
  checkmate::assert_number(variance_volatility, lower = .Machine$double.eps, finite = TRUE)
  checkmate::assert_number(correlation_eq_var, lower = -1, upper = 1, finite = TRUE)
  checkmate::assert_number(correlation_eq_rate, lower = -1, upper = 1, finite = TRUE)
  checkmate::assert_integerish(integration_points, lower = 1, len = 1)

  integration_points <- as.integer(integration_points)
  lambda <- mean_reversion
  eta <- short_rate_vol
  kappa <- variance_mean_reversion
  theta_v <- variance_long_term
  gamma_v <- variance_volatility
  rho_xv <- correlation_eq_var
  rho_xr <- correlation_eq_rate

  compute_raw_cf <- function(arg_u) {
    iu <- 1i * arg_u
    u_sq <- arg_u^2

    d1 <- sqrt((kappa - gamma_v * rho_xv * iu)^2 + (u_sq + iu) * gamma_v^2)
    g <- (kappa - gamma_v * rho_xv * iu - d1) / (kappa - gamma_v * rho_xv * iu + d1)
    exp_neg_d1T <- exp(-d1 * maturity)
    d_term <- (1 - exp_neg_d1T) / (gamma_v^2 * (1 - g * exp_neg_d1T)) * (kappa - gamma_v * rho_xv * iu - d1)
    c_term <- (iu - 1) / lambda * (1 - exp(-lambda * maturity))

    theta_integrand <- function(z) {
      (1 - exp(-lambda * z)) * theta_fun(pmax(maturity - z, 0))
    }
    i1_adj <- (iu - 1) * trapezoidal_integral(theta_integrand, 0, maturity, integration_points)

    log_term <- log((1 - g * exp_neg_d1T) / (1 - g))
    i2 <- (maturity / gamma_v^2) * (kappa - gamma_v * rho_xv * iu - d1) - (2 / gamma_v^2) * log_term

    i3 <- ((1i + arg_u)^2) / (2 * lambda^3) * (3 + exp(-2 * lambda * maturity) - 4 * exp(-lambda * maturity) - 2 * lambda * maturity)

    mean_sqrt <- function(t) {
      mean_cir_sqrt(t, kappa = kappa, v0 = initial_variance, v_bar = theta_v, sigma = gamma_v)
    }

    i4 <- -(1 / lambda) * (iu + u_sq) * trapezoidal_integral(
      function(z) mean_sqrt(pmax(maturity - z, 0)) * (1 - exp(-lambda * z)),
      0,
      maturity,
      integration_points
    )

    a_term <- i1_adj + kappa * theta_v * i2 + 0.5 * eta^2 * i3 + eta * rho_xr * i4

    exp(a_term + c_term * initial_rate + d_term * initial_variance)
  }

  raw_values <- compute_raw_cf(u)
  if (any(u == 0)) {
    zero_val <- raw_values[which(u == 0)[1]]
  } else {
    zero_val <- compute_raw_cf(0)
  }

  raw_values / zero_val
}

#' Price H1-HW Options via COS Method
#'
#' @inheritParams price_bshw_option_cos
#' @export
price_h1_hw_option_cos <- function(process_spec,
                                   strikes,
                                   maturity,
                                   n_terms = 512L,
                                   truncation = 10,
                                   integration_points = 2000L) {
  checkmate::assert_class(process_spec, "h1_hw_spec")
  checkmate::assert_numeric(strikes, lower = .Machine$double.eps, any.missing = FALSE, finite = TRUE)
  checkmate::assert_number(maturity, lower = 0, finite = TRUE)
  checkmate::assert_integerish(n_terms, lower = 1, len = 1)
  checkmate::assert_number(truncation, lower = 0, finite = TRUE)
  checkmate::assert_integerish(integration_points, lower = 1, len = 1)

  state <- process_spec$method$engine_state
  discount_factor <- state$discount_fun(maturity)

  cf <- function(u) {
    h1_hw_characteristic_function(
      u = u,
      maturity = maturity,
      theta_fun = state$theta_fun,
      initial_rate = state$initial_rate,
      mean_reversion = state$mean_reversion,
      short_rate_vol = state$short_rate_vol,
      initial_variance = state$heston$initial_variance,
      variance_mean_reversion = state$heston$mean_reversion,
      variance_long_term = state$heston$long_term,
      variance_volatility = state$heston$volatility,
      correlation_eq_var = state$correlation_eq_var,
      correlation_eq_rate = state$correlation_eq_rate,
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

# ---------------------------------------------------------------------------
# Schoebel-Zhu Hull-White specification
# ---------------------------------------------------------------------------

szhw_c_term <- function(u, tau, lambda) {
  (1i * u - 1) / lambda * (1 - exp(-lambda * tau))
}

szhw_d_term <- function(u, tau, kappa, rho_xsigma, gamma) {
  a0 <- -0.5 * u * (1i + u)
  a1 <- 2 * (gamma * rho_xsigma * 1i * u - kappa)
  a2 <- 2 * gamma^2
  d_val <- sqrt(a1^2 - 4 * a0 * a2)
  g_val <- (-a1 - d_val) / (-a1 + d_val)
  exp_neg_dt <- exp(-d_val * tau)
  (-a1 - d_val) / (2 * a2 * (1 - g_val * exp_neg_dt)) * (1 - exp_neg_dt)
}

szhw_e_term <- function(u,
                        tau,
                        lambda,
                        gamma,
                        rho_xsigma,
                        rho_rsigma,
                        rho_xr,
                        eta,
                        kappa,
                        sigma_bar) {
  a0 <- -0.5 * u * (1i + u)
  a1 <- 2 * (gamma * rho_xsigma * 1i * u - kappa)
  a2 <- 2 * gamma^2
  d_val <- sqrt(a1^2 - 4 * a0 * a2)
  g_val <- (-a1 - d_val) / (-a1 + d_val)
  c1 <- gamma * rho_xsigma * 1i * u - kappa - 0.5 * (a1 + d_val)

  denom <- 1 - g_val * exp(-d_val * tau)
  exp_c1_tau <- exp(c1 * tau)

  f1 <- (1 - exp(-c1 * tau)) / c1 + (exp(-(c1 + d_val) * tau) - 1) / (c1 + d_val)
  f2 <- (1 - exp(-c1 * tau)) / c1 + (exp(-(c1 + lambda) * tau) - 1) / (c1 + lambda)
  f3 <- (exp(-(c1 + d_val) * tau) - 1) / (c1 + d_val) + (1 - exp(-(c1 + d_val + lambda) * tau)) / (c1 + d_val + lambda)
  f4 <- 1 / c1 - 1 / (c1 + d_val) - 1 / (c1 + lambda) + 1 / (c1 + d_val + lambda)
  f5 <- exp(-(c1 + d_val + lambda) * tau) * (
    exp(lambda * tau) * (1 / (c1 + d_val) - exp(d_val * tau) / c1) +
      exp(d_val * tau) / (c1 + lambda) - 1 / (c1 + d_val + lambda)
  )

  a2_inv <- 1 / a2
  term1 <- kappa * sigma_bar * a2_inv * (-a1 - d_val) * f1
  term2 <- eta * rho_xr * 1i * u * (1i * u - 1) / lambda * (f2 + g_val * f3)
  term3 <- -rho_rsigma * eta * gamma * a2_inv / lambda * (a1 + d_val) * (1i * u - 1) * (f4 + f5)

  exp_c1_tau / denom * (term1 + term2 + term3)
}

szhw_a_term <- function(u,
                        tau,
                        lambda,
                        eta,
                        gamma,
                        rho_xsigma,
                        rho_rsigma,
                        rho_xr,
                        kappa,
                        sigma_bar,
                        integration_points) {
  a0 <- -0.5 * u * (1i + u)
  a1 <- 2 * (gamma * rho_xsigma * 1i * u - kappa)
  a2 <- 2 * gamma^2
  d_val <- sqrt(a1^2 - 4 * a0 * a2)
  g_val <- (-a1 - d_val) / (-a1 + d_val)
  exp_neg_dt <- exp(-d_val * tau)

  f6 <- eta^2 / (4 * lambda^3) * (1i + u)^2 * (3 + exp(-2 * lambda * tau) - 4 * exp(-lambda * tau) - 2 * lambda * tau)
  a1_base <- 0.25 * ((-a1 - d_val) * tau - 2 * log((1 - g_val * exp_neg_dt) / (1 - g_val))) + f6

  grid <- seq(0, tau, length.out = integration_points + 1)
  e_vals <- szhw_e_term(u, grid, lambda, gamma, rho_xsigma, rho_rsigma, rho_xr, eta, kappa, sigma_bar)
  c_vals <- szhw_c_term(u, grid, lambda)
  integrand <- (kappa * sigma_bar + 0.5 * gamma^2 * e_vals + gamma * eta * rho_rsigma * c_vals) * e_vals
  integral <- trapz_complex(grid, integrand)

  a1_base + integral
}

#' Schoebel-Zhu Hull-White Specification
#'
#' @inheritParams h1_hw_spec
#' @param initial_volatility Initial instantaneous volatility.
#' @param vol_mean_reversion Volatility mean reversion (kappa).
#' @param vol_long_term Long-term volatility level.
#' @param vol_volatility Volatility of volatility (gamma).
#' @param correlation_eq_vol Correlation between equity and volatility drivers.
#' @param correlation_rate_vol Correlation between rate and volatility drivers.
#' @param correlation_eq_rate Correlation between equity and rate drivers.
#'
#' @return A `szhw_spec` inheriting from `hybrid_spec`.
#' @export
szhw_spec <- function(spot,
                      short_rate,
                      initial_volatility,
                      vol_mean_reversion,
                      vol_long_term,
                      vol_volatility,
                      correlation_eq_vol,
                      correlation_rate_vol,
                      correlation_eq_rate) {
  checkmate::assert_number(spot, lower = .Machine$double.eps, finite = TRUE)
  checkmate::assert_class(short_rate, "short_rate_spec")
  state <- require_hull_white(short_rate)
  checkmate::assert_number(initial_volatility, lower = 0, finite = TRUE)
  checkmate::assert_number(vol_mean_reversion, lower = .Machine$double.eps, finite = TRUE)
  checkmate::assert_number(vol_long_term, lower = 0, finite = TRUE)
  checkmate::assert_number(vol_volatility, lower = .Machine$double.eps, finite = TRUE)
  checkmate::assert_number(correlation_eq_vol, lower = -1, upper = 1, finite = TRUE)
  checkmate::assert_number(correlation_rate_vol, lower = -1, upper = 1, finite = TRUE)
  checkmate::assert_number(correlation_eq_rate, lower = -1, upper = 1, finite = TRUE)

  spec <- new_hybrid_spec(
    class = "szhw_spec",
    args = list(
      spot = spot,
      short_rate = short_rate,
      initial_volatility = initial_volatility,
      vol_mean_reversion = vol_mean_reversion,
      vol_long_term = vol_long_term,
      vol_volatility = vol_volatility,
      correlation_eq_vol = correlation_eq_vol,
      correlation_rate_vol = correlation_rate_vol,
      correlation_eq_rate = correlation_eq_rate
    ),
    process_type = "szhw"
  )

  spec <- set_process_metadata(
    spec,
    model_variant = "szhw_hull_white",
    state_variables = c("equity_price", "short_rate", "volatility"),
    engines = list(price = "cos", simulate = "monte_carlo")
  )

  spec$method$engine_state <- list(
    discount_fun = state$discount_fun,
    initial_rate = short_rate$args$initial_rate,
    mean_reversion = short_rate$args$mean_reversion,
    theta_fun = state$theta_fun_vec,
    short_rate_vol = short_rate$args$volatility,
    ou = list(
      initial = initial_volatility,
      mean_reversion = vol_mean_reversion,
      long_term = vol_long_term,
      volatility = vol_volatility
    ),
    correlation_eq_vol = correlation_eq_vol,
    correlation_rate_vol = correlation_rate_vol,
    correlation_eq_rate = correlation_eq_rate
  )

  spec
}

#' @export
simulate_paths.szhw_spec <- function(process_spec,
                                     n_paths,
                                     n_steps,
                                     maturity,
                                     seed = sample.int(.Machine$integer.max, 1),
                                     ...) {
  checkmate::assert_class(process_spec, "szhw_spec")
  checkmate::assert_int(n_paths, lower = 1)
  checkmate::assert_int(n_steps, lower = 1)
  checkmate::assert_number(maturity, lower = 0)
  checkmate::assert_int(seed, lower = 0)

  if (maturity == 0) {
    return(tibble::tibble(
      path_id = seq_len(n_paths),
      time = 0,
      spot = rep(process_spec$spot, n_paths),
      short_rate = rep(process_spec$short_rate$args$initial_rate, n_paths),
      volatility = rep(process_spec$initial_volatility, n_paths),
      money_market = rep(1, n_paths),
      discount_factor = rep(1, n_paths)
    ))
  }

  state <- process_spec$method$engine_state
  dt <- maturity / n_steps
  sqrt_dt <- sqrt(dt)
  time_grid <- seq(0, maturity, length.out = n_steps + 1)

  lambda <- state$mean_reversion
  eta <- state$short_rate_vol
  theta_fun <- state$theta_fun

  kappa <- state$ou$mean_reversion
  sigma_bar <- state$ou$long_term
  gamma <- state$ou$volatility

  rho_xsigma <- state$correlation_eq_vol
  rho_rsigma <- state$correlation_rate_vol
  rho_xr <- state$correlation_eq_rate

  corr_matrix <- matrix(
    c(
      1, rho_xsigma, rho_xr,
      rho_xsigma, 1, rho_rsigma,
      rho_xr, rho_rsigma, 1
    ),
    nrow = 3,
    byrow = TRUE
  )
  if (any(abs(eigen(corr_matrix, symmetric = TRUE, only.values = TRUE)$values) < 1e-10)) {
    rlang::abort("Correlation matrix must be positive definite")
  }
  chol_factor <- chol(corr_matrix)

  random_terms <- with_random_seed(
    seed,
    list(
      x = generate_standardized_normals(n_paths, n_steps),
      sigma = generate_standardized_normals(n_paths, n_steps),
      r = generate_standardized_normals(n_paths, n_steps)
    )
  )

  correlated_steps <- purrr::map(
    seq_len(n_steps),
    function(step_idx) {
      base_step <- cbind(
        random_terms$x[, step_idx],
        random_terms$sigma[, step_idx],
        random_terms$r[, step_idx]
      )
      base_step %*% t(chol_factor)
    }
  )

  initial_state <- list(
    log_spot = rep(log(process_spec$spot), n_paths),
    sigma = rep(state$ou$initial, n_paths),
    rate = rep(state$initial_rate, n_paths),
    money_market = rep(1, n_paths)
  )

  path_states <- purrr::accumulate(
    seq_len(n_steps),
    function(prev, step_idx) {
      increments <- correlated_steps[[step_idx]]
      dWx <- sqrt_dt * increments[, 1]
      dWsigma <- sqrt_dt * increments[, 2]
      dWr <- sqrt_dt * increments[, 3]

      sigma_prev <- prev$sigma
      rate_prev <- prev$rate
      log_prev <- prev$log_spot
      mm_prev <- prev$money_market
      t_prev <- time_grid[step_idx]

      sigma_next <- sigma_prev + kappa * (sigma_bar - sigma_prev) * dt + gamma * dWsigma

      theta_t <- theta_fun(t_prev)
      rate_next <- rate_prev + lambda * (theta_t - rate_prev) * dt + eta * dWr

      money_market_next <- mm_prev * exp(0.5 * (rate_prev + rate_next) * dt)

      log_next_raw <- log_prev + (rate_prev - 0.5 * sigma_prev^2) * dt + sigma_prev * dWx
      adjustment <- process_spec$spot / mean(exp(log_next_raw) / money_market_next)
      log_next <- log_next_raw + log(adjustment)

      list(
        log_spot = log_next,
        sigma = sigma_next,
        rate = rate_next,
        money_market = money_market_next
      )
    },
    .init = initial_state
  )

  log_spot_paths <- do.call(cbind, purrr::map(path_states, "log_spot"))
  sigma_paths <- do.call(cbind, purrr::map(path_states, "sigma"))
  rate_paths <- do.call(cbind, purrr::map(path_states, "rate"))
  money_market_paths <- do.call(cbind, purrr::map(path_states, "money_market"))

  discount_paths <- purrr::accumulate(
    seq_len(n_steps),
    function(prev, step_idx) {
      prev * exp(-0.5 * (rate_paths[, step_idx] + rate_paths[, step_idx + 1]) * dt)
    },
    .init = rep(1, n_paths)
  )
  discount_paths <- do.call(cbind, discount_paths)

  tibble::tibble(
    path_id = rep(seq_len(n_paths), each = length(time_grid)),
    time = rep(time_grid, times = n_paths),
    spot = as.vector(t(exp(log_spot_paths))),
    short_rate = as.vector(t(rate_paths)),
    volatility = as.vector(t(sigma_paths)),
    money_market = as.vector(t(money_market_paths)),
    discount_factor = as.vector(t(discount_paths))
  )
}

#' @export
print.szhw_spec <- function(x, ...) {
  cli::cli_h2("Schoebel-Zhu Hull-White Specification")
  cli::cli_dl(c(
    "Spot" = cli::col_cyan(format(x$spot, digits = 6)),
    "Initial volatility" = cli::col_cyan(format(x$initial_volatility, digits = 6)),
    "Volatility mean reversion" = cli::col_cyan(format(x$vol_mean_reversion, digits = 6)),
    "Volatility long-term level" = cli::col_cyan(format(x$vol_long_term, digits = 6)),
    "Volatility of volatility" = cli::col_cyan(format(x$vol_volatility, digits = 6)),
    "rho_xsigma" = cli::col_cyan(format(x$correlation_eq_vol, digits = 6)),
    "rho_rsigma" = cli::col_cyan(format(x$correlation_rate_vol, digits = 6)),
    "R_{xr}" = cli::col_cyan(format(x$correlation_eq_rate, digits = 6))
  ))
  cli::cli_text("")
  cli::cli_text("{.strong Lineage:} {format_spec_lineage(x)}")
  invisible(x)
}

#' Schoebel-Zhu Hull-White Characteristic Function
#'
#' @inheritParams h1_hw_characteristic_function
#' @param discount_fun Discount function of the short-rate model.
#' @param initial_volatility Initial level of the OU volatility factor.
#' @param vol_mean_reversion OU mean reversion.
#' @param vol_long_term OU long-run level.
#' @param vol_volatility OU volatility of volatility.
#' @param correlation_eq_vol Correlation between equity and volatility drivers.
#' @param correlation_rate_vol Correlation between short rate and volatility drivers.
#' @param correlation_eq_rate Correlation between equity and short rate drivers.
#' @export
szhw_characteristic_function <- function(u,
                                         maturity,
                                         discount_fun,
                                         initial_rate,
                                         mean_reversion,
                                         short_rate_vol,
                                         initial_volatility,
                                         vol_mean_reversion,
                                         vol_long_term,
                                         vol_volatility,
                                         correlation_eq_vol,
                                         correlation_rate_vol,
                                         correlation_eq_rate,
                                         integration_points = 2000L) {
  if (maturity <= 0) {
    return(rep(1 + 0i, length(u)))
  }

  lambda <- mean_reversion
  eta <- short_rate_vol
  kappa <- vol_mean_reversion
  sigma_bar <- vol_long_term
  gamma <- vol_volatility
  rho_xsigma <- correlation_eq_vol
  rho_rsigma <- correlation_rate_vol
  rho_xr <- correlation_eq_rate

  compute_single <- function(ui) {
    d_val <- szhw_d_term(ui, maturity, kappa, rho_xsigma, gamma)
    e_val <- szhw_e_term(ui, maturity, lambda, gamma, rho_xsigma, rho_rsigma, rho_xr, eta, kappa, sigma_bar)
    a_val <- szhw_a_term(ui, maturity, lambda, eta, gamma, rho_xsigma, rho_rsigma, rho_xr, kappa, sigma_bar, integration_points)

    hlp <- eta^2 / (2 * lambda^2) * (
      maturity + (2 / lambda) * (exp(-lambda * maturity) - 1) -
        (1 / (2 * lambda)) * (exp(-2 * lambda * maturity) - 1)
    )
    correction <- (1i * ui - 1) * (log(1 / discount_fun(maturity)) + hlp) + szhw_c_term(ui, maturity, lambda) * initial_rate

    exp((initial_volatility^2) * d_val + initial_volatility * e_val + a_val + correction)
  }

  raw_values <- vapply(u, compute_single, complex(1))
  if (any(u == 0)) {
    zero_val <- raw_values[which(u == 0)[1]]
  } else {
    zero_val <- compute_single(0)
  }

  raw_values / zero_val
}

#' Price SZHW Options via COS Method
#'
#' @inheritParams price_h1_hw_option_cos
#' @param process_spec A `szhw_spec`.
#' @export
price_szhw_option_cos <- function(process_spec,
                                  strikes,
                                  maturity,
                                  n_terms = 512L,
                                  truncation = 10,
                                  integration_points = 2000L) {
  checkmate::assert_class(process_spec, "szhw_spec")
  checkmate::assert_numeric(strikes, lower = .Machine$double.eps, any.missing = FALSE, finite = TRUE)
  checkmate::assert_number(maturity, lower = 0, finite = TRUE)
  checkmate::assert_integerish(n_terms, lower = 1, len = 1)
  checkmate::assert_number(truncation, lower = 0, finite = TRUE)
  checkmate::assert_integerish(integration_points, lower = 1, len = 1)

  state <- process_spec$method$engine_state
  discount_factor <- state$discount_fun(maturity)

  cf <- function(u) {
    szhw_characteristic_function(
      u = u,
      maturity = maturity,
      discount_fun = state$discount_fun,
      initial_rate = state$initial_rate,
      mean_reversion = state$mean_reversion,
      short_rate_vol = state$short_rate_vol,
      initial_volatility = state$ou$initial,
      vol_mean_reversion = state$ou$mean_reversion,
      vol_long_term = state$ou$long_term,
      vol_volatility = state$ou$volatility,
      correlation_eq_vol = state$correlation_eq_vol,
      correlation_rate_vol = state$correlation_rate_vol,
      correlation_eq_rate = state$correlation_eq_rate,
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

#' SZHW Diversification Value via Monte Carlo
#'
#' Replicates the diversification payoff from the `SZHW_Diversification.py`
#' lecture by simulating correlated stock, volatility, and short-rate paths
#' under the Schoebel-Zhu Hull-White dynamics and discounting the terminal
#' payoff back to today.
#'
#' @param process_spec An object created by [szhw_spec()].
#' @param omega Numeric vector of diversification weights between the equity
#'   and bond legs.
#' @param maturity Option maturity.
#' @param settlement Forward start settlement date used for the bond leg.
#' @param n_paths Number of Monte Carlo paths.
#' @param n_steps Number of time steps for the Euler scheme.
#' @param seed Random seed for reproducibility. Defaults to a random draw.
#'
#' @return Tibble with `omega` and the present value of the diversification
#'   payoff.
#'
#' @export
szhw_diversification_value <- function(process_spec,
                                       omega,
                                       maturity,
                                       settlement,
                                       n_paths = 5000L,
                                       n_steps = 360L,
                                       seed = sample.int(.Machine$integer.max, 1)) {
  checkmate::assert_class(process_spec, "szhw_spec")
  checkmate::assert_numeric(omega, any.missing = FALSE, finite = TRUE)
  checkmate::assert_number(maturity, lower = 0)
  checkmate::assert_number(settlement, lower = maturity)
  checkmate::assert_int(n_paths, lower = 1)
  checkmate::assert_int(n_steps, lower = 1)
  checkmate::assert_int(seed, lower = 0)

  mc_paths <- simulate_paths(
    process_spec = process_spec,
    n_paths = n_paths,
    n_steps = n_steps,
    maturity = maturity,
    seed = seed
  )

  terminal_slice <- mc_paths |>
    dplyr::filter(abs(time - maturity) < 1e-10)

  short_rate_spec <- process_spec$args$short_rate
  short_rate_state <- short_rate_state(short_rate_spec)
  discount_settlement <- short_rate_state$discount_fun(settlement)

  zero_coupon_T_T1 <- price_zcb(
    object = short_rate_spec,
    maturities = rep(settlement, nrow(terminal_slice)),
    valuation_time = maturity,
    short_rate = terminal_slice$short_rate
  )

  payoff_component <- function(weight) {
    bond_leg <- zero_coupon_T_T1 / discount_settlement
    equity_leg <- terminal_slice$spot / process_spec$spot
    payoff <- weight * equity_leg + (1 - weight) * bond_leg
    mean(pmax(payoff, 0) / terminal_slice$money_market)
  }

  tibble::tibble(
    omega = omega,
    present_value = vapply(omega, payoff_component, numeric(1))
  )
}

# ---------------------------------------------------------------------------
# FX Heston-Hull-White specification
# ---------------------------------------------------------------------------

fx_c_term <- function(u, tau, kappa, gamma, rho_xv) {
  a0 <- -0.5 * u * (1i + u)
  a1 <- 2 * (gamma * rho_xv * 1i * u - kappa)
  a2 <- 2 * gamma^2
  d_val <- sqrt(a1^2 - 4 * a0 * a2)
  g_val <- (kappa - gamma * rho_xv * 1i * u - d_val) / (kappa - gamma * rho_xv * 1i * u + d_val)
  (1 - exp(-d_val * tau)) / (gamma^2 * (1 - g_val * exp(-d_val * tau))) * (kappa - gamma * rho_xv * 1i * u - d_val)
}

fx_bd <- function(t, T, lambda) {
  (exp(-lambda * (T - t)) - 1) / lambda
}

fx_characteristic_component <- function(u,
                                        maturity,
                                        domestic_state,
                                        foreign_state,
                                        heston_state,
                                        correlations,
                                        integration_points) {
  vapply(
    u,
    function(ui) {
      lambda_d <- domestic_state$mean_reversion
      lambda_f <- foreign_state$mean_reversion
      eta_d <- domestic_state$volatility
      eta_f <- foreign_state$volatility

      kappa <- heston_state$mean_reversion
      theta_v <- heston_state$long_term
      gamma <- heston_state$volatility
      v0 <- heston_state$initial_variance

      rho_xv <- correlations$eq_var
      rho_xrd <- correlations$eq_domestic
      rho_xrf <- correlations$eq_foreign
      rho_vrd <- correlations$var_domestic
      rho_vrf <- correlations$var_foreign
      rho_rdrf <- correlations$domestic_foreign

      c_val <- fx_c_term(ui, maturity, kappa, gamma, rho_xv)

      z <- seq(0, maturity, length.out = integration_points + 1)

      g_vals <- mean_cir_sqrt(maturity - z, kappa, v0, theta_v, gamma)
      bd_vals <- fx_bd(maturity - z, maturity, lambda_d)
      bf_vals <- fx_bd(maturity - z, maturity, lambda_f)

      temp1 <- kappa * theta_v + rho_vrd * gamma * eta_d * g_vals * bd_vals
      temp2 <- -rho_vrd * gamma * eta_d * g_vals * bd_vals * 1i * ui
      temp3 <- rho_vrf * gamma * eta_f * g_vals * bf_vals * 1i * ui
      c_z <- fx_c_term(ui, z, kappa, gamma, rho_xv)
      int1 <- trapz_complex(z, (temp1 + temp2 + temp3) * c_z)

      zeta_vals <- (
        (rho_xrd * eta_d * bd_vals - rho_xrf * eta_f * bf_vals) * g_vals +
          rho_rdrf * eta_d * eta_f * bd_vals * bf_vals -
          0.5 * (eta_d^2 * bd_vals^2 + eta_f^2 * bf_vals^2)
      )
      int2 <- (ui^2 + 1i * ui) * trapz_complex(z, zeta_vals)

      exp(int1 + int2 + v0 * c_val)
    },
    complex(1)
  )
}

#' FX Heston-Hull-White Specification
#'
#' @param spot_fx Spot FX rate (domestic per foreign unit).
#' @param domestic_short_rate Hull-White domestic short-rate spec.
#' @param foreign_short_rate Hull-White foreign short-rate spec.
#' @inheritParams h1_hw_spec
#' @param correlation_eq_domestic Correlation between equity and domestic rate drivers.
#' @param correlation_eq_foreign Correlation between equity and foreign rate drivers.
#' @param correlation_var_domestic Correlation between variance and domestic rate drivers.
#' @param correlation_var_foreign Correlation between variance and foreign rate drivers.
#' @param correlation_domestic_foreign Correlation between domestic and foreign rate drivers.
#'
#' @return A `fx_h1_hw_spec` inheriting from `hybrid_spec`.
#' @export
fx_h1_hw_spec <- function(spot_fx,
                          domestic_short_rate,
                          foreign_short_rate,
                          initial_variance,
                          variance_mean_reversion,
                          variance_long_term,
                          variance_volatility,
                          correlation_eq_var,
                          correlation_eq_domestic,
                          correlation_eq_foreign,
                          correlation_var_domestic,
                          correlation_var_foreign,
                          correlation_domestic_foreign) {
  checkmate::assert_number(spot_fx, lower = .Machine$double.eps, finite = TRUE)
  checkmate::assert_class(domestic_short_rate, "short_rate_spec")
  checkmate::assert_class(foreign_short_rate, "short_rate_spec")
  d_state <- require_hull_white(domestic_short_rate)
  f_state <- require_hull_white(foreign_short_rate)

  checkmate::assert_number(initial_variance, lower = 0, finite = TRUE)
  checkmate::assert_number(variance_mean_reversion, lower = .Machine$double.eps, finite = TRUE)
  checkmate::assert_number(variance_long_term, lower = 0, finite = TRUE)
  checkmate::assert_number(variance_volatility, lower = .Machine$double.eps, finite = TRUE)
  checkmate::assert_number(correlation_eq_var, lower = -1, upper = 1, finite = TRUE)
  checkmate::assert_number(correlation_eq_domestic, lower = -1, upper = 1, finite = TRUE)
  checkmate::assert_number(correlation_eq_foreign, lower = -1, upper = 1, finite = TRUE)
  checkmate::assert_number(correlation_var_domestic, lower = -1, upper = 1, finite = TRUE)
  checkmate::assert_number(correlation_var_foreign, lower = -1, upper = 1, finite = TRUE)
  checkmate::assert_number(correlation_domestic_foreign, lower = -1, upper = 1, finite = TRUE)

  spec <- new_hybrid_spec(
    class = "fx_h1_hw_spec",
    args = list(
      spot_fx = spot_fx,
      domestic_short_rate = domestic_short_rate,
      foreign_short_rate = foreign_short_rate,
      initial_variance = initial_variance,
      variance_mean_reversion = variance_mean_reversion,
      variance_long_term = variance_long_term,
      variance_volatility = variance_volatility,
      correlation_eq_var = correlation_eq_var,
      correlation_eq_domestic = correlation_eq_domestic,
      correlation_eq_foreign = correlation_eq_foreign,
      correlation_var_domestic = correlation_var_domestic,
      correlation_var_foreign = correlation_var_foreign,
      correlation_domestic_foreign = correlation_domestic_foreign
    ),
    process_type = "fx_h1_hw"
  )

  spec <- set_process_metadata(
    spec,
    model_variant = "fx_heston_hull_white",
    state_variables = c("fx_rate", "domestic_rate", "foreign_rate", "variance"),
    engines = list(price = "cos")
  )

  spec$method$engine_state <- list(
    domestic = list(
      discount_fun = d_state$discount_fun,
      theta_fun = d_state$theta_fun_vec,
      initial_rate = domestic_short_rate$args$initial_rate,
      mean_reversion = domestic_short_rate$args$mean_reversion,
      volatility = domestic_short_rate$args$volatility
    ),
    foreign = list(
      discount_fun = f_state$discount_fun,
      theta_fun = f_state$theta_fun_vec,
      initial_rate = foreign_short_rate$args$initial_rate,
      mean_reversion = foreign_short_rate$args$mean_reversion,
      volatility = foreign_short_rate$args$volatility
    ),
    heston = list(
      initial_variance = initial_variance,
      mean_reversion = variance_mean_reversion,
      long_term = variance_long_term,
      volatility = variance_volatility
    ),
    correlations = list(
      eq_var = correlation_eq_var,
      eq_domestic = correlation_eq_domestic,
      eq_foreign = correlation_eq_foreign,
      var_domestic = correlation_var_domestic,
      var_foreign = correlation_var_foreign,
      domestic_foreign = correlation_domestic_foreign
    )
  )

  spec
}

#' @export
print.fx_h1_hw_spec <- function(x, ...) {
  cli::cli_h2("FX Heston-Hull-White Specification")
  cli::cli_dl(c(
    "Spot FX" = cli::col_cyan(format(x$spot_fx, digits = 6)),
    "Variance (v0)" = cli::col_cyan(format(x$initial_variance, digits = 6)),
    "Variance mean reversion" = cli::col_cyan(format(x$variance_mean_reversion, digits = 6)),
    "Variance long-term level" = cli::col_cyan(format(x$variance_long_term, digits = 6)),
    "Variance volatility" = cli::col_cyan(format(x$variance_volatility, digits = 6))
  ))
  cli::cli_text("")
  cli::cli_text("{.strong Lineage:} {format_spec_lineage(x)}")
  invisible(x)
}

#' FX Heston-Hull-White Characteristic Function
#'
#' @param u Numeric integration vector.
#' @param maturity Time to maturity.
#' @param state Internal state returned by [fx_h1_hw_spec()].
#' @param integration_points Integer number of trapezoidal steps.
#' @return Complex vector of characteristic function evaluations.
#' @export
fx_h1_hw_characteristic_function <- function(u,
                                             maturity,
                                             state,
                                             integration_points = 1500L) {
  if (maturity <= 0) {
    return(rep(1 + 0i, length(u)))
  }

  compute_single <- function(ui) {
    fx_characteristic_component(
      u = ui,
      maturity = maturity,
      domestic_state = state$domestic,
      foreign_state = state$foreign,
      heston_state = state$heston,
      correlations = state$correlations,
      integration_points = integration_points
    )
  }

  raw_values <- vapply(u, compute_single, complex(1))
  if (any(u == 0)) {
    zero_val <- raw_values[which(u == 0)[1]]
  } else {
    zero_val <- compute_single(0)
  }

  raw_values / zero_val
}

#' Price FX H1-HW Options via COS Method
#'
#' @param process_spec An `fx_h1_hw_spec` specification.
#' @inheritParams price_h1_hw_option_cos
#' @param integration_points Integer number of points for deterministic integrals.
#'
#' @return Tibble with present value call and put prices (domestic currency).
#' @export
price_fx_h1_hw_option_cos <- function(process_spec,
                                      strikes,
                                      maturity,
                                      n_terms = 512L,
                                      truncation = 8,
                                      integration_points = 1500L) {
  checkmate::assert_class(process_spec, "fx_h1_hw_spec")
  checkmate::assert_numeric(strikes, lower = .Machine$double.eps, any.missing = FALSE, finite = TRUE)
  checkmate::assert_number(maturity, lower = 0, finite = TRUE)
  checkmate::assert_integerish(n_terms, lower = 1, len = 1)
  checkmate::assert_number(truncation, lower = 0, finite = TRUE)
  checkmate::assert_integerish(integration_points, lower = 1, len = 1)

  state <- process_spec$method$engine_state
  discount_d <- state$domestic$discount_fun(maturity)
  discount_f <- state$foreign$discount_fun(maturity)
  forward <- process_spec$spot_fx * discount_f / discount_d

  cf <- function(u) {
    fx_h1_hw_characteristic_function(
      u = u,
      maturity = maturity,
      state = state,
      integration_points = integration_points
    )
  }

  call_prices <- cos_call_put_price_stoch_ir(
    cf = cf,
    option_type = "call",
    spot = forward,
    maturity = maturity,
    strikes = strikes,
    discount_factor = discount_d,
    n_terms = n_terms,
    truncation = truncation
  )

  put_prices <- cos_call_put_price_stoch_ir(
    cf = cf,
    option_type = "put",
    spot = forward,
    maturity = maturity,
    strikes = strikes,
    discount_factor = discount_d,
    n_terms = n_terms,
    truncation = truncation
  )

  tibble::tibble(
    strike = strikes,
    call_price = call_prices$price,
    put_price = put_prices$price
  )
}
