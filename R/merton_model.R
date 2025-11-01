#' Create Merton Jump-Diffusion Specification
#'
#' Defines the parameters for the Merton (1976) jump-diffusion model where the
#' logarithmic price follows a diffusion with normally distributed jumps.
#'
#' @param initial_price Numeric. Initial asset price \eqn{S_0 > 0}.
#' @param risk_free_rate Numeric. Continuously compounded risk-free rate \eqn{r}.
#' @param volatility Numeric. Diffusive volatility \eqn{\sigma > 0}.
#' @param jump_intensity Numeric. Jump arrival intensity \eqn{\lambda > 0}.
#' @param jump_mean Numeric. Mean of the logarithmic jump size \eqn{\mu_J}.
#' @param jump_sd Numeric. Standard deviation of the logarithmic jump size \eqn{\sigma_J > 0}.
#'
#' @return An object of class `merton_spec` inheriting from `process_spec`.
#'
#' @examples
#' spec <- merton_spec(
#'   initial_price = 100,
#'   risk_free_rate = 0.03,
#'   volatility = 0.2,
#'   jump_intensity = 1.0,
#'   jump_mean = -0.1,
#'   jump_sd = 0.2
#' )
#'
#' spec |>
#'   simulate_paths(n_paths = 100, n_steps = 252, maturity = 1)
#'
#' @export
merton_spec <- function(initial_price,
                        risk_free_rate,
                        volatility,
                        jump_intensity,
                        jump_mean,
                        jump_sd) {
  checkmate::assert_number(initial_price, lower = .Machine$double.eps, finite = TRUE)
  checkmate::assert_number(risk_free_rate, finite = TRUE)
  checkmate::assert_number(volatility, lower = 0, finite = TRUE)
  checkmate::assert_number(jump_intensity, lower = 0, finite = TRUE)
  checkmate::assert_number(jump_mean, finite = TRUE)
  checkmate::assert_number(jump_sd, lower = 0, finite = TRUE)

  structure(
    list(
      initial_price = initial_price,
      risk_free_rate = risk_free_rate,
      volatility = volatility,
      jump_intensity = jump_intensity,
      jump_mean = jump_mean,
      jump_sd = jump_sd
    ),
    class = c("merton_spec", "process_spec")
  )
}

#' @export
print.merton_spec <- function(x, ...) {
  cli::cli_h2("Merton Jump-Diffusion Specification")
  cli::cli_dl(c(
    "Initial price (S0)" = cli::col_cyan("{format(x$initial_price, digits = 6)}"),
    "Risk-free rate (r)" = cli::col_magenta("{format(x$risk_free_rate, digits = 6)}"),
    "Diffusion volatility (\u03c3)" = cli::col_blue("{format(x$volatility, digits = 6)}"),
    "Jump intensity (\u03bb)" = cli::col_green("{format(x$jump_intensity, digits = 6)}"),
    "Jump mean (\u03bc_J)" = cli::col_yellow("{format(x$jump_mean, digits = 6)}"),
    "Jump SD (\u03c3_J)" = cli::col_yellow("{format(x$jump_sd, digits = 6)}")
  ))
  cli::cli_text("")
  cli::cli_alert_info("Use {.fn simulate_paths} to generate jump-diffusion paths")
  invisible(x)
}

#' Simulate Merton Jump-Diffusion Paths
#'
#' Generates Monte Carlo sample paths for the Merton jump-diffusion model using
#' an Euler discretisation with normally distributed jump sizes.
#'
#' @param process_spec A `merton_spec` object.
#' @inheritParams simulate_paths
#'
#' @return A tibble with columns `path_id`, `time`, and `stock_price`.
#'
#' @export
simulate_paths.merton_spec <- function(process_spec,
                                       n_paths,
                                       n_steps,
                                       maturity,
                                       seed = 123,
                                       ...) {
  checkmate::assert_integerish(n_paths, lower = 1, len = 1)
  checkmate::assert_integerish(n_steps, lower = 1, len = 1)
  checkmate::assert_number(maturity, lower = 0, finite = TRUE)
  checkmate::assert_integerish(seed, len = 1)

  n_paths <- as.integer(n_paths)
  n_steps <- as.integer(n_steps)
  seed <- as.integer(seed)

  set.seed(seed)

  dt <- maturity / n_steps
  sqrt_dt <- sqrt(dt)
  time_grid <- seq(0, maturity, length.out = n_steps + 1)

  S0 <- process_spec$initial_price
  r <- process_spec$risk_free_rate
  sigma <- process_spec$volatility
  lambda <- process_spec$jump_intensity
  mu_j <- process_spec$jump_mean
  sigma_j <- process_spec$jump_sd

  kappa <- exp(mu_j + 0.5 * sigma_j^2) - 1

  diffusion_normals <- generate_standardized_normals(n_paths, n_steps)
  diffusion_increments <- sigma * sqrt_dt * diffusion_normals
  drift_increments <- (r - 0.5 * sigma^2 - lambda * kappa) * dt

  jump_counts <- matrix(
    stats::rpois(n_paths * n_steps, lambda * dt),
    nrow = n_paths,
    ncol = n_steps
  )
  jump_totals <- jump_counts * mu_j + sqrt(jump_counts) * sigma_j * matrix(
    stats::rnorm(n_paths * n_steps),
    nrow = n_paths,
    ncol = n_steps
  )

  log_increments <- drift_increments + diffusion_increments + jump_totals
  log_price_matrix <- compute_cumulative_paths(
    initial_value = log(S0),
    increments = log_increments
  )
  stock_matrix <- exp(log_price_matrix)

  paths_tidy <- purrr::map(
    seq_len(n_paths),
    purrr::in_parallel(
      \(path_idx) {
        tibble::tibble(
          path_id = path_idx,
          time = time_grid,
          stock_price = stock_matrix[path_idx, ]
        )
      },
      time_grid = time_grid,
      stock_matrix = stock_matrix
    )
  ) |> purrr::list_rbind()

  attr(paths_tidy, "process_type") <- "merton"
  attr(paths_tidy, "intensity") <- lambda
  attr(paths_tidy, "spec") <- process_spec

  paths_tidy
}

#' Characteristic Function of the Merton Jump-Diffusion Model
#'
#' Computes the characteristic function of the log-return under the risk-neutral
#' measure for the Merton jump-diffusion model.
#'
#' @param process_spec A `merton_spec` object.
#' @param maturity Numeric. Time to maturity in years.
#' @param risk_free_rate Optional override for the risk-free rate.
#' @param dividend_yield Optional override for the continuous dividend yield.
#'
#' @return A function accepting a numeric vector `u` and returning the complex
#'   characteristic function values.
#' @export
merton_characteristic_function <- function(process_spec,
                                           maturity,
                                           risk_free_rate = process_spec$risk_free_rate,
                                           dividend_yield = 0) {
  checkmate::assert_class(process_spec, "merton_spec")
  checkmate::assert_number(maturity, lower = 0, finite = TRUE)
  checkmate::assert_number(risk_free_rate, finite = TRUE)
  checkmate::assert_number(dividend_yield, finite = TRUE)
  sigma <- process_spec$volatility
  lambda <- process_spec$jump_intensity
  mu_j <- process_spec$jump_mean
  sigma_j <- process_spec$jump_sd
  kappa <- exp(mu_j + 0.5 * sigma_j^2) - 1
  drift <- (risk_free_rate - dividend_yield - 0.5 * sigma^2 - lambda * kappa) * maturity
  function(u) {
    u <- as.complex(u)
    i <- 1i
    exp(i * u * drift - 0.5 * sigma^2 * maturity * u^2 +
      lambda * maturity * (exp(i * u * mu_j - 0.5 * sigma_j^2 * u^2) - 1))
  }
}
