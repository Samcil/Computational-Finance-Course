#' Create CIR Process Specification
#'
#' Creates a specification for the Cox-Ingersoll-Ross (CIR) mean-reverting
#' square-root diffusion process.
#'
#' @param initial_value Numeric. Initial variance/rate \eqn{v_0}.
#' @param mean_reversion Numeric. Mean reversion speed \eqn{\kappa > 0}.
#' @param long_term_mean Numeric. Long-term equilibrium level \eqn{\theta > 0}.
#' @param volatility Numeric. Volatility of volatility \eqn{\sigma > 0}.
#'
#' @details
#' ## Mathematical Formulation
#'
#' The CIR process follows:
#' \deqn{dv(t) = \kappa(\theta - v(t))dt + \sigma\sqrt{v(t)}dW(t)}
#'
#' Where:
#' \itemize{
#'   \item \eqn{\kappa}: Mean reversion speed
#'   \item \eqn{\theta}: Long-term mean (equilibrium)
#'   \item \eqn{\sigma}: Volatility of volatility
#'   \item \eqn{v(t)}: Variance/rate at time t
#' }
#'
#' **Feller Condition**: \eqn{2\kappa\theta \geq \sigma^2} ensures \eqn{v(t) > 0}.
#'
#' @return A cir_spec object containing the process parameters.
#'
#' @examples
#' # Interest rate model
#' spec <- cir_spec(0.05, 3.0, 0.04, 0.2)
#'
#' # Pipeline workflow
#' cir_spec(0.05, 3.0, 0.04, 0.2) |>
#'   simulate_paths(50, 252, 1.0) |>
#'   plot_paths()
#'
#' @export
cir_spec <- function(initial_value, mean_reversion, long_term_mean, volatility) {
  # Validation
  checkmate::assert_number(initial_value, lower = 0, finite = TRUE)
  checkmate::assert_number(mean_reversion, lower = 0, finite = TRUE)
  checkmate::assert_number(long_term_mean, lower = 0, finite = TRUE)
  checkmate::assert_number(volatility, lower = 0, finite = TRUE)

  # Feller condition check
  feller <- 2 * mean_reversion * long_term_mean
  if (feller < volatility^2) {
    cli::cli_alert_warning(
      "Feller condition not met: 2\u03ba\u03b8 = {round(feller, 4)} < \u03c3\u00b2 = {round(volatility^2, 4)}"
    )
    cli::cli_alert_info("Process may reach zero with positive probability")
  }

  structure(
    list(
      initial_value = initial_value,
      mean_reversion = mean_reversion,
      long_term_mean = long_term_mean,
      volatility = volatility
    ),
    class = c("cir_spec", "process_spec")
  )
}

#' @export
print.cir_spec <- function(x, ...) {
  cli::cli_h2("CIR Process Specification")
  cli::cli_text("Process: dv(t) = \u03ba(\u03b8 - v(t))dt + \u03c3\u221av(t)dW(t)")
  cli::cli_text("")
  cli::cli_dl(c(
    "Initial value (v\u2080)" = cli::col_cyan("{x$initial_value}"),
    "Mean reversion (\u03ba)" = cli::col_green("{x$mean_reversion}"),
    "Long-term mean (\u03b8)" = cli::col_blue("{x$long_term_mean}"),
    "Volatility (\u03c3)" = cli::col_magenta("{x$volatility}")
  ))
  cli::cli_text("")
  feller <- 2 * x$mean_reversion * x$long_term_mean
  if (feller >= x$volatility^2) {
    cli::cli_alert_success("Feller condition satisfied: 2\u03ba\u03b8 \u2265 \u03c3\u00b2")
  } else {
    cli::cli_alert_warning("Feller condition NOT satisfied")
  }
  cli::cli_text("")
  cli::cli_alert_info("Use {.fn simulate_paths} to generate sample paths")
  invisible(x)
}

#' Simulate CIR Process Paths
#'
#' @param process_spec A cir_spec object.
#' @param n_paths Integer. Number of paths to simulate.
#' @param n_steps Integer. Number of time steps.
#' @param maturity Numeric. Time to maturity in years.
#' @param seed Integer. Random seed for reproducibility.
#'
#' @return A tibble with columns: path_id, time, variance.
#'
#' @export
simulate_paths.cir_spec <- function(process_spec, n_paths, n_steps, maturity, seed = 123, ...) {
  # Validation
  checkmate::assert_integerish(n_paths, lower = 1, len = 1)
  checkmate::assert_integerish(n_steps, lower = 1, len = 1)
  checkmate::assert_number(maturity, lower = 0, finite = TRUE)
  checkmate::assert_integerish(seed, len = 1)

  # Coercion
  n_paths <- as.integer(n_paths)
  n_steps <- as.integer(n_steps)
  seed <- as.integer(seed)

  # Set seed
  set.seed(seed)

  # Time grid
  dt <- maturity / n_steps
  time_grid <- seq(0, maturity, length.out = n_steps + 1)

  # Extract parameters
  v0 <- process_spec$initial_value
  kappa <- process_spec$mean_reversion
  theta <- process_spec$long_term_mean
  sigma <- process_spec$volatility

  # Generate Brownian increments (scaled by sqrt(dt))
  dW <- matrix(
    stats::rnorm(n_paths * n_steps),
    nrow = n_paths,
    ncol = n_steps
  ) * sqrt(dt)

  # Functional Euler-Maruyama with reflection at zero
  variance_history <- purrr::accumulate(
    .x = seq_len(n_steps),
    .init = rep(v0, n_paths),
    .f = \(state, step_idx) {
      v_prev <- pmax(state, 0)
      drift <- kappa * (theta - v_prev) * dt
      diffusion <- sigma * sqrt(v_prev) * dW[, step_idx]
      pmax(v_prev + drift + diffusion, 0)
    }
  )

  variance_matrix <- do.call(cbind, variance_history) |> as.matrix()

  paths_tidy <- matrix_to_tidy(
    value_matrix = variance_matrix,
    time_grid = time_grid,
    value_name = "variance"
  )

  # Add process type attribute
  attr(paths_tidy, "process_type") <- "cir"
  attr(paths_tidy, "long_term_mean") <- theta

  paths_tidy
}
