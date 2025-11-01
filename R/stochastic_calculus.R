#' Simulate Brownian Motion and Associated Integrals
#'
#' Generates Brownian motion paths alongside the Riemann integral \eqn{\int_0^t W(s)\,ds}
#' and the It\u00f4 integral \eqn{\int_0^t g(W(s))\,dW(s)} for a supplied integrand.
#' The routine matches the lecture examples while returning a tidy tibble.
#'
#' @param maturity Numeric time horizon \eqn{T > 0}.
#' @param n_paths Integer number of Monte Carlo paths.
#' @param n_steps Integer number of time steps.
#' @param integrand Function applied to the Brownian motion within the It\u00f4
#'   integral. Defaults to the identity function.
#' @param seed Optional integer seed for reproducibility.
#'
#' @return A tibble with columns `path_id`, `time`, `brownian`,
#'   `riemann_integral`, `ito_integral`, and `quadratic_variation`.
#' @export
simulate_stochastic_integrals <- function(maturity,
                                          n_paths,
                                          n_steps,
                                          integrand = base::identity,
                                          seed = NULL) {
  checkmate::assert_number(maturity, lower = .Machine$double.eps, finite = TRUE)
  checkmate::assert_integerish(n_paths, lower = 1, len = 1)
  checkmate::assert_integerish(n_steps, lower = 1, len = 1)
  checkmate::assert_function(integrand)
  if (!is.null(seed)) {
    checkmate::assert_integerish(seed, len = 1)
    set.seed(as.integer(seed))
  }

  n_paths <- as.integer(n_paths)
  n_steps <- as.integer(n_steps)
  dt <- maturity / n_steps

  dW <- sqrt(dt) * matrix(stats::rnorm(n_paths * n_steps), nrow = n_paths, ncol = n_steps)
  brownian_matrix <- compute_cumulative_paths(0, dW)

  current_brownian <- brownian_matrix[, seq_len(n_steps), drop = FALSE]
  integrand_vals <- integrand(as.vector(current_brownian))
  checkmate::assert_numeric(integrand_vals, len = length(current_brownian), any.missing = FALSE)
  integrand_matrix <- matrix(integrand_vals, nrow = n_paths, ncol = n_steps)

  ito_increments <- integrand_matrix * dW
  ito_matrix <- compute_cumulative_paths(0, ito_increments)

  riemann_increments <- current_brownian * dt
  riemann_matrix <- compute_cumulative_paths(0, riemann_increments)

  quadratic_increments <- dW^2
  quadratic_matrix <- compute_cumulative_paths(0, quadratic_increments)

  time_grid <- seq(0, maturity, length.out = n_steps + 1)

  tibble::tibble(
    path_id = rep(seq_len(n_paths), each = n_steps + 1),
    time = rep(time_grid, times = n_paths),
    brownian = as.vector(t(brownian_matrix)),
    riemann_integral = as.vector(t(riemann_matrix)),
    ito_integral = as.vector(t(ito_matrix)),
    quadratic_variation = as.vector(t(quadratic_matrix))
  )
}

#' Quadratic Variation Diagnostics for Brownian Motion
#'
#' Estimates the first two moments of squared Brownian increments across a grid
#' of discretisations, mirroring the lecture diagnostic.
#'
#' @param maturity Numeric time horizon \eqn{T > 0}.
#' @param n_paths Integer number of Monte Carlo paths.
#' @param step_grid Integer vector of step counts.
#' @param seed Optional integer seed for reproducibility.
#'
#' @return A tibble with columns `n_steps`, `dt`, `mean_sq_increment`,
#'   `var_sq_increment`, `expected_mean`, and `expected_variance`.
#' @export
quadratic_variation_diagnostics <- function(maturity,
                                            n_paths,
                                            step_grid,
                                            seed = NULL) {
  checkmate::assert_number(maturity, lower = .Machine$double.eps, finite = TRUE)
  checkmate::assert_integerish(n_paths, lower = 1, len = 1)
  checkmate::assert_integerish(step_grid, lower = 1, any.missing = FALSE)
  if (!is.null(seed)) {
    checkmate::assert_integerish(seed, len = 1)
  }

  n_paths <- as.integer(n_paths)
  step_grid <- as.integer(step_grid)

  purrr::map2_dfr(
    .x = step_grid,
    .y = seq_along(step_grid),
    .f = \(n_steps, idx) {
      if (!is.null(seed)) {
        set.seed(as.integer(seed) + idx - 1)
      }
      dt <- maturity / n_steps
      dW <- sqrt(dt) * matrix(stats::rnorm(n_paths * n_steps), nrow = n_paths, ncol = n_steps)
      sq_increments <- dW^2
      mean_sq <- mean(sq_increments)
      var_sq <- stats::var(as.vector(sq_increments))

      tibble::tibble(
        n_steps = n_steps,
        dt = dt,
        mean_sq_increment = mean_sq,
        var_sq_increment = var_sq,
        expected_mean = dt,
        expected_variance = 2 * dt^2
      )
    }
  )
}
