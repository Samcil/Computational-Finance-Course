#' Deterministic Quadrature Utilities
#'
#' Implements trapezoidal and Simpson's rule approximations for one-dimensional
#' integrals. The helpers mirror the Python lecture utilities used to benchmark
#' Monte Carlo integrators by providing deterministic reference values.
#'
#' @param g Function to integrate. Must accept a numeric vector and return a
#'   numeric vector of the same length.
#' @param lower Numeric lower bound.
#' @param upper Numeric upper bound with `upper > lower`.
#' @param n_steps Integer number of subintervals.
#'
#' @return Numeric approximation to the integral $\int_{lower}^{upper} g(x) dx$.
#' @export
#' @examples
#' trapezoidal_integral(sin, 0, pi, n_steps = 200)
trapezoidal_integral <- function(g, lower, upper, n_steps) {
  checkmate::assert_function(g)
  checkmate::assert_number(lower, finite = TRUE)
  checkmate::assert_number(upper, finite = TRUE)
  if (upper <= lower) {
    rlang::abort("`upper` must be strictly greater than `lower`.")
  }
  checkmate::assert_integerish(n_steps, lower = 1, len = 1)

  n_steps <- as.integer(n_steps)
  h <- (upper - lower) / n_steps
  grid <- seq(lower, upper, length.out = n_steps + 1)
  values <- g(grid)
  checkmate::assert_numeric(values, len = n_steps + 1, any.missing = FALSE)

  h * (0.5 * values[1] + sum(values[2:n_steps]) + 0.5 * values[n_steps + 1])
}

#' @rdname trapezoidal_integral
#' @export
#' @examples
#' simpson_integral(exp, 0, 1, n_steps = 200)
simpson_integral <- function(g, lower, upper, n_steps) {
  checkmate::assert_function(g)
  checkmate::assert_number(lower, finite = TRUE)
  checkmate::assert_number(upper, finite = TRUE)
  if (upper <= lower) {
    rlang::abort("`upper` must be strictly greater than `lower`.")
  }
  checkmate::assert_integerish(n_steps, lower = 2, len = 1)

  n_steps <- as.integer(n_steps)
  if (n_steps %% 2 == 1) {
    rlang::abort("Simpson's rule requires an even number of steps.")
  }

  h <- (upper - lower) / n_steps
  grid <- seq(lower, upper, length.out = n_steps + 1)
  values <- g(grid)
  checkmate::assert_numeric(values, len = n_steps + 1, any.missing = FALSE)

  odd_idx <- seq(2, n_steps, by = 2)
  even_idx <- seq(3, n_steps - 1, by = 2)
  h / 3 * (values[1] + values[n_steps + 1] + 4 * sum(values[odd_idx]) + 2 * sum(values[even_idx]))
}

#' Benchmark Deterministic Integrator Accuracy
#'
#' Runs trapezoidal and Simpson's rule over a set of step counts and compares
#' the estimates against a supplied ground truth value.
#'
#' @param g Function to integrate.
#' @param lower,upper Integration bounds.
#' @param step_grid Integer vector of step counts to evaluate.
#' @param exact Numeric ground truth for the integral.
#'
#' @return Tibble with columns `method`, `n_steps`, `estimate`, `abs_error`, and
#'   `log10_error`.
#' @export
#' @examples
#' benchmark_deterministic_integrators(
#'   g = sin,
#'   lower = 0,
#'   upper = pi,
#'   step_grid = c(20, 40, 80, 160),
#'   exact = 2
#' )
benchmark_deterministic_integrators <- function(g,
                                                lower,
                                                upper,
                                                step_grid,
                                                exact) {
  checkmate::assert_function(g)
  checkmate::assert_number(lower, finite = TRUE)
  checkmate::assert_number(upper, finite = TRUE)
  if (upper <= lower) {
    rlang::abort("`upper` must be strictly greater than `lower`.")
  }
  checkmate::assert_integerish(step_grid, lower = 2, any.missing = FALSE)
  checkmate::assert_number(exact, finite = TRUE)

  purrr::map_dfr(
    .x = unique(as.integer(step_grid)),
    .f = function(n) {
      trap_est <- trapezoidal_integral(g, lower, upper, n)
      simp_est <- if (n %% 2 == 0) {
        simpson_integral(g, lower, upper, n)
      } else {
        NA_real_
      }

      tibble::tibble(
        method = c("trapezoidal", "simpson"),
        n_steps = n,
        estimate = c(trap_est, simp_est)
      )
    }
  ) |>
    dplyr::filter(!is.na(estimate)) |>
    dplyr::mutate(
      abs_error = abs(estimate - exact),
      log10_error = log10(abs_error)
    )
}
