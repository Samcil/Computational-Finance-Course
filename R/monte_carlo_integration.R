#' Monte Carlo Integration via Hit-or-Miss Sampling
#'
#' Estimates a definite integral using hit-or-miss Monte Carlo sampling over a
#' rectangular bounding box. The routine follows the lecture example but returns
#' tidy summaries and standard errors.
#'
#' @param g Function representing the integrand. Must accept and return numeric
#'   vectors.
#' @param lower Numeric lower bound of the integration interval.
#' @param upper Numeric upper bound of the integration interval.
#' @param lower_bound Numeric lower bound of the sampling rectangle in the
#'   vertical direction.
#' @param upper_bound Numeric upper bound of the sampling rectangle in the
#'   vertical direction.
#' @param n_samples Integer number of Monte Carlo samples.
#' @param seed Optional integer seed for reproducibility.
#'
#' @return A tibble with columns `method`, `estimate`, `std_error`, and
#'   `n_samples`.
#' @export
mc_integrate_hit_or_miss <- function(g,
                                     lower,
                                     upper,
                                     lower_bound,
                                     upper_bound,
                                     n_samples = 10000L,
                                     seed = NULL) {
  checkmate::assert_function(g)
  checkmate::assert_number(lower, finite = TRUE)
  checkmate::assert_number(upper, finite = TRUE)
  if (upper <= lower) {
    rlang::abort("`upper` must be strictly greater than `lower`.")
  }
  checkmate::assert_number(lower_bound, finite = TRUE)
  checkmate::assert_number(upper_bound, finite = TRUE)
  if (upper_bound <= lower_bound) {
    rlang::abort("`upper_bound` must be strictly greater than `lower_bound`.")
  }
  checkmate::assert_integerish(n_samples, lower = 1, len = 1)
  if (!is.null(seed)) {
    checkmate::assert_integerish(seed, len = 1)
    set.seed(as.integer(seed))
  }

  n_samples <- as.integer(n_samples)
  width <- upper - lower
  height <- upper_bound - lower_bound

  x_samples <- stats::runif(n_samples, lower, upper)
  y_samples <- stats::runif(n_samples, lower_bound, upper_bound)

  g_values <- g(x_samples)
  checkmate::assert_numeric(g_values, len = n_samples, any.missing = FALSE)

  if (any(g_values < lower_bound - sqrt(.Machine$double.eps)) ||
      any(g_values > upper_bound + sqrt(.Machine$double.eps))) {
    rlang::warn("Integrand values exceed the specified bounding box; estimates may be biased.")
  }

  hits <- g_values > y_samples
  hit_rate <- mean(hits)
  estimate <- width * height * hit_rate
  variance <- width^2 * height^2 * hit_rate * (1 - hit_rate) / n_samples
  std_error <- sqrt(variance)

  tibble::tibble(
    method = "hit-or-miss",
    estimate = estimate,
    std_error = std_error,
    n_samples = n_samples
  )
}

#' Monte Carlo Integration via Sample Mean
#'
#' Estimates a definite integral using the sample mean Monte Carlo method,
#' drawing uniform samples over the integration domain.
#'
#' @inheritParams mc_integrate_hit_or_miss
#' @param lower_bound,upper_bound Ignored for this method.
#'
#' @return A tibble with columns `method`, `estimate`, `std_error`, and
#'   `n_samples`.
#' @export
mc_integrate_sample_mean <- function(g,
                                     lower,
                                     upper,
                                     n_samples = 10000L,
                                     seed = NULL) {
  checkmate::assert_function(g)
  checkmate::assert_number(lower, finite = TRUE)
  checkmate::assert_number(upper, finite = TRUE)
  if (upper <= lower) {
    rlang::abort("`upper` must be strictly greater than `lower`.")
  }
  checkmate::assert_integerish(n_samples, lower = 1, len = 1)
  if (!is.null(seed)) {
    checkmate::assert_integerish(seed, len = 1)
    set.seed(as.integer(seed))
  }

  n_samples <- as.integer(n_samples)
  width <- upper - lower

  x_samples <- stats::runif(n_samples, lower, upper)
  g_values <- g(x_samples)
  checkmate::assert_numeric(g_values, len = n_samples, any.missing = FALSE)

  mean_val <- mean(g_values)
  sd_val <- stats::sd(g_values)
  estimate <- width * mean_val
  std_error <- width * sd_val / sqrt(n_samples)

  tibble::tibble(
    method = "sample-mean",
    estimate = estimate,
    std_error = std_error,
    n_samples = n_samples
  )
}
