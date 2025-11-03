#' Convexity-Corrected Forward Rates
#'
#' Compute convexity-adjusted forward rates and their money-market present
#' values for FRA-style cash flows under Gaussian short-rate specifications.
#' The implementation mirrors the lecture script `ConvexityCorrection.py` while
#' reusing deterministic integrals already exposed by [`short_rate_spec`]
#' utilities.
#'
#' @param spec A [`short_rate_spec`] created by [short_rate_spec()].
#' @param start Numeric vector of forward accrual start times (in years).
#' @param end Numeric vector of forward accrual end times (in years). Each end
#'   point must exceed the corresponding start.
#'
#' @return A tibble with columns
#'   \describe{
#'     \item{start}{Start of the forward accrual period.}
#'     \item{end}{End of the forward accrual period.}
#'     \item{tenor}{Accrual length `end - start`.}
#'     \item{discount_factor}{Discount factor `P(0, start)` implied by the
#'       specification.}
#'     \item{forward_rate}{Forward `L(0; start, end)` ignoring convexity.}
#'     \item{adjusted_forward}{Convexity-adjusted forward implied by Hull-White
#'       analytics.}
#'     \item{convexity_adjustment}{Difference between adjusted and naive
#'       forwards.}
#'     \item{present_value}{Money-market discounted expectation of the forward
#'       cash flow.}
#'   }
#'
#' @examples
#' curve <- tibble::tibble(tenor = 0:5, discount_factor = exp(-0.02 * tenor))
#' spec <- short_rate_spec(
#'   model = "hull_white",
#'   volatility = 0.01,
#'   mean_reversion = 0.1,
#'   curve = curve
#' )
#' convexity_correction_forward_rate(spec, start = 1, end = 2)
#'
#' @export
convexity_correction_forward_rate <- function(spec, start, end) {
  checkmate::assert_class(spec, "short_rate_spec")
  checkmate::assert_numeric(start, lower = 0, finite = TRUE, any.missing = FALSE)
  checkmate::assert_numeric(end, lower = 0, finite = TRUE, any.missing = FALSE)

  if (length(start) == 0 || length(end) == 0) {
    rlang::abort("`start` and `end` must contain at least one value")
  }

  size <- max(length(start), length(end))
  start <- rep(start, length.out = size)
  end <- rep(end, length.out = size)

  if (any(end <= start + 1e-12)) {
    rlang::abort("Each `end` must be greater than the corresponding `start`")
  }

  args <- spec$args
  sigma <- args$volatility
  a <- args$mean_reversion
  initial_rate <- args$initial_rate

  state <- short_rate_state(spec)
  theta_fun <- state$theta_fun
  theta_fun_vec <- state$theta_fun_vec

  tenor <- end - start

  discount_start <- price_zcb(spec, start)
  discount_end <- price_zcb(spec, end)
  forward_naive <- (discount_start / discount_end - 1) / tenor

  if (sigma < 1e-12) {
    adjusted_forward <- forward_naive
    convexity_adj <- rep(0, length(tenor))
    present_value <- discount_start * adjusted_forward
    return(build_convexity_output(start, end, tenor, discount_start, forward_naive, adjusted_forward, convexity_adj, present_value))
  }

  theta_integral_start <- purrr::map_dbl(start, ~ theta_integral(theta_fun, a, 0, .x))
  int_theta_start <- purrr::map_dbl(start, ~ safe_theta_integral(theta_fun_vec, .x))

  idx <- seq_along(start)
  short_rate_mean <- purrr::map_dbl(idx, ~ expected_short_rate(
    mean_reversion = a,
    initial_rate = initial_rate,
    theta_int = theta_integral_start[.x],
    integral_theta = int_theta_start[.x],
    time = start[.x]
  ))

  integral_mean <- purrr::map_dbl(idx, ~ expected_integral_mean(
    mean_reversion = a,
    initial_rate = initial_rate,
    theta_int = theta_integral_start[.x],
    time = start[.x]
  ))

  short_rate_variance <- purrr::map_dbl(start, ~ short_rate_variance_gaussian(a, sigma, .x))
  integral_variance <- purrr::map_dbl(start, ~ sigma^2 * b_squared_integral(a, 0, .x))
  covariance_ir <- purrr::map_dbl(start, ~ covariance_integral_short_rate(a, sigma, .x))

  theta_segment <- purrr::map2_dbl(start, end, ~ theta_integral(theta_fun, a, .x, .y))
  variance_segment <- purrr::map2_dbl(start, end, ~ b_squared_integral(a, .x, .y))
  b_values <- purrr::map_dbl(tenor, ~ b_factor(a, .x))

  prefactor <- exp(theta_segment - 0.5 * sigma^2 * variance_segment)

  y_mean <- -integral_mean + b_values * short_rate_mean
  y_variance <- integral_variance + (b_values^2) * short_rate_variance - 2 * b_values * covariance_ir
  y_variance <- pmax(y_variance, 0)

  pv_high <- prefactor * exp(y_mean + 0.5 * y_variance)
  pv_low <- exp(-integral_mean + 0.5 * integral_variance)

  present_value <- (pv_high - pv_low) / tenor
  adjusted_forward <- present_value / discount_start
  convexity_adj <- adjusted_forward - forward_naive

  build_convexity_output(start, end, tenor, discount_start, forward_naive, adjusted_forward, convexity_adj, present_value)
}

build_convexity_output <- function(start,
                                   end,
                                   tenor,
                                   discount_start,
                                   forward_naive,
                                   adjusted_forward,
                                   convexity_adj,
                                   present_value) {
  tibble::tibble(
    start = start,
    end = end,
    tenor = tenor,
    discount_factor = discount_start,
    forward_rate = forward_naive,
    adjusted_forward = adjusted_forward,
    convexity_adjustment = convexity_adj,
    present_value = present_value
  )
}

safe_theta_integral <- function(theta_fun_vec, time) {
  if (time <= 1e-12) {
    return(0)
  }
  stats::integrate(theta_fun_vec, lower = 0, upper = time, rel.tol = 1e-6, subdivisions = 200)$value
}

expected_short_rate <- function(mean_reversion,
                                initial_rate,
                                theta_int,
                                integral_theta,
                                time) {
  if (time <= 1e-12) {
    return(initial_rate)
  }

  if (mean_reversion > 1e-12) {
    initial_rate * exp(-mean_reversion * time) + integral_theta - mean_reversion * theta_int
  } else {
    initial_rate + integral_theta
  }
}

expected_integral_mean <- function(mean_reversion,
                                   initial_rate,
                                   theta_int,
                                   time) {
  if (time <= 1e-12) {
    return(0)
  }

  if (mean_reversion > 1e-12) {
    initial_rate * (1 - exp(-mean_reversion * time)) / mean_reversion + theta_int
  } else {
    initial_rate * time + theta_int
  }
}

short_rate_variance_gaussian <- function(mean_reversion, sigma, time) {
  if (time <= 1e-12 || sigma <= 0) {
    return(0)
  }

  if (mean_reversion > 1e-12) {
    sigma^2 * (1 - exp(-2 * mean_reversion * time)) / (2 * mean_reversion)
  } else {
    sigma^2 * time
  }
}

covariance_integral_short_rate <- function(mean_reversion, sigma, time) {
  if (time <= 1e-12 || sigma <= 0) {
    return(0)
  }

  if (mean_reversion > 1e-12) {
    sigma^2 * (1 - exp(-mean_reversion * time))^2 / (2 * mean_reversion^2)
  } else {
    sigma^2 * time^2 / 2
  }
}
