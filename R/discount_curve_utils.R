#' Discount Curve Utility Functions
#'
#' Shared helpers for constructing discount-factor splines and their
#' derivatives. These utilities are used by short-rate specifications and
#' term-structure modules to avoid duplicating spline logic.
#'
#' @keywords internal
make_log_discount_spline <- function(curve_data) {
  stats::splinefun(
    x = curve_data$tenor,
    y = log(curve_data$discount_factor),
    method = "natural"
  )
}

#' @keywords internal
make_discount_function <- function(log_discount_fun) {
  function(t) {
    checkmate::assert_numeric(t, lower = 0, finite = TRUE, any.missing = FALSE)
    t <- pmax(t, 0)
    exp(log_discount_fun(t))
  }
}

#' @keywords internal
make_forward_function <- function(log_discount_fun) {
  function(t) {
    checkmate::assert_numeric(t, lower = 0, finite = TRUE, any.missing = FALSE)
    t <- pmax(t, 0)
    -log_discount_fun(t, deriv = 1)
  }
}

#' @keywords internal
make_forward_derivative_function <- function(log_discount_fun) {
  function(t) {
    checkmate::assert_numeric(t, lower = 0, finite = TRUE, any.missing = FALSE)
    t <- pmax(t, 0)
    -log_discount_fun(t, deriv = 2)
  }
}
