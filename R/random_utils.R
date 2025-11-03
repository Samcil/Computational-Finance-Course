#' Random Number Utilities
#'
#' Internal helpers for controlling random number generation consistently
#' across simulation engines.
#'
#' @keywords internal
with_random_seed <- function(seed, code) {
  checkmate::assert_integerish(seed, len = 1, any.missing = FALSE)
  withr::with_seed(as.integer(seed), code)
}
