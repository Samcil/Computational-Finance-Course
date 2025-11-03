#' Augment a Fitted Object
#'
#' Lightweight alternative to broom's `augment()` tailored for this package.
#'
#' @param x A fitted object.
#' @param ... Additional arguments passed to methods.
#'
#' @return Typically a tibble with augmented results.
#'
#' @export
augment <- function(x, ...) {
  UseMethod("augment")
}


#' @export
augment.default <- function(x, ...) {
  rlang::abort(
    message = "No augment() method exists for this object",
    class = "augment_method_not_implemented"
  )
}
