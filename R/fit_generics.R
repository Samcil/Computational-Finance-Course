#' Fit a Model Specification
#'
#' Generic S3 function mirroring tidymodels semantics. Specifications created
#' within this package (e.g., [term_structure_spec()]) provide dedicated
#' methods.
#'
#' @param object Model specification object.
#' @param ... Additional arguments passed to methods.
#'
#' @return The result of the dispatched method.
#'
#' @examples
#' \dontrun{
#' fit(gbm_spec(initial_value = 100, drift = 0.05, volatility = 0.2), data)
#' }
#'
#' @export
fit <- function(object, ...) {
  UseMethod("fit")
}


#' @export
fit.default <- function(object, ...) {
  rlang::abort(
    message = "No fit() method exists for this object",
    class = "fit_method_not_implemented"
  )
}
