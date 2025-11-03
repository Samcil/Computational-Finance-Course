#' Price Fixed-Income Instruments
#'
#' Generic pricing dispatcher. Specific model objects provide concrete methods.
#'
#' @param object Model object.
#' @param ... Forwarded to methods.
#'
#' @return Model-specific value.
#'
#' @export
price_zcb <- function(object, ...) {
  UseMethod("price_zcb")
}


#' @export
price_zcb.default <- function(object, ...) {
  rlang::abort(
    message = "No price_zcb() method exists for this object",
    class = "price_zcb_method_not_implemented"
  )
}
