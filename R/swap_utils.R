#' Validate Swap Cashflow Schedules
#'
#' Internal utility to verify and normalise swap cashflow schedules used by
#' pricing, exposure, and calibration routines. Ensures required columns are
#' present, coerces numeric types, and orders rows by payment time.
#'
#' @keywords internal
validate_swap_schedule <- function(schedule) {
  checkmate::assert_data_frame(schedule, any.missing = FALSE)

  required_cols <- c("start", "end", "pay_time", "accrual_fraction", "notional")
  missing_cols <- setdiff(required_cols, names(schedule))
  if (length(missing_cols) > 0) {
    rlang::abort(
      message = "`schedule` is missing required columns",
      class = "swap_schedule_missing_columns",
      columns = missing_cols
    )
  }

  schedule <- tibble::as_tibble(schedule) |>
    dplyr::mutate(
      start = as.numeric(.data$start),
      end = as.numeric(.data$end),
      pay_time = as.numeric(.data$pay_time),
      accrual_fraction = as.numeric(.data$accrual_fraction),
      notional = as.numeric(.data$notional)
    ) |>
    dplyr::arrange(.data$pay_time)

  if (any(!is.finite(schedule$start)) || any(!is.finite(schedule$end))) {
    rlang::abort("`start` and `end` columns must contain finite values")
  }

  if (any(schedule$end <= schedule$start)) {
    rlang::abort("Each cashflow must satisfy `end` > `start`")
  }

  if (any(!is.finite(schedule$pay_time))) {
    rlang::abort("`pay_time` must contain finite values")
  }

  if (any(schedule$accrual_fraction < 0) || any(!is.finite(schedule$accrual_fraction))) {
    rlang::abort("`accrual_fraction` must be finite and non-negative")
  }

  if (any(!is.finite(schedule$notional))) {
    rlang::abort("`notional` must contain finite values")
  }

  schedule
}
