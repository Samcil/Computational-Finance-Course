#' Term Structure Sensitivity Utilities
#'
#' Helpers that convert calibrated term-structure fits into reusable pricing
#' inputs and bump-and-revalue sensitivity diagnostics. These functions mirror
#' the Python `YieldCurveBuildGreeks.py` reference by providing discount
#' interpolators, swap valuation, and quote-level delta estimation.
#'
#' @name term_structure_sensitivity_utils
#' @keywords internal
NULL

#' Build a Discount Function from a Term Structure Fit
#'
#' Construct a discount factor interpolation from a calibrated
#' `term_structure_fit`. The resulting function can be supplied to pricing
#' helpers that expect a continuous discount curve.
#'
#' @param fit A `term_structure_fit` produced by [fit()].
#' @param curve Optional character identifier selecting the curve to extract
#'   when the calibration contains multiple curves.
#' @param method Interpolation scheme. Either "spline" (default) for a natural
#'   cubic spline in log-discount space or "linear" for piecewise-linear
#'   interpolation of discount factors.
#'
#' @return A function accepting numeric maturities (in years) and returning
#'   discount factors.
#'
#' @examples
#' spec <- term_structure_spec()
#' quotes <- tibble::tibble(tenor = c(0.5, 1, 2), quote = c(0.02, 0.021, 0.024))
#' fit_obj <- fit(spec, quotes)
#' df_fun <- term_structure_discount_function(fit_obj)
#' df_fun(c(0.5, 1))
#'
#' @export
term_structure_discount_function <- function(fit,
                                             curve = NULL,
                                             method = c("spline", "linear")) {
  method <- rlang::arg_match(method)
  checkmate::assert_class(fit, "term_structure_fit")

  calibration <- augment(fit)
  if (!"curve" %in% names(calibration)) {
    default_curve <- fit$spec$args$curve_type
    if (is.null(default_curve)) {
      default_curve <- "curve"
    }
    calibration$curve <- default_curve
  }

  if (is.null(curve)) {
    unique_curves <- unique(calibration$curve)
    if (length(unique_curves) != 1) {
      rlang::abort(
        message = "`curve` must be supplied when the calibration includes multiple curves.",
        class = "term_structure_multiple_curves"
      )
    }
    curve <- unique_curves
  } else {
    checkmate::assert_string(curve, min.chars = 1)
  }

  curve_slice <- calibration[calibration$curve == curve, , drop = FALSE]
  if (nrow(curve_slice) == 0) {
    rlang::abort(
      message = paste0("No calibrated data found for curve '", curve, "'."),
      class = "term_structure_unknown_curve"
    )
  }

  curve_slice <- curve_slice |>
    dplyr::select(dplyr::all_of(c("tenor", "discount_factor"))) |>
    dplyr::mutate(
      tenor = as.numeric(.data$tenor),
      discount_factor = as.numeric(.data$discount_factor)
    ) |>
    dplyr::arrange(.data$tenor) |>
    dplyr::distinct(dplyr::across(dplyr::all_of("tenor")), .keep_all = TRUE)

  if (!any(abs(curve_slice$tenor) < 1e-12)) {
    curve_slice <- dplyr::bind_rows(
      tibble::tibble(tenor = 0, discount_factor = 1),
      curve_slice
    ) |>
      dplyr::arrange(.data$tenor)
  } else {
    zero_idx <- which.min(abs(curve_slice$tenor))
    curve_slice$discount_factor[zero_idx] <- 1
  }

  if (nrow(curve_slice) < 2) {
    rlang::abort("At least two tenor pillars are required to build a discount function.")
  }

  if (method == "spline") {
    log_discount_fun <- make_log_discount_spline(curve_slice)
    discount_fun <- make_discount_function(log_discount_fun)
  } else {
    approx_fun <- stats::approxfun(
      x = curve_slice$tenor,
      y = curve_slice$discount_factor,
      method = "linear",
      rule = 2
    )
    discount_fun <- function(t) {
      checkmate::assert_numeric(t, lower = 0, finite = TRUE, any.missing = FALSE)
      approx_fun(pmax(t, 0))
    }
  }

  function(t) {
    discount_fun(t)
  }
}

#' Price a Swap Using a Calibrated Term Structure
#'
#' Evaluate the present value and DV01 of a fixed-for-floating interest rate
#' swap using discount factors implied by a `term_structure_fit`.
#'
#' @inheritParams price_swap
#' @param fit A `term_structure_fit` resulting from [fit()].
#' @param discount_curve Optional identifier selecting the discount curve used
#'   for cash-flow present value calculations. Required when multiple curves
#'   are present in the calibration.
#' @param projection_curve Optional identifier selecting the projection curve
#'   used to derive forward rates. Defaults to `discount_curve` when omitted.
#'
#' @return Tibble with metrics `pv` and `dv01`.
#'
#' @examples
#' spec <- term_structure_spec()
#' quotes <- tibble::tibble(tenor = c(0.5, 1, 2), quote = c(0.02, 0.021, 0.024))
#' fit_obj <- fit(spec, quotes)
#' schedule <- swap_cashflow_schedule(0, 2, frequency = 2, notional = 1e6)
#' price_swap_term_structure(fit_obj, schedule, fixed_rate = 0.025)
#'
#' @export
price_swap_term_structure <- function(fit,
                                      schedule,
                                      fixed_rate,
                                      type = c("payer", "receiver"),
                                      discount_curve = NULL,
                                      projection_curve = NULL) {
  type <- rlang::arg_match(type)
  checkmate::assert_class(fit, "term_structure_fit")
  schedule <- validate_swap_schedule(schedule)
  checkmate::assert_number(fixed_rate, finite = TRUE)

  discount_fun <- term_structure_discount_function(fit, curve = discount_curve)

  if (is.null(projection_curve)) {
    projection_curve <- discount_curve
  }

  projection_fun <- term_structure_discount_function(fit, curve = projection_curve)

  df_pay <- unname(discount_fun(schedule$pay_time))
  accrual <- schedule$accrual_fraction

  proj_start <- unname(projection_fun(schedule$start))
  proj_end <- unname(projection_fun(schedule$end))

  forward_rates <- numeric(length(accrual))
  active_idx <- abs(accrual) > 1e-12
  if (any(active_idx)) {
    ratio <- proj_start[active_idx] / proj_end[active_idx]
    forward_rates[active_idx] <- (ratio - 1) / accrual[active_idx]
  }

  float_leg <- sum(schedule$notional * accrual * forward_rates * df_pay)
  fixed_leg <- sum(schedule$notional * fixed_rate * accrual * df_pay)

  pv <- if (type == "payer") {
    float_leg - fixed_leg
  } else {
    fixed_leg - float_leg
  }

  dv01 <- sum(schedule$notional * accrual * df_pay) * 1e-4

  tibble::tibble(
    metric = c("pv", "dv01"),
    value = c(pv, dv01)
  )
}

#' Quote-Level Sensitivities via Bump-and-Revalue
#'
#' Compute finite-difference sensitivities of a valuation functional with
#' respect to the market quotes used in term-structure calibration. The
#' procedure mirrors the Python `YieldCurveBuildGreeks.py` workflow by
#' recalibrating the curve after bumping each quote individually.
#'
#' @param spec A `term_structure_spec` configured with the desired calibration
#'   engine.
#' @param market_data Tidy data frame of market quotes supplied to [fit()].
#' @param valuation_fun Function accepting a `term_structure_fit` and returning
#'   a numeric scalar valuation (for example swap PV).
#' @param bump_size Numeric bump applied to each quote. Defaults to `1e-4`.
#' @param bump_type Either "additive" (default) or "relative". Relative bumps
#'   scale `bump_size` by the absolute quote magnitude, falling back to
#'   `bump_size` when the quote is zero.
#' @param bump_column Column within `market_data` that is perturbed. Defaults to
#'   "quote".
#'
#' @return Tibble containing the original quote metadata alongside the bumped
#'   valuation and sensitivity.
#'
#' @export
term_structure_quote_sensitivities <- function(spec,
                                               market_data,
                                               valuation_fun,
                                               bump_size = 1e-4,
                                               bump_type = c("additive", "relative"),
                                               bump_column = "quote") {
  bump_type <- rlang::arg_match(bump_type)
  checkmate::assert_class(spec, "term_structure_spec")
  checkmate::assert_function(valuation_fun)
  checkmate::assert_number(bump_size, lower = 0, finite = TRUE)

  data_tbl <- tibble::as_tibble(market_data)
  if (!bump_column %in% names(data_tbl)) {
    rlang::abort(
      message = paste0("Column '", bump_column, "' not found in market data."),
      class = "term_structure_missing_bump_column"
    )
  }

  quotes_vec <- as.numeric(data_tbl[[bump_column]])
  checkmate::assert_numeric(quotes_vec, any.missing = FALSE, finite = TRUE)
  data_tbl[[bump_column]] <- quotes_vec

  base_fit <- fit(spec, data_tbl)
  base_value <- valuation_fun(base_fit)
  checkmate::assert_number(base_value, finite = TRUE)

  idx_seq <- seq_len(nrow(data_tbl))
  sensitivity_tbl <- purrr::map_dfr(idx_seq, \(idx) {
    bumped_tbl <- data_tbl
    quote_val <- bumped_tbl[[bump_column]][idx]
    bump_amount <- if (identical(bump_type, "relative")) {
      adjusted <- abs(quote_val) * bump_size
      if (adjusted == 0) bump_size else adjusted
    } else {
      bump_size
    }

    bumped_tbl[[bump_column]][idx] <- quote_val + bump_amount

    bumped_fit <- fit(spec, bumped_tbl)
    bumped_value <- valuation_fun(bumped_fit)

    tibble::tibble(
      data_index = idx,
      curve = if ("curve" %in% names(data_tbl)) as.character(data_tbl$curve[idx]) else NA_character_,
      tenor = data_tbl$tenor[idx],
      instrument = if ("instrument" %in% names(data_tbl)) as.character(data_tbl$instrument[idx]) else NA_character_,
      quote = quote_val,
      bump = bump_amount,
      bumped_value = bumped_value,
      sensitivity = (bumped_value - base_value) / bump_amount
    )
  })

  sensitivity_tbl |>
    dplyr::mutate(
      base_value = base_value,
      bump_type = bump_type
    )
}
