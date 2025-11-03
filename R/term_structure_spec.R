#' Term Structure Specification
#'
#' Construct a specification for deterministic term-structure calibration that
#' follows tidymodels and hardhat design conventions. The specification stores
#' curve-level metadata, expected quote types, and engine settings used during
#' calibration.
#'
#' @param curve_type Character. The underlying curve being calibrated. Options
#'   are `"ois"` for overnight indexed swap discounting or `"libor"` for
#'   forward rate curves. Defaults to `"ois"`.
#' @param quote_type Character. Market quotes supplied in the calibration data.
#'   Supported values are `"zero_rate"` (continuously compounded) and
#'   `"discount_factor"`. Defaults to `"zero_rate"`.
#' @param engine Character. Engine identifier. Supported values are
#'   `"direct"` (closed-form conversion of zero/discount quotes),
#'   `"bootstrap"` (piecewise construction from a mix of zero/discount and par
#'   swap quotes), `"newton"` (nonlinear solver matching instrument
#'   residuals), `"multi_curve_newton"` (joint discount/projection curve
#'   calibration), and `"treasury_bootstrap"` (coupon bond bootstrapping for
#'   treasury securities).
#' @param engine_options Named list with engine-specific configuration. Common
#'   entries include `instrument_col` (defaults to "instrument"),
#'   `fixed_leg_frequency` (number of fixed coupons per year, defaults to 1),
#'   and `max_iter` (maximum iterations for the Newton engine).
#'
#' @return A `term_structure_spec` object that can be passed to
#'   [fit()] along with market quotes.
#'
#' @examples
#' spec <- term_structure_spec()
#' spec
#'
#' @export
term_structure_spec <- function(curve_type = c("ois", "libor"),
                                quote_type = c("zero_rate", "discount_factor"),
                                engine = "direct",
                                engine_options = list()) {
  curve_type <- rlang::arg_match(curve_type)
  quote_type <- rlang::arg_match(quote_type)

  checkmate::assert_list(engine_options, names = "unique", null.ok = FALSE)

  spec <- new_model_spec(
    class = "term_structure_spec",
    args = list(
      curve_type = curve_type,
      quote_type = quote_type
    ),
    mode = "calibration"
  )

  spec <- set_spec_metadata(
    spec,
    spec_type = "term_structure",
    engines = list(
      fit = c("direct", "bootstrap", "newton", "multi_curve_newton", "treasury_bootstrap"),
      predict = "linear_interp"
    ),
    notes = "Deterministic curve calibration with linear interpolation of zero rates"
  )

  engine_options <- rlang::list2(!!!engine_options)
  spec <- set_engine(spec, engine = engine, !!!engine_options)
  spec
}


#' @export
print.term_structure_spec <- function(x, ...) {
  cli::cli_h2("Term Structure Specification")
  engine_label <- if (is.null(x$method$engine)) "<unset>" else x$method$engine
  cli::cli_dl(c(
    "Curve type" = cli::col_cyan(x$args$curve_type),
    "Quote type" = cli::col_cyan(x$args$quote_type),
    "Engine" = cli::col_cyan(engine_label)
  ))
  cli::cli_alert_info("Use fit() with market quotes to calibrate the curve")
  invisible(x)
}


#' Fit Term Structure Specification
#'
#' Calibrate a deterministic discount curve from market quotes. The function
#' validates and tidies the incoming data using hardhat blueprints, then
#' delegates the calibration to the configured engine.
#'
#' @param object A `term_structure_spec` created by [term_structure_spec()].
#' @param data A data frame or tibble containing market quotes with at least the
#'   columns `tenor` (in years) and `quote` (value consistent with `quote_type`).
#' @param ... Reserved for future engine-specific arguments.
#'
#' @return An object of class `term_structure_fit`.
#'
#' @examples
#' spec <- term_structure_spec()
#' quotes <- tibble::tibble(
#'   tenor = c(0.5, 1, 2, 3),
#'   quote = c(0.03, 0.032, 0.035, 0.037)
#' )
#' fit(spec, quotes)
#'
#' @export
fit.term_structure_spec <- function(object, data, ...) {
  rlang::check_dots_empty()
  checkmate::assert_data_frame(data, any.missing = FALSE, min.rows = 1)

  required_cols <- c("tenor", "quote")
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    rlang::abort(
      message = "`data` is missing required columns",
      class = "term_structure_missing_columns",
      missing = missing_cols
    )
  }

  blueprint <- hardhat::default_xy_blueprint(
    intercept = FALSE
  )

  forged <- hardhat::mold(
    x = tibble::as_tibble(data),
    y = NULL,
    blueprint = blueprint
  )

  engine <- object$method$engine
  if (is.null(engine)) {
    engine <- "direct"
  }

  engine_args <- object$eng_args
  if (is.null(engine_args)) {
    engine_args <- list()
  }

  calibration <- switch(engine,
    direct = calibrate_term_structure_direct(
      predictors = forged$predictors,
      quote_type = object$args$quote_type,
      engine_args = engine_args,
      curve_type = object$args$curve_type
    ),
    bootstrap = calibrate_term_structure_bootstrap(
      predictors = forged$predictors,
      quote_type = object$args$quote_type,
      engine_args = engine_args,
      curve_type = object$args$curve_type
    ),
    newton = calibrate_term_structure_newton(
      predictors = forged$predictors,
      quote_type = object$args$quote_type,
      engine_args = engine_args,
      curve_type = object$args$curve_type
    ),
    multi_curve_newton = calibrate_term_structure_multi_curve_newton(
      predictors = forged$predictors,
      quote_type = object$args$quote_type,
      engine_args = engine_args,
      default_curve_type = object$args$curve_type
    ),
    treasury_bootstrap = calibrate_term_structure_treasury(
      predictors = forged$predictors,
      quote_type = object$args$quote_type,
      engine_args = engine_args,
      curve_type = object$args$curve_type
    ),
    rlang::abort(
      message = paste0("Engine '", engine, "' is not implemented for term_structure_spec"),
      class = "term_structure_unsupported_engine",
      engine = engine
    )
  )

  fit <- list(
    spec = object,
    blueprint = forged$blueprint,
    calibration = calibration
  )
  class(fit) <- "term_structure_fit"
  fit
}


#' @export
print.term_structure_fit <- function(x, ...) {
  cli::cli_h2("Term Structure Fit")
  print(x$spec)
  cli::cli_rule()
  cli::cli_text("Calibrated tenors: {.value {paste(format(x$calibration$tenor, digits = 4), collapse = ', ')}}")
  invisible(x)
}


#' Augment Term Structure Calibration
#'
#' Returns the calibrated discount factors, zero rates, and forward rates in a
#' tidy tibble.
#'
#' @param x A `term_structure_fit` produced by [fit()].
#' @param ... Unused.
#'
#' @return Tibble containing the calibrated term structure.
#'
#' @export
augment.term_structure_fit <- function(x, ...) {
  rlang::check_dots_empty()
  x$calibration
}


#' Predict from a Term Structure Fit
#'
#' @param object A `term_structure_fit`.
#' @param new_data Optional data frame with a `tenor` column specifying the
#'   maturities (in years) at which to evaluate the curve. If omitted, the
#'   calibration tenors are returned.
#' @param type Character. Prediction type: `"discount"`, `"zero"`, or
#'   `"forward"`.
#' @param ... Unused.
#'
#' @return Tibble with columns `tenor` and `.pred` containing the requested
#'   quantity.
#'
#' @export
predict.term_structure_fit <- function(object,
                                       new_data = NULL,
                                       type = c("discount", "zero", "forward"),
                                       ...) {
  rlang::check_dots_empty()
  type <- rlang::arg_match(type)
  calibration <- object$calibration

  if (!"curve" %in% names(calibration)) {
    curve_fallback <- object$spec$args$curve_type
    if (is.null(curve_fallback)) {
      curve_fallback <- "curve"
    }
    calibration$curve <- curve_fallback
  }

  curves <- unique(calibration$curve)
  multi_curve <- length(curves) > 1

  interpolate_curve <- function(curve_data, tenor_values) {
    if (length(tenor_values) == 0) {
      return(tibble::tibble(
        tenor = numeric(0),
        discount_factor = numeric(0),
        zero_rate = numeric(0),
        forward_rate = numeric(0)
      ))
    }

    zero_interp <- stats::approx(
      x = curve_data$tenor,
      y = curve_data$zero_rate,
      xout = tenor_values,
      rule = 2,
      ties = "ordered"
    )$y

    discount_interp <- exp(-zero_interp * tenor_values)

    forward_interp <- if (nrow(curve_data) > 1) {
      stats::approx(
        x = curve_data$tenor,
        y = curve_data$forward_rate,
        xout = tenor_values,
        rule = 2,
        ties = "ordered"
      )$y
    } else {
      rep(curve_data$forward_rate, length(tenor_values))
    }

    tibble::tibble(
      tenor = tenor_values,
      discount_factor = discount_interp,
      zero_rate = zero_interp,
      forward_rate = forward_interp
    )
  }

  evaluated <- if (is.null(new_data)) {
    calibration
  } else {
    checkmate::assert_data_frame(new_data, any.missing = FALSE, min.rows = 1)
    if (!"tenor" %in% names(new_data)) {
      rlang::abort(
        message = "`new_data` must contain a `tenor` column",
        class = "term_structure_missing_tenor"
      )
    }
    tenor <- as.numeric(new_data$tenor)
    if (any(!is.finite(tenor)) || any(tenor < 0)) {
      rlang::abort("`tenor` values must be non-negative and finite")
    }

    if (multi_curve) {
      if (!"curve" %in% names(new_data)) {
        rlang::abort(
          message = "`new_data` must include a `curve` column for multi-curve predictions",
          class = "term_structure_missing_curve"
        )
      }

      new_data <- tibble::as_tibble(new_data)
      new_data |>
        dplyr::group_split(.data$curve, .keep = TRUE) |>
        purrr::map_dfr(
          \(curve_df) {
            curve_name <- curve_df$curve[[1]]
            curve_data <- calibration[calibration$curve == curve_name, , drop = FALSE]
            if (nrow(curve_data) == 0) {
              rlang::abort(
                message = paste0("No calibrated curve data found for `", curve_name, "`."),
                class = "term_structure_unknown_curve",
                curve = curve_name
              )
            }
            interpolated <- interpolate_curve(curve_data, as.numeric(curve_df$tenor))
            interpolated$curve <- curve_name
            interpolated
          }
        )
    } else {
      interpolated <- interpolate_curve(calibration, tenor)
      interpolated$curve <- calibration$curve[1]
      interpolated
    }
  }

  result <- switch(type,
    discount = tibble::tibble(
      curve = evaluated$curve,
      tenor = evaluated$tenor,
      .pred = evaluated$discount_factor
    ),
    zero = tibble::tibble(
      curve = evaluated$curve,
      tenor = evaluated$tenor,
      .pred = evaluated$zero_rate
    ),
    forward = tibble::tibble(
      curve = evaluated$curve,
      tenor = evaluated$tenor,
      .pred = evaluated$forward_rate
    )
  )

  if (!multi_curve) {
    result$curve <- NULL
  }

  result
}


#' @export
set_engine.term_structure_spec <- function(object, engine, ...) {
  checkmate::assert_choice(engine, choices = c("direct", "bootstrap", "newton", "multi_curve_newton", "treasury_bootstrap"))
  eng_args <- rlang::list2(...)
  set_engine_base(object, engine = engine, eng_args = eng_args)
}
