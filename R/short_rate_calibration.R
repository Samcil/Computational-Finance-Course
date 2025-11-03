#' Calibrate Short-Rate Volatility Using Caplet Quotes
#'
#' Estimate the volatility parameter of a Gaussian short-rate specification by
#' matching analytic caplet prices to supplied market quotes. The calibration is
#' performed via a one-dimensional search that minimises the weighted squared
#' error between model and market prices.
#'
#' @param spec A [`short_rate_spec`] object that provides the discount curve and
#'   mean-reversion input for analytic pricing.
#' @param caplet_quotes Data frame or tibble with at least the columns `reset`,
#'   `payment`, `strike`, and `quote`. Optional columns include `notional`,
#'   `accrual`, and `weight`.
#' @param lower Numeric scalar. Lower bound for the volatility search interval.
#' @param upper Numeric scalar. Upper bound for the volatility search interval.
#'   Defaults to a multiple of the current specification volatility when not
#'   supplied.
#' @param tol Numeric scalar. Optimisation tolerance passed to
#'   [stats::optimize()].
#'
#' @return An object of class `short_rate_calibration_result` containing the
#'   calibrated specification, parameter estimates, and fitted quote details.
#' @export
calibrate_short_rate_volatility <- function(spec,
                                            caplet_quotes,
                                            lower = 1e-6,
                                            upper = NULL,
                                            tol = 1e-8) {
  checkmate::assert_class(spec, "short_rate_spec")
  args <- spec$args
  engine <- spec$method$engine
  if (is.null(engine)) {
    engine <- "monte_carlo"
  }
  engine_options <- spec$eng_args
  if (is.null(engine_options)) {
    engine_options <- list()
  }
  args <- spec$args
  engine <- spec$method$engine
  if (is.null(engine)) {
    engine <- "monte_carlo"
  }
  engine_options <- spec$eng_args
  if (is.null(engine_options)) {
    engine_options <- list()
  }
  quotes <- tibble::as_tibble(caplet_quotes)
  required_cols <- c("reset", "payment", "strike", "quote")
  missing_cols <- setdiff(required_cols, names(quotes))
  if (length(missing_cols) > 0) {
    rlang::abort(
      message = "`caplet_quotes` is missing required columns",
      class = "caplet_quotes_missing_columns",
      columns = missing_cols
    )
  }

  if (nrow(quotes) == 0) {
    rlang::abort("`caplet_quotes` must contain at least one row")
  }

  quotes$notional <- if ("notional" %in% names(quotes)) quotes$notional else rep(1, nrow(quotes))
  quotes$accrual <- if ("accrual" %in% names(quotes)) quotes$accrual else quotes$payment - quotes$reset
  quotes$weight <- if ("weight" %in% names(quotes)) quotes$weight else rep(1, nrow(quotes))

  checkmate::assert_numeric(quotes$notional, finite = TRUE, any.missing = FALSE)
  checkmate::assert_numeric(quotes$accrual, lower = 0, finite = TRUE, any.missing = FALSE)
  checkmate::assert_numeric(quotes$weight, lower = 0, finite = TRUE, any.missing = FALSE)

  current_vol <- args$volatility
  if (is.null(upper)) {
    upper <- max(current_vol * 5, lower * 10, 1)
  }
  if (upper <= lower) {
    rlang::abort("`upper` must be greater than `lower`")
  }

  objective <- function(vol) {
    candidate <- short_rate_spec(
      model = args$model,
      volatility = vol,
      mean_reversion = args$mean_reversion,
      curve = args$curve,
      initial_rate = args$initial_rate,
      engine = engine,
      engine_options = engine_options
    )
    model_prices <- caplet_price_vector(candidate, quotes)
    errors <- model_prices - quotes$quote
    sum(quotes$weight * errors^2)
  }

  opt <- stats::optimize(objective, interval = c(lower, upper), tol = tol)
  calibrated_vol <- opt$minimum

  calibrated_spec <- short_rate_spec(
    model = args$model,
    volatility = calibrated_vol,
    mean_reversion = args$mean_reversion,
    curve = args$curve,
    initial_rate = args$initial_rate,
    engine = engine,
    engine_options = engine_options
  )

  fitted_prices <- caplet_price_vector(calibrated_spec, quotes)
  fitted_quotes <- quotes
  fitted_quotes$model_price <- fitted_prices
  fitted_quotes$error <- fitted_prices - quotes$quote

  estimates <- tibble::tibble(
    parameter = "volatility",
    estimate = calibrated_vol,
    objective = opt$objective,
    n_quotes = nrow(quotes)
  )

  result <- list(
    spec = calibrated_spec,
    estimates = estimates,
    fitted = fitted_quotes
  )
  class(result) <- "short_rate_calibration_result"
  result
}


caplet_price_vector <- function(spec, quotes) {
  purrr::pmap_dbl(
    list(
      quotes$reset,
      quotes$payment,
      quotes$strike,
      quotes$notional,
      quotes$accrual
    ),
    function(reset, payment, strike, notional, accrual) {
      price_caplet(
        spec = spec,
        reset = reset,
        payment = payment,
        strike = strike,
        notional = notional,
        accrual = accrual
      )$value
    }
  )
}


#' Calibrate Short-Rate Volatility Using Swaption Quotes
#'
#' Estimate the volatility parameter of a Gaussian short-rate specification by
#' matching analytic swaption prices (via Jamshidian's trick) to supplied market
#' quotes. The calibration minimises the weighted squared pricing errors.
#'
#' @param spec A [`short_rate_spec`] object.
#' @param swaption_quotes Tibble or data frame with columns `schedule` (a list
#'   of tibbles built by [swap_cashflow_schedule()]), `fixed_rate`, `type`
#'   ("payer" or "receiver"), and `quote`. Optional column `weight` supplies
#'   non-negative calibration weights.
#' @param lower Numeric scalar. Lower bound for the volatility search interval.
#' @param upper Numeric scalar. Upper bound for the search interval. Defaults to
#'   a multiple of the current specification volatility when not supplied.
#' @param tol Numeric scalar. Optimisation tolerance passed to
#'   [stats::optimize()].
#'
#' @return An object of class `short_rate_calibration_result` containing the
#'   calibrated specification, parameter estimates, and fitted quote details.
#' @export
calibrate_short_rate_swaption_volatility <- function(spec,
                                                     swaption_quotes,
                                                     lower = 1e-6,
                                                     upper = NULL,
                                                     tol = 1e-8) {
  checkmate::assert_class(spec, "short_rate_spec")
  args <- spec$args
  engine <- spec$method$engine
  if (is.null(engine)) {
    engine <- "monte_carlo"
  }
  engine_options <- spec$eng_args
  if (is.null(engine_options)) {
    engine_options <- list()
  }
  quotes <- tibble::as_tibble(swaption_quotes)
  required_cols <- c("schedule", "fixed_rate", "type", "quote")
  missing_cols <- setdiff(required_cols, names(quotes))
  if (length(missing_cols) > 0) {
    rlang::abort(
      message = "`swaption_quotes` is missing required columns",
      class = "swaption_quotes_missing_columns",
      columns = missing_cols
    )
  }

  if (nrow(quotes) == 0) {
    rlang::abort("`swaption_quotes` must contain at least one row")
  }

  if (!is.list(quotes$schedule)) {
    rlang::abort("`swaption_quotes$schedule` must be a list column of tibbles")
  }

  required_schedule_cols <- c("start", "end", "pay_time", "accrual_fraction", "notional")
  invalid_schedule <- purrr::map_lgl(quotes$schedule, function(x) {
    !checkmate::test_data_frame(x, any.missing = FALSE) ||
      length(setdiff(required_schedule_cols, names(x))) > 0
  })
  if (any(invalid_schedule)) {
    rlang::abort(
      message = "Each schedule must be a tibble from `swap_cashflow_schedule()`",
      class = "swaption_quotes_invalid_schedule",
      rows = which(invalid_schedule)
    )
  }

  quotes$type <- tolower(as.character(quotes$type))
  checkmate::assert_character(quotes$type, any.missing = FALSE)
  invalid_type <- !quotes$type %in% c("payer", "receiver")
  if (any(invalid_type)) {
    rlang::abort(
      message = "`swaption_quotes$type` must be 'payer' or 'receiver'",
      class = "swaption_quotes_invalid_type",
      rows = which(invalid_type)
    )
  }

  quotes$weight <- if ("weight" %in% names(quotes)) quotes$weight else rep(1, nrow(quotes))
  checkmate::assert_numeric(quotes$weight, lower = 0, finite = TRUE, any.missing = FALSE)
  checkmate::assert_numeric(quotes$fixed_rate, finite = TRUE, any.missing = FALSE)
  checkmate::assert_numeric(quotes$quote, lower = 0, finite = TRUE, any.missing = FALSE)

  current_vol <- args$volatility
  if (is.null(upper)) {
    upper <- max(current_vol * 5, lower * 10, 1)
  }
  if (upper <= lower) {
    rlang::abort("`upper` must be greater than `lower`")
  }

  objective <- function(vol) {
    candidate <- short_rate_spec(
      model = args$model,
      volatility = vol,
      mean_reversion = args$mean_reversion,
      curve = args$curve,
      initial_rate = args$initial_rate,
      engine = engine,
      engine_options = engine_options
    )
    model_prices <- swaption_price_vector(candidate, quotes)
    errors <- model_prices - quotes$quote
    sum(quotes$weight * errors^2)
  }

  opt <- stats::optimize(objective, interval = c(lower, upper), tol = tol)
  calibrated_vol <- opt$minimum

  calibrated_spec <- short_rate_spec(
    model = args$model,
    volatility = calibrated_vol,
    mean_reversion = args$mean_reversion,
    curve = args$curve,
    initial_rate = args$initial_rate,
    engine = engine,
    engine_options = engine_options
  )

  fitted_prices <- swaption_price_vector(calibrated_spec, quotes)
  fitted_quotes <- quotes
  fitted_quotes$model_price <- fitted_prices
  fitted_quotes$error <- fitted_prices - quotes$quote

  estimates <- tibble::tibble(
    parameter = "volatility",
    estimate = calibrated_vol,
    objective = opt$objective,
    n_quotes = nrow(quotes)
  )

  result <- list(
    spec = calibrated_spec,
    estimates = estimates,
    fitted = fitted_quotes
  )
  class(result) <- "short_rate_calibration_result"
  result
}


swaption_price_vector <- function(spec, quotes) {
  purrr::pmap_dbl(
    list(
      quotes$schedule,
      quotes$fixed_rate,
      quotes$type
    ),
    function(schedule, fixed_rate, type) {
      price_swaption(
        spec = spec,
        schedule = schedule,
        fixed_rate = fixed_rate,
        type = type
      )["value"][[1]]
    }
  )
}


#' Fit a Short-Rate Specification to Market Quotes
#'
#' Calibrate the volatility parameter of a Gaussian short-rate specification
#' using either caplet or swaption market quotes while following hardhat
#' blueprint conventions.
#'
#' @param object A [`short_rate_spec`] object created by [short_rate_spec()].
#' @param data Tibble or data frame containing the calibration instruments.
#'   Expected columns depend on `method`.
#' @param method Character. Calibration routine to apply. Supported values are
#'   `"caplet_volatility"` and `"swaption_volatility"`.
#' @param ... Additional arguments passed to the underlying calibration helper
#'   (e.g., `lower`, `upper`, `tol`).
#'
#' @return An object of class `short_rate_fit` containing the calibrated
#'   specification and fitted quote diagnostics.
#' @export
fit.short_rate_spec <- function(object,
                                data,
                                method = c("caplet_volatility", "swaption_volatility"),
                                ...) {
  method <- rlang::arg_match(method)
  dots <- rlang::list2(...)

  checkmate::assert_data_frame(data, any.missing = FALSE, min.rows = 1)
  quotes_tbl <- tibble::as_tibble(data)

  required_cols <- switch(method,
    caplet_volatility = c("reset", "payment", "strike", "quote"),
    swaption_volatility = c("schedule", "fixed_rate", "type", "quote")
  )
  missing_cols <- setdiff(required_cols, names(quotes_tbl))
  if (length(missing_cols) > 0) {
    rlang::abort(
      message = "`data` is missing required columns",
      class = "short_rate_calibration_missing_columns",
      columns = missing_cols
    )
  }

  blueprint <- hardhat::default_xy_blueprint(
    intercept = FALSE,
    composition = "tibble"
  )

  forged <- hardhat::mold(
    x = quotes_tbl,
    y = NULL,
    blueprint = blueprint
  )

  cleaned_quotes <- forged$predictors

  calibration_result <- switch(method,
    caplet_volatility = rlang::exec(
      calibrate_short_rate_volatility,
      spec = object,
      caplet_quotes = cleaned_quotes,
      !!!dots
    ),
    swaption_volatility = rlang::exec(
      calibrate_short_rate_swaption_volatility,
      spec = object,
      swaption_quotes = cleaned_quotes,
      !!!dots
    )
  )

  fit <- list(
    spec = calibration_result$spec,
    calibration = calibration_result,
    method = method,
    blueprint = forged$blueprint,
    original_spec = object
  )
  class(fit) <- "short_rate_fit"
  fit
}


#' @export
print.short_rate_fit <- function(x, ...) {
  cli::cli_h2("Short-Rate Fit")
  cli::cli_text("Calibration method: {cli::col_cyan(x$method)}")
  cli::cli_rule()
  print(x$spec)
  cli::cli_rule()
  cli::cli_dl(c(
    "Objective" = format(x$calibration$estimates$objective, digits = 6),
    "Quotes fitted" = x$calibration$estimates$n_quotes
  ))
  invisible(x)
}


#' @export
augment.short_rate_fit <- function(x, ...) {
  rlang::check_dots_empty()
  tibble::as_tibble(x$calibration$fitted)
}


#' @export
print.short_rate_calibration_result <- function(x, ...) {
  cli::cli_h2("Short-Rate Calibration")
  cli::cli_text(
    "Calibrated volatility: {format(x$estimates$estimate, digits = 6)}"
  )
  cli::cli_text(
    "Objective value: {format(x$estimates$objective, digits = 6)}"
  )
  cli::cli_text(
    "Caplet quotes fitted: {x$estimates$n_quotes}"
  )
  invisible(x)
}
