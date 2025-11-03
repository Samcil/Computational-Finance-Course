#' Exposure Specification for Interest-Rate Portfolios
#'
#' Assemble a portfolio of interest-rate swaps and capture netting set
#' assignments together with reporting preferences for exposure analytics. The
#' resulting specification plugs into [simulate_exposure()] which reuses the
#' short-rate Monte Carlo engines to produce swap-level, netting-set, and
#' portfolio summaries in line with the lecture scripts
#' `Exposures_HW_Netting.py` and `ConvexityCorrection.py`.
#'
#' @param trades List of trades. Each trade must be a list containing
#'   `schedule` (output of [swap_cashflow_schedule()]), `fixed_rate`, `type`
#'   (`"payer"` or `"receiver"`), optional `trade_id`, and optional
#'   `netting_set` (defaults to the trade identifier).
#' @param quantiles Numeric vector of tail probabilities in `(0, 1)` at which
#'   potential future exposure (PFE) is reported.
#' @param discount Logical. When `TRUE` (default) the summary reports discounted
#'   positive exposure in addition to undiscounted metrics.
#'
#' @return An object of class `exposure_spec` that can be passed to
#'   [simulate_exposure()].
#'
#' @export
exposure_spec <- function(trades,
                          quantiles = 0.95,
                          discount = TRUE) {
  trade_list <- validate_exposure_trades(trades)
  checkmate::assert_numeric(quantiles, lower = 0, upper = 1, any.missing = FALSE)
  quantiles <- unique(sort(quantiles))
  checkmate::assert_flag(discount)

  spec <- new_model_spec(
    class = "exposure_spec",
    args = list(
      trades = trade_list,
      quantiles = quantiles,
      discount = discount
    ),
    mode = "exposure_analysis"
  )

  netting_sets <- sort(unique(purrr::map_chr(trade_list, "netting_set")))
  spec <- set_spec_metadata(
    spec,
    trades = purrr::map_chr(trade_list, "trade_id"),
    netting_sets = netting_sets
  )
  spec
}

#' Simulate Portfolio Exposure
#'
#' Run Monte Carlo exposure analytics for a portfolio created via
#' [exposure_spec()]. The engine simulates short-rate paths once and revalues all
#' trades per path, producing tidy summaries for trades, netting sets, and the
#' aggregate portfolio. Netting aggregation follows the lecture workflow in
#' `Exposures_HW_Netting.py` by applying positive-part logic *after* summing
#' trades within each set.
#'
#' @param spec Exposure specification created by [exposure_spec()].
#' @param process_spec Short-rate specification produced by [short_rate_spec()].
#' @param n_paths Integer. Number of Monte Carlo paths.
#' @param n_steps Integer. Number of time steps.
#' @param maturity Numeric. Horizon in years. Defaults to the maximum payment
#'   date across all trades.
#' @param seed Integer seed passed to [simulate_paths()].
#' @param keep_paths Logical. When `TRUE`, path-level exposures are included in
#'   the result for downstream VaR calculations.
#'
#' @return An object of class `portfolio_exposure_result` containing portfolio,
#'   netting-set, and trade summaries plus optional path-level detail.
#'
#' @export
simulate_exposure <- function(spec,
                              process_spec,
                              ...) {
  UseMethod("simulate_exposure")
}

#' @export
simulate_exposure.exposure_spec <- function(spec,
                                            process_spec,
                                            n_paths,
                                            n_steps,
                                            maturity = NULL,
                                            seed = 123,
                                            keep_paths = FALSE,
                                            ...) {
  rlang::check_dots_empty(...)
  checkmate::assert_class(process_spec, "short_rate_spec")
  checkmate::assert_int(n_paths, lower = 1)
  checkmate::assert_int(n_steps, lower = 1)
  checkmate::assert_number(seed, finite = TRUE)
  checkmate::assert_flag(keep_paths)

  trades <- spec$args$trades
  quantiles <- spec$args$quantiles

  final_payments <- purrr::map_dbl(trades, ~ max(.x$schedule$pay_time))
  if (is.null(maturity)) {
    maturity <- max(final_payments)
  }
  if (maturity < max(final_payments) - 1e-10) {
    rlang::abort("`maturity` must cover the latest payment date in `trades`")
  }

  path_data <- simulate_paths(
    process_spec = process_spec,
    n_paths = n_paths,
    n_steps = n_steps,
    maturity = maturity,
    seed = seed
  )

  trade_paths <- purrr::map_dfr(
    trades,
    compute_trade_paths,
    path_data = path_data,
    process_spec = process_spec
  )

  netting_paths <- aggregate_exposure_paths(trade_paths, c("netting_set"))
  portfolio_paths <- aggregate_exposure_paths(trade_paths, character())

  trade_summary <- summarise_exposure_paths(trade_paths, group_cols = c("trade_id"), quantiles = quantiles, discount = spec$args$discount)
  netting_summary <- summarise_exposure_paths(netting_paths, group_cols = c("netting_set"), quantiles = quantiles, discount = spec$args$discount)
  portfolio_summary <- summarise_exposure_paths(portfolio_paths, group_cols = character(), quantiles = quantiles, discount = spec$args$discount)

  result <- list(
    portfolio = portfolio_summary,
    netting = netting_summary,
    trades = trade_summary
  )

  if (keep_paths) {
    result$paths <- list(
      trades = trade_paths,
      netting = netting_paths,
      portfolio = portfolio_paths
    )
  }

  class(result) <- "portfolio_exposure_result"
  result
}

#' @export
print.portfolio_exposure_result <- function(x, ...) {
  cli::cli_h2("Portfolio Exposure Summary")
  cli::cli_text("Portfolio level metrics:")
  print(x$portfolio, ...)
  cli::cli_rule()
  cli::cli_text("Netting sets:")
  print(x$netting, ...)
  cli::cli_rule()
  cli::cli_text("Trades:")
  print(x$trades, ...)
  invisible(x)
}

#' Risk-Measure Specification for VaR and ES
#'
#' Capture portfolio composition and reporting thresholds for Value-at-Risk and
#' Expected Shortfall analytics. Monte Carlo VaR reuses
#' [simulate_exposure()] to obtain the distribution of future portfolio values
#' under a short-rate model, matching the `MonteCarloVaR.py` script. Historical
#' VaR consumes shocked discount curves aligned with `HistoricalVaR_Calculation.py`.
#'
#' @param trades List of trades as accepted by [exposure_spec()].
#' @param horizon Numeric scalar. Horizon in years used for VaR reporting.
#' @param alpha Numeric vector of tail probabilities (e.g., `0.05` for 95% VaR).
#' @param method Character. Either `"monte_carlo"` or `"historical"`.
#'
#' @return An object of class `risk_measure_spec` for use with
#'   [simulate_risk_measures()].
#'
#' @export
risk_measure_spec <- function(trades,
                              horizon,
                              alpha = 0.05,
                              method = c("monte_carlo", "historical")) {
  method <- rlang::arg_match(method)
  trade_list <- validate_exposure_trades(trades)
  checkmate::assert_number(horizon, lower = 0, finite = TRUE)
  checkmate::assert_numeric(alpha, lower = 0, upper = 1, any.missing = FALSE)
  alpha <- unique(sort(alpha))

  spec <- new_model_spec(
    class = "risk_measure_spec",
    args = list(
      trades = trade_list,
      horizon = horizon,
      alpha = alpha,
      method = method
    ),
    mode = "risk_measure"
  )

  spec <- set_spec_metadata(
    spec,
    trades = purrr::map_chr(trade_list, "trade_id"),
    method = method
  )
  spec
}

#' Simulate VaR and ES for a Portfolio
#'
#' @param spec Risk specification created by [risk_measure_spec()].
#' @param process_spec Short-rate specification. Required when
#'   `method = "monte_carlo"`.
#' @param ... Additional arguments passed to the underlying engine. For Monte
#'   Carlo this includes `n_paths`, `n_steps`, `maturity`, and `seed`. For
#'   historical VaR pass `curve_scenarios` (list of discount-factor tibbles) and
#'   an optional `base_curve`.
#'
#' @return Tibble containing VaR and ES estimates for each requested `alpha`.
#'
#' @export
simulate_risk_measures <- function(spec,
                                   process_spec = NULL,
                                   ...) {
  UseMethod("simulate_risk_measures")
}

#' @export
simulate_risk_measures.risk_measure_spec <- function(spec,
                                                     process_spec = NULL,
                                                     ...) {
  args <- rlang::list2(...)
  method <- spec$args$method

  if (identical(method, "monte_carlo")) {
    required <- c("n_paths", "n_steps")
    present <- names(args)
    missing <- setdiff(required, present)
    if (length(missing) > 0) {
      rlang::abort(paste("Missing arguments for Monte Carlo VaR:", paste(missing, collapse = ", ")))
    }
    res <- run_monte_carlo_var(
      spec = spec,
      process_spec = process_spec,
      n_paths = args$n_paths,
      n_steps = args$n_steps,
      maturity = args$maturity,
      seed = args$seed
    )
  } else {
    res <- run_historical_var(
      spec = spec,
      curve_scenarios = args$curve_scenarios,
      base_curve = args$base_curve
    )
  }
  res
}

#' @export
print.risk_measure_spec <- function(x, ...) {
  cli::cli_h2("Risk Measure Specification")
  cli::cli_dl(c(
    "Method" = x$args$method,
    "Horizon" = format(x$args$horizon),
    "Tail probs" = paste(format(x$args$alpha), collapse = ", ")
  ))
  cli::cli_text("Trades: {length(x$args$trades)}")
  invisible(x)
}

run_monte_carlo_var <- function(spec,
                                process_spec,
                                n_paths,
                                n_steps,
                                maturity = NULL,
                                seed = 123) {
  if (is.null(process_spec)) {
    rlang::abort("`process_spec` must be supplied for Monte Carlo VaR")
  }

  exposure_spec_obj <- exposure_spec(
    trades = spec$args$trades,
    quantiles = spec$args$alpha,
    discount = TRUE
  )

  exposure_res <- simulate_exposure(
    spec = exposure_spec_obj,
    process_spec = process_spec,
    n_paths = n_paths,
    n_steps = n_steps,
    maturity = maturity,
    seed = seed,
    keep_paths = TRUE
  )

  horizon_time <- select_horizon_time(exposure_res$paths$portfolio$time, spec$args$horizon)
  horizon_slice <- exposure_res$paths$portfolio |>
    dplyr::filter(abs(.data$time - horizon_time) < 1e-10)

  values <- horizon_slice$value
  if (length(values) == 0) {
    rlang::abort("Failed to extract horizon valuations for VaR computation")
  }

  compute_var_es(values, spec$args$alpha)
}

run_historical_var <- function(spec,
                               curve_scenarios,
                               base_curve = NULL) {
  if (is.null(curve_scenarios) || length(curve_scenarios) == 0) {
    rlang::abort("Provide `curve_scenarios` for historical VaR")
  }

  scenario_values <- purrr::map_dbl(
    curve_scenarios,
    ~ evaluate_portfolio_with_curve(spec$args$trades, .x)
  )

  baseline <- if (is.null(base_curve)) {
    evaluate_portfolio_with_curve(spec$args$trades, curve_scenarios[[1]])
  } else {
    evaluate_portfolio_with_curve(spec$args$trades, base_curve)
  }

  compute_var_es(scenario_values - baseline, spec$args$alpha)
}

compute_var_es <- function(values, alpha) {
  losses <- -values
  var_values <- unname(stats::quantile(losses, probs = alpha, names = FALSE))
  es_values <- purrr::map_dbl(alpha, function(prob) {
    var_cutoff <- stats::quantile(losses, probs = prob, names = FALSE)
    mean(losses[losses >= var_cutoff])
  })

  dplyr::bind_rows(
    tibble::tibble(metric = "VaR", alpha = alpha, value = var_values),
    tibble::tibble(metric = "ES", alpha = alpha, value = es_values)
  )
}

select_horizon_time <- function(time_grid, horizon) {
  time_grid[which.min(abs(time_grid - horizon))]
}

compute_trade_paths <- function(trade, path_data, process_spec) {
  schedule <- trade$schedule
  path_data |>
    dplyr::group_by(.data$path_id) |>
    dplyr::group_modify(~ compute_swap_exposure_path(
      .x,
      spec = process_spec,
      schedule = schedule,
      fixed_rate = trade$fixed_rate,
      type = trade$type
    )) |>
    dplyr::ungroup() |>
    dplyr::transmute(
      path_id = .data$path_id,
      time = .data$time,
      discount_factor = .data$discount_factor,
      value = .data$swap_value,
      positive_exposure = .data$positive_exposure,
      discounted_positive_exposure = .data$discounted_positive_exposure,
      trade_id = trade$trade_id,
      netting_set = trade$netting_set
    )
}

aggregate_exposure_paths <- function(paths, extra_group) {
  group_cols <- c("path_id", "time", extra_group)
  paths |>
    dplyr::group_by(dplyr::across(all_of(group_cols))) |>
    dplyr::summarise(
      discount_factor = dplyr::first(.data$discount_factor),
      value = sum(.data$value),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      positive_exposure = pmax(.data$value, 0),
      discounted_positive_exposure = .data$positive_exposure * .data$discount_factor
    )
}

summarise_exposure_paths <- function(paths, group_cols, quantiles, discount) {
  quantiles <- quantiles[quantiles > 0 & quantiles < 1]
  quantile_labels <- purrr::map_chr(quantiles, format_quantile_label)
  quantile_exprs <- build_quantile_expressions(quantiles, quantile_labels)

  summary_tbl <- paths |>
    dplyr::group_by(dplyr::across(all_of(c(group_cols, "time")))) |>
    dplyr::summarise(
      expected_exposure = mean(.data$positive_exposure),
      expected_discounted_exposure = mean(.data$discounted_positive_exposure),
      !!!quantile_exprs,
      .groups = "drop"
    ) |>
    dplyr::arrange(dplyr::across(all_of(group_cols)), .data$time)

  if (!discount) {
    summary_tbl <- dplyr::select(summary_tbl, -"expected_discounted_exposure")
  }

  summary_tbl
}

validate_exposure_trades <- function(trades) {
  checkmate::assert_list(trades, min.len = 1)
  purrr::imap(trades, normalise_trade)
}

normalise_trade <- function(trade, idx) {
  schedule <- trade$schedule
  if (is.null(schedule)) {
    rlang::abort("Each trade must supply a `schedule`")
  }
  schedule <- validate_swap_schedule(schedule)
  type <- trade$type
  type <- rlang::arg_match(type, c("payer", "receiver"))
  checkmate::assert_number(trade$fixed_rate, finite = TRUE)

  list(
    trade_id = trade$trade_id %||% paste0("trade_", idx),
    schedule = schedule,
    fixed_rate = trade$fixed_rate,
    type = type,
    netting_set = trade$netting_set %||% trade$trade_id %||% paste0("trade_", idx)
  )
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

evaluate_portfolio_with_curve <- function(trades, curve_tbl) {
  discount_fun <- build_discount_function_from_curve(curve_tbl)
  values <- purrr::map_dbl(
    trades,
    ~ price_swap_with_discount_fun(
      schedule = .x$schedule,
      fixed_rate = .x$fixed_rate,
      type = .x$type,
      discount_fun = discount_fun
    )
  )
  sum(values)
}

build_discount_function_from_curve <- function(curve_tbl) {
  checkmate::assert_data_frame(curve_tbl, any.missing = FALSE)
  required <- c("tenor", "discount_factor")
  missing <- setdiff(required, names(curve_tbl))
  if (length(missing) > 0) {
    rlang::abort(paste("Curve is missing columns:", paste(missing, collapse = ", ")))
  }
  curve_tbl <- tibble::as_tibble(curve_tbl) |>
    dplyr::arrange(.data$tenor) |>
    dplyr::distinct(.data$tenor, .keep_all = TRUE)
  if (!any(abs(curve_tbl$tenor) < 1e-12)) {
    curve_tbl <- dplyr::bind_rows(
      tibble::tibble(tenor = 0, discount_factor = 1),
      curve_tbl
    ) |>
      dplyr::arrange(.data$tenor)
  }
  log_fun <- make_log_discount_spline(curve_tbl)
  make_discount_function(log_fun)
}

price_swap_with_discount_fun <- function(schedule, fixed_rate, type, discount_fun) {
  df_pay <- discount_fun(schedule$pay_time)
  df_start <- discount_fun(schedule$start)
  df_end <- discount_fun(schedule$end)

  float_leg <- sum(schedule$notional * (df_start - df_end))
  fixed_leg <- sum(schedule$notional * fixed_rate * schedule$accrual_fraction * df_pay)

  if (type == "payer") {
    float_leg - fixed_leg
  } else {
    fixed_leg - float_leg
  }
}
