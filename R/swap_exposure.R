#' Simulate Interest-Rate Swap Exposure Profiles
#'
#' Generate Monte Carlo exposure profiles for a fixed-for-floating interest rate
#' swap under a Gaussian short-rate specification. The routine evaluates pathwise
#' swap values, positive exposures, and discounted exposures that are suitable
#' for CVA-style analytics while returning tidy summaries of expected exposure
#' (EE) and potential future exposure (PFE) at supplied confidence levels.
#'
#' @param spec A [`short_rate_spec`] constructed via [short_rate_spec()].
#' @param schedule Tibble produced by [swap_cashflow_schedule()] describing the
#'   swap cash-flow schedule.
#' @param fixed_rate Numeric scalar giving the fixed leg coupon in decimal form.
#' @param type Character string indicating the swap direction. Either
#'   `"payer"` (pay fixed, receive float) or `"receiver"`.
#' @param n_paths Integer. Number of Monte Carlo paths.
#' @param n_steps Integer. Number of simulation steps over the horizon.
#' @param maturity Numeric scalar. Simulation horizon in years. Defaults to the
#'   maximum payment time in `schedule`.
#' @param quantiles Numeric vector of probabilities in `(0, 1)` at which PFEs
#'   are reported. Defaults to `0.95`.
#' @param seed Integer seed for reproducibility.
#' @param keep_paths Logical. When `TRUE`, the returned list contains the full
#'   path-level exposures in addition to the aggregated summary.
#'
#' @return A list with elements:
#'   - `summary`: tibble containing expected exposure (undiscounted), expected
#'     discounted exposure, and PFE columns for each requested quantile.
#'   - `paths`: tibble of pathwise results (present when `keep_paths = TRUE`) with
#'     columns `path_id`, `time`, `short_rate`, `discount_factor`,
#'     `swap_value`, `positive_exposure`, and `discounted_positive_exposure`.
#'
#' @examples
#' curve <- tibble::tibble(tenor = 0:5, discount_factor = exp(-0.02 * tenor))
#' spec <- short_rate_spec(model = "ho_lee", volatility = 0.01, curve = curve)
#' schedule <- swap_cashflow_schedule(start = 0, end = 5, frequency = 2, notional = 1e6)
#' exposures <- simulate_swap_exposure(
#'   spec = spec,
#'   schedule = schedule,
#'   fixed_rate = 0.02,
#'   type = "payer",
#'   n_paths = 1000,
#'   n_steps = 200,
#'   quantiles = c(0.9, 0.95)
#' )
#'
#' @export
simulate_swap_exposure <- function(spec,
                                   schedule,
                                   fixed_rate,
                                   type = c("payer", "receiver"),
                                   n_paths,
                                   n_steps,
                                   maturity = max(schedule[["pay_time"]]),
                                   quantiles = 0.95,
                                   seed = 123,
                                   keep_paths = FALSE) {
  type <- rlang::arg_match(type)
  checkmate::assert_class(spec, "short_rate_spec")
  schedule <- validate_swap_schedule(schedule)
  checkmate::assert_number(fixed_rate, finite = TRUE)
  checkmate::assert_int(n_paths, lower = 1)
  checkmate::assert_int(n_steps, lower = 1)
  checkmate::assert_number(maturity, lower = 0)
  checkmate::assert_numeric(quantiles, lower = 0, upper = 1, any.missing = FALSE)
  checkmate::assert_number(seed, finite = TRUE)
  checkmate::assert_flag(keep_paths)

  if (maturity < max(schedule[["pay_time"]]) - 1e-10) {
    rlang::abort("`maturity` must cover the final payment date in `schedule`")
  }

  path_data <- simulate_paths(
    process_spec = spec,
    n_paths = n_paths,
    n_steps = n_steps,
    maturity = maturity,
    seed = seed
  )

  path_exposures <- path_data |>
    dplyr::group_by(.data$path_id) |>
    dplyr::group_modify(~ compute_swap_exposure_path(
      .x,
      spec = spec,
      schedule = schedule,
      fixed_rate = fixed_rate,
      type = type
    )) |>
    dplyr::ungroup()

  summary <- summarise_swap_exposure(path_exposures, quantiles = quantiles)

  result <- list(summary = summary)
  if (keep_paths) {
    result$paths <- path_exposures
  }

  class(result) <- "swap_exposure_result"
  result
}


compute_swap_exposure_path <- function(path_tbl,
                                       spec,
                                       schedule,
                                       fixed_rate,
                                       type) {
  times <- path_tbl$time
  rates <- path_tbl$short_rate
  discount <- path_tbl$discount_factor
  values <- compute_swap_values(
    spec = spec,
    schedule = schedule,
    fixed_rate = fixed_rate,
    type = type,
    times = times,
    rates = rates
  )

  positive <- pmax(values, 0)
  discounted_positive <- positive * discount

  path_tbl |>
    dplyr::mutate(
      swap_value = values,
      positive_exposure = positive,
      discounted_positive_exposure = discounted_positive
    )
}


compute_swap_values <- function(spec,
                                schedule,
                                fixed_rate,
                                type,
                                times,
                                rates,
                                tol = 1e-10) {
  pay_times <- schedule[["pay_time"]]
  start_times <- schedule[["start"]]
  end_times <- schedule[["end"]]
  accruals <- schedule[["accrual_fraction"]]
  notionals <- schedule[["notional"]]

  evaluate_swap_value <- function(t_val, r_val) {
    future_idx <- pay_times > t_val + tol
    if (!any(future_idx)) {
      return(0)
    }

    future_count <- sum(future_idx)
    rate_vec <- rep(r_val, future_count)

    end_df <- price_zcb(
      spec,
      maturities = end_times[future_idx],
      valuation_time = t_val,
      short_rate = rate_vec
    )

    start_df <- rep(1, future_count)
    future_start_idx <- start_times[future_idx] > t_val + tol
    if (any(future_start_idx)) {
      start_df[future_start_idx] <- price_zcb(
        spec,
        maturities = start_times[future_idx][future_start_idx],
        valuation_time = t_val,
        short_rate = rep(r_val, sum(future_start_idx))
      )
    }

    future_notionals <- notionals[future_idx]
    future_accruals <- accruals[future_idx]

    float_leg <- sum(future_notionals * (start_df - end_df))
    fixed_leg <- sum(future_notionals * fixed_rate * future_accruals * end_df)

    if (type == "payer") {
      float_leg - fixed_leg
    } else {
      fixed_leg - float_leg
    }
  }

  purrr::map2_dbl(times, rates, evaluate_swap_value)
}


summarise_swap_exposure <- function(path_exposures,
                                    quantiles) {
  quantiles <- quantiles[quantiles > 0 & quantiles < 1]
  quantile_labels <- purrr::map_chr(quantiles, format_quantile_label)

  quantile_exprs <- build_quantile_expressions(quantiles, quantile_labels)

  summary_tbl <- path_exposures |>
    dplyr::group_by(.data$time) |>
    dplyr::summarise(
      expected_exposure = mean(.data$positive_exposure),
      expected_discounted_exposure = mean(.data$discounted_positive_exposure),
      !!!quantile_exprs,
      .groups = "drop"
    )

  summary_tbl
}


build_quantile_expressions <- function(quantiles, labels) {
  if (length(quantiles) == 0) {
    return(list())
  }

  exprs <- purrr::imap(
    quantiles,
    \(prob, idx) rlang::expr(stats::quantile(.data$positive_exposure, probs = !!prob, names = FALSE))
  )
  stats::setNames(exprs, labels)
}


format_quantile_label <- function(q) {
  q_percent <- q * 100
  if (abs(q_percent - round(q_percent)) < 1e-8) {
    q_str <- format(round(q_percent))
  } else {
    q_str <- gsub("\\.", "p", format(q_percent, trim = TRUE))
  }
  paste0("pfe_", q_str)
}


#' @export
print.swap_exposure_result <- function(x, ...) {
  cli::cli_h2("Swap Exposure Summary")
  cli::cli_text("Metrics computed for {nrow(x$summary)} time points.")
  cli::cli_rule()
  print(x$summary, ...) # nocov
  invisible(x)
}
