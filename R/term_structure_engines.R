#' Term Structure Engine Implementations
#'
#' Internal calibration engines supporting `term_structure_spec()`. Each engine
#' accepts predictor data prepared by hardhat blueprints and returns a tidy
#' tibble containing discount, zero, and forward rates.
#'
#' @keywords internal
#' @noRd

tenor_key <- function(x) {
  sprintf("%.8f", x)
}

resolve_engine_option <- function(engine_args, name, default) {
  value <- engine_args[[name]]
  if (is.null(value)) default else value
}

normalize_instrument_type <- function(values, default_type) {
  values <- as.character(values)
  empty_idx <- is.na(values) | !nzchar(trimws(values))
  values[empty_idx] <- default_type
  canonical_input <- tolower(trimws(values))

  mapping <- c(
    discount_factor = "discount_factor",
    discount = "discount_factor",
    df = "discount_factor",
    zero_rate = "zero_rate",
    zero = "zero_rate",
    spot = "zero_rate",
    par_swap = "par_swap",
    swap = "par_swap",
    swap_rate = "par_swap",
    par = "par_swap",
    treasury = "treasury_price",
    treasury_price = "treasury_price",
    treasury_bond = "treasury_price",
    bond = "treasury_price"
  )

  canonical <- mapping[canonical_input]

  if (any(is.na(canonical))) {
    bad <- unique(values[is.na(canonical)])
    rlang::abort(
      message = "Unsupported instrument types supplied to term structure calibration.",
      class = "term_structure_unknown_instrument",
      instrument = bad
    )
  }

  canonical
}

prepare_term_structure_inputs <- function(predictors,
                                          quote_type,
                                          engine_args,
                                          require_unique = FALSE,
                                          curve_type = NULL,
                                          allow_multiple_curves = FALSE) {
  tibble_data <- tibble::as_tibble(predictors)

  if (!all(c("tenor", "quote") %in% names(tibble_data))) {
    rlang::abort("Predictors must contain `tenor` and `quote` columns")
  }

  tibble_data <- tibble_data |>
    dplyr::mutate(
      tenor = as.numeric(.data$tenor),
      quote = as.numeric(.data$quote)
    ) |>
    dplyr::arrange(.data$tenor)

  if (any(!is.finite(tibble_data$tenor)) || any(tibble_data$tenor < 0)) {
    rlang::abort("Tenors must be finite and non-negative")
  }
  if (any(!is.finite(tibble_data$quote))) {
    rlang::abort("Quotes must be finite numeric values")
  }

  instrument_col <- resolve_engine_option(engine_args, "instrument_col", "instrument")
  if (instrument_col %in% names(tibble_data)) {
    instrument_raw <- tibble_data[[instrument_col]]
  } else {
    instrument_raw <- rep(NA_character_, nrow(tibble_data))
  }

  instrument <- normalize_instrument_type(instrument_raw, default_type = quote_type)
  tibble_data$instrument <- instrument

  curve_col <- resolve_engine_option(engine_args, "curve_col", NULL)

  if (allow_multiple_curves) {
    if (is.null(curve_col) || !curve_col %in% names(tibble_data)) {
      rlang::abort(
        message = "Multi-curve calibration requires a curve identifier column.",
        class = "term_structure_missing_curve_column"
      )
    }
    tibble_data$curve <- as.character(tibble_data[[curve_col]])
    if (any(!nzchar(trimws(tibble_data$curve)))) {
      rlang::abort("Curve identifiers must be non-empty for multi-curve calibration.")
    }
  } else {
    if (!is.null(curve_col) && curve_col %in% names(tibble_data)) {
      curve_values <- unique(as.character(tibble_data[[curve_col]]))
      if (length(curve_values) > 1) {
        rlang::abort(
          message = "Single-curve engines received multiple curve identifiers.",
          class = "term_structure_unexpected_curves",
          curves = curve_values
        )
      }
      curve_type <- curve_values
    }
    if (is.null(curve_type)) {
      curve_type <- quote_type
    }
    tibble_data$curve <- rep_len(curve_type, nrow(tibble_data))
  }

  if (require_unique) {
    if (allow_multiple_curves) {
      dup_mask <- duplicated(tibble_data[, c("curve", "tenor")])
    } else {
      dup_mask <- duplicated(tibble_data$tenor)
    }
    if (any(dup_mask)) {
      rlang::abort("Term structure calibration currently requires unique tenor pillars.")
    }
  }

  tibble_data <- tibble_data |>
    dplyr::arrange(.data$curve, .data$tenor)

  list(data = tibble_data, instrument_col = instrument_col, curve_col = curve_col)
}

build_term_structure_output <- function(curve, tenor, discount, quote_reference) {
  stopifnot(length(curve) == length(tenor), length(tenor) == length(discount))

  order_idx <- order(curve, tenor)
  curve <- curve[order_idx]
  tenor <- tenor[order_idx]
  discount <- discount[order_idx]

  zero_rate <- ifelse(tenor > 0, -log(discount) / tenor, 0)

  if (length(discount) > 1) {
    forward_raw <- -diff(log(discount)) / diff(tenor)
    forward_rate <- c(forward_raw, utils::tail(forward_raw, 1))
  } else {
    forward_rate <- zero_rate
  }

  quote_lookup <- quote_reference |>
    dplyr::arrange(.data$curve, .data$tenor) |>
    dplyr::distinct(.data$curve, .data$tenor, .keep_all = TRUE) |>
    dplyr::select(.data$curve, .data$tenor, input_quote = .data$quote)

  key_lookup <- paste(quote_lookup$curve, tenor_key(quote_lookup$tenor))
  keys <- paste(curve, tenor_key(tenor))
  input_quote <- quote_lookup$input_quote[match(keys, key_lookup)]

  tibble::tibble(
    curve = curve,
    tenor = tenor,
    input_quote = input_quote,
    discount_factor = discount,
    zero_rate = zero_rate,
    forward_rate = forward_rate
  )
}

calibrate_term_structure_direct <- function(predictors,
                                            quote_type,
                                            engine_args = list(),
                                            curve_type = NULL) {
  prepared <- prepare_term_structure_inputs(
    predictors = predictors,
    quote_type = quote_type,
    engine_args = engine_args,
    require_unique = TRUE,
    curve_type = curve_type,
    allow_multiple_curves = FALSE
  )

  data <- prepared$data
  if (!all(data$instrument %in% c("zero_rate", "discount_factor"))) {
    rlang::abort("Direct engine supports only zero rate or discount factor quotes.")
  }

  discount <- ifelse(
    data$instrument == "discount_factor",
    data$quote,
    exp(-data$quote * data$tenor)
  )

  if (any(!is.finite(discount)) || any(discount <= 0)) {
    rlang::abort("Direct engine produced invalid discount factors.")
  }

  build_term_structure_output(
    curve = data$curve,
    tenor = data$tenor,
    discount = discount,
    quote_reference = data
  )
}

calibrate_term_structure_bootstrap <- function(predictors,
                                               quote_type,
                                               engine_args = list(),
                                               curve_type = NULL) {
  prepared <- prepare_term_structure_inputs(
    predictors = predictors,
    quote_type = quote_type,
    engine_args = engine_args,
    require_unique = TRUE,
    curve_type = curve_type,
    allow_multiple_curves = FALSE
  )

  data <- prepared$data
  frequency <- resolve_engine_option(engine_args, "fixed_leg_frequency", 1)
  checkmate::assert_number(frequency, lower = 1)
  step <- 1 / frequency

  bootstrap_rows <- split(data, seq_len(nrow(data)))

  bootstrap_step <- function(state, row_df) {
    map <- state$map
    tenor_val <- row_df$tenor[[1]]
    quote_val <- row_df$quote[[1]]
    instrument <- row_df$instrument[[1]]

    discount <- switch(
      instrument,
      discount_factor = quote_val,
      zero_rate = exp(-quote_val * tenor_val),
      {
        n_payments <- as.integer(round(tenor_val / step))
        if (abs(n_payments * step - tenor_val) > 1e-8) {
          rlang::abort(
            message = paste0("Tenor ", tenor_val, " is not an integer multiple of the coupon step for bootstrap engine."),
            class = "term_structure_invalid_coupon_grid",
            tenor = tenor_val
          )
        }

        payment_times <- seq_len(n_payments) * step
        accrual <- diff(c(0, payment_times))

        prior_count <- length(payment_times) - 1
        sum_prior <- if (prior_count > 0) {
          prior_times <- payment_times[seq_len(prior_count)]
          prior_keys <- tenor_key(prior_times)
          existing_keys <- names(map)
          if (is.null(existing_keys)) {
            existing_keys <- character(0)
          }
          missing <- prior_times[!prior_keys %in% existing_keys]
          if (length(missing) > 0) {
            rlang::abort(
              message = "Bootstrap engine requires discount factors for all earlier coupon dates.",
              class = "term_structure_missing_pillars",
              missing = missing
            )
          }
          prior_discounts <- purrr::map_dbl(prior_keys, \(key) map[[key]])
          sum(accrual[seq_len(prior_count)] * prior_discounts)
        } else {
          0
        }

        alpha_n <- accrual[length(accrual)]
        numerator <- 1 - quote_val * sum_prior
        denominator <- 1 + quote_val * alpha_n
        if (denominator <= 0) {
          rlang::abort("Par swap equation produced non-positive denominator during bootstrap calibration.")
        }

        numerator / denominator
      }
    )

    if (!is.finite(discount) || discount <= 0) {
      rlang::abort("Bootstrap engine produced invalid discount factor.")
    }

    map[[tenor_key(tenor_val)]] <- discount

    list(
      map = map,
      last_discount = discount
    )
  }

  states <- purrr::accumulate(
    .x = bootstrap_rows,
    .init = list(map = list(), last_discount = NA_real_),
    .f = bootstrap_step
  )[-1]

  discount_values <- purrr::map_dbl(states, "last_discount") |> unname()

  build_term_structure_output(
    curve = data$curve,
    tenor = data$tenor,
    discount = discount_values,
    quote_reference = data
  )
}

calibrate_term_structure_treasury <- function(predictors,
                                              quote_type,
                                              engine_args = list(),
                                              curve_type = NULL) {
  prepared <- prepare_term_structure_inputs(
    predictors = predictors,
    quote_type = quote_type,
    engine_args = engine_args,
    require_unique = TRUE,
    curve_type = curve_type,
    allow_multiple_curves = FALSE
  )

  data <- prepared$data

  if (!"coupon_rate" %in% names(data)) {
    rlang::abort(
      message = "Treasury calibration requires a `coupon_rate` column.",
      class = "term_structure_missing_coupon_rate"
    )
  }

  face_value <- resolve_engine_option(engine_args, "face_value", 100)
  checkmate::assert_number(face_value, lower = .Machine$double.eps)

  default_frequency <- resolve_engine_option(engine_args, "coupon_frequency", 2)
  checkmate::assert_number(default_frequency, lower = 1)

  has_coupon_frequency <- "coupon_frequency" %in% names(data)

  treasury_rows <- split(data, seq_len(nrow(data)))

  treasury_step <- function(state, row_df) {
    map <- state$map
    tenor_val <- row_df$tenor[[1]]
    quote_val <- row_df$quote[[1]]
    instrument <- row_df$instrument[[1]]

    discount <- switch(
      instrument,
      discount_factor = quote_val,
      zero_rate = exp(-quote_val * tenor_val),
      treasury_price = {
        coupon_rate <- as.numeric(row_df$coupon_rate[[1]])
        if (!is.finite(coupon_rate) || coupon_rate < 0) {
          rlang::abort(
            message = "Treasury instruments require non-negative coupon rates.",
            class = "term_structure_invalid_coupon_rate",
            coupon_rate = coupon_rate
          )
        }

        frequency_row <- default_frequency
        if (has_coupon_frequency) {
          freq_candidate <- as.numeric(row_df$coupon_frequency[[1]])
          if (is.finite(freq_candidate) && freq_candidate > 0) {
            frequency_row <- freq_candidate
          }
        }
        checkmate::assert_number(frequency_row, lower = 1)

        step <- 1 / frequency_row
        n_payments <- as.integer(round(tenor_val / step))
        if (abs(n_payments * step - tenor_val) > 1e-8) {
          rlang::abort(
            message = paste0("Tenor ", tenor_val, " is not aligned with the coupon frequency."),
            class = "term_structure_invalid_coupon_grid",
            tenor = tenor_val
          )
        }

        payment_times <- seq_len(n_payments) * step
        coupon_per_period <- coupon_rate / frequency_row

        price_per_unit <- quote_val / face_value

        prior_contribution <- if (n_payments > 1) {
          prior_times <- payment_times[-length(payment_times)]
          prior_keys <- tenor_key(prior_times)
          existing_keys <- names(map)
          if (is.null(existing_keys)) {
            existing_keys <- character(0)
          }
          missing <- prior_times[!prior_keys %in% existing_keys]
          if (length(missing) > 0) {
            rlang::abort(
              message = "Treasury bootstrap requires discount factors for earlier coupon dates.",
              class = "term_structure_missing_pillars",
              missing = missing
            )
          }
          prior_discounts <- purrr::map_dbl(prior_keys, \(key) map[[key]])
          sum(coupon_per_period * prior_discounts)
        } else {
          0
        }

        last_coupon <- coupon_per_period + 1
        (price_per_unit - prior_contribution) / last_coupon
      },
      {
        rlang::abort(
          message = paste0("Instrument type '", instrument, "' is not supported by the treasury engine."),
          class = "term_structure_unknown_instrument",
          instrument = instrument
        )
      }
    )

    if (!is.finite(discount) || discount <= 0) {
      rlang::abort("Treasury calibration produced invalid discount factors.")
    }

    map[[tenor_key(tenor_val)]] <- discount

    list(
      map = map,
      last_discount = discount
    )
  }

  treasury_states <- purrr::accumulate(
    .x = treasury_rows,
    .init = list(map = list(), last_discount = NA_real_),
    .f = treasury_step
  )[-1]

  discount_values <- purrr::map_dbl(treasury_states, "last_discount") |> unname()

  build_term_structure_output(
    curve = data$curve,
    tenor = data$tenor,
    discount = discount_values,
    quote_reference = data
  )
}

calibrate_term_structure_newton <- function(predictors,
                                            quote_type,
                                            engine_args = list(),
                                            curve_type = NULL) {
  prepared <- prepare_term_structure_inputs(
    predictors = predictors,
    quote_type = quote_type,
    engine_args = engine_args,
    require_unique = TRUE,
    curve_type = curve_type,
    allow_multiple_curves = FALSE
  )

  data <- prepared$data
  frequency <- resolve_engine_option(engine_args, "fixed_leg_frequency", 1)
  checkmate::assert_number(frequency, lower = 1)
  step <- 1 / frequency

  max_iter <- resolve_engine_option(engine_args, "max_iter", 500)
  tolerance <- resolve_engine_option(engine_args, "tolerance", 1e-8)
  checkmate::assert_integerish(max_iter, lower = 1, len = 1)
  checkmate::assert_number(tolerance, lower = 0)

  unique_tenors <- data$tenor
  tenor_keys <- tenor_key(unique_tenors)
  tenor_index <- stats::setNames(seq_along(unique_tenors), tenor_keys)

  par_rows <- which(data$instrument == "par_swap")
  if (length(par_rows) > 0) {
    purrr::walk(
      par_rows,
      \(row) {
        maturity <- data$tenor[row]
        n_payments <- as.integer(round(maturity / step))
        if (abs(n_payments * step - maturity) > 1e-8) {
          rlang::abort(
            message = paste0("Tenor ", maturity, " is not an integer multiple of the coupon step for Newton engine."),
            class = "term_structure_invalid_coupon_grid",
            tenor = maturity
          )
        }
        payment_times <- seq_len(n_payments) * step
        missing_idx <- !(tenor_key(payment_times) %in% names(tenor_index))
        if (any(missing_idx)) {
          rlang::abort(
            message = "Newton engine requires discount pillars for all coupon dates.",
            class = "term_structure_missing_pillars",
            missing = payment_times[missing_idx]
          )
        }
      }
    )
  }

  initial_discount <- purrr::map_dbl(
    unique_tenors,
    \(tenor_val) {
      row_idx <- which(data$tenor == tenor_val)[1]
      quote_val <- data$quote[row_idx]
      instrument <- data$instrument[row_idx]
      disc <- if (instrument == "discount_factor") {
        quote_val
      } else {
        exp(-quote_val * tenor_val)
      }
      if (!is.finite(disc) || disc <= 0) {
        exp(-0.02 * tenor_val)
      } else {
        disc
      }
    }
  )

  residual_function <- function(log_discount) {
    discounts <- exp(log_discount)

    purrr::map_dbl(
      seq_len(nrow(data)),
      \(row_idx) {
        tenor_val <- data$tenor[row_idx]
        quote_val <- data$quote[row_idx]
        instrument <- data$instrument[row_idx]
        discount_idx <- tenor_index[[tenor_key(tenor_val)]]
        discount <- discounts[discount_idx]

        if (instrument == "discount_factor") {
          return(discount - quote_val)
        }

        if (instrument == "zero_rate") {
          zero_val <- if (tenor_val > 0) -log(discount) / tenor_val else 0
          return(zero_val - quote_val)
        }

        n_payments <- as.integer(round(tenor_val / step))
        payment_times <- seq_len(n_payments) * step
        payment_keys <- tenor_key(payment_times)
        payment_indices <- unname(tenor_index[payment_keys])
        discount_payments <- discounts[payment_indices]
        accrual <- diff(c(0, payment_times))
        fixed_leg <- sum(accrual * discount_payments)
        float_leg <- 1 - discount_payments[length(discount_payments)]
        model_rate <- float_leg / fixed_leg
        model_rate - quote_val
      }
    )
  }

  objective <- function(log_discount) {
    res <- residual_function(log_discount)
    sum(res * res)
  }

  optimisation <- stats::optim(
    par = log(initial_discount),
    fn = objective,
    method = "BFGS",
    control = list(maxit = max_iter, reltol = tolerance)
  )

  if (optimisation$convergence != 0) {
    rlang::abort("Newton term structure calibration did not converge.")
  }

  calibrated_discount <- exp(optimisation$par)

  build_term_structure_output(
    curve = rep_len(data$curve[1], length(unique_tenors)),
    tenor = unique_tenors,
    discount = calibrated_discount,
    quote_reference = data
  )
}

calibrate_term_structure_multi_curve_newton <- function(predictors,
                                                        quote_type,
                                                        engine_args = list(),
                                                        default_curve_type = NULL) {
  prepared <- prepare_term_structure_inputs(
    predictors = predictors,
    quote_type = quote_type,
    engine_args = engine_args,
    require_unique = TRUE,
    curve_type = default_curve_type,
    allow_multiple_curves = TRUE
  )

  data <- prepared$data

  default_discount_curve <- default_curve_type
  if (is.null(default_discount_curve)) {
    default_discount_curve <- unique(data$curve)[1]
  }

  discount_curve_id <- resolve_engine_option(
    engine_args,
    "discount_curve_id",
    default_discount_curve
  )

  if (!discount_curve_id %in% data$curve) {
    rlang::abort(
      message = paste0("Discount curve '", discount_curve_id, "' not present in calibration data."),
      class = "term_structure_missing_discount_curve",
      curve = discount_curve_id
    )
  }

  frequency <- resolve_engine_option(engine_args, "fixed_leg_frequency", 1)
  checkmate::assert_number(frequency, lower = 1)
  accrual_step <- 1 / frequency

  max_iter <- resolve_engine_option(engine_args, "max_iter", 600)
  tolerance <- resolve_engine_option(engine_args, "tolerance", 1e-8)
  checkmate::assert_integerish(max_iter, lower = 1, len = 1)
  checkmate::assert_number(tolerance, lower = 0)

  unique_points <- data |>
    dplyr::distinct(.data$curve, .data$tenor) |>
    dplyr::arrange(.data$curve, .data$tenor)

  key_vector <- paste(unique_points$curve, tenor_key(unique_points$tenor))
  tenor_index <- stats::setNames(seq_len(nrow(unique_points)), key_vector)

  initial_discount <- purrr::map2_dbl(
    unique_points$curve,
    unique_points$tenor,
    \(curve_id, tenor_val) {
      rows <- which(data$curve == curve_id & abs(data$tenor - tenor_val) < 1e-12)
      disc <- NA_real_
      if (length(rows) > 0) {
        row_df <- data[rows, , drop = FALSE]
        discount_idx <- which(row_df$instrument == "discount_factor")
        if (length(discount_idx) > 0) {
          disc <- row_df$quote[discount_idx[[1]]]
        } else {
          zero_idx <- which(row_df$instrument == "zero_rate")
          if (length(zero_idx) > 0) {
            zero_quote <- row_df$quote[zero_idx[[1]]]
            disc <- exp(-zero_quote * tenor_val)
          }
        }
      }
      if (!is.finite(disc) || is.na(disc) || disc <= 0) {
        exp(-0.02 * tenor_val)
      } else {
        disc
      }
    }
  )

  get_discount <- function(log_discount, curve_id, tenor_val) {
    if (abs(tenor_val) < 1e-12) {
      return(1)
    }
    key <- paste(curve_id, tenor_key(tenor_val))
    idx <- tenor_index[[key]]
    if (is.null(idx)) {
      rlang::abort(
        message = paste0("Missing pillar for curve '", curve_id, "' at tenor ", tenor_val),
        class = "term_structure_missing_pillars",
        curve = curve_id,
        tenor = tenor_val
      )
    }
    exp(log_discount[idx])
  }

  residual_function <- function(log_discount) {
    purrr::map_dbl(
      seq_len(nrow(data)),
      \(row_idx) {
        curve_id <- data$curve[row_idx]
        tenor_val <- data$tenor[row_idx]
        quote_val <- data$quote[row_idx]
        instrument <- data$instrument[row_idx]

        key <- paste(curve_id, tenor_key(tenor_val))
        idx <- tenor_index[[key]]
        discount_val <- exp(log_discount[idx])

        if (instrument == "discount_factor") {
          return(discount_val - quote_val)
        }

        if (instrument == "zero_rate") {
          zero_val <- if (tenor_val > 0) -log(discount_val) / tenor_val else 0
          return(zero_val - quote_val)
        }

        if (instrument != "par_swap") {
          rlang::abort(
            message = paste0("Unsupported instrument '", instrument, "' for multi-curve calibration."),
            class = "term_structure_unknown_instrument",
            instrument = instrument
          )
        }

        n_payments <- as.integer(round(tenor_val / accrual_step))
        if (abs(n_payments * accrual_step - tenor_val) > 1e-8) {
          rlang::abort(
            message = paste0("Tenor ", tenor_val, " is not aligned with the coupon frequency."),
            class = "term_structure_invalid_coupon_grid",
            tenor = tenor_val
          )
        }

        payment_times <- seq_len(n_payments) * accrual_step
        accrual <- diff(c(0, payment_times))

        discount_df <- purrr::map_dbl(
          payment_times,
          \(t_j) get_discount(log_discount, discount_curve_id, t_j)
        )
        prev_times <- c(0, utils::head(payment_times, -1))
        prev_proj <- purrr::map_dbl(
          prev_times,
          \(t_prev) get_discount(log_discount, curve_id, t_prev)
        )
        curr_proj <- purrr::map_dbl(
          payment_times,
          \(t_j) get_discount(log_discount, curve_id, t_j)
        )
        forward_rate <- (prev_proj / curr_proj - 1) / accrual

        discount_leg <- sum(accrual * discount_df)
        floating_leg <- sum(accrual * discount_df * forward_rate)
        floating_leg - quote_val * discount_leg
      }
    )
  }

  objective <- function(log_discount) {
    res <- residual_function(log_discount)
    sum(res * res)
  }

  optimisation <- stats::optim(
    par = log(initial_discount),
    fn = objective,
    method = "BFGS",
    control = list(maxit = max_iter, reltol = tolerance)
  )

  if (optimisation$convergence != 0) {
    rlang::abort("Multi-curve Newton calibration did not converge.")
  }

  calibrated_discount <- exp(optimisation$par)

  build_term_structure_output(
    curve = unique_points$curve,
    tenor = unique_points$tenor,
    discount = calibrated_discount,
    quote_reference = data
  )
}
