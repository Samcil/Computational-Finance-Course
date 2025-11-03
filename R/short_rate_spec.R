#' Short-Rate Model Specification
#'
#' Construct a specification for one-factor Gaussian short-rate models with a
#' focus on the Ho-Lee setup. The specification stores the volatility parameter

#' together with a discount curve that is used to infer the instantaneous
#' forward and drift adjustments required for arbitrage-free simulations.

#' @param model Character. One of "ho_lee" or "hull_white".
#' @param volatility Numeric. Instantaneous volatility of the short rate.
#' @param mean_reversion Numeric. Mean-reversion speed (required for Hull-White,
#'   ignored for Ho-Lee).
#' @param curve Either a `term_structure_fit` produced by [fit()], or a data
#'   frame/tibble with columns `tenor` (in years) and `discount_factor`.
#' @param initial_rate Optional numeric override for the initial short rate. If
#'   omitted, the rate is inferred from the supplied discount curve via the
#'   instantaneous forward rate at time zero.
#' @param engine Character identifier for the simulation engine. Defaults to
#'   "monte_carlo".
#' @param engine_options Named list of engine-specific options.
#'
#' @return An object of class `short_rate_spec` (inheriting from
#'   `model_spec`) that can be passed to [simulate_paths()].
#'
#' @examples
#' curve <- tibble::tibble(
#'   tenor = c(0, 1, 2, 3),
#'   discount_factor = exp(-0.02 * tenor)
#' )
#'
#' spec <- short_rate_spec(
#'   model = "ho_lee",
#'   volatility = 0.01,
#'   curve = curve
#' )
#'
#' @export
short_rate_spec <- function(model = c("ho_lee", "hull_white"),
                            volatility,
                            mean_reversion = NULL,
                            curve = NULL,
                            initial_rate = NULL,
                            engine = "monte_carlo",
                            engine_options = list()) {
  model <- rlang::arg_match(model)
  checkmate::assert_number(volatility, lower = 0, finite = TRUE)

  if (identical(model, "hull_white")) {
    if (is.null(mean_reversion)) {
      rlang::abort("`mean_reversion` must be supplied for the Hull-White model")
    }
    checkmate::assert_number(mean_reversion, lower = 0, finite = TRUE)
    if (mean_reversion <= 0) {
      rlang::abort("`mean_reversion` must be strictly positive for Hull-White")
    }
  } else {
    mean_reversion <- 0
  }

  checkmate::assert_list(engine_options, names = "unique", null.ok = FALSE)

  curve_data <- prepare_discount_curve(curve, initial_rate)
  log_discount_fun <- make_log_discount_spline(curve_data)
  discount_fun <- make_discount_function(log_discount_fun)
  forward_fun <- make_forward_function(log_discount_fun)
  forward_derivative_fun <- make_forward_derivative_function(log_discount_fun)

  if (is.null(initial_rate)) {
    initial_rate <- forward_fun(1e-4)
  } else {
    checkmate::assert_number(initial_rate, finite = TRUE)
  }

  theta_raw <- make_theta_function(
    forward_fun = forward_fun,
    forward_derivative_fun = forward_derivative_fun,
    volatility = volatility,
    mean_reversion = mean_reversion,
    model = model
  )

  discount_fun_vec <- Vectorize(discount_fun)
  forward_fun_vec <- Vectorize(forward_fun)
  forward_derivative_fun_vec <- Vectorize(forward_derivative_fun)
  theta_fun_vec <- Vectorize(theta_raw)

  spec <- new_process_spec(
    class = "short_rate_spec",
    args = list(
      model = model,
      volatility = volatility,
      mean_reversion = mean_reversion,
      curve = curve_data,
      initial_rate = initial_rate
    ),
    process_type = "short_rate",
    inheritance = "ou_family_spec"
  )

  spec <- set_process_metadata(
    spec,
    model_variant = model,
    curve_inputs = c("tenor", "discount_factor"),
    state_variables = c("short_rate", "discount_factor"),
    engines = list(
      simulate = "monte_carlo",
      price = "analytic_zcb"
    )
  )

  engine_state <- list(
    discount_fun = discount_fun,
    discount_fun_vec = discount_fun_vec,
    forward_fun = forward_fun,
    forward_fun_vec = forward_fun_vec,
    forward_derivative_fun = forward_derivative_fun,
    forward_derivative_fun_vec = forward_derivative_fun_vec,
    theta_fun = theta_raw,
    theta_fun_vec = theta_fun_vec
  )

  engine_options <- rlang::list2(!!!engine_options)
  spec <- set_engine(spec, engine = engine, !!!engine_options)
  spec$method$engine_state <- engine_state
  spec
}
#' @export
print.short_rate_spec <- function(x, ...) {
  cli::cli_h2("Short-Rate Specification")
  engine_label <- if (is.null(x$method$engine)) "<unset>" else x$method$engine
  args <- x$args
  cli::cli_dl(c(
    "Model" = cli::col_cyan(args$model),
    "Volatility" = cli::col_cyan(format(args$volatility, digits = 4)),
    "Mean reversion" = cli::col_cyan(format(args$mean_reversion, digits = 4)),
    "Engine" = cli::col_cyan(engine_label)
  ))
  cli::cli_text("{.strong Lineage:} {format_spec_lineage(x)}")
  cli::cli_alert_info("Call simulate_paths() to generate Monte Carlo scenarios")
  invisible(x)
}


#' @export
simulate_paths.short_rate_spec <- function(process_spec,
                                           n_paths,
                                           n_steps,
                                           maturity,
                                           seed = 123,
                                           ...) {
  rlang::check_dots_empty()

  engine <- process_spec$method$engine
  if (is.null(engine)) {
    engine <- "monte_carlo"
  }

  engine_args <- process_spec$eng_args
  if (is.null(engine_args)) {
    engine_args <- list()
  }

  simulation <- switch(engine,
    monte_carlo = rlang::exec(
      simulate_gaussian_short_rate_paths,
      process_spec = process_spec,
      n_paths = n_paths,
      n_steps = n_steps,
      maturity = maturity,
      seed = seed,
      !!!engine_args
    ),
    rlang::abort(
      message = paste0("Engine '", engine, "' is not implemented for simulate_paths.short_rate_spec"),
      class = "short_rate_unsupported_engine",
      engine = engine
    )
  )

  tibble::tibble(
    path_id = rep(seq_len(n_paths), each = length(simulation$time)),
    time = rep(simulation$time, times = n_paths),
    short_rate = as.vector(t(simulation$rates)),
    discount_factor = as.vector(t(simulation$discounts))
  )
}


#' @export
price_zcb.short_rate_spec <- function(object,
                                      maturities,
                                      valuation_time = 0,
                                      short_rate = NULL,
                                      ...) {
  rlang::check_dots_empty()

  checkmate::assert_numeric(maturities, lower = 0, finite = TRUE, any.missing = FALSE)
  checkmate::assert_number(valuation_time, lower = 0, finite = TRUE)
  maturities <- as.numeric(maturities)

  args <- object$args
  state <- short_rate_state(object)

  if (any(maturities < valuation_time - 1e-12)) {
    rlang::abort("`maturities` must be greater than or equal to `valuation_time`")
  }

  if (is.null(short_rate)) {
    if (valuation_time > 0) {
      rlang::abort("`short_rate` must be supplied when valuation_time > 0")
    }
    short_rate <- rep(args$initial_rate, length(maturities))
  }

  checkmate::assert_numeric(short_rate, finite = TRUE, any.missing = FALSE)

  if (length(maturities) == 1L && length(short_rate) > 1L) {
    maturities <- rep(maturities, length(short_rate))
  } else if (length(short_rate) == 1L && length(maturities) > 1L) {
    short_rate <- rep(short_rate, length(maturities))
  }

  if (length(short_rate) != length(maturities)) {
    rlang::abort("`short_rate` must have length 1 or match `maturities`")
  }

  a <- args$mean_reversion
  sigma <- args$volatility
  theta_fun <- state$theta_fun

  purrr::map2_dbl(
    maturities,
    short_rate,
    function(maturity, rate_t) {
      if (abs(maturity - valuation_time) < 1e-12) {
        return(1)
      }

      delta <- maturity - valuation_time
      b_tT <- b_factor(a, delta)
      theta_int <- theta_integral(theta_fun, a, valuation_time, maturity)
      variance_term <- sigma^2 * b_squared_integral(a, valuation_time, maturity)

      exp(-rate_t * b_tT - theta_int + 0.5 * variance_term)
    }
  )
}


prepare_discount_curve <- function(curve, initial_rate) {
  if (inherits(curve, "term_structure_fit")) {
    curve_tbl <- augment(curve) |>
      dplyr::select(tenor, discount_factor)
  } else if (is.null(curve)) {
    if (is.null(initial_rate)) {
      rlang::abort("Provide either `curve` data or `initial_rate`")
    }
    curve_tbl <- tibble::tibble(
      tenor = c(0, 30),
      discount_factor = exp(-initial_rate * c(0, 30))
    )
  } else {
    checkmate::assert_data_frame(curve, any.missing = FALSE, min.rows = 1)
    if (!all(c("tenor", "discount_factor") %in% names(curve))) {
      rlang::abort("`curve` must contain `tenor` and `discount_factor`")
    }
    curve_tbl <- tibble::as_tibble(curve[, c("tenor", "discount_factor")])
  }

  curve_tbl <- curve_tbl |>
    dplyr::mutate(
      tenor = as.numeric(tenor),
      discount_factor = as.numeric(discount_factor)
    ) |>
    dplyr::arrange(tenor) |>
    dplyr::distinct(tenor, .keep_all = TRUE)

  if (any(!is.finite(curve_tbl$tenor)) || any(curve_tbl$tenor < 0)) {
    rlang::abort("All tenors must be finite and non-negative")
  }

  if (any(!is.finite(curve_tbl$discount_factor)) || any(curve_tbl$discount_factor <= 0)) {
    rlang::abort("Discount factors must be positive and finite")
  }

  if (!any(abs(curve_tbl$tenor) < 1e-12)) {
    curve_tbl <- dplyr::bind_rows(
      tibble::tibble(tenor = 0, discount_factor = 1),
      curve_tbl
    ) |>
      dplyr::arrange(tenor)
  } else {
    zero_index <- which.min(abs(curve_tbl$tenor))
    curve_tbl$discount_factor[zero_index] <- 1
  }

  curve_tbl
}


make_theta_function <- function(forward_fun,
                                forward_derivative_fun,
                                volatility,
                                mean_reversion,
                                model) {
  sigma_sq <- volatility^2
  if (identical(model, "ho_lee")) {
    function(t) {
      forward_derivative_fun(t) + sigma_sq * t
    }
  } else {
    function(t) {
      forward_fun(t) + forward_derivative_fun(t) / mean_reversion +
        (sigma_sq / (2 * mean_reversion^2)) * (1 - exp(-2 * mean_reversion * t))
    }
  }
}


b_factor <- function(mean_reversion, delta) {
  if (mean_reversion > 0) {
    (1 - exp(-mean_reversion * delta)) / mean_reversion
  } else {
    delta
  }
}


theta_integral <- function(theta_fun, mean_reversion, lower, upper) {
  if (abs(upper - lower) < 1e-12) {
    return(0)
  }

  integrand <- function(u) {
    delta <- upper - u
    b_val <- if (mean_reversion > 0) {
      (1 - exp(-mean_reversion * delta)) / mean_reversion
    } else {
      delta
    }
    theta_fun(u) * b_val
  }

  stats::integrate(
    integrand,
    lower = lower,
    upper = upper,
    rel.tol = 1e-6,
    subdivisions = 200
  )$value
}


b_squared_integral <- function(mean_reversion, lower, upper) {
  delta <- upper - lower
  if (mean_reversion > 0) {
    term1 <- delta / (mean_reversion^2)
    term2 <- 2 * (1 - exp(-mean_reversion * delta)) / (mean_reversion^3)
    term3 <- (1 - exp(-2 * mean_reversion * delta)) / (2 * mean_reversion^3)
    term1 - term2 + term3
  } else {
    (delta^3) / 3
  }
}


#' @export
set_engine.short_rate_spec <- function(object, engine, ...) {
  checkmate::assert_choice(engine, choices = c("monte_carlo"))
  eng_args <- rlang::list2(...)
  set_engine_base(object, engine = engine, eng_args = eng_args)
}


#' Retrieve Engine State for Short-Rate Specifications
#'
#' Expose the engine state produced by [short_rate_spec()] for downstream
#' utilities that need direct access to discount functions or model metadata.
#'
#' @param spec A short-rate specification created by [short_rate_spec()].
#'
#' @return List containing engine-specific state, including discount
#'   functions.
#' @export
short_rate_state <- function(spec) {
  state <- spec$method$engine_state
  if (is.null(state)) {
    rlang::abort("Short-rate specification has no engine state; construct via short_rate_spec() or set_engine().")
  }
  state
}
