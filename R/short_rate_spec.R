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
short_rate_spec <- function(model = c("ho_lee", "hull_white", "g2pp"),
                            volatility,
                            mean_reversion = NULL,
                            second_volatility = NULL,
                            second_mean_reversion = NULL,
                            correlation = 0,
                            curve = NULL,
                            initial_rate = NULL,
                            engine = "monte_carlo",
                            engine_options = list()) {
  model <- rlang::arg_match(model)
  checkmate::assert_number(volatility, lower = 0, finite = TRUE)

  if (identical(model, "ho_lee")) {
    mean_reversion <- 0
    second_volatility <- 0
    second_mean_reversion <- 0
    correlation <- 0
  } else if (identical(model, "hull_white")) {
    if (is.null(mean_reversion)) {
      rlang::abort("`mean_reversion` must be supplied for the Hull-White model")
    }
    checkmate::assert_number(mean_reversion, lower = 0, finite = TRUE)
    if (mean_reversion <= 0) {
      rlang::abort("`mean_reversion` must be strictly positive for Hull-White")
    }
    second_volatility <- 0
    second_mean_reversion <- 0
    correlation <- 0
  } else {
    if (is.null(mean_reversion) || is.null(second_mean_reversion)) {
      rlang::abort("`mean_reversion` and `second_mean_reversion` must be supplied for the G2++ model")
    }
    if (is.null(second_volatility)) {
      rlang::abort("`second_volatility` must be supplied for the G2++ model")
    }
    checkmate::assert_number(mean_reversion, lower = .Machine$double.eps, finite = TRUE)
    checkmate::assert_number(second_mean_reversion, lower = .Machine$double.eps, finite = TRUE)
    checkmate::assert_number(second_volatility, lower = 0, finite = TRUE)
    checkmate::assert_number(correlation, lower = -1, upper = 1, finite = TRUE)
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

  discount_fun_vec <- Vectorize(discount_fun)
  forward_fun_vec <- Vectorize(forward_fun)
  forward_derivative_fun_vec <- Vectorize(forward_derivative_fun)

  theta_raw <- NULL
  theta_fun_vec <- NULL
  phi_fun <- NULL
  phi_fun_vec <- NULL

  if (identical(model, "g2pp")) {
    phi_raw <- make_phi_function_g2pp(
      forward_fun = forward_fun,
      volatility1 = volatility,
      volatility2 = second_volatility,
      mean_reversion1 = mean_reversion,
      mean_reversion2 = second_mean_reversion,
      correlation = correlation
    )
    phi_fun <- phi_raw
    phi_fun_vec <- Vectorize(phi_raw)
  } else {
    theta_raw <- make_theta_function(
      forward_fun = forward_fun,
      forward_derivative_fun = forward_derivative_fun,
      volatility = volatility,
      mean_reversion = mean_reversion,
      model = model
    )
    theta_fun_vec <- Vectorize(theta_raw)
  }

  spec <- new_process_spec(
    class = "short_rate_spec",
    args = list(
      model = model,
      volatility = volatility,
      mean_reversion = mean_reversion,
      second_volatility = second_volatility,
      second_mean_reversion = second_mean_reversion,
      correlation = correlation,
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
    theta_fun_vec = theta_fun_vec,
    phi_fun = phi_fun,
    phi_fun_vec = phi_fun_vec
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

  args <- process_spec$args
  engine <- process_spec$method$engine
  if (is.null(engine)) {
    engine <- "monte_carlo"
  }

  engine_args <- process_spec$eng_args
  if (is.null(engine_args)) {
    engine_args <- list()
  }

  simulation <- switch(engine,
    monte_carlo = if (identical(args$model, "g2pp")) {
      rlang::exec(
        simulate_g2pp_paths,
        process_spec = process_spec,
        n_paths = n_paths,
        n_steps = n_steps,
        maturity = maturity,
        seed = seed,
        !!!engine_args
      )
    } else {
      rlang::exec(
        simulate_gaussian_short_rate_paths,
        process_spec = process_spec,
        n_paths = n_paths,
        n_steps = n_steps,
        maturity = maturity,
        seed = seed,
        !!!engine_args
      )
    },
    rlang::abort(
      message = paste0("Engine '", engine, "' is not implemented for simulate_paths.short_rate_spec"),
      class = "short_rate_unsupported_engine",
      engine = engine
    )
  )

  result_tbl <- tibble::tibble(
    path_id = rep(seq_len(n_paths), each = length(simulation$time)),
    time = rep(simulation$time, times = n_paths),
    short_rate = as.vector(t(simulation$rates)),
    discount_factor = as.vector(t(simulation$discounts))
  )

  if (!is.null(simulation$factor1)) {
    result_tbl <- result_tbl |>
      dplyr::mutate(
        factor_1 = as.vector(t(simulation$factor1)),
        factor_2 = as.vector(t(simulation$factor2))
      )
  }

  result_tbl
}


#' @export
price_zcb.short_rate_spec <- function(object,
                                      maturities,
                                      valuation_time = 0,
                                      short_rate = NULL,
                                      factor_state = NULL,
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

  if (identical(args$model, "g2pp")) {
    parsed_factors <- normalize_g2pp_factor_state(factor_state, length(maturities))

    purrr::map2_dbl(
      seq_along(maturities),
      maturities,
      function(idx, maturity) {
        if (abs(maturity - valuation_time) < 1e-12) {
          return(1)
        }

        g2pp_discount_point(
          discount_fun = state$discount_fun_vec,
          valuation_time = valuation_time,
          maturity = maturity,
          x_state = parsed_factors[idx, 1],
          y_state = parsed_factors[idx, 2],
          lambda1 = args$mean_reversion,
          lambda2 = args$second_mean_reversion,
          eta1 = args$volatility,
          eta2 = args$second_volatility,
          rho = args$correlation
        )
      }
    )
  } else {
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


make_phi_function_g2pp <- function(forward_fun,
                                   volatility1,
                                   volatility2,
                                   mean_reversion1,
                                   mean_reversion2,
                                   correlation) {
  function(t) {
    kappa1 <- g2pp_kappa(mean_reversion1, t)
    kappa2 <- g2pp_kappa(mean_reversion2, t)

    forward_fun(t) +
      0.5 * volatility1^2 * kappa1^2 +
      0.5 * volatility2^2 * kappa2^2 +
      correlation * volatility1 * volatility2 * kappa1 * kappa2
  }
}


g2pp_kappa <- function(lambda, delta) {
  if (lambda <= 1e-8) {
    delta
  } else {
    -expm1(-lambda * delta) / lambda
  }
}


g2pp_V_term <- function(lambda1,
                         lambda2,
                         eta1,
                         eta2,
                         rho,
                         lower,
                         upper) {
  delta <- upper - lower

  v1 <- if (lambda1 <= 1e-8) {
    eta1^2 * delta^3 / 3
  } else {
    exp1 <- exp(-lambda1 * delta)
    exp2 <- exp(-2 * lambda1 * delta)
    eta1^2 / (lambda1^2) * (
      delta +
        (2 / lambda1) * exp1 -
        (exp2 / (2 * lambda1)) -
        (3 / (2 * lambda1))
    )
  }

  v2 <- if (lambda2 <= 1e-8) {
    eta2^2 * delta^3 / 3
  } else {
    exp1 <- exp(-lambda2 * delta)
    exp2 <- exp(-2 * lambda2 * delta)
    eta2^2 / (lambda2^2) * (
      delta +
        (2 / lambda2) * exp1 -
        (exp2 / (2 * lambda2)) -
        (3 / (2 * lambda2))
    )
  }

  cross <- function(lambda_a, lambda_b) {
    if (abs(lambda_a + lambda_b) <= 1e-8) {
      delta^2 / 2
    } else {
      delta +
        expm1(-lambda_a * delta) / lambda_a +
        expm1(-lambda_b * delta) / lambda_b -
        expm1(-(lambda_a + lambda_b) * delta) / (lambda_a + lambda_b)
    }
  }

  cross_term <- 2 * rho * eta1 * eta2 / (lambda1 * lambda2) * cross(lambda1, lambda2)

  v1 + v2 + cross_term
}


g2pp_discount_point <- function(discount_fun,
                                 valuation_time,
                                 maturity,
                                 x_state,
                                 y_state,
                                 lambda1,
                                 lambda2,
                                 eta1,
                                 eta2,
                                 rho) {
  delta <- maturity - valuation_time
  b1 <- g2pp_kappa(lambda1, delta)
  b2 <- g2pp_kappa(lambda2, delta)

  v_tT <- g2pp_V_term(lambda1, lambda2, eta1, eta2, rho, valuation_time, maturity)
  v_0T <- g2pp_V_term(lambda1, lambda2, eta1, eta2, rho, 0, maturity)
  v_0t <- g2pp_V_term(lambda1, lambda2, eta1, eta2, rho, 0, valuation_time)

  int_phi <- -log(discount_fun(maturity) / discount_fun(valuation_time) * exp(-0.5 * (v_0T - v_0t)))

  exp(-int_phi - b1 * x_state - b2 * y_state + 0.5 * v_tT)
}


normalize_g2pp_factor_state <- function(factor_state, n) {
  if (is.null(factor_state)) {
    return(matrix(0, nrow = n, ncol = 2))
  }

  if (is.data.frame(factor_state)) {
    factor_state <- as.matrix(factor_state)
  }

  if (is.list(factor_state) && !is.matrix(factor_state)) {
    if (all(c("factor1", "factor2") %in% names(factor_state))) {
      factor_state <- cbind(factor_state$factor1, factor_state$factor2)
    } else if (all(c("x", "y") %in% names(factor_state))) {
      factor_state <- cbind(factor_state$x, factor_state$y)
    } else {
      rlang::abort("`factor_state` list must contain `factor1`/`factor2` or `x`/`y` components")
    }
  }

  if (is.numeric(factor_state) && !is.matrix(factor_state)) {
    if (length(factor_state) == 2) {
      factor_state <- matrix(rep(factor_state, times = n), nrow = n, byrow = TRUE)
    } else if (length(factor_state) == n * 2) {
      factor_state <- matrix(factor_state, ncol = 2, byrow = TRUE)
    } else {
      rlang::abort("Numeric `factor_state` must have length 2 or 2 * length(maturities)")
    }
  }

  factor_matrix <- as.matrix(factor_state)

  if (nrow(factor_matrix) == 1L && n > 1L) {
    factor_matrix <- matrix(rep(factor_matrix, each = n), nrow = n, byrow = TRUE)
  }

  if (nrow(factor_matrix) != n || ncol(factor_matrix) < 2) {
    rlang::abort("`factor_state` must provide two factors for each maturity")
  }

  factor_matrix[, 1:2, drop = FALSE]
}
