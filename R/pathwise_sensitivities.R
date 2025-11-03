#' Pathwise Delta Estimator for Geometric Brownian Motion
#'
#' Computes the Monte Carlo pathwise estimator of the option delta for
#' geometric Brownian motion simulations. Supports both call and put European
#' options and accepts multiple strikes.
#'
#' @param paths Tibble of simulated paths produced by `simulate_paths()` for a
#'   GBM specification containing columns `path_id`, `time`, and `stock_price`.
#' @param strikes Numeric vector of strike prices.
#' @param maturity Numeric time to maturity in years. Defaults to the maximum
#'   `time` present in `paths`.
#' @param risk_free_rate Numeric continuously compounded risk-free rate.
#' @param option_type Character string, one of "call" or "put".
#' @param initial_price Optional numeric initial stock price. When omitted the
#'   value is retrieved from the path specification attribute.
#'
#' @return Tibble with columns `strike` and `delta_pathwise`.
#' @details
#'   Strikes are processed sequentially with `purrr::map()`. To parallelise the
#'   workload, wrap calls in your preferred parallel purrr backend or shard the
#'   strike vector across workers.
#' @examples
#' spec <- gbm_spec(initial_value = 100, drift = 0.05, volatility = 0.2)
#' paths <- simulate_paths(spec, n_paths = 5000, n_steps = 128, maturity = 1, seed = 42)
#' pathwise_delta(paths, strikes = c(90, 100, 110), risk_free_rate = 0.05)
#'
#' @export
pathwise_delta <- function(paths,
                           strikes,
                           maturity = max(paths$time, na.rm = TRUE),
                           risk_free_rate,
                           option_type = c("call", "put"),
                           initial_price = NULL) {
  option_type <- rlang::arg_match(option_type)
  checkmate::assert_data_frame(paths)
  checkmate::assert_numeric(strikes, lower = .Machine$double.eps, any.missing = FALSE, finite = TRUE)
  checkmate::assert_number(maturity, lower = .Machine$double.eps, finite = TRUE)
  checkmate::assert_number(risk_free_rate, finite = TRUE)

  if (is.null(initial_price)) {
    spec <- attr(paths, "spec")
    if (!is.null(spec) && !is.null(spec$initial_value)) {
      initial_price <- spec$initial_value
    } else if (!is.null(spec) && !is.null(spec$initial_price)) {
      initial_price <- spec$initial_price
    }
  }
  checkmate::assert_number(initial_price, lower = .Machine$double.eps, finite = TRUE)

  terminal_time <- maturity
  terminal_prices <- paths |>
    dplyr::filter(.data$time >= terminal_time - .Machine$double.eps) |>
    dplyr::pull(.data$stock_price)

  discount_factor <- exp(-risk_free_rate * maturity)
  sign_multiplier <- if (option_type == "call") 1 else -1

  delta_vals <- pathwise_delta_cpp(
    terminal_prices = terminal_prices,
    strikes = strikes,
    initial_price = initial_price,
    discount_factor = discount_factor,
    sign_multiplier = sign_multiplier,
    is_call = identical(option_type, "call")
  )

  tibble::tibble(
    strike = strikes,
    delta_pathwise = delta_vals
  )
}

#' Pathwise Vega Estimator for Geometric Brownian Motion
#'
#' Computes the Monte Carlo pathwise estimator of the option vega for
#' geometric Brownian motion simulations.
#'
#' @inheritParams pathwise_delta
#' @param volatility Numeric volatility used in the simulation.
#'
#' @return Tibble with columns `strike` and `vega_pathwise`.
#' @details
#'   Strikes are evaluated sequentially. Use your preferred parallel backend to
#'   distribute the work if lower latency is required.
#' @examples
#' spec <- gbm_spec(initial_value = 100, drift = 0.05, volatility = 0.2)
#' paths <- simulate_paths(spec, n_paths = 5000, n_steps = 128, maturity = 1, seed = 42)
#' pathwise_vega(paths, strikes = 100, risk_free_rate = 0.05, volatility = 0.2)
#'
#' @export
pathwise_vega <- function(paths,
                          strikes,
                          maturity = max(paths$time, na.rm = TRUE),
                          risk_free_rate,
                          option_type = c("call", "put"),
                          volatility,
                          initial_price = NULL) {
  option_type <- rlang::arg_match(option_type)
  checkmate::assert_data_frame(paths)
  checkmate::assert_numeric(strikes, lower = .Machine$double.eps, any.missing = FALSE, finite = TRUE)
  checkmate::assert_number(maturity, lower = .Machine$double.eps, finite = TRUE)
  checkmate::assert_number(risk_free_rate, finite = TRUE)
  checkmate::assert_number(volatility, lower = .Machine$double.eps, finite = TRUE)

  if (is.null(initial_price)) {
    spec <- attr(paths, "spec")
    if (!is.null(spec) && !is.null(spec$initial_value)) {
      initial_price <- spec$initial_value
    } else if (!is.null(spec) && !is.null(spec$initial_price)) {
      initial_price <- spec$initial_price
    }
  }
  checkmate::assert_number(initial_price, lower = .Machine$double.eps, finite = TRUE)

  terminal_time <- maturity
  terminal_prices <- paths |>
    dplyr::filter(.data$time >= terminal_time - .Machine$double.eps) |>
    dplyr::pull(.data$stock_price)

  discount_factor <- exp(-risk_free_rate * maturity)
  sign_multiplier <- if (option_type == "call") 1 else -1
  log_term <- log(terminal_prices / initial_price)
  adjustment <- log_term - (risk_free_rate + 0.5 * volatility^2) * maturity

  vega_vals <- pathwise_vega_cpp(
    terminal_prices = terminal_prices,
    adjustment = adjustment,
    strikes = strikes,
    volatility = volatility,
    discount_factor = discount_factor,
    sign_multiplier = sign_multiplier,
    is_call = identical(option_type, "call")
  )

  tibble::tibble(
    strike = strikes,
    vega_pathwise = vega_vals
  )
}

#' Compare Pathwise and Finite Difference Greek Estimators
#'
#' Generates Monte Carlo estimates for delta and vega using both the pathwise
#' method and central finite differences, returning a tidy comparison together
#' with analytical Black-Scholes benchmarks.
#'
#' @param process_spec A `gbm_spec` object describing the process dynamics.
#' @param strikes Numeric vector of strike prices.
#' @param maturity Numeric option maturity in years.
#' @param n_paths Integer number of Monte Carlo paths.
#' @param n_steps Integer number of time steps per path.
#' @param risk_free_rate Numeric risk-free rate. Defaults to the drift stored
#'   on the specification.
#' @param volatility Numeric volatility. Defaults to the volatility stored on
#'   the specification.
#' @param option_type Character string, one of "call" or "put".
#' @param bump_spot Numeric bump size applied to the initial price when computing
#'   finite-difference delta. Default is 0.1.
#' @param bump_vol Numeric bump size applied to the volatility when computing
#'   finite-difference vega. Default is 1e-3.
#' @param seed Optional integer random seed ensuring common random numbers across
#'   the simulations.
#'
#' @return Tibble containing columns `strike`, `delta_pathwise`, `delta_fd`,
#'   `delta_analytic`, `vega_pathwise`, `vega_fd`, and `vega_analytic`.
#' @details
#'   Strike-wise computations execute sequentially. For parallel workflows,
#'   invoke external purrr or future backends to distribute the strike grid.
#' @examples
#' spec <- gbm_spec(initial_value = 100, drift = 0.05, volatility = 0.2)
#' compare_with_finite_diff(spec, strikes = 100, maturity = 1, n_paths = 10000, n_steps = 128)
#'
#' @export
compare_with_finite_diff <- function(process_spec,
                                     strikes,
                                     maturity,
                                     n_paths,
                                     n_steps,
                                     risk_free_rate = process_spec$drift,
                                     volatility = process_spec$volatility,
                                     option_type = c("call", "put"),
                                     bump_spot = 0.1,
                                     bump_vol = 1e-3,
                                     seed = NULL) {
  option_type <- rlang::arg_match(option_type)
  checkmate::assert_class(process_spec, "gbm_spec")
  checkmate::assert_numeric(strikes, lower = .Machine$double.eps, any.missing = FALSE, finite = TRUE)
  checkmate::assert_number(maturity, lower = .Machine$double.eps, finite = TRUE)
  checkmate::assert_int(n_paths, lower = 1)
  checkmate::assert_int(n_steps, lower = 1)
  checkmate::assert_number(risk_free_rate, finite = TRUE)
  checkmate::assert_number(volatility, lower = .Machine$double.eps, finite = TRUE)
  checkmate::assert_number(bump_spot, lower = 0, finite = TRUE)
  checkmate::assert_number(bump_vol, lower = 0, finite = TRUE)
  if (!is.null(seed)) {
    checkmate::assert_int(seed)
  }

  run_simulation <- function(spec) {
    if (is.null(seed)) {
      simulate_paths(
        process_spec = spec,
        n_paths = n_paths,
        n_steps = n_steps,
        maturity = maturity
      )
    } else {
      simulate_paths(
        process_spec = spec,
        n_paths = n_paths,
        n_steps = n_steps,
        maturity = maturity,
        seed = seed
      )
    }
  }

  initial_price <- process_spec$initial_value
  base_paths <- run_simulation(process_spec)

  delta_pw <- pathwise_delta(
    base_paths,
    strikes = strikes,
    maturity = maturity,
    risk_free_rate = risk_free_rate,
    option_type = option_type,
    initial_price = initial_price
  )

  vega_pw <- pathwise_vega(
    base_paths,
    strikes = strikes,
    maturity = maturity,
    risk_free_rate = risk_free_rate,
    option_type = option_type,
    volatility = volatility,
    initial_price = initial_price
  )

  mc_price <- function(path_tbl, strike, option_type, maturity, risk_free_rate) {
    terminal_prices <- path_tbl |>
      dplyr::filter(.data$time >= maturity - .Machine$double.eps) |>
      dplyr::pull(.data$stock_price)
    payoff <- if (option_type == "call") {
      pmax(terminal_prices - strike, 0)
    } else {
      pmax(strike - terminal_prices, 0)
    }
    exp(-risk_free_rate * maturity) * mean(payoff)
  }

  gbm_engine <- simulate_paths.gbm_spec

  scenario_prices <- function(init_price, vol_level) {
    spec <- gbm_spec(
      initial_value = init_price,
      drift = risk_free_rate,
      volatility = vol_level
    )

    paths_tbl <- if (is.null(seed)) {
      gbm_engine(
        process_spec = spec,
        n_paths = n_paths,
        n_steps = n_steps,
        maturity = maturity
      )
    } else {
      gbm_engine(
        process_spec = spec,
        n_paths = n_paths,
        n_steps = n_steps,
        maturity = maturity,
        seed = seed
      )
    }

    purrr::map_dbl(
      strikes,
      \(strike) mc_price(paths_tbl, strike, option_type, maturity, risk_free_rate)
    )
  }

  if (bump_spot > 0) {
    if (initial_price - bump_spot <= 0) {
      rlang::abort("Initial price minus bump_spot must remain positive")
    }

    spot_scenarios <- list(
      list(initial = initial_price + bump_spot, volatility = volatility),
      list(initial = initial_price - bump_spot, volatility = volatility)
    )

    spot_prices <- purrr::map(
      spot_scenarios,
      \(scenario) {
        scenario_prices(scenario$initial, scenario$volatility)
      }
    )

    price_spot_up <- spot_prices[[1]]
    price_spot_down <- spot_prices[[2]]
    delta_fd_vec <- (price_spot_up - price_spot_down) / (2 * bump_spot)
  } else {
    delta_fd_vec <- rep(NA_real_, length(strikes))
  }

  if (bump_vol > 0) {
    if (volatility - bump_vol <= 0) {
      rlang::abort("Volatility minus bump_vol must remain positive")
    }

    vol_scenarios <- list(
      list(initial = initial_price, volatility = volatility + bump_vol),
      list(initial = initial_price, volatility = volatility - bump_vol)
    )

    vol_prices <- purrr::map(
      vol_scenarios,
      \(scenario) {
        scenario_prices(scenario$initial, scenario$volatility)
      }
    )

    price_vol_up <- vol_prices[[1]]
    price_vol_down <- vol_prices[[2]]
    vega_fd_vec <- (price_vol_up - price_vol_down) / (2 * bump_vol)
  } else {
    vega_fd_vec <- rep(NA_real_, length(strikes))
  }

  analytic_greeks <- purrr::map(
    strikes,
    \(strike) {
      bs_spec <- black_scholes_spec(
        option_type = option_type,
        strike = strike,
        maturity = maturity,
        risk_free_rate = risk_free_rate
      )
      price_options.black_scholes_spec(bs_spec, spot = initial_price, volatility = volatility) |>
        dplyr::mutate(strike = strike)
    }
  ) |> purrr::list_rbind()

  delta_pw |>
    dplyr::left_join(vega_pw, by = "strike") |>
    dplyr::mutate(
      delta_fd = delta_fd_vec,
      vega_fd = vega_fd_vec
    ) |>
    dplyr::left_join(
      analytic_greeks |>
        dplyr::select(strike, delta_analytic = delta, vega_analytic = vega),
      by = "strike"
    )
}
