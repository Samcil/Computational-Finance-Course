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
#'   Work is distributed with `purrr::in_parallel()`, so the calculations run
#'   sequentially by default and switch to parallel execution when callers start
#'   mirai daemons (e.g. `mirai::daemons(6)` before the call and
#'   `mirai::daemons(0)` afterwards).
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

  purrr::map(
    strikes,
    purrr::in_parallel(
      \(strike) {
        indicator <- if (option_type == "call") {
          terminal_prices > strike
        } else {
          terminal_prices < strike
        }
        estimator <- discount_factor * sign_multiplier * mean((terminal_prices / initial_price) * indicator)
        tibble::tibble(strike = strike, delta_pathwise = estimator)
      },
      terminal_prices = terminal_prices,
      discount_factor = discount_factor,
      sign_multiplier = sign_multiplier,
      initial_price = initial_price,
      option_type = option_type
    )
  ) |> purrr::list_rbind()
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
#'   Work is assigned through `purrr::in_parallel()`, remaining sequential unless
#'   mirai daemons are running. Start workers with `mirai::daemons(n)` before
#'   calling and shut them down with `mirai::daemons(0)` afterwards to leverage
#'   multiple cores.
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

  purrr::map(
    strikes,
    purrr::in_parallel(
      \(strike) {
        indicator <- if (option_type == "call") {
          terminal_prices > strike
        } else {
          terminal_prices < strike
        }
        estimator <- discount_factor * sign_multiplier * mean((terminal_prices / volatility) * adjustment * indicator)
        tibble::tibble(strike = strike, vega_pathwise = estimator)
      },
      terminal_prices = terminal_prices,
      discount_factor = discount_factor,
      sign_multiplier = sign_multiplier,
      volatility = volatility,
      adjustment = adjustment,
      option_type = option_type
    )
  ) |> purrr::list_rbind()
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
#'   Strike-wise computations run through `purrr::in_parallel()`. Without active
#'   mirai daemons the routine stays sequential; start workers with
#'   `mirai::daemons(n)` beforehand and stop them with `mirai::daemons(0)` when
#'   finished to parallelise the workflow.
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

  simulate_gbm <- function(initial, vol) {
    gbm_spec(
      initial_value = initial,
      drift = risk_free_rate,
      volatility = vol
    ) |>
      run_simulation()
  }

  if (bump_spot > 0) {
    if (initial_price - bump_spot <= 0) {
      rlang::abort("Initial price minus bump_spot must remain positive")
    }
    paths_up <- simulate_gbm(initial_price + bump_spot, volatility)
    paths_down <- simulate_gbm(initial_price - bump_spot, volatility)
    price_spot_up <- purrr::map_dbl(
      strikes,
      purrr::in_parallel(
        \(strike) {
          mc_price(paths_up, strike, option_type, maturity, risk_free_rate)
        },
        paths_up = paths_up,
        option_type = option_type,
        maturity = maturity,
        risk_free_rate = risk_free_rate,
        mc_price = mc_price
      )
    )
    price_spot_down <- purrr::map_dbl(
      strikes,
      purrr::in_parallel(
        \(strike) {
          mc_price(paths_down, strike, option_type, maturity, risk_free_rate)
        },
        paths_down = paths_down,
        option_type = option_type,
        maturity = maturity,
        risk_free_rate = risk_free_rate,
        mc_price = mc_price
      )
    )
    delta_fd_vec <- (price_spot_up - price_spot_down) / (2 * bump_spot)
  } else {
    delta_fd_vec <- rep(NA_real_, length(strikes))
  }

  if (bump_vol > 0) {
    if (volatility - bump_vol <= 0) {
      rlang::abort("Volatility minus bump_vol must remain positive")
    }
    paths_vol_up <- simulate_gbm(initial_price, volatility + bump_vol)
    paths_vol_down <- simulate_gbm(initial_price, volatility - bump_vol)
    price_vol_up <- purrr::map_dbl(
      strikes,
      purrr::in_parallel(
        \(strike) {
          mc_price(paths_vol_up, strike, option_type, maturity, risk_free_rate)
        },
        paths_vol_up = paths_vol_up,
        option_type = option_type,
        maturity = maturity,
        risk_free_rate = risk_free_rate,
        mc_price = mc_price
      )
    )
    price_vol_down <- purrr::map_dbl(
      strikes,
      purrr::in_parallel(
        \(strike) {
          mc_price(paths_vol_down, strike, option_type, maturity, risk_free_rate)
        },
        paths_vol_down = paths_vol_down,
        option_type = option_type,
        maturity = maturity,
        risk_free_rate = risk_free_rate,
        mc_price = mc_price
      )
    )
    vega_fd_vec <- (price_vol_up - price_vol_down) / (2 * bump_vol)
  } else {
    vega_fd_vec <- rep(NA_real_, length(strikes))
  }

  price_options_fn <- price_options.black_scholes_spec

  analytic_greeks <- purrr::map(
    strikes,
    purrr::in_parallel(
      \(strike) {
        bs_spec <- black_scholes_spec(
          option_type = option_type,
          strike = strike,
          maturity = maturity,
          risk_free_rate = risk_free_rate
        )
        price_options_fn(bs_spec, spot = initial_price, volatility = volatility) |>
          dplyr::mutate(strike = strike)
      },
      option_type = option_type,
      maturity = maturity,
      risk_free_rate = risk_free_rate,
      initial_price = initial_price,
      volatility = volatility,
      black_scholes_spec = black_scholes_spec,
      price_options_fn = price_options.black_scholes_spec
    )
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
