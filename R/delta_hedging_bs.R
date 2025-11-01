#' Simulate Delta Hedging Strategy under Black-Scholes
#'
#' Simulates discrete-time delta hedging of a European option under the
#' Black-Scholes model. The function returns the time series of hedge ratios,
#' option and portfolio values, and the terminal hedging error.
#'
#' @param spot_paths Tibble of simulated underlying paths produced by
#'   [simulate_paths()] with columns `path_id`, `time`, and `stock_price`.
#' @param option_spec A `black_scholes_spec` describing the option contract to
#'   be hedged.
#' @param volatility Numeric. Assumed (constant) volatility used in the hedge.
#' @param hedge_steps Integer. Number of discrete hedging rebalancing steps
#'   between time 0 and maturity. Must align with the provided path grid.
#' @param transaction_cost Numeric proportional transaction cost applied on each
#'   trade. Default is 0 (no costs).
#' @param initial_capital Optional numeric specifying the initial capital in the
#'   hedging portfolio. Defaults to the Black-Scholes option price at time zero.
#'
#' @return A tibble containing `path_id`, `time`, `stock_price`, `option_value`,
#'   `delta`, `bond_position`, `portfolio_value`, and terminal `hedging_error`.
#'
#' @export
simulate_delta_hedge_bs <- function(spot_paths,
                                    option_spec,
                                    volatility,
                                    hedge_steps,
                                    transaction_cost = 0,
                                    initial_capital = NULL) {
  checkmate::assert_data_frame(spot_paths)
  required_cols <- c("path_id", "time", "stock_price")
  checkmate::assert_subset(required_cols, names(spot_paths))
  checkmate::assert_class(option_spec, "black_scholes_spec")
  checkmate::assert_number(volatility, lower = .Machine$double.eps, finite = TRUE)
  checkmate::assert_int(hedge_steps, lower = 1)
  checkmate::assert_number(transaction_cost, lower = 0, finite = TRUE)

  spot_paths <- spot_paths |>
    dplyr::mutate(
      path_id = as.integer(.data$path_id)
    )

  maturity <- option_spec$maturity
  risk_free_rate <- option_spec$risk_free_rate
  dividend_yield <- option_spec$dividend_yield

  time_grid <- spot_paths |>
    dplyr::distinct(.data$time) |>
    dplyr::arrange(.data$time) |>
    dplyr::pull(.data$time)

  if (!dplyr::near(dplyr::first(time_grid), 0, tol = 1e-10) ||
    !dplyr::near(dplyr::last(time_grid), maturity, tol = 1e-8)) {
    rlang::abort("Path times must start at 0 and end at option maturity")
  }

  rebalancing_times <- seq(0, maturity, length.out = hedge_steps + 1)
  rebalancing_on_grid <- purrr::map_lgl(
    rebalancing_times,
    ~ any(abs(.x - time_grid) <= 1e-10)
  )
  if (!all(rebalancing_on_grid)) {
    rlang::abort("Hedging times must coincide with the provided path grid")
  }

  if (is.null(initial_capital)) {
    initial_spot <- spot_paths |>
      dplyr::filter(.data$time == 0) |>
      dplyr::slice_head(n = 1) |>
      dplyr::pull(.data$stock_price)
    initial_capital <- price_options(
      option_spec,
      spot = initial_spot,
      volatility = volatility
    ) |>
      dplyr::pull(.data$price)
  }

  spot_paths |>
    dplyr::group_by(.data$path_id) |>
    dplyr::arrange(.data$time, .by_group = TRUE) |>
    dplyr::group_modify(
      ~ simulate_delta_hedge_path(
        path_tbl = .x,
        option_spec = option_spec,
        volatility = volatility,
        rebalancing_times = rebalancing_times,
        initial_capital = initial_capital,
        risk_free_rate = risk_free_rate,
        dividend_yield = dividend_yield,
        transaction_cost = transaction_cost
      )
    ) |>
    dplyr::ungroup()
}

simulate_delta_hedge_path <- function(path_tbl,
                                      option_spec,
                                      volatility,
                                      rebalancing_times,
                                      initial_capital,
                                      risk_free_rate,
                                      dividend_yield,
                                      transaction_cost) {
  times <- path_tbl$time
  stock_prices <- path_tbl$stock_price
  n_steps <- length(times)

  option_values <- numeric(n_steps)
  theoretical_delta <- numeric(n_steps)
  held_delta <- numeric(n_steps)
  bond_position <- numeric(n_steps)
  portfolio_value <- numeric(n_steps)
  hedging_error <- rep(NA_real_, n_steps)

  cash_account <- initial_capital
  current_delta <- 0

  for (idx in seq_len(n_steps)) {
    if (idx > 1) {
      dt <- times[idx] - times[idx - 1]
      cash_account <- cash_account * exp(risk_free_rate * dt)
      if (!dplyr::near(dividend_yield, 0, tol = 1e-12)) {
        cash_account <- cash_account + current_delta * stock_prices[idx - 1] * dividend_yield * dt
      }
    }

    tau <- max(option_spec$maturity - times[idx], 0)
    metrics <- bs_option_metrics(
      spot = stock_prices[idx],
      remaining_time = tau,
      option_spec = option_spec,
      volatility = volatility
    )

    option_values[idx] <- metrics$price
    theoretical_delta[idx] <- metrics$delta

    rebalance_now <- idx < n_steps && any(abs(times[idx] - rebalancing_times) <= 1e-10)
    if (idx == 1 || rebalance_now) {
      delta_target <- theoretical_delta[idx]
      trade_size <- delta_target - current_delta
      if (!dplyr::near(trade_size, 0, tol = 1e-12)) {
        trade_notional <- trade_size * stock_prices[idx]
        cash_account <- cash_account - trade_notional - transaction_cost * abs(trade_notional)
        current_delta <- delta_target
      }
    }

    held_delta[idx] <- current_delta
    bond_position[idx] <- cash_account
    portfolio_value[idx] <- cash_account + current_delta * stock_prices[idx]
  }

  payoff <- bs_option_payoff(stock_prices[length(stock_prices)], option_spec)
  hedging_error[n_steps] <- portfolio_value[n_steps] - payoff

  path_tbl |>
    dplyr::mutate(
      option_value = option_values,
      delta = held_delta,
      bond_position = bond_position,
      portfolio_value = portfolio_value,
      hedging_error = hedging_error
    )
}

bs_option_metrics <- function(spot, remaining_time, option_spec, volatility) {
  if (remaining_time <= 1e-10) {
    payoff <- bs_option_payoff(spot, option_spec)
    delta <- bs_terminal_delta(spot, option_spec)
    return(list(price = payoff, delta = delta))
  }

  strike <- option_spec$strike
  r <- option_spec$risk_free_rate
  q <- option_spec$dividend_yield
  option_type <- option_spec$option_type

  sqrt_time <- sqrt(remaining_time)
  forward_factor <- exp(-q * remaining_time)
  discount_factor <- exp(-r * remaining_time)

  d1 <- (
    log(spot / strike) + (r - q + 0.5 * volatility^2) * remaining_time
  ) / (volatility * sqrt_time)
  d2 <- d1 - volatility * sqrt_time

  if (option_type == "call") {
    price <- spot * forward_factor * stats::pnorm(d1) -
      strike * discount_factor * stats::pnorm(d2)
    delta <- forward_factor * stats::pnorm(d1)
  } else {
    price <- strike * discount_factor * stats::pnorm(-d2) -
      spot * forward_factor * stats::pnorm(-d1)
    delta <- forward_factor * (stats::pnorm(d1) - 1)
  }

  list(price = price, delta = delta)
}

bs_option_payoff <- function(spot, option_spec) {
  strike <- option_spec$strike
  if (option_spec$option_type == "call") {
    return(max(spot - strike, 0))
  }
  max(strike - spot, 0)
}

bs_terminal_delta <- function(spot, option_spec) {
  strike <- option_spec$strike
  if (option_spec$option_type == "call") {
    if (spot > strike) {
      return(1)
    }
    if (spot < strike) {
      return(0)
    }
    return(0.5)
  }
  if (spot < strike) {
    return(-1)
  }
  if (spot > strike) {
    return(0)
  }
  -0.5
}
