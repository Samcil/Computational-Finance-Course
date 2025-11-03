#' Compute Black-Scholes Implied Volatility for Calls
#'
#' Estimates implied volatilities for European call options by inverting the
#' Black-Scholes pricing formula. The function accepts a tidy tibble of observed
#' option prices and returns the same tibble augmented with an
#' `implied_volatility` column.
#'
#' @param option_data A data frame containing at least the columns `strike` and
#'   `price` representing observed call prices.
#' @param spot Numeric scalar spot price \eqn{S_0 > 0}.
#' @param maturity Numeric time to maturity in years.
#' @param risk_free_rate Numeric continuously compounded risk-free rate.
#' @param dividend_yield Numeric continuous dividend yield. Default is 0.
#' @param lower Numeric lower bound for the implied volatility search interval.
#'   Default is `1e-4`.
#' @param upper Numeric upper bound for the implied volatility search interval.
#'   Default is `3`.
#' @param tol Numeric tolerance supplied to `uniroot()`. Default is `1e-6`.
#' @param max_iter Integer maximum number of iterations for `uniroot()`.
#'
#' @return A tibble with the original columns plus `implied_volatility`.
#' @examples
#' strikes <- c(90, 100, 110)
#' true_vols <- c(0.18, 0.2, 0.22)
#' prices <- purrr::map2_dbl(
#'   strikes,
#'   true_vols,
#'   ~ price_options(
#'       black_scholes_spec(
#'         option_type = "call",
#'         strike = .x,
#'         maturity = 1,
#'         risk_free_rate = 0.02
#'       ),
#'       spot = 100,
#'       volatility = .y
#'     )$price
#' )
#'
#' option_tbl <- tibble::tibble(strike = strikes, price = prices)
#' implied_volatility_call(
#'   option_tbl,
#'   spot = 100,
#'   maturity = 1,
#'   risk_free_rate = 0.02
#' )
#'
#' @export
implied_volatility_call <- function(option_data,
                                    spot,
                                    maturity,
                                    risk_free_rate,
                                    dividend_yield = 0,
                                    lower = 1e-4,
                                    upper = 3,
                                    tol = 1e-6,
                                    max_iter = 100) {
  compute_implied_volatility(
    option_data = option_data,
    option_type = "call",
    spot = spot,
    maturity = maturity,
    risk_free_rate = risk_free_rate,
    dividend_yield = dividend_yield,
    lower = lower,
    upper = upper,
    tol = tol,
    max_iter = max_iter
  )
}

#' Compute Black-Scholes Implied Volatility for Puts
#'
#' Estimates implied volatilities for European put options by inverting the
#' Black-Scholes pricing formula.
#'
#' @inheritParams implied_volatility_call
#'
#' @return A tibble with the original columns plus `implied_volatility`.
#' @examples
#' strikes <- c(90, 100, 110)
#' true_vols <- c(0.18, 0.2, 0.22)
#' prices <- purrr::map2_dbl(
#'   strikes,
#'   true_vols,
#'   ~ price_options(
#'       black_scholes_spec(
#'         option_type = "put",
#'         strike = .x,
#'         maturity = 1,
#'         risk_free_rate = 0.02
#'       ),
#'       spot = 100,
#'       volatility = .y
#'     )$price
#' )
#'
#' option_tbl <- tibble::tibble(strike = strikes, price = prices)
#' implied_volatility_put(
#'   option_tbl,
#'   spot = 100,
#'   maturity = 1,
#'   risk_free_rate = 0.02
#' )
#'
#' @export
implied_volatility_put <- function(option_data,
                                   spot,
                                   maturity,
                                   risk_free_rate,
                                   dividend_yield = 0,
                                   lower = 1e-4,
                                   upper = 3,
                                   tol = 1e-6,
                                   max_iter = 100) {
  compute_implied_volatility(
    option_data = option_data,
    option_type = "put",
    spot = spot,
    maturity = maturity,
    risk_free_rate = risk_free_rate,
    dividend_yield = dividend_yield,
    lower = lower,
    upper = upper,
    tol = tol,
    max_iter = max_iter
  )
}

#' Plot Implied Volatility Smile
#'
#' Produces a scatter and optional line plot of implied volatility against
#' strike, supporting faceted or coloured series via an optional column name.
#'
#' @param volatility_data A data frame returned by
#'   `implied_volatility_call()` or `implied_volatility_put()` containing at
#'   least `strike` and `implied_volatility`.
#' @param colour_var Optional string naming a column in `volatility_data` used
#'   to colour the smile (e.g., maturity or option type).
#' @param line Logical; if `TRUE` (default) a line is drawn through the points.
#' @param point_size Numeric size of the plotted points. Default is 2.
#'
#' @return A `ggplot` object visualising the volatility smile.
#' @export
plot_volatility_smile <- function(volatility_data,
                                  colour_var = NULL,
                                  line = TRUE,
                                  point_size = 2) {
  checkmate::assert_data_frame(volatility_data)
  required_cols <- c("strike", "implied_volatility")
  missing_cols <- setdiff(required_cols, names(volatility_data))
  if (length(missing_cols) > 0) {
    rlang::abort("volatility_data must contain columns: 'strike', 'implied_volatility'")
  }

  checkmate::assert_flag(line)
  checkmate::assert_number(point_size, lower = 0, finite = TRUE)

  plot_tbl <- tibble::as_tibble(volatility_data)

  mapping <- if (is.null(colour_var)) {
    ggplot2::aes(x = strike, y = implied_volatility)
  } else {
    checkmate::assert_choice(colour_var, names(plot_tbl))
    colour_sym <- rlang::sym(colour_var)
    ggplot2::aes(x = strike, y = implied_volatility, colour = !!colour_sym)
  }

  smile_plot <- ggplot2::ggplot(plot_tbl, mapping) +
    ggplot2::geom_point(size = point_size)

  if (line) {
    smile_plot <- smile_plot + ggplot2::geom_line()
  }

  smile_plot +
    ggplot2::labs(
      x = "Strike",
      y = "Implied volatility",
      colour = if (!is.null(colour_var)) colour_var else NULL,
      title = "Implied Volatility Smile"
    ) +
    ggplot2::theme_minimal()
}

compute_implied_volatility <- function(option_data,
                                       option_type,
                                       spot,
                                       maturity,
                                       risk_free_rate,
                                       dividend_yield,
                                       lower,
                                       upper,
                                       tol,
                                       max_iter) {
  option_tbl <- tibble::as_tibble(option_data)
  required_cols <- c("strike", "price")
  missing_cols <- setdiff(required_cols, names(option_tbl))
  if (length(missing_cols) > 0) {
    rlang::abort("option_data must contain columns: 'strike', 'price'")
  }

  checkmate::assert_number(spot, lower = .Machine$double.eps, finite = TRUE)
  checkmate::assert_number(maturity, lower = .Machine$double.eps, finite = TRUE)
  checkmate::assert_number(risk_free_rate, finite = TRUE)
  checkmate::assert_number(dividend_yield, finite = TRUE)
  checkmate::assert_numeric(option_tbl$strike, lower = .Machine$double.eps, any.missing = FALSE, finite = TRUE)
  checkmate::assert_numeric(option_tbl$price, lower = 0, any.missing = FALSE, finite = TRUE)
  checkmate::assert_number(lower, lower = 0, finite = TRUE)
  checkmate::assert_number(upper, lower = lower, finite = TRUE)
  checkmate::assert_number(tol, lower = 0, finite = TRUE)
  checkmate::assert_integerish(max_iter, lower = 1, len = 1)

  implied <- purrr::map2_dbl(
    option_tbl$price,
    option_tbl$strike,
    \(target_price, strike) {
      solve_implied_volatility(
        option_type = option_type,
        target_price = target_price,
        strike = strike,
        spot = spot,
        maturity = maturity,
        risk_free_rate = risk_free_rate,
        dividend_yield = dividend_yield,
        lower = lower,
        upper = upper,
        tol = tol,
        max_iter = max_iter,
        price_options_fn = price_options.black_scholes_spec,
        black_scholes_spec_fn = black_scholes_spec
      )
    }
  )

  dplyr::mutate(option_tbl, implied_volatility = implied)
}

solve_implied_volatility <- function(option_type,
                                     target_price,
                                     strike,
                                     spot,
                                     maturity,
                                     risk_free_rate,
                                     dividend_yield,
                                     lower,
                                     upper,
                                     tol,
                                     max_iter,
                                     price_options_fn = price_options.black_scholes_spec,
                                     black_scholes_spec_fn = black_scholes_spec) {
  forward_factor <- exp(-dividend_yield * maturity)
  discount_factor <- exp(-risk_free_rate * maturity)

  intrinsic <- if (option_type == "call") {
    max(spot * forward_factor - strike * discount_factor, 0)
  } else {
    max(strike * discount_factor - spot * forward_factor, 0)
  }
  upper_bound <- if (option_type == "call") {
    spot * forward_factor
  } else {
    strike * discount_factor
  }

  if (is.na(target_price) || target_price < intrinsic - 1e-8 || target_price > upper_bound + 1e-8) {
    return(NA_real_)
  }

  bs_spec <- black_scholes_spec_fn(
    option_type = option_type,
    strike = strike,
    maturity = maturity,
    risk_free_rate = risk_free_rate,
    dividend_yield = dividend_yield
  )

  pricing_difference <- function(vol) {
    price_options_fn(bs_spec, spot = spot, volatility = vol)$price - target_price
  }

  lower_eval <- pricing_difference(lower)
  upper_eval <- pricing_difference(upper)

  if (abs(lower_eval) < tol) {
    return(lower)
  }
  if (abs(upper_eval) < tol) {
    return(upper)
  }

  if (lower_eval * upper_eval > 0) {
    expansion <- c(2, 5, 10)
    expansion_results <- purrr::map(
      expansion,
      \(mult) {
        candidate <- upper * mult
        value <- pricing_difference(candidate)
        tibble::tibble(
          candidate = candidate,
          value = value,
          valid = lower_eval * value <= 0
        )
      }
    ) |> dplyr::bind_rows()

    valid_candidate <- expansion_results |>
      dplyr::filter(valid) |>
      dplyr::slice_head(n = 1)

    if (nrow(valid_candidate) > 0) {
      upper <- valid_candidate$candidate
      upper_eval <- valid_candidate$value
    }
  }

  if (lower_eval * upper_eval > 0) {
    return(NA_real_)
  }

  tryCatch(
    stats::uniroot(pricing_difference, lower = lower, upper = upper, tol = tol, maxiter = max_iter)$root,
    error = function(...) NA_real_
  )
}
