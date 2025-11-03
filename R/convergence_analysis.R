#' Euler Convergence Study for Geometric Brownian Motion
#'
#' Estimates the weak convergence rate of the Euler-Maruyama discretisation for
#' geometric Brownian motion by comparing simulated terminal prices with the
#' exact solution using common random numbers.
#'
#' @param initial_price Numeric. Initial price \eqn{S_0 > 0}.
#' @param drift Numeric. Drift coefficient \eqn{\mu}.
#' @param volatility Numeric. Volatility coefficient \eqn{\sigma > 0}.
#' @param maturity Numeric. Time horizon in years.
#' @param n_paths Integer. Number of Monte Carlo paths.
#' @param step_grid Integer vector of time step counts to analyse.
#' @param seed Optional integer random seed for reproducibility.
#'
#' @return A tibble with columns `scheme`, `n_steps`, `dt`, and `rmse`.
#' @export
euler_convergence_study <- function(initial_price,
                                    drift,
                                    volatility,
                                    maturity,
                                    n_paths,
                                    step_grid,
                                    seed = NULL) {
  checkmate::assert_number(initial_price, lower = .Machine$double.eps, finite = TRUE)
  checkmate::assert_number(drift, finite = TRUE)
  checkmate::assert_number(volatility, lower = 0, finite = TRUE)
  checkmate::assert_number(maturity, lower = .Machine$double.eps, finite = TRUE)
  checkmate::assert_integerish(n_paths, lower = 1, len = 1)
  checkmate::assert_integerish(step_grid, lower = 1, any.missing = FALSE)
  if (!is.null(seed)) {
    checkmate::assert_integerish(seed, len = 1)
  }

  purrr::map2(
    .x = as.integer(step_grid),
    .y = seq_along(step_grid),
    .f = \(n_steps, idx) {
      if (!is.null(seed)) {
        set.seed(seed + idx - 1)
      }
      dt <- maturity / n_steps
      dW <- sqrt(dt) * matrix(stats::rnorm(n_paths * n_steps), nrow = n_paths, ncol = n_steps)
      paths <- matrix(initial_price, nrow = n_paths, ncol = n_steps + 1)

      paths <- purrr::reduce(
        .x = seq_len(n_steps),
        .init = paths,
        .f = \(state, step) {
          current_state <- state[, step]
          increment <- drift * current_state * dt + volatility * current_state * dW[, step]
          state[, step + 1] <- current_state + increment
          state
        }
      )

      terminal_euler <- paths[, n_steps + 1]
      brownian_terminal <- rowSums(dW)
      terminal_exact <- initial_price * exp((drift - 0.5 * volatility^2) * maturity + volatility * brownian_terminal)
      rmse <- sqrt(mean((terminal_euler - terminal_exact)^2))

      tibble::tibble(
        scheme = "euler",
        n_steps = n_steps,
        dt = dt,
        rmse = rmse
      )
    }
  ) |> dplyr::bind_rows()
}

#' Milstein Convergence Study for Geometric Brownian Motion
#'
#' Estimates the convergence of the Milstein discretisation by benchmarking
#' against the exact geometric Brownian motion solution using common random
#' numbers.
#'
#' @inheritParams euler_convergence_study
#'
#' @return A tibble with columns `scheme`, `n_steps`, `dt`, and `rmse`.
#' @export
milstein_convergence_study <- function(initial_price,
                                       drift,
                                       volatility,
                                       maturity,
                                       n_paths,
                                       step_grid,
                                       seed = NULL) {
  checkmate::assert_number(initial_price, lower = .Machine$double.eps, finite = TRUE)
  checkmate::assert_number(drift, finite = TRUE)
  checkmate::assert_number(volatility, lower = 0, finite = TRUE)
  checkmate::assert_number(maturity, lower = .Machine$double.eps, finite = TRUE)
  checkmate::assert_integerish(n_paths, lower = 1, len = 1)
  checkmate::assert_integerish(step_grid, lower = 1, any.missing = FALSE)
  if (!is.null(seed)) {
    checkmate::assert_integerish(seed, len = 1)
  }

  purrr::map2(
    .x = as.integer(step_grid),
    .y = seq_along(step_grid),
    .f = \(n_steps, idx) {
      if (!is.null(seed)) {
        set.seed(seed + idx - 1)
      }
      dt <- maturity / n_steps
      dW <- sqrt(dt) * matrix(stats::rnorm(n_paths * n_steps), nrow = n_paths, ncol = n_steps)
      paths <- matrix(initial_price, nrow = n_paths, ncol = n_steps + 1)

      paths <- purrr::reduce(
        .x = seq_len(n_steps),
        .init = paths,
        .f = \(state, step) {
          current_state <- state[, step]
          dW_step <- dW[, step]
          increment <- drift * current_state * dt + volatility * current_state * dW_step +
            0.5 * volatility^2 * current_state * (dW_step^2 - dt)
          state[, step + 1] <- current_state + increment
          state
        }
      )

      terminal_milstein <- paths[, n_steps + 1]
      brownian_terminal <- rowSums(dW)
      terminal_exact <- initial_price * exp((drift - 0.5 * volatility^2) * maturity + volatility * brownian_terminal)
      rmse <- sqrt(mean((terminal_milstein - terminal_exact)^2))

      tibble::tibble(
        scheme = "milstein",
        n_steps = n_steps,
        dt = dt,
        rmse = rmse
      )
    }
  ) |> dplyr::bind_rows()
}

#' Plot Convergence Rates
#'
#' Produces a log-log plot of root mean square errors against time step sizes for
#' Euler and Milstein convergence studies, optionally overlaying a reference
#' slope.
#'
#' @param convergence_results Tibble returned by `euler_convergence_study()` or
#'   `milstein_convergence_study()` (or their row-bound combination).
#' @param reference_rate Numeric slope of the reference line in log-log space.
#'   Default is 1.
#' @param show_reference Logical. If `TRUE`, adds a dashed reference line.
#'
#' @return A `ggplot` object.
#' @export
plot_convergence_rates <- function(convergence_results,
                                   reference_rate = 1,
                                   show_reference = TRUE) {
  checkmate::assert_data_frame(convergence_results)
  required_cols <- c("scheme", "dt", "rmse")
  missing_cols <- setdiff(required_cols, names(convergence_results))
  if (length(missing_cols) > 0) {
    rlang::abort("convergence_results must contain columns: scheme, dt, rmse")
  }
  checkmate::assert_number(reference_rate, finite = TRUE)
  checkmate::assert_flag(show_reference)

  plot_data <- convergence_results |>
    dplyr::mutate(
      scheme = factor(scheme)
    )

  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = dt, y = rmse, colour = scheme)) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::scale_x_log10() +
    ggplot2::scale_y_log10() +
    ggplot2::labs(
      x = "Time step (dt)",
      y = "RMSE",
      colour = "Scheme"
    ) +
    ggplot2::theme_minimal()

  if (show_reference) {
    dt_range <- range(plot_data$dt)
    anchor_idx <- which.min(plot_data$dt)
    anchor_rmse <- plot_data$rmse[anchor_idx]
    anchor_dt <- plot_data$dt[anchor_idx]
    reference_const <- anchor_rmse / (anchor_dt^reference_rate)
    reference_df <- tibble::tibble(
      dt = dt_range,
      rmse = reference_const * (dt_range^reference_rate)
    )

    p <- p + ggplot2::geom_line(
      data = reference_df,
      ggplot2::aes(x = dt, y = rmse),
      inherit.aes = FALSE,
      colour = "grey30",
      linetype = "dashed"
    )
  }

  p
}
