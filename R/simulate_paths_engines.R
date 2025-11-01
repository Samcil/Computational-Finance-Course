#' Simulation Engines for Stochastic Processes
#'
#' Internal implementation functions (engines) for different stochastic processes.
#' These functions handle the computational complexity and are called via dispatch
#' from the user-facing simulate_paths() function.
#'
#' Following hardhat principles, these engines are separated from the user interface.
#'
#' @name simulation_engines
#' @keywords internal
NULL


#' Simulate GBM Paths (Engine Implementation)
#'
#' @keywords internal
#' @export
simulate_paths.gbm_spec <- function(process_spec,
                                    n_paths,
                                    n_steps,
                                    maturity,
                                    seed = 123,
                                    ...) {
  # Set seed
  set.seed(seed)

  # Extract parameters
  S0 <- process_spec$initial_value
  r <- process_spec$drift
  sigma <- process_spec$volatility

  # Time step
  dt <- maturity / n_steps

  # Generate standardized random samples
  Z <- generate_standardized_normals(n_paths, n_steps)

  # Compute log-price increments
  log_increments <- (r - 0.5 * sigma^2) * dt + sigma * sqrt(dt) * Z

  # Cumulative sum along time dimension using pure functional approach
  log_prices <- compute_cumulative_paths(
    initial_value = log(S0),
    increments = log_increments
  )

  # Convert to stock prices
  stock_prices <- exp(log_prices)

  # Create time grid
  time_grid <- seq(0, maturity, length.out = n_steps + 1)

  # Convert to tidy tibble using modern purrr
  paths_tidy <- matrix_to_tidy(
    value_matrix = stock_prices,
    time_grid = time_grid,
    value_name = "stock_price"
  )

  # Add process type attribute
  attr(paths_tidy, "process_type") <- "gbm"
  attr(paths_tidy, "spec") <- process_spec

  paths_tidy
}


#' Simulate ABM Paths (Engine Implementation)
#'
#' @keywords internal
#' @export
simulate_paths.abm_spec <- function(process_spec,
                                    n_paths,
                                    n_steps,
                                    maturity,
                                    seed = 123,
                                    ...) {
  # Set seed
  set.seed(seed)

  # Extract parameters
  X0 <- process_spec$initial_value
  mu <- process_spec$drift
  sigma <- process_spec$volatility

  # Time step
  dt <- maturity / n_steps

  # Generate standardized random samples
  Z <- generate_standardized_normals(n_paths, n_steps)

  # Compute increments
  increments <- mu * dt + sigma * sqrt(dt) * Z

  # Cumulative sum along time dimension
  paths <- compute_cumulative_paths(
    initial_value = X0,
    increments = increments
  )

  # Create time grid
  time_grid <- seq(0, maturity, length.out = n_steps + 1)

  # Convert to tidy tibble
  paths_tidy <- matrix_to_tidy(
    value_matrix = paths,
    time_grid = time_grid,
    value_name = "value"
  )

  # Add process type attribute
  attr(paths_tidy, "process_type") <- "abm"
  attr(paths_tidy, "spec") <- process_spec

  paths_tidy
}


#' Generate Standardized Normal Random Variables
#'
#' Generates a matrix of random normal samples and standardizes each time step
#' to have mean 0 and variance 1.
#'
#' @param n_paths Integer. Number of paths.
#' @param n_steps Integer. Number of time steps.
#'
#' @return Matrix of size (n_paths x n_steps) with standardized normal samples.
#'
#' @keywords internal
generate_standardized_normals <- function(n_paths, n_steps) {
  # Generate and standardize samples column-wise using vectorised scale
  Z <- matrix(
    stats::rnorm(n_paths * n_steps),
    nrow = n_paths,
    ncol = n_steps
  )

  if (n_paths > 1) {
    Z <- scale(Z, center = TRUE, scale = TRUE)
  }

  Z
}


#' Compute Cumulative Paths from Increments
#'
#' Efficiently computes cumulative sums along time dimension for multiple paths
#' using functional programming.
#'
#' @param initial_value Numeric. Initial value (scalar or vector of length n_paths).
#' @param increments Matrix. Increments of size (n_paths x n_steps).
#'
#' @return Matrix of cumulative values including initial value (n_paths x n_steps+1).
#'
#' @keywords internal
compute_cumulative_paths <- function(initial_value, increments) {
  n_paths <- nrow(increments)
  n_steps <- ncol(increments)

  X0 <- if (length(initial_value) == 1) {
    rep(initial_value, n_paths)
  } else {
    initial_value
  }

  path_history <- purrr::accumulate(
    .x = seq_len(n_steps),
    .init = X0,
    .f = \(state, step_idx) state + increments[, step_idx]
  )

  do.call(cbind, path_history) |> as.matrix()
}


#' Convert Matrix to Tidy Tibble
#'
#' Converts a matrix of paths to tidy (long format) tibble.
#' Uses modern purrr::list_rbind instead of deprecated map_dfr.
#'
#' @param value_matrix Matrix. Values of size (n_paths x n_timesteps).
#' @param time_grid Numeric vector. Time points.
#' @param value_name Character. Name for the value column.
#'
#' @return Tibble in long format with path_id, time, and value column.
#'
#' @keywords internal
matrix_to_tidy <- function(value_matrix, time_grid, value_name = "value") {
  value_matrix <- as.matrix(value_matrix)
  n_paths <- nrow(value_matrix)

  # Create list of tibbles, one per path
  # Using anonymous function \(x) syntax (R >= 4.1.0)
  path_tibbles <- purrr::map(
    seq_len(n_paths),
    \(i) {
      path_tbl <- tibble::tibble(
        path_id = i,
        time = time_grid,
        value = value_matrix[i, ]
      )
      stats::setNames(path_tbl, c("path_id", "time", value_name))
    }
  )

  # Use modern list_rbind instead of map_dfr
  purrr::list_rbind(path_tibbles)
}
