#' Create Poisson Process Specification
#'
#' Specification for standard and compensated Poisson process simulation.
#'
#' @param intensity Numeric. Jump intensity \eqn{\lambda} (average jumps per unit time).
#' @param initial_value Numeric. Initial count \eqn{N_0}. Default is 0.
#'
#' @return An S3 object of class `poisson_spec` containing the process parameters.
#'
#' @details
#' ## Mathematical Formulation
#'
#' The Poisson process \eqn{N(t)} counts the number of jumps by time \eqn{t}
#' with intensity \eqn{\lambda}:
#' \deqn{P(N(t) = k) = \frac{(\lambda t)^k e^{-\lambda t}}{k!}}
#'
#' Properties:
#' \itemize{
#'   \item \eqn{E[N(t)] = \lambda t}
#'   \item \eqn{Var[N(t)] = \lambda t}
#'   \item Jump times follow exponential distribution with rate \eqn{\lambda}
#'   \item Independent increments
#'   \item Stationary increments
#' }
#'
#' The **compensated Poisson process** is a martingale:
#' \deqn{\tilde{N}(t) = N(t) - \lambda t}
#'
#' with properties:
#' \itemize{
#'   \item \eqn{E[\tilde{N}(t)] = 0} (zero mean)
#'   \item \eqn{Var[\tilde{N}(t)] = \lambda t}
#'   \item Martingale property: \eqn{E[\tilde{N}(t) | \mathcal{F}_s] = \tilde{N}(s)} for \eqn{s < t}
#' }
#'
#' ## Simulation Returns
#'
#' The simulation returns a tibble with columns:
#' \itemize{
#'   \item `path_id`: Path identifier
#'   \item `time`: Time grid
#'   \item `count`: Standard Poisson process \eqn{N(t)}
#'   \item `compensated_count`: Compensated process \eqn{\tilde{N}(t) = N(t) - \lambda t}
#' }
#'
#' @examples
#' # Create Poisson process specification
#' spec <- poisson_spec(intensity = 1.0, initial_value = 0)
#'
#' # Simulate paths
#' paths <- simulate_paths(spec, n_paths = 25, n_steps = 500, maturity = 30)
#'
#' # Plot standard process
#' plot_paths(paths, plot_type = "standard")
#'
#' # Plot compensated (martingale) process
#' plot_paths(paths, plot_type = "compensated")
#'
#' # Pipeline workflow
#' poisson_spec(2.0) |>
#'   simulate_paths(50, 1000, 20) |>
#'   plot_paths(n_paths_plot = 10)
#'
#' @export
poisson_spec <- function(intensity, initial_value = 0) {
  # Validation
  checkmate::assert_number(intensity, lower = 0, finite = TRUE)
  checkmate::assert_number(initial_value, finite = TRUE)
  
  # Create specification object
  spec <- list(
    intensity = intensity,
    initial_value = initial_value
  )
  
  class(spec) <- c("poisson_spec", "process_spec")
  spec
}


#' Print Method for Poisson Process Specification
#'
#' @param x A `poisson_spec` object.
#' @param ... Additional arguments (ignored).
#'
#' @return Invisibly returns the input object.
#' @export
print.poisson_spec <- function(x, ...) {
  cli::cli_h2("Poisson Process Specification")
  cli::cli_text("Process: {.emph N(t) ~ Poisson(\u03bbt)} with jumps")
  cli::cli_text("")
  cli::cli_dl(c(
    "Initial value (N\u2080)" = cli::col_cyan("{x$initial_value}"),
    "Intensity (\u03bb)" = cli::col_green("{x$intensity}")
  ))
  cli::cli_text("")
  cli::cli_text("{.strong Process types:} standard, compensated")
  cli::cli_text("{.emph Compensated:} \u00d1(t) = N(t) - \u03bbt (martingale)")
  cli::cli_text("")
  cli::cli_alert_info("Use {.fn simulate_paths} to generate sample paths")
  invisible(x)
}


#' Simulate Poisson Process Paths
#'
#' S3 method for simulating standard and compensated Poisson process paths.
#'
#' @param process_spec A `poisson_spec` object.
#' @param n_paths Integer. Number of paths to simulate.
#' @param n_steps Integer. Number of time steps.
#' @param maturity Numeric. Time horizon.
#' @param seed Integer. Random seed for reproducibility. Default is 123.
#' @param ... Additional arguments (ignored).
#'
#' @return A tibble with columns:
#'   \item{path_id}{Path identifier}
#'   \item{time}{Time grid}
#'   \item{count}{Standard Poisson process N(t)}
#'   \item{compensated_count}{Compensated process Ñ(t) = N(t) - λt}
#'
#' @details
#' Simulates paths using the independent increments property:
#' \deqn{\Delta N_i \sim Poisson(\lambda \Delta t)}
#'
#' where \eqn{\Delta t} is the time step size.
#'
#' @export
simulate_paths.poisson_spec <- function(process_spec, n_paths, n_steps, maturity,
                                         seed = 123, ...) {
  # Validation
  checkmate::assert_integerish(n_paths, lower = 1, len = 1, any.missing = FALSE)
  checkmate::assert_integerish(n_steps, lower = 1, len = 1, any.missing = FALSE)
  checkmate::assert_number(maturity, lower = 0, finite = TRUE)
  checkmate::assert_integerish(seed, len = 1, any.missing = FALSE)
  
  # Coercion
  n_paths <- as.integer(n_paths)
  n_steps <- as.integer(n_steps)
  seed <- as.integer(seed)
  
  # Extract parameters
  intensity <- process_spec$intensity
  initial_value <- process_spec$initial_value
  
  # Set seed for reproducibility
  set.seed(seed)
  
  # Time step
  dt <- maturity / n_steps
  
  # Time grid
  time_grid <- seq(0, maturity, length.out = n_steps + 1)
  
  # Generate Poisson increments for each time step
  # Each increment is Poisson(λ * dt)
  increments <- purrr::map(
    seq_len(n_steps),
    \(step) stats::rpois(n_paths, intensity * dt)
  )
  
  # Compute cumulative counts (standard Poisson process)
  # N(t) = N(0) + sum of increments
  counts_standard <- purrr::accumulate(
    increments,
    \(prev_count, increment) prev_count + increment,
    .init = rep(initial_value, n_paths)
  )
  
  # Convert list to matrix: rows = paths, cols = time points
  counts_standard_matrix <- do.call(cbind, counts_standard)
  
  # Compensated Poisson process: Ñ(t) = N(t) - λt
  # Subtract λ * time from each column
  compensation <- purrr::map_dbl(seq_along(time_grid) - 1, \(step) intensity * dt * step)
  
  counts_compensated_matrix <- sweep(
    counts_standard_matrix,
    MARGIN = 2,
    STATS = compensation,
    FUN = `-`
  )
  
  # Convert to tidy format
  paths_tidy <- purrr::map_dfr(
    seq_len(n_paths),
    \(path_idx) {
      tibble::tibble(
        path_id = path_idx,
        time = time_grid,
        count = counts_standard_matrix[path_idx, ],
        compensated_count = counts_compensated_matrix[path_idx, ]
      )
    }
  )
  
  # Add process type attribute
  attr(paths_tidy, "process_type") <- "poisson"
  attr(paths_tidy, "intensity") <- intensity
  
  paths_tidy
}
