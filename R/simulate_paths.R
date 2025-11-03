#' Simulate Stochastic Process Paths
#'
#' User-facing interface for simulating various stochastic processes following
#' tidymodels/hardhat design principles. This function provides a consistent
#' interface while hiding implementation complexity.
#'
#' @param process_spec A model specification object (inherits from
#'   `model_spec`) created by one of:
#'   \itemize{
#'     \item \code{gbm_spec()} for Geometric Brownian Motion
#'     \item \code{abm_spec()} for Arithmetic Brownian Motion
#'     \item \code{poisson_spec()} for Poisson processes
#'     \item \code{cir_spec()} for CIR processes
#'     \item \code{correlated_bm_spec()} for correlated Brownian motion
#'     \item \code{heston_spec()} for Heston stochastic volatility model
#'     \item \code{merton_spec()} for Merton jump-diffusion
#'     \item \code{short_rate_spec()} for short-rate models (e.g., Ho-Lee)
#'   }
#' @param n_paths Integer. Number of Monte Carlo paths to simulate.
#' @param n_steps Integer. Number of time steps for discretization.
#' @param maturity Numeric. Time to maturity in years.
#' @param seed Integer. Random seed for reproducibility. Default is 123.
#' @param ... Process-specific options forwarded to the underlying simulation
#'   engine (e.g. schemes or variance reduction controls).
#'
#' @return A tibble in tidy format with columns:
#'   \describe{
#'     \item{path_id}{Path identifier (1 to n_paths)}
#'     \item{time}{Time point}
#'     \item{value}{Simulated process value}
#'   }
#'   Additional columns depend on the specific process.
#'
#' @details
#' This function implements the "spec-fit-predict" pattern from tidymodels:
#' 1. Create specification with process parameters (e.g., \code{gbm_spec()})
#' 2. Simulate paths using \code{simulate_paths()}
#' 3. Extract or visualize results
#'
#' Following hardhat principles, this provides a simple, consistent interface
#' for users while isolating complexity in engine implementations.
#'
#' @examples
#' # Geometric Brownian Motion
#' gbm <- gbm_spec(
#'   initial_value = 100,
#'   drift = 0.05,
#'   volatility = 0.2
#' )
#'
#' paths <- simulate_paths(
#'   process_spec = gbm,
#'   n_paths = 100,
#'   n_steps = 252,
#'   maturity = 1.0
#' )
#'
#' # Pipeline-friendly usage
#' gbm_spec(initial_value = 100, drift = 0.05, volatility = 0.2) |>
#'   simulate_paths(n_paths = 50, n_steps = 252, maturity = 1.0) |>
#'   plot_paths()
#'
#' @export
simulate_paths <- function(process_spec,
                           n_paths,
                           n_steps,
                           maturity,
                           seed = 123,
                           ...) {
  # Input validation using checkmate
  # Use assert_integerish to allow numeric values coercible to integers (e.g., 100 instead of 100L)
  if (!inherits(process_spec, "model_spec")) {
    rlang::abort("`process_spec` must be created by a supported specification constructor.")
  }
  checkmate::assert_integerish(n_paths, lower = 1, len = 1, any.missing = FALSE)
  checkmate::assert_integerish(n_steps, lower = 1, len = 1, any.missing = FALSE)
  checkmate::assert_number(maturity, lower = 0, finite = TRUE)
  checkmate::assert_integerish(seed, len = 1, any.missing = FALSE)

  # Coerce to integer for internal use
  n_paths <- as.integer(n_paths)
  n_steps <- as.integer(n_steps)
  seed <- as.integer(seed)

  # Dispatch to appropriate engine based on process type
  UseMethod("simulate_paths", process_spec)
}
