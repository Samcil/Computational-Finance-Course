#' Simulate Stochastic Process Paths
#'
#' User-facing interface for simulating various stochastic processes following
#' tidymodels/hardhat design principles. This function provides a consistent
#' interface while hiding implementation complexity.
#'
#' @param process_spec A process specification object created by one of:
#'   \itemize{
#'     \item \code{gbm_spec()} for Geometric Brownian Motion
#'     \item \code{abm_spec()} for Arithmetic Brownian Motion
#'     \item \code{poisson_spec()} for Poisson processes (future)
#'     \item \code{cir_spec()} for CIR processes (future)
#'     \item \code{heston_spec()} for Heston model (future)
#'   }
#' @param n_paths Integer. Number of Monte Carlo paths to simulate.
#' @param n_steps Integer. Number of time steps for discretization.
#' @param maturity Numeric. Time to maturity in years.
#' @param seed Integer. Random seed for reproducibility. Default is 123.
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
                           seed = 123) {
  
  # Input validation using checkmate
  checkmate::assert_class(process_spec, "process_spec")
  checkmate::assert_int(n_paths, lower = 1)
  checkmate::assert_int(n_steps, lower = 1)
  checkmate::assert_number(maturity, lower = 0, finite = TRUE)
  checkmate::assert_int(seed)
  
  # Dispatch to appropriate engine based on process type
  UseMethod("simulate_paths", process_spec)
}


#' Create Geometric Brownian Motion Specification
#'
#' Specification for GBM process under risk-neutral measure:
#'   dS(t) = r*S(t)*dt + sigma*S(t)*dW(t)
#'
#' @param initial_value Numeric. Initial value S(0).
#' @param drift Numeric. Drift parameter (risk-free rate in risk-neutral measure).
#' @param volatility Numeric. Volatility parameter (annualized).
#'
#' @return A gbm_spec object (inherits from process_spec)
#'
#' @examples
#' spec <- gbm_spec(initial_value = 100, drift = 0.05, volatility = 0.2)
#'
#' @export
gbm_spec <- function(initial_value, drift, volatility) {
  
  # Validation
  checkmate::assert_number(initial_value, lower = 0, finite = TRUE)
  checkmate::assert_number(drift, finite = TRUE)
  checkmate::assert_number(volatility, lower = 0, finite = TRUE)
  
  # Create specification object
  spec <- list(
    initial_value = initial_value,
    drift = drift,
    volatility = volatility,
    process_type = "gbm"
  )
  
  # Set class hierarchy for dispatch
  class(spec) <- c("gbm_spec", "process_spec", "list")
  
  spec
}


#' @export
print.gbm_spec <- function(x, ...) {
  cat("Geometric Brownian Motion Specification\n")
  cat("────────────────────────────────────────\n")
  cat(sprintf("Initial value:  %g\n", x$initial_value))
  cat(sprintf("Drift:          %g\n", x$drift))
  cat(sprintf("Volatility:     %g\n", x$volatility))
  cat("\nUse simulate_paths() to generate sample paths.\n")
  invisible(x)
}


#' Create Arithmetic Brownian Motion Specification
#'
#' Specification for ABM process:
#'   dX(t) = mu*dt + sigma*dW(t)
#'
#' @param initial_value Numeric. Initial value X(0).
#' @param drift Numeric. Drift parameter.
#' @param volatility Numeric. Volatility parameter (annualized).
#'
#' @return An abm_spec object (inherits from process_spec)
#'
#' @examples
#' spec <- abm_spec(initial_value = 0, drift = 0.03, volatility = 0.2)
#'
#' @export
abm_spec <- function(initial_value, drift, volatility) {
  
  # Validation
  checkmate::assert_number(initial_value, finite = TRUE)
  checkmate::assert_number(drift, finite = TRUE)
  checkmate::assert_number(volatility, lower = 0, finite = TRUE)
  
  # Create specification object
  spec <- list(
    initial_value = initial_value,
    drift = drift,
    volatility = volatility,
    process_type = "abm"
  )
  
  # Set class hierarchy for dispatch
  class(spec) <- c("abm_spec", "process_spec", "list")
  
  spec
}


#' @export
print.abm_spec <- function(x, ...) {
  cat("Arithmetic Brownian Motion Specification\n")
  cat("─────────────────────────────────────────\n")
  cat(sprintf("Initial value:  %g\n", x$initial_value))
  cat(sprintf("Drift:          %g\n", x$drift))
  cat(sprintf("Volatility:     %g\n", x$volatility))
  cat("\nUse simulate_paths() to generate sample paths.\n")
  invisible(x)
}
