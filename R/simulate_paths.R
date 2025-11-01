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
  # Use assert_integerish to allow numeric values coercible to integers (e.g., 100 instead of 100L)
  checkmate::assert_class(process_spec, "process_spec")
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


#' Create Geometric Brownian Motion Specification
#'
#' Specification for GBM process under risk-neutral measure.
#'
#' @param initial_value Numeric. Initial value \eqn{S_0}.
#' @param drift Numeric. Drift parameter \eqn{\mu} (risk-free rate in risk-neutral measure).
#' @param volatility Numeric. Volatility parameter \eqn{\sigma} (annualized).
#' 
#' @details
#' ## Mathematical Formulation
#' 
#' The Geometric Brownian Motion follows the stochastic differential equation:
#' \deqn{dS(t) = \mu S(t) dt + \sigma S(t) dW(t)}
#' 
#' Where:
#' \itemize{
#'   \item \eqn{S(t)} = Stock price at time \eqn{t}
#'   \item \eqn{\mu} = Drift parameter (risk-free rate under risk-neutral measure)
#'   \item \eqn{\sigma} = Volatility parameter (annualized standard deviation)
#'   \item \eqn{W(t)} = Standard Brownian motion
#' }
#' 
#' The analytical solution is:
#' \deqn{S(t) = S_0 \exp\left[(\mu - \frac{\sigma^2}{2})t + \sigma W(t)\right]}
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
  cli::cli_h2("Geometric Brownian Motion Specification")
  cli::cli_text("Process: {.emph dS(t) = \u03bc S(t) dt + \u03c3 S(t) dW(t)}")
  cli::cli_text("")
  cli::cli_dl(c(
    "Initial value (S\u2080)" = cli::col_cyan("{x$initial_value}"),
    "Drift (\u03bc)" = cli::col_green("{x$drift}"),
    "Volatility (\u03c3)" = cli::col_blue("{x$volatility}")
  ))
  cli::cli_text("")
  cli::cli_alert_info("Use {.fn simulate_paths} to generate sample paths")
  invisible(x)
}


#' Create Arithmetic Brownian Motion Specification
#'
#' Specification for ABM process (also known as Brownian motion with drift).
#'
#' @param initial_value Numeric. Initial value \eqn{X_0}.
#' @param drift Numeric. Drift parameter \eqn{\mu}.
#' @param volatility Numeric. Volatility parameter \eqn{\sigma} (annualized).
#' 
#' @details
#' ## Mathematical Formulation
#' 
#' The Arithmetic Brownian Motion follows the stochastic differential equation:
#' \deqn{dX(t) = \mu dt + \sigma dW(t)}
#' 
#' Where:
#' \itemize{
#'   \item \eqn{X(t)} = Process value at time \eqn{t}
#'   \item \eqn{\mu} = Drift parameter (deterministic trend)
#'   \item \eqn{\sigma} = Volatility parameter (annualized standard deviation)
#'   \item \eqn{W(t)} = Standard Brownian motion
#' }
#' 
#' The analytical solution is:
#' \deqn{X(t) = X_0 + \mu t + \sigma W(t)}
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
  cli::cli_h2("Arithmetic Brownian Motion Specification")
  cli::cli_text("Process: {.emph dX(t) = \u03bc dt + \u03c3 dW(t)}")
  cli::cli_text("")
  cli::cli_dl(c(
    "Initial value (X\u2080)" = cli::col_cyan("{x$initial_value}"),
    "Drift (\u03bc)" = cli::col_green("{x$drift}"),
    "Volatility (\u03c3)" = cli::col_blue("{x$volatility}")
  ))
  cli::cli_text("")
  cli::cli_alert_info("Use {.fn simulate_paths} to generate sample paths")
  invisible(x)
}
