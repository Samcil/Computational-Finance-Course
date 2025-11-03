#' Brownian Motion Specifications
#'
#' Provides specification constructors for basic Brownian motion processes
#' following tidymodels-style "spec" objects. These are separated from the
#' simulation engines to keep user-facing interfaces cohesive.
#'
#' @name bm_specs
NULL


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
#' @return A gbm_spec object (inherits from model_spec)
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

  spec <- new_diffusion_spec(
    class = "gbm_spec",
    args = list(
      initial_value = initial_value,
      drift = drift,
      volatility = volatility
    ),
    process_type = "gbm",
    inheritance = "abm_spec"
  )

  spec <- set_process_metadata(
    spec,
    parent_spec = "abm_spec",
    state_transform = "exp",
    notes = "Log-state Euler scheme shares implementation with ABM",
    engines = list(simulate = "euler")
  )

  spec
}


#' @export
print.gbm_spec <- function(x, ...) {
  cli::cli_h2("Geometric Brownian Motion Specification")
  cli::cli_text("Process: {.emph dS(t) = \u03bc S(t) dt + \u03c3 S(t) dW(t)}")
  cli::cli_text("")
  cli::cli_dl(c(
    "Initial value (S\u2080)" = cli::col_cyan("{format(x$initial_value, digits = 6)}"),
    "Drift (\u03bc)" = cli::col_green("{format(x$drift, digits = 6)}"),
    "Volatility (\u03c3)" = cli::col_blue("{format(x$volatility, digits = 6)}")
  ))
  cli::cli_text("")
  cli::cli_text("{.strong Lineage:} {format_spec_lineage(x)}")
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
#' @return An abm_spec object (inherits from model_spec)
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

  spec <- new_diffusion_spec(
    class = "abm_spec",
    args = list(
      initial_value = initial_value,
      drift = drift,
      volatility = volatility
    ),
    process_type = "abm"
  )

  spec <- set_process_metadata(
    spec,
    state_transform = "identity",
    notes = "Baseline diffusion used as parent for GBM",
    engines = list(simulate = "euler")
  )

  spec
}


#' @export
print.abm_spec <- function(x, ...) {
  cli::cli_h2("Arithmetic Brownian Motion Specification")
  cli::cli_text("Process: {.emph dX(t) = \u03bc dt + \u03c3 dW(t)}")
  cli::cli_text("")
  cli::cli_dl(c(
    "Initial value (X\u2080)" = cli::col_cyan("{format(x$initial_value, digits = 6)}"),
    "Drift (\u03bc)" = cli::col_green("{format(x$drift, digits = 6)}"),
    "Volatility (\u03c3)" = cli::col_blue("{format(x$volatility, digits = 6)}")
  ))
  cli::cli_text("")
  cli::cli_text("{.strong Lineage:} {format_spec_lineage(x)}")
  cli::cli_text("")
  cli::cli_alert_info("Use {.fn simulate_paths} to generate sample paths")
  invisible(x)
}
