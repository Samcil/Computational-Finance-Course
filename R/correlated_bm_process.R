#' Create Correlated Brownian Motion Specification
#'
#' Defines a multivariate Brownian motion with drift and covariance structure
#' suitable for simulating correlated asset paths. The specification follows the
#' tidymodels/hardhat design pattern and integrates with `simulate_paths()`.
#'
#' @param initial_values Numeric vector. Initial values for each component.
#' @param drift Numeric vector. Drift parameters \eqn{\boldsymbol{\mu}} for each
#'   component (annualised).
#' @param covariance_matrix Numeric matrix. Positive semi-definite covariance
#'   matrix \eqn{\Sigma} describing instantaneous covariances.
#' @param component_names Character vector. Optional names for the components.
#'   Defaults to `paste0("component_", seq_along(initial_values))` when `NULL`.
#'
#' @details
#' ## Mathematical Formulation
#'
#' A correlated Brownian motion \eqn{\mathbf{X}(t)} with drift satisfies
#' \deqn{d\mathbf{X}(t) = \boldsymbol{\mu}\, dt + L\, d\mathbf{W}(t)}
#' where \eqn{L} is a Cholesky factor of \eqn{\Sigma} (i.e. \eqn{LL^\top = \Sigma})
#' and \eqn{\mathbf{W}(t)} is a vector of independent Brownian motions.
#'
#' @return An object of class `correlated_bm_spec` inheriting from
#'   `process_spec`.
#'
#' @examples
#' spec <- correlated_bm_spec(
#'   initial_values = c(100, 95),
#'   drift = c(0.05, 0.04),
#'   covariance_matrix = matrix(c(0.04, 0.03, 0.03, 0.09), nrow = 2),
#'   component_names = c("Asset_A", "Asset_B")
#' )
#'
#' spec |>
#'   simulate_paths(n_paths = 100, n_steps = 252, maturity = 1) |>
#'   plot_paths()
#'
#' @export
correlated_bm_spec <- function(initial_values,
                               drift,
                               covariance_matrix,
                               component_names = NULL) {
  checkmate::assert_numeric(initial_values, finite = TRUE, any.missing = FALSE)
  checkmate::assert_numeric(drift, len = length(initial_values), finite = TRUE, any.missing = FALSE)
  checkmate::assert_matrix(
    covariance_matrix,
    mode = "numeric",
    any.missing = FALSE,
    nrows = length(initial_values),
    ncols = length(initial_values)
  )
  checkmate::assert_true(isTRUE(all.equal(covariance_matrix, t(covariance_matrix))))
  eigenvalues <- eigen(covariance_matrix, symmetric = TRUE, only.values = TRUE)$values
  checkmate::assert_true(all(eigenvalues >= -sqrt(.Machine$double.eps)))

  if (is.null(component_names)) {
    component_names <- paste0("component_", seq_along(initial_values))
  }
  checkmate::assert_character(component_names, len = length(initial_values), any.missing = FALSE)

  structure(
    list(
      initial_values = initial_values,
      drift = drift,
      covariance_matrix = covariance_matrix,
      component_names = component_names
    ),
    class = c("correlated_bm_spec", "process_spec")
  )
}

#' @export
print.correlated_bm_spec <- function(x, ...) {
  cli::cli_h2("Correlated Brownian Motion Specification")
  cli::cli_text("Process: dX(t) = \u03bc dt + L dW(t)")
  cli::cli_text("")
  cli::cli_dl(c(
    "Dimensions" = cli::col_cyan(length(x$initial_values)),
    "Component names" = cli::col_green(paste(x$component_names, collapse = ", "))
  ))
  cli::cli_text("")
  cli::cli_alert_info("Use {.fn simulate_paths} to generate correlated paths")
  invisible(x)
}

#' Simulate Correlated Brownian Motion Paths
#'
#' Generates multivariate Brownian motion sample paths using an
#' Euler-Maruyama discretisation with shared randomness to enforce the desired
#' covariance structure.
#'
#' @param process_spec A `correlated_bm_spec` object.
#' @inheritParams simulate_paths
#'
#' @return Tibble with columns `path_id`, `time`, `component`, and `value`.
#'
#' @export
simulate_paths.correlated_bm_spec <- function(process_spec,
                                              n_paths,
                                              n_steps,
                                              maturity,
                                              seed = 123,
                                              ...) {
  checkmate::assert_integerish(n_paths, lower = 1, len = 1, any.missing = FALSE)
  checkmate::assert_integerish(n_steps, lower = 1, len = 1, any.missing = FALSE)
  checkmate::assert_number(maturity, lower = 0, finite = TRUE)
  checkmate::assert_integerish(seed, len = 1, any.missing = FALSE)

  n_paths <- as.integer(n_paths)
  n_steps <- as.integer(n_steps)
  seed <- as.integer(seed)

  set.seed(seed)

  dt <- maturity / n_steps
  time_grid <- seq(0, maturity, length.out = n_steps + 1)
  dimension <- length(process_spec$initial_values)
  component_names <- process_spec$component_names
  initial_state <- stats::setNames(process_spec$initial_values, component_names)
  drift_step <- process_spec$drift * dt
  chol_factor <- chol(process_spec$covariance_matrix)

  standard_normals <- matrix(
    stats::rnorm(n_paths * n_steps * dimension),
    ncol = dimension
  )
  diffusion_flat <- (standard_normals %*% t(chol_factor)) * sqrt(dt)
  increments_flat <- diffusion_flat +
    matrix(
      rep(drift_step, each = n_paths * n_steps),
      ncol = dimension,
      byrow = FALSE
    )
  increments_array <- array(
    data = increments_flat,
    dim = c(n_paths, n_steps, dimension)
  )
  path_increments <- purrr::array_tree(increments_array, margin = 1)

  path_tibbles <- purrr::imap(
    path_increments,
    \(increments_matrix, path_idx) {
      state_history <- purrr::accumulate(
        .x = seq_len(n_steps),
        .init = initial_state,
        .f = \(state, step_idx) {
          next_state <- state + increments_matrix[step_idx, ]
          stats::setNames(next_state, component_names)
        }
      )
      state_matrix <- do.call(rbind, state_history)
      tibble::as_tibble(state_matrix, .name_repair = "minimal") |>
        dplyr::mutate(
          path_id = path_idx,
          time = time_grid,
          .before = 1
        ) |>
        tidyr::pivot_longer(
          cols = dplyr::all_of(component_names),
          names_to = "component",
          values_to = "value"
        )
    }
  )

  paths_tidy <- purrr::list_rbind(path_tibbles)
  attr(paths_tidy, "process_type") <- "correlated_bm"
  attr(paths_tidy, "component_names") <- component_names
  attr(paths_tidy, "spec") <- process_spec

  paths_tidy
}
