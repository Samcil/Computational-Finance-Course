#' Engine Registry for Process Specifications
#'
#' Returns a tibble describing the available engines for key CompFinanceR
#' specifications along with their class lineage and process types. The registry
#' is derived from lightweight prototype specifications built with canonical
#' parameter sets to avoid heavy computations.
#'
#' @return Tidy tibble with columns:
#'   * `spec_class` – concrete specification class name
#'   * `process_type` – short descriptor used by simulation outputs
#'   * `lineage` – ordered class lineage illustrating inheritance
#'   * `simulate_engines` – comma separated list of registered simulation engines
#'   * `price_engines` – comma separated list of registered pricing engines
#'
#' @export
engine_registry <- function() {
  constructors <- list(
    gbm_spec = function() gbm_spec(initial_value = 100, drift = 0.05, volatility = 0.2),
    abm_spec = function() abm_spec(initial_value = 0, drift = 0.05, volatility = 0.2),
    poisson_spec = function() poisson_spec(intensity = 1, initial_value = 0),
    cir_spec = function() cir_spec(initial_value = 0.05, mean_reversion = 1.5, long_term_mean = 0.04, volatility = 0.2),
    correlated_bm_spec = function() {
      correlated_bm_spec(
        initial_values = c(100, 95),
        drift = c(0.05, 0.04),
        covariance_matrix = matrix(c(0.04, 0.02, 0.02, 0.09), nrow = 2)
      )
    },
    merton_spec = function() {
      merton_spec(
        initial_price = 100,
        risk_free_rate = 0.02,
        volatility = 0.2,
        jump_intensity = 0.8,
        jump_mean = -0.1,
        jump_sd = 0.3
      )
    },
    heston_spec = function() {
      heston_spec(
        initial_price = 100,
        initial_variance = 0.04,
        risk_free_rate = 0.02,
        dividend_yield = 0,
        mean_reversion = 1.5,
        long_term_variance = 0.04,
        vol_of_vol = 0.5,
        correlation = -0.7
      )
    },
    bshw_spec = function() {
      curve <- tibble::tibble(
        tenor = c(0, 1, 2, 5, 10),
        discount_factor = exp(-0.02 * tenor)
      )
      short_rate <- short_rate_spec(
        model = "hull_white",
        volatility = 0.01,
        mean_reversion = 0.1,
        curve = curve
      )
      bshw_spec(
        spot = 100,
        short_rate = short_rate,
        equity_vol = 0.2,
        correlation = 0.3
      )
    },
    short_rate_spec = function() {
      curve <- tibble::tibble(
        tenor = c(0, 1, 2, 3),
        discount_factor = exp(-0.02 * tenor)
      )
      short_rate_spec(
        model = "ho_lee",
        volatility = 0.01,
        curve = curve
      )
    },
    term_structure_spec = function() term_structure_spec(),
    mortgage_annuity_spec = function() {
      mortgage_annuity_spec(
        principal = 250000,
        coupon_rate = 0.03,
        maturity = 30,
        frequency = 12
      )
    },
    amortizing_swap_spec = function() {
      mortgage <- mortgage_annuity_spec(
        principal = 250000,
        coupon_rate = 0.03,
        maturity = 30,
        frequency = 12
      )
      amortizing_swap_spec(
        mortgage = mortgage,
        fixed_rate = mortgage$args$coupon_rate,
        type = "receiver"
      )
    }
  )

  collapse_or_na <- function(values) {
    if (is.null(values) || length(values) == 0) {
      NA_character_
    } else {
      paste(values, collapse = ", ")
    }
  }

  purrr::imap_dfr(constructors, function(builder, name) {
    spec <- builder()
    metadata <- spec$metadata
    if (is.null(metadata)) {
      metadata <- list()
    }
    engines <- metadata$engines
    if (is.null(engines)) {
      engines <- list()
    }

    spec_type <- if (!is.null(spec$process_type)) {
      spec$process_type
    } else if (!is.null(metadata$spec_type)) {
      metadata$spec_type
    } else if (!is.null(spec$mode)) {
      spec$mode
    } else {
      NA_character_
    }

    lineage <- if (inherits(spec, "process_spec")) {
      format_spec_lineage(spec)
    } else {
      paste(class(spec), collapse = " -> ")
    }

    tibble::tibble(
      spec_class = class(spec)[1],
      process_type = spec_type,
      lineage = lineage,
      simulate_engines = collapse_or_na(engines$simulate),
      price_engines = collapse_or_na(engines$price),
      fit_engines = collapse_or_na(engines$fit),
      predict_engines = collapse_or_na(engines$predict)
    )
  }) |>
    dplyr::distinct()
}
