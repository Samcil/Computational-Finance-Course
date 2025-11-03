bond_option_price_short_rate <- function(spec,
                                         option,
                                         strike,
                                         expiry,
                                         maturity) {
  option <- rlang::arg_match(option, c("call", "put"))
  checkmate::assert_class(spec, "short_rate_spec")
  checkmate::assert_number(strike, lower = 0, finite = TRUE)
  checkmate::assert_number(expiry, lower = 0, finite = TRUE)
  checkmate::assert_number(maturity, lower = expiry, finite = TRUE)

  state <- short_rate_state(spec)
  discount_fun <- state$discount_fun_vec
  discount0_expiry <- discount_fun(expiry)
  discount0_maturity <- discount_fun(maturity)

  forward_bond <- discount0_maturity / discount0_expiry

  args <- spec$args
  sigma <- args$volatility
  a <- args$mean_reversion

  if (sigma < 1e-10) {
    intrinsic <- if (option == "call") {
      max(forward_bond - strike, 0)
    } else {
      max(strike - forward_bond, 0)
    }
    return(discount0_expiry * intrinsic)
  }

  b_val <- b_factor(a, maturity - expiry)
  chi_term <- if (a > 0) {
    (1 - exp(-2 * a * expiry)) / (2 * a)
  } else {
    expiry
  }

  sigma_p_sq <- sigma^2 * b_val^2 * chi_term
  sigma_p <- sqrt(sigma_p_sq)

  if (sigma_p < 1e-10) {
    intrinsic <- if (option == "call") {
      max(forward_bond - strike, 0)
    } else {
      max(strike - forward_bond, 0)
    }
    return(discount0_expiry * intrinsic)
  }

  log_ratio <- log(discount0_maturity / (strike * discount0_expiry))
  d1 <- (log_ratio + 0.5 * sigma_p_sq) / sigma_p
  d2 <- d1 - sigma_p

  call_value <- discount0_maturity * stats::pnorm(d1) - strike * discount0_expiry * stats::pnorm(d2)

  if (option == "call") {
    call_value
  } else {
    call_value - discount0_maturity + strike * discount0_expiry
  }
}


find_short_rate_root <- function(fn,
                                 lower = -1,
                                 upper = 1,
                                 tol = 1e-8) {
  brackets <- list(
    c(lower, upper),
    c(lower * 5, upper * 5),
    c(-10, 10),
    c(-25, 25)
  )

  bracket <- purrr::detect(
    brackets,
    \(br) {
      f_low <- fn(br[1])
      f_high <- fn(br[2])
      is.finite(f_low) && is.finite(f_high) && f_low * f_high <= 0
    }
  )

  if (is.null(bracket)) {
    rlang::abort("Unable to bracket Jamshidian root for the supplied swaption")
  }

  stats::uniroot(fn, interval = bracket, tol = tol)$root
}


#' Price European Swaptions Under Short-Rate Models
#'
#' Evaluate a payer or receiver swaption by applying Jamshidian's trick to the
#' fixed-leg cash flows of the referenced swap.
#'
#' @param spec A [`short_rate_spec`] object.
#' @param schedule Tibble produced by [swap_cashflow_schedule()] describing the
#'   underlying swap.
#' @param fixed_rate Numeric scalar. Fixed coupon rate in decimal form.
#' @param type Character. Either `"payer"` or `"receiver"`.
#'
#' @return A tibble with columns `metric = "price"` and `value`.
#' @export
price_swaption <- function(spec,
                           schedule,
                           fixed_rate,
                           type = c("payer", "receiver")) {
  type <- rlang::arg_match(type)
  checkmate::assert_class(spec, "short_rate_spec")
  schedule <- validate_swap_schedule(schedule)
  checkmate::assert_number(fixed_rate, finite = TRUE)
  exercise <- min(schedule$start)

  swap_value <- price_swap(spec, schedule, fixed_rate = fixed_rate, type = type) |>
    dplyr::filter(.data$metric == "pv") |>
    dplyr::pull(.data$value)

  args <- spec$args

  if (exercise <= 1e-10 || args$volatility < 1e-10) {
    option_value <- max(swap_value, 0)
    return(tibble::tibble(metric = "price", value = option_value))
  }

  future_rows <- schedule$pay_time > exercise + 1e-10
  if (!any(future_rows)) {
    option_value <- max(swap_value, 0)
    return(tibble::tibble(metric = "price", value = option_value))
  }

  pay_times <- schedule$pay_time[future_rows]
  accruals <- schedule$accrual_fraction[future_rows]
  notionals <- schedule$notional[future_rows]

  sign <- if (type == "payer") 1 else -1
  coeff0 <- sign * notionals[1]
  coeffs <- -sign * fixed_rate * accruals * notionals
  coeffs[length(coeffs)] <- coeffs[length(coeffs)] - sign * notionals[length(notionals)]

  root_fn <- function(r) {
    bond_prices <- price_zcb(
      spec,
      maturities = pay_times,
      valuation_time = exercise,
      short_rate = rep(r, length(pay_times))
    )
    coeff0 + sum(coeffs * bond_prices)
  }

  r_star <- find_short_rate_root(root_fn)

  strike_bonds <- price_zcb(
    spec,
    maturities = pay_times,
    valuation_time = exercise,
    short_rate = rep(r_star, length(pay_times))
  )

  contributions <- purrr::pmap_dbl(
    list(weight = coeffs, strike = strike_bonds, maturity = pay_times),
    \(weight, strike, maturity) {
      if (abs(weight) < 1e-12) {
        return(0)
      }
      option_type <- if (weight > 0) "call" else "put"
      option_price <- bond_option_price_short_rate(
        spec = spec,
        option = option_type,
        strike = strike,
        expiry = exercise,
        maturity = maturity
      )
      abs(weight) * option_price
    }
  )

  option_value <- max(sum(contributions), 0)
  tibble::tibble(metric = "price", value = option_value)
}
