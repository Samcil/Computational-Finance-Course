#' Construct an Annuity Mortgage Schedule
#'
#' Build an amortisation schedule for a level-payment (annuity) mortgage using
#' the specified principal, coupon rate, maturity, and payment frequency. The
#' routine normalises the schedule into a tidy tibble that shares the core
#' columns required by swap pricing utilities while adding mortgage-specific
#' breakdowns of interest and principal components.
#'
#' @param principal Numeric scalar describing the initial mortgage balance.
#' @param coupon_rate Numeric scalar giving the nominal annual mortgage rate.
#' @param maturity Numeric scalar for the mortgage maturity in years.
#' @param frequency Integer number of payments per year (for example `12` for
#'   monthly).
#' @param daycount Optional numeric value for the accrual fraction per period.
#'   Defaults to `1 / frequency` when omitted.
#' @param payment Optional numeric override for the level payment each period.
#'   When `NULL` the payment is solved to exactly amortise the mortgage.
#' @param prepayment Optional numeric scalar or vector describing the constant
#'   or time-varying conditional prepayment rate (CPR). Accepts a scalar (used
#'   for every period), a vector of length equal to the number of periods, or a
#'   vector of length `periods + 1` mirroring the Python reference
#'   implementation where the first element corresponds to time zero.
#'
#' @return Tibble with columns `period`, `start`, `end`, `pay_time`,
#'   `accrual_fraction`, `notional`, `payment`, `interest_component`,
#'   `principal_component`, `outstanding_start`, `outstanding_end`, and
#'   `period_rate`.
#'
#' @examples
#' schedule <- mortgage_annuity_schedule(
#'   principal = 250000,
#'   coupon_rate = 0.03,
#'   maturity = 30,
#'   frequency = 12
#' )
#'
#' @export
mortgage_annuity_schedule <- function(principal,
                                      coupon_rate,
                                      maturity,
                                      frequency = 12,
                                      daycount = NULL,
                                      payment = NULL,
                                      prepayment = 0) {
  checkmate::assert_number(principal, lower = 0, finite = TRUE)
  checkmate::assert_number(coupon_rate, finite = TRUE)
  checkmate::assert_number(maturity, lower = 0, finite = TRUE)
  checkmate::assert_int(frequency, lower = 1)
  if (!is.null(daycount)) {
    checkmate::assert_number(daycount, lower = 0, finite = TRUE)
  }
  if (!is.null(payment)) {
    checkmate::assert_number(payment, lower = 0, finite = TRUE)
  }
  checkmate::assert_numeric(prepayment, any.missing = FALSE, finite = TRUE, lower = 0)

  step <- 1 / frequency
  periods_raw <- maturity / step
  if (!checkmate::test_number(periods_raw, lower = 1, finite = TRUE) ||
    abs(periods_raw - round(periods_raw)) > 1e-8) {
    rlang::abort("`maturity` must be an integer multiple of `1 / frequency`")
  }
  n_periods <- as.integer(round(periods_raw))
  if (n_periods < 1) {
    rlang::abort("Mortgage schedule requires at least one payment period")
  }

  accrual_fraction <- if (is.null(daycount)) {
    rep(step, n_periods)
  } else {
    rep(daycount, n_periods)
  }

  base_schedule <- swap_cashflow_schedule(
    start = 0,
    end = maturity,
    frequency = frequency,
    notional = rep(1, n_periods),
    daycount = accrual_fraction
  )

  period_rate <- coupon_rate / frequency

  if (is.null(payment)) {
    payment_value <- if (abs(period_rate) < 1e-12) {
      principal / n_periods
    } else {
      principal * period_rate / (1 - (1 + period_rate)^(-n_periods))
    }
  } else {
    payment_value <- payment
  }

  if (length(prepayment) == 1L) {
    prepayment_vec <- rep(prepayment, n_periods)
  } else if (length(prepayment) == n_periods) {
    prepayment_vec <- as.numeric(prepayment)
  } else if (length(prepayment) == n_periods + 1L) {
    prepayment_vec <- as.numeric(prepayment[-1])
  } else {
    rlang::abort("`prepayment` must be scalar or have length equal to the number of periods (or periods + 1)")
  }

  compute_payment <- function(current_balance, remaining_periods) {
    if (is.null(payment)) {
      if (abs(period_rate) < 1e-12) {
        current_balance / remaining_periods
      } else {
        current_balance * period_rate / (1 - (1 + period_rate)^(-remaining_periods))
      }
    } else {
      payment_value
    }
  }

  evolve_period <- function(state, idx) {
    balance <- state$balance
    remaining_periods <- n_periods - (idx - 1)
    payment_current <- compute_payment(balance, remaining_periods)
    interest_current <- balance * period_rate
    principal_current <- payment_current - interest_current
    if (principal_current < -1e-8) {
      rlang::abort("Level payment does not amortise the mortgage under the supplied coupon rate")
    }

    scheduled_outstanding <- balance - principal_current
    prepayment_current <- prepayment_vec[[idx]] * scheduled_outstanding

    new_balance <- scheduled_outstanding - prepayment_current
    if (new_balance < 0 && new_balance > -1e-6) {
      new_balance <- 0
    }

    record <- tibble::tibble(
      outstanding_start = balance,
      payment = payment_current,
      interest_component = interest_current,
      principal_component = principal_current,
      prepayment_component = prepayment_current,
      outstanding_end = new_balance
    )

    list(
      balance = new_balance,
      records = append(state$records, list(record))
    )
  }

  initial_state <- list(balance = principal, records = list())
  states <- purrr::accumulate(seq_len(n_periods), evolve_period, .init = initial_state)
  final_state <- states[[length(states)]]
  schedule_details <- purrr::list_rbind(final_state$records)

  total_principal <- sum(schedule_details$principal_component + schedule_details$prepayment_component)
  diff_principal <- principal - total_principal
  if (abs(diff_principal) > max(1e-6, principal * 1e-8)) {
    schedule_details$principal_component[[n_periods]] <- schedule_details$principal_component[[n_periods]] + diff_principal
    schedule_details$payment[[n_periods]] <- schedule_details$interest_component[[n_periods]] + schedule_details$principal_component[[n_periods]]
    schedule_details$outstanding_end[[n_periods]] <- max(
      schedule_details$outstanding_start[[n_periods]] - schedule_details$principal_component[[n_periods]] - schedule_details$prepayment_component[[n_periods]],
      0
    )
  }

  base_schedule |>
    dplyr::mutate(
      notional = schedule_details$outstanding_start,
      payment = schedule_details$payment,
      interest_component = schedule_details$interest_component,
      principal_component = schedule_details$principal_component,
      prepayment_component = schedule_details$prepayment_component,
      outstanding_start = schedule_details$outstanding_start,
      outstanding_end = schedule_details$outstanding_end,
      period_rate = rep(period_rate, n_periods),
      cpr = prepayment_vec
    )
}


#' Construct a Bullet Mortgage Schedule
#'
#' Generate an amortisation schedule for a bullet mortgage where interest is
#' paid periodically and the outstanding principal is settled at maturity. The
#' routine mirrors the lecture `BulletMortgage.py` script while reusing the
#' swap schedule utilities so mortgage and swap analytics operate on consistent
#' columns.
#'
#' @inheritParams mortgage_annuity_schedule
#'
#' @return Tibble containing the bullet mortgage schedule with the same column
#'   set as [mortgage_annuity_schedule()].
#'
#' @export
mortgage_bullet_schedule <- function(principal,
                                     coupon_rate,
                                     maturity,
                                     frequency = 12,
                                     daycount = NULL,
                                     prepayment = 0) {
  checkmate::assert_number(principal, lower = 0, finite = TRUE)
  checkmate::assert_number(coupon_rate, finite = TRUE)
  checkmate::assert_number(maturity, lower = 0, finite = TRUE)
  checkmate::assert_int(frequency, lower = 1)
  if (!is.null(daycount)) {
    checkmate::assert_number(daycount, lower = 0, finite = TRUE)
  }
  checkmate::assert_numeric(prepayment, any.missing = FALSE, finite = TRUE, lower = 0)

  step <- 1 / frequency
  periods_raw <- maturity / step
  if (!checkmate::test_number(periods_raw, lower = 1, finite = TRUE) ||
    abs(periods_raw - round(periods_raw)) > 1e-8) {
    rlang::abort("`maturity` must be an integer multiple of `1 / frequency`")
  }
  n_periods <- as.integer(round(periods_raw))

  accrual_fraction <- if (is.null(daycount)) {
    rep(step, n_periods)
  } else {
    rep(daycount, n_periods)
  }

  base_schedule <- swap_cashflow_schedule(
    start = 0,
    end = maturity,
    frequency = frequency,
    notional = rep(1, n_periods),
    daycount = accrual_fraction
  )

  period_rate <- coupon_rate / frequency

  if (length(prepayment) == 1L) {
    prepayment_vec <- rep(prepayment, n_periods)
  } else if (length(prepayment) == n_periods) {
    prepayment_vec <- as.numeric(prepayment)
  } else if (length(prepayment) == n_periods + 1L) {
    prepayment_vec <- as.numeric(prepayment[-1])
  } else {
    rlang::abort("`prepayment` must be scalar or match the number of periods (or periods + 1)")
  }

  evolve_period <- function(state, idx) {
    balance <- state$balance
    interest_current <- balance * period_rate
    principal_current <- if (idx == n_periods) balance else 0
    scheduled_outstanding <- balance - principal_current
    prepayment_current <- prepayment_vec[[idx]] * scheduled_outstanding

    new_balance <- scheduled_outstanding - prepayment_current
    if (idx == n_periods && abs(new_balance) > 1e-8) {
      new_balance <- 0
    }

    record <- tibble::tibble(
      outstanding_start = balance,
      payment = interest_current + principal_current + prepayment_current,
      interest_component = interest_current,
      principal_component = principal_current,
      prepayment_component = prepayment_current,
      outstanding_end = new_balance
    )

    list(
      balance = new_balance,
      records = append(state$records, list(record))
    )
  }

  initial_state <- list(balance = principal, records = list())
  states <- purrr::accumulate(seq_len(n_periods), evolve_period, .init = initial_state)
  schedule_details <- purrr::list_rbind(states[[length(states)]]$records)

  base_schedule |>
    dplyr::mutate(
      notional = schedule_details$outstanding_start,
      payment = schedule_details$payment,
      interest_component = schedule_details$interest_component,
      principal_component = schedule_details$principal_component,
      prepayment_component = schedule_details$prepayment_component,
      outstanding_start = schedule_details$outstanding_start,
      outstanding_end = schedule_details$outstanding_end,
      period_rate = rep(period_rate, n_periods),
      cpr = prepayment_vec
    )
}


#' Prepayment Incentive Function
#'
#' Compute the conditional prepayment rate (CPR) incentive using a smooth
#' logistic function of the mortgage--swap spread, mirroring the
#' `PrepaymentFunction.py` lecture script. The curve shifts smoothly between the
#' baseline CPR and a higher incentive level as the spread crosses the chosen
#' threshold.
#'
#' @param rate_spread Numeric vector of mortgage--swap spreads expressed in
#'   decimal form.
#' @param base Baseline CPR (decimal) when the spread is far below the
#'   threshold.
#' @param amplitude Maximum uplift above the baseline CPR as the spread widens.
#' @param slope Steepness of the logistic transition.
#' @param threshold Spread level where the incentive is halfway between `base`
#'   and `base + amplitude`.
#'
#' @return Numeric vector with the incentive CPR values.
#'
#' @export
mortgage_prepayment_incentive <- function(rate_spread,
                                          base = 0.04,
                                          amplitude = 0.1,
                                          slope = 115,
                                          threshold = 0.02) {
  checkmate::assert_numeric(rate_spread, any.missing = FALSE, finite = TRUE)
  checkmate::assert_number(base, lower = 0)
  checkmate::assert_number(amplitude, lower = 0)
  checkmate::assert_number(slope, finite = TRUE)
  checkmate::assert_number(threshold, finite = TRUE)

  base + amplitude / (1 + exp(slope * (threshold - rate_spread)))
}


#' Mortgage Annuity Specification
#'
#' Create a specification wrapping the annuity mortgage schedule so that pricing
#' and swap utilities can consume a consistent contract object. The schedule is
#' stored on the specification for reuse across valuation and simulation
#' routines.
#'
#' @inheritParams mortgage_annuity_schedule
#' @param engine Character scalar identifying the engine. Currently only
#'   `'deterministic'` is supported.
#' @param engine_options Named list of engine-specific options. Reserved for
#'   future extensions.
#'
#' @return Object of class `mortgage_annuity_spec` inheriting from `model_spec`.
#'
#' @export
mortgage_annuity_spec <- function(principal,
                                  coupon_rate,
                                  maturity,
                                  frequency = 12,
                                  daycount = NULL,
                                  payment = NULL,
                                  prepayment = 0,
                                  engine = "deterministic",
                                  engine_options = list()) {
  checkmate::assert_list(engine_options, names = "unique", null.ok = FALSE)

  schedule <- mortgage_annuity_schedule(
    principal = principal,
    coupon_rate = coupon_rate,
    maturity = maturity,
    frequency = frequency,
    daycount = daycount,
    payment = payment,
    prepayment = prepayment
  )

  spec <- new_model_spec(
    class = "mortgage_annuity_spec",
    args = list(
      principal = principal,
      coupon_rate = coupon_rate,
      maturity = maturity,
      frequency = frequency,
      payment = schedule$payment[[1]],
      prepayment = prepayment
    ),
    mode = "valuation"
  )

  spec$schedule <- schedule
  spec <- set_spec_metadata(
    spec,
    spec_type = "mortgage",
    engines = list(price = "deterministic"),
    schedule_columns = names(schedule)
  )

  engine_options <- rlang::list2(!!!engine_options)
  spec <- set_engine(spec, engine = engine, !!!engine_options)
  spec
}


#' @export
print.mortgage_annuity_spec <- function(x, ...) {
  cli::cli_h2("Mortgage Annuity Specification")
  engine_label <- if (is.null(x$method$engine)) "<unset>" else x$method$engine
  cli::cli_dl(c(
    "Principal" = cli::col_cyan(format(x$args$principal, digits = 6)),
    "Coupon" = cli::col_cyan(format(x$args$coupon_rate, digits = 4)),
    "Maturity" = cli::col_cyan(format(x$args$maturity, digits = 4)),
    "Frequency" = cli::col_cyan(x$args$frequency),
    "Engine" = cli::col_cyan(engine_label)
  ))
  cli::cli_alert_info("Use mortgage_schedule() to inspect the amortisation table")
  invisible(x)
}


#' @export
set_engine.mortgage_annuity_spec <- function(object, engine, ...) {
  checkmate::assert_choice(engine, choices = c("deterministic"))
  eng_args <- rlang::list2(...)
  set_engine_base(object, engine = engine, eng_args = eng_args)
}


#' Extract Mortgage Schedule
#'
#' Return the amortisation schedule stored on a mortgage specification.
#'
#' @param object A `mortgage_annuity_spec` produced by
#'   [mortgage_annuity_spec()].
#'
#' @return Tibble containing the mortgage schedule.
#'
#' @export
mortgage_schedule <- function(object) {
  checkmate::assert_true(inherits(object, "mortgage_annuity_spec"))
  object$schedule
}


#' Price an Annuity Mortgage
#'
#' Compute the present value of the mortgage from the perspective of the lender
#' or borrower by discounting the scheduled payments with a supplied short-rate
#' specification.
#'
#' @param mortgage A `mortgage_annuity_spec` created by
#'   [mortgage_annuity_spec()].
#' @param discount_spec A `short_rate_spec` providing discount factors via
#'   [price_zcb()].
#' @param type Character scalar. `"lender"` reports the present value from the
#'   lender's perspective (positive when payments cover funding at the supplied
#'   curve). `"borrower"` flips the sign.
#'
#' @return Tibble with metrics `pv`, `pv_principal`, and `pv_interest`.
#'
#' @export
price_mortgage <- function(mortgage,
                           discount_spec,
                           type = c("lender", "borrower")) {
  type <- rlang::arg_match(type)
  checkmate::assert_true(inherits(mortgage, "mortgage_annuity_spec"))
  checkmate::assert_class(discount_spec, "short_rate_spec")

  schedule <- mortgage_schedule(mortgage)
  state <- short_rate_state(discount_spec)
  discount_fun <- state$discount_fun_vec
  discount_factors <- discount_fun(schedule$pay_time)

  pv_principal <- sum((schedule$principal_component + schedule$prepayment_component) * discount_factors)
  pv_interest <- sum(schedule$interest_component * discount_factors)
  initial_principal <- mortgage$args$principal
  pv_total <- pv_principal + pv_interest - initial_principal

  if (type == "borrower") {
    pv_total <- -pv_total
    pv_principal <- -pv_principal
    pv_interest <- -pv_interest
  }

  tibble::tibble(
    metric = c("pv", "pv_principal", "pv_interest"),
    value = c(pv_total, pv_principal, pv_interest)
  )
}


#' Construct an Amortizing Swap Schedule
#'
#' Transform a mortgage amortisation schedule into a swap cash-flow schedule by
#' mapping the outstanding balance to the notional profile and retaining the
#' principal repayments as additional metadata.
#'
#' @param schedule Either a `mortgage_annuity_spec` or a tibble produced by
#'   [mortgage_annuity_schedule()].
#'
#' @return Tibble compatible with [swap_cashflow_schedule()] containing the
#'   amortising notional path.
#'
#' @export
amortizing_swap_schedule <- function(schedule) {
  if (inherits(schedule, "mortgage_annuity_spec")) {
    schedule <- mortgage_schedule(schedule)
  }

  checkmate::assert_data_frame(schedule, any.missing = FALSE, min.rows = 1)
  required_cols <- c(
    "period", "start", "end", "pay_time", "accrual_fraction",
    "outstanding_start", "principal_component", "prepayment_component"
  )
  missing_cols <- setdiff(required_cols, names(schedule))
  if (length(missing_cols) > 0) {
    rlang::abort(
      message = "Mortgage schedule is missing required columns",
      class = "amortizing_schedule_missing_columns",
      columns = missing_cols
    )
  }

  tibble::tibble(
    period = schedule$period,
    start = schedule$start,
    end = schedule$end,
    pay_time = schedule$pay_time,
    accrual_fraction = schedule$accrual_fraction,
    notional = schedule$outstanding_start,
    principal_payment = schedule$principal_component,
    prepayment_payment = schedule$prepayment_component,
    payment = schedule$payment,
    outstanding_end = schedule$outstanding_end
  )
}


#' Amortizing Swap Specification
#'
#' Create a specification that represents a fixed-for-floating swap whose
#' notional follows the amortisation pattern of a supplied mortgage schedule.
#'
#' @param mortgage Either a `mortgage_annuity_spec` or a tibble produced by
#'   [mortgage_annuity_schedule()].
#' @param fixed_rate Numeric scalar giving the fixed rate applied to the
#'   amortising notional. When omitted and `mortgage` is a specification the
#'   mortgage coupon rate is used.
#' @param type Character scalar. Either `"payer"` (pay fixed, receive float) or
#'   `"receiver"`.
#' @param engine Character scalar identifying the engine. Defaults to
#'   `'deterministic'`.
#' @param engine_options Named list of engine-specific parameters.
#'
#' @return Object of class `amortizing_swap_spec` inheriting from `model_spec`.
#'
#' @export
amortizing_swap_spec <- function(mortgage,
                                 fixed_rate = NULL,
                                 type = c("payer", "receiver"),
                                 engine = "deterministic",
                                 engine_options = list()) {
  schedule <- amortizing_swap_schedule(mortgage)
  checkmate::assert_list(engine_options, names = "unique", null.ok = FALSE)

  if (is.null(fixed_rate)) {
    if (inherits(mortgage, "mortgage_annuity_spec")) {
      fixed_rate <- mortgage$args$coupon_rate
    } else {
      rlang::abort("`fixed_rate` must be supplied when `mortgage` is not a specification")
    }
  }

  checkmate::assert_number(fixed_rate, finite = TRUE)
  swap_type <- rlang::arg_match(type)

  spec <- new_model_spec(
    class = "amortizing_swap_spec",
    args = list(
      fixed_rate = fixed_rate,
      type = swap_type
    ),
    mode = "valuation"
  )

  spec$schedule <- schedule
  spec <- set_spec_metadata(
    spec,
    spec_type = "amortizing_swap",
    engines = list(price = "deterministic", exposure = "monte_carlo"),
    schedule_columns = names(schedule)
  )

  engine_options <- rlang::list2(!!!engine_options)
  spec <- set_engine(spec, engine = engine, !!!engine_options)
  spec
}


#' @export
print.amortizing_swap_spec <- function(x, ...) {
  cli::cli_h2("Amortizing Swap Specification")
  engine_label <- if (is.null(x$method$engine)) "<unset>" else x$method$engine
  cli::cli_dl(c(
    "Fixed rate" = cli::col_cyan(format(x$args$fixed_rate, digits = 4)),
    "Type" = cli::col_cyan(x$args$type),
    "Engine" = cli::col_cyan(engine_label)
  ))
  cli::cli_alert_info("Use price_amortizing_swap() or simulate_swap_exposure() with the stored schedule")
  invisible(x)
}


#' @export
set_engine.amortizing_swap_spec <- function(object, engine, ...) {
  checkmate::assert_choice(engine, choices = c("deterministic"))
  eng_args <- rlang::list2(...)
  set_engine_base(object, engine = engine, eng_args = eng_args)
}


#' Price an Amortizing Swap
#'
#' Delegate pricing of an amortizing swap specification to [price_swap()],
#' allowing optional overrides of the fixed rate or swap direction.
#'
#' @param amortizing_swap An `amortizing_swap_spec` created by
#'   [amortizing_swap_spec()].
#' @param discount_spec A `short_rate_spec` supplying discount factors.
#' @param fixed_rate Optional numeric override for the fixed rate.
#' @param type Optional character override for the payer/receiver direction.
#'
#' @return Tibble identical to [price_swap()].
#'
#' @export
price_amortizing_swap <- function(amortizing_swap,
                                  discount_spec,
                                  fixed_rate = NULL,
                                  type = NULL) {
  checkmate::assert_true(inherits(amortizing_swap, "amortizing_swap_spec"))
  checkmate::assert_class(discount_spec, "short_rate_spec")

  schedule <- amortizing_swap$schedule
  rate <- if (is.null(fixed_rate)) {
    amortizing_swap$args$fixed_rate
  } else {
    checkmate::assert_number(fixed_rate, finite = TRUE)
    fixed_rate
  }

  swap_type <- if (is.null(type)) {
    amortizing_swap$args$type
  } else {
    rlang::arg_match(type, c("payer", "receiver"))
  }

  price_swap(
    spec = discount_spec,
    schedule = schedule,
    fixed_rate = rate,
    type = swap_type
  )
}
