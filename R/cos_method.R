#' Fourier Cosine (COS) Method for European Options
#'
#' Implements the Fourier Cosine (COS) expansion method of Fang and Oosterlee
#' (2008) for pricing European call and put options under models with known
#' characteristic functions. The implementation follows the specification-first
#' design of the package and returns tidy tibbles for convenient downstream use.
#'
#' @param cf Function returning the characteristic function evaluated at complex
#'   arguments. The function must be vectorised over the input `u` and return a
#'   complex vector.
#' @param option_type Character string, either "call" or "put".
#' @param spot Numeric scalar. Current underlying price.
#' @param risk_free_rate Numeric scalar. Continuously compounded risk-free rate.
#' @param maturity Numeric scalar. Time to maturity in years.
#' @param strikes Numeric vector of strike prices.
#' @param n_terms Integer. Number of COS expansion terms. Default is 256.
#' @param truncation Numeric scalar controlling the integration domain size.
#'   Default is 8. Ignored when `lower_bound` and `upper_bound` are supplied.
#' @param lower_bound Optional numeric scalar specifying the lower limit of the
#'   truncated integration domain.
#' @param upper_bound Optional numeric scalar specifying the upper limit of the
#'   truncated integration domain.
#'
#' @return A tibble with columns `strike` and `price` containing COS prices.
#'
#' @references Fang, F., & Oosterlee, C. W. (2008). A novel pricing method for
#'   European options based on Fourier-cosine series expansions. SIAM Journal on
#'   Scientific Computing, 31(2), 826-848.
#'
#' @export
cos_call_put_price <- function(cf,
                               option_type = c("call", "put"),
                               spot,
                               risk_free_rate,
                               maturity,
                               strikes,
                               n_terms = 256L,
                               truncation = 8,
                               lower_bound = NULL,
                               upper_bound = NULL) {
  option_type <- rlang::arg_match(option_type)
  checkmate::assert_function(cf)
  checkmate::assert_number(spot, lower = .Machine$double.eps, finite = TRUE)
  checkmate::assert_number(risk_free_rate, finite = TRUE)
  checkmate::assert_number(maturity, lower = 0, finite = TRUE)
  checkmate::assert_numeric(strikes, lower = .Machine$double.eps, any.missing = FALSE, finite = TRUE)
  checkmate::assert_integerish(n_terms, lower = 1, len = 1)
  if (is.null(lower_bound) && is.null(upper_bound)) {
    checkmate::assert_number(truncation, lower = 0, finite = TRUE)
  } else {
    checkmate::assert_number(lower_bound, finite = TRUE)
    checkmate::assert_number(upper_bound, finite = TRUE)
    if (lower_bound >= upper_bound) {
      rlang::abort("lower_bound must be strictly less than upper_bound")
    }
  }

  n_terms <- as.integer(n_terms)
  strikes <- as.numeric(strikes)

  if (!is.null(lower_bound) && !is.null(upper_bound)) {
    a <- lower_bound
    b <- upper_bound
  } else {
    a <- -truncation * sqrt(maturity)
    b <- truncation * sqrt(maturity)
  }
  k <- seq_len(n_terms) - 1
  u <- k * pi / (b - a)

  coefficients <- cos_coefficients(option_type, a, b, k)
  characteristic_values <- cf(u)
  characteristic_values[1] <- characteristic_values[1] * 0.5
  weights <- coefficients * characteristic_values

  x0 <- log(spot / strikes)
  exponential_matrix <- exp(1i * outer(x0 - a, u))
  prices <- exp(-risk_free_rate * maturity) * strikes * Re(exponential_matrix %*% weights)

  tibble::tibble(
    strike = strikes,
    price = as.numeric(prices)
  )
}

#' Recover Probability Densities via COS Expansion
#'
#' Applies the COS method to reconstruct the probability density function of a
#' random variable given its characteristic function. The implementation follows
#' Fang and Oosterlee (2008) and provides a tidy tibble for the evaluated grid.
#'
#' @param cf Function returning characteristic function values.
#' @param x Numeric vector of evaluation points.
#' @param maturity Numeric scalar representing the time horizon of the random
#'   variable (used to scale the truncation interval).
#' @param n_terms Integer. Number of cosine expansion terms. Default is 256.
#' @param truncation Numeric scalar controlling the truncation domain. Default is
#'   8.
#'
#' @return A tibble with columns `x` and `density`.
#' @export
cos_density_recovery <- function(cf,
                                 x,
                                 maturity,
                                 n_terms = 256L,
                                 truncation = 8) {
  checkmate::assert_function(cf)
  checkmate::assert_numeric(x, any.missing = FALSE, finite = TRUE)
  checkmate::assert_number(maturity, lower = 0, finite = TRUE)
  checkmate::assert_integerish(n_terms, lower = 1, len = 1)
  checkmate::assert_number(truncation, lower = 0, finite = TRUE)

  n_terms <- as.integer(n_terms)
  x <- as.numeric(x)

  a <- -truncation * sqrt(maturity)
  b <- truncation * sqrt(maturity)
  k <- seq_len(n_terms) - 1
  u <- k * pi / (b - a)

  coefficients <- 2 / (b - a) * Re(cf(u) * exp(-1i * u * a))
  coefficients[1] <- coefficients[1] * 0.5

  cosine_matrix <- cos(outer(x - a, u))
  densities <- cosine_matrix %*% coefficients

  tibble::tibble(
    x = x,
    density = as.numeric(densities)
  )
}

#' COS Method Coefficients for Payoff Functions
#'
#' Internal helper computing payoff-dependent cosine coefficients for call and
#' put options.
#'
#' @param option_type Character string, one of "call" or "put".
#' @param a Numeric scalar. Lower truncation bound.
#' @param b Numeric scalar. Upper truncation bound.
#' @param k Integer vector of cosine indices.
#'
#' @return Numeric vector of coefficients.
#' @keywords internal
cos_coefficients <- function(option_type, a, b, k) {
  option_type <- rlang::arg_match(option_type, c("call", "put"))
  width <- b - a
  checkmate::assert_number(width, lower = .Machine$double.eps, finite = TRUE)
  k <- as.numeric(k)

  if (option_type == "call") {
    c <- 0
    d <- b
    if (b <= 0) {
      return(rep(0, length(k)))
    }
    coef <- chi_psi_functions(a, b, c, d, k)
    2 / width * (coef$chi - coef$psi)
  } else {
    c <- a
    d <- 0
    coef <- chi_psi_functions(a, b, c, d, k)
    2 / width * (-coef$chi + coef$psi)
  }
}

#' Compute Chi and Psi Terms for COS Expansion
#'
#' Internal helper matching the notation of Fang and Oosterlee (2008).
#'
#' @param a Numeric scalar. Lower truncation bound.
#' @param b Numeric scalar. Upper truncation bound.
#' @param c Numeric scalar defining payoff interval lower limit.
#' @param d Numeric scalar defining payoff interval upper limit.
#' @param k Numeric vector of cosine indices.
#'
#' @return List with elements `chi` and `psi`.
#' @keywords internal
chi_psi_functions <- function(a, b, c, d, k) {
  width <- b - a
  k_pi_over_width <- k * pi / width

  psi <- sin(k_pi_over_width * (d - a)) - sin(k_pi_over_width * (c - a))
  non_zero <- k != 0
  psi[non_zero] <- psi[non_zero] * width / (k[non_zero] * pi)
  psi[!non_zero] <- d - c

  chi <- cos(k_pi_over_width * (d - a)) * exp(d) - cos(k_pi_over_width * (c - a)) * exp(c)
  chi <- chi + k_pi_over_width * (sin(k_pi_over_width * (d - a)) * exp(d) - sin(k_pi_over_width * (c - a)) * exp(c))
  chi <- chi / (1 + k_pi_over_width^2)

  list(chi = chi, psi = psi)
}
