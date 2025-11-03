#' Recover Probability Density via FFT
#'
#' Implements a Fourier inversion routine based on the Fast Fourier Transform
#' (FFT) to recover a probability density from its characteristic function.
#' The implementation mirrors the lecture material approach but returns a tidy
#' tibble suitable for further analysis.
#'
#' @param cf Function returning characteristic function evaluations. Must accept
#'   a numeric vector `u` and return a complex vector of the same length.
#' @param x Numeric vector of evaluation points.
#' @param n Integer number of FFT grid points. Defaults to `8192L`.
#' @param u_max Numeric upper bound of the frequency domain grid. Defaults to 20.
#'
#' @return A tibble with columns `x` and `density` containing the reconstructed
#'   density values.
#' @export
fft_density_recovery <- function(cf,
                                 x,
                                 n = 8192L,
                                 u_max = 20) {
  checkmate::assert_function(cf)
  checkmate::assert_numeric(x, any.missing = FALSE, finite = TRUE)
  checkmate::assert_integerish(n, lower = 1, len = 1)
  checkmate::assert_number(u_max, lower = 0, finite = TRUE)

  n <- as.integer(n)
  u_max <- as.numeric(u_max)
  x <- as.numeric(x)

  du <- u_max / n
  u <- du * (seq_len(n) - 1)

  x_min <- min(x)
  dx <- 2 * pi / (n * du)
  x_grid <- x_min + dx * (seq_len(n) - 1)

  x_max <- max(x)
  grid_max <- x_grid[length(x_grid)]
  if (x_max > grid_max) {
    rlang::abort("`x` extends beyond the FFT grid. Increase `n` or `u_max`.")
  }

  cf_u <- cf(u)
  if (!is.complex(cf_u)) {
    cf_u <- as.complex(cf_u)
  }
  if (length(cf_u) != length(u)) {
    rlang::abort("Characteristic function must return values matching input length.")
  }

  phase_shift <- exp(-1i * x_min * u)
  phi <- phase_shift * cf_u

  gamma_1 <- exp(-1i * x_grid * u[1]) * cf(u[1])
  gamma_2 <- exp(-1i * x_grid * u[length(u)]) * cf(u[length(u)])
  phi_boundary <- 0.5 * (gamma_1 + gamma_2)

  fft_values <- stats::fft(phi)
  densities <- du / pi * Re(fft_values - phi_boundary)

  interpolator <- stats::splinefun(x_grid, densities, method = "fmm")

  tibble::tibble(
    x = x,
    density = interpolator(x)
  )
}
