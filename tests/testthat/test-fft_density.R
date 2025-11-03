test_that("fft_density_recovery matches normal pdf", {
  mu <- 0
  sigma <- 1
  cf <- function(u) exp(1i * mu * u - 0.5 * sigma^2 * u^2)
  x <- seq(-4, 4, length.out = 25)

  result <- fft_density_recovery(cf, x, n = 2^12, u_max = 25)
  expected <- tibble::tibble(x = x, density = dnorm(x, mean = mu, sd = sigma))

  expect_equal(result$density, expected$density, tolerance = 5e-3)
})
