test_that("q_measure_paths produces consistent P and Q simulations", {
  spec <- gbm_spec(initial_value = 100, drift = 0.08, volatility = 0.2)
  risk_free <- 0.03
  dividend <- 0.01
  maturity <- 1
  paths <- q_measure_paths(
    spec = spec,
    risk_free_rate = risk_free,
    n_paths = 5000,
    n_steps = 252,
    maturity = maturity,
    dividend_yield = dividend,
    seed = 314
  )

  expect_s3_class(paths, c("measure_paths", "tbl_df"))
  expect_true(all(c("measure", "path_id", "time", "stock_price", "radon_nikodym") %in% names(paths)))

  lambda_expected <- (spec$drift - (risk_free - dividend)) / spec$volatility
  expect_equal(attr(paths, "market_price_of_risk"), lambda_expected, tolerance = 1e-8)
  expect_equal(attr(paths, "risk_free_rate"), risk_free)
  expect_equal(attr(paths, "dividend_yield"), dividend)

  p_terminal <- dplyr::filter(paths, measure == "P", time == maturity)
  q_terminal <- dplyr::filter(paths, measure == "Q", time == maturity)

  expect_true(all(p_terminal$radon_nikodym > 0))
  radon_start <- dplyr::filter(paths, measure == "P", time == 0)$radon_nikodym
  expect_true(all(abs(radon_start - 1) < 1e-8))

  theoretical <- spec$initial_value * exp((risk_free - dividend) * maturity)
  expect_equal(mean(q_terminal$stock_price), theoretical, tolerance = 0.5)
  expect_equal(mean(p_terminal$stock_price * p_terminal$radon_nikodym), theoretical, tolerance = 0.5)

  q_weights <- dplyr::filter(paths, measure == "Q", time == maturity)$radon_nikodym
  expect_true(all(abs(q_weights - 1) < 1e-12))
})

test_that("q_measure_paths respects explicit market price of risk", {
  spec <- gbm_spec(initial_value = 50, drift = 0.07, volatility = 0.25)
  lambda <- 0.4
  res <- q_measure_paths(
    spec = spec,
    risk_free_rate = 0.02,
    n_paths = 1000,
    n_steps = 64,
    maturity = 0.5,
    market_price_of_risk = lambda,
    seed = 11
  )

  expect_equal(attr(res, "market_price_of_risk"), lambda)
})

test_that("short_rate_measure_paths returns consistent discount expectations", {
  curve <- tibble::tibble(tenor = c(0, 1, 2, 3), discount_factor = exp(-0.02 * tenor))
  spec <- short_rate_spec(
    model = "ho_lee",
    volatility = 0.01,
    curve = curve
  )

  maturity <- 2
  paths <- short_rate_measure_paths(
    spec = spec,
    n_paths = 4000,
    n_steps = 200,
    maturity = maturity,
    market_price_of_risk = 0.4,
    seed = 2718
  )

  expect_s3_class(paths, c("short_rate_measure_paths", "measure_paths", "tbl_df"))
  expect_true(all(c("measure", "path_id", "time", "short_rate", "discount_factor", "radon_nikodym") %in% names(paths)))

  radon_start <- dplyr::filter(paths, measure == "P", time == 0)$radon_nikodym
  expect_true(all(abs(radon_start - 1) < 1e-8))
  expect_true(all(dplyr::filter(paths, measure == "P")$radon_nikodym > 0))

  q_terminal <- dplyr::filter(paths, measure == "Q", time == maturity)
  p_terminal <- dplyr::filter(paths, measure == "P", time == maturity)

  target_discount <- short_rate_state(spec)$discount_fun(maturity)
  expect_equal(mean(q_terminal$discount_factor), target_discount, tolerance = 1.5e-2)
  expect_equal(mean(p_terminal$discount_factor * p_terminal$radon_nikodym), target_discount, tolerance = 1.5e-2)

  lambda_grid <- attr(paths, "market_price_of_risk")
  expect_s3_class(lambda_grid, "tbl_df")
  expect_equal(lambda_grid$lambda[1], 0.4, tolerance = 1e-8)
})

test_that("short_rate_measure_paths accepts functional risk price", {
  curve <- tibble::tibble(tenor = c(0, 0.5, 1.5), discount_factor = exp(-0.025 * tenor))
  spec <- short_rate_spec(model = "ho_lee", volatility = 0.012, curve = curve)
  lambda_fun <- function(t) 0.1 + 0.05 * t

  res <- short_rate_measure_paths(
    spec = spec,
    n_paths = 1000,
    n_steps = 60,
    maturity = 1.5,
    market_price_of_risk = lambda_fun,
    seed = 99
  )

  grid <- attr(res, "market_price_of_risk")
  expect_true(nrow(grid) > 0)
  expect_true(any(abs(grid$lambda - (0.1 + 0.05 * grid$time)) < 1e-8))
})
