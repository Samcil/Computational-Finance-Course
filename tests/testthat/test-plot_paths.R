test_that("plot_paths creates ggplot object", {
  spec <- gbm_spec(100, 0.05, 0.2)
  paths <- simulate_paths(spec, 10, 50, 1.0)
  
  p <- plot_paths(paths)
  
  expect_s3_class(p, "ggplot")
})


test_that("plot_paths works with filtered paths", {
  spec <- gbm_spec(100, 0.05, 0.2)
  paths <- simulate_paths(spec, 20, 50, 1.0)
  
  p <- plot_paths(paths, n_paths_plot = 5)
  
  expect_s3_class(p, "ggplot")
})


test_that("plot_paths works with different themes", {
  spec <- abm_spec(0, 0.03, 0.15)
  paths <- simulate_paths(spec, 5, 50, 1.0)
  
  p1 <- plot_paths(paths, theme = "minimal")
  p2 <- plot_paths(paths, theme = "classic")
  p3 <- plot_paths(paths, theme = "bw")
  
  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
  expect_s3_class(p3, "ggplot")
})


test_that("plot_paths validates inputs", {
  spec <- gbm_spec(100, 0.05, 0.2)
  paths <- simulate_paths(spec, 10, 50, 1.0)
  
  expect_error(
    plot_paths(paths, n_paths_plot = 0),
    "n_paths_plot"
  )
  
  expect_error(
    plot_paths(paths, alpha = 1.5),
    "alpha"
  )
  
  expect_error(
    plot_paths("not a data frame"),
    "paths_data"
  )
})


test_that("demo_paths works", {
  spec <- gbm_spec(100, 0.05, 0.2)
  result <- demo_paths(spec, n_paths = 10)
  
  expect_type(result, "list")
  expect_true(all(c("paths", "plot") %in% names(result)))
  expect_s3_class(result$paths, "tbl_df")
  expect_s3_class(result$plot, "ggplot")
})


test_that("pipeline with plotting works", {
  # Full pipeline test
  result <- gbm_spec(100, 0.05, 0.3) |>
    simulate_paths(n_paths = 15, n_steps = 100, maturity = 1.0) |>
    plot_paths(n_paths_plot = 10, alpha = 0.4)
  
  expect_s3_class(result, "ggplot")
})
