test_that("engine_registry captures baseline specs", {
  registry <- engine_registry()

  expect_s3_class(registry, "tbl_df")
  expect_true("gbm_spec" %in% registry$spec_class)

  gbm_entry <- registry[registry$spec_class == "gbm_spec", , drop = FALSE]
  expect_equal(nrow(gbm_entry), 1)
  expect_equal(gbm_entry$simulate_engines, "euler")
  expect_match(gbm_entry$lineage, "gbm_spec")

  heston_entry <- registry[registry$spec_class == "heston_spec", , drop = FALSE]
  expect_equal(heston_entry$price_engines, "cos")

  short_rate_entry <- registry[registry$spec_class == "short_rate_spec", , drop = FALSE]
  expect_equal(short_rate_entry$process_type, "short_rate")
  expect_equal(short_rate_entry$simulate_engines, "monte_carlo")
  expect_equal(short_rate_entry$price_engines, "analytic_zcb")

  term_structure_entry <- registry[registry$spec_class == "term_structure_spec", , drop = FALSE]
  expect_equal(term_structure_entry$process_type, "term_structure")
  expect_equal(term_structure_entry$fit_engines, "direct, bootstrap, newton, multi_curve_newton, treasury_bootstrap")
  expect_equal(term_structure_entry$predict_engines, "linear_interp")
})
