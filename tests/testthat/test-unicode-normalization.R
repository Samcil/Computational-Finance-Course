test_that("normalize_to_ascii converts Greek letters and superscripts", {
  sample <- c(
    "X ~ N(\u03bc, \u03c3\u00b2)",
    "Characteristic function \u03c6(u)"
  )

  normalized <- normalize_to_ascii(sample)

  expect_equal(normalized[1], "X ~ N(mu, sigma^2)")
  expect_equal(normalized[2], "Characteristic function phi(u)")
  expect_false(any(grepl("[^\\x00-\\x7F]", normalized, perl = TRUE)))
})

test_that("normalize_to_ascii handles subscripts and macrons", {
  sample <- "\u03ba\u0304(T\u2081) = 4\u03ba v\u0304 / \u03b3\u00b2"

  normalized <- normalize_to_ascii(sample)

  expect_equal(normalized, "kappa_bar(T_1) = 4kappa v_bar / gamma^2")
  expect_false(grepl("[^\\x00-\\x7F]", normalized, perl = TRUE))
})
