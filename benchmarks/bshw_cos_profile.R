#!/usr/bin/env Rscript

"%||%" <- function(lhs, rhs) {
  if (is.null(lhs) || is.na(lhs) || lhs == "") {
    rhs
  } else {
    lhs
  }
}

resolve_repo_root <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_prefix <- "--file="
  script_path <- args[grepl(file_prefix, args, fixed = TRUE)]
  if (length(script_path) == 0) {
    return(getwd())
  }
  script_dir <- dirname(normalizePath(sub(file_prefix, "", script_path)))
  normalizePath(file.path(script_dir, ".."))
}

with_dir <- function(path, code) {
  old <- getwd()
  on.exit(setwd(old), add = TRUE)
  setwd(path)
  force(code)
}

default_iter <- as.integer(Sys.getenv("BSHW_PROFILE_ITERATIONS") %||% "25")
default_terms <- as.integer(Sys.getenv("BSHW_PROFILE_N_TERMS") %||% "512")
default_trunc <- as.numeric(Sys.getenv("BSHW_PROFILE_TRUNCATION") %||% "10")
default_steps <- as.integer(Sys.getenv("BSHW_PROFILE_INTEGRATION_STEPS") %||% "2000")

target_repo <- resolve_repo_root()
with_dir(target_repo, {
  suppressPackageStartupMessages({
    if (!requireNamespace("devtools", quietly = TRUE)) {
      stop("Package 'devtools' is required to load the local package context.")
    }
    devtools::load_all(quiet = TRUE)
  })

  curve <- tibble::tibble(
    tenor = seq(0, 20, by = 1),
    discount_factor = exp(-0.05 * tenor)
  )

  short_rate <- short_rate_spec(
    model = "hull_white",
    volatility = 0.05,
    mean_reversion = 0.1,
    curve = curve
  )

  spec <- bshw_spec(
    spot = 100,
    short_rate = short_rate,
    equity_vol = 0.2,
    correlation = 0.3
  )

  strikes <- seq(60, 140, length.out = 5)
  profile_path <- file.path("benchmarks", "bshw_cos_profile.out")

  message("Profiling price_bshw_option_cos() with ", default_iter, " iterations...")
  Rprof(profile_path, interval = 0.001)
  purrr::walk(
    seq_len(default_iter),
    ~ price_bshw_option_cos(
      spec,
      strikes = strikes,
      maturity = 5,
      n_terms = default_terms,
      truncation = default_trunc,
      integration_points = default_steps
    )
  )
  Rprof(NULL)

  profile_summary <- summaryRprof(profile_path)
  by_total <- profile_summary$by.total
  if (is.null(by_total) || nrow(by_total) == 0) {
    stop("Profiler did not capture any samples.")
  }

  cat("\nTop 10 functions by total time (seconds):\n")
  print(head(by_total, 10))

  cat("\nSample interval:", profile_summary$sample.interval, "seconds\n")
  cat("Total sampled time:", profile_summary$sampling.time, "seconds\n")
  cat("Profile written to:", normalizePath(profile_path), "\n")
})
