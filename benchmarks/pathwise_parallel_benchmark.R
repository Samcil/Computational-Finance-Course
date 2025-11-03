if (!requireNamespace("mirai", quietly = TRUE)) {
  stop("mirai package is required for benchmarking.")
}
if (!requireNamespace("pkgload", quietly = TRUE)) {
  stop("pkgload package is required for benchmarking.")
}

pkgload::load_all(helpers = FALSE, quiet = TRUE)
library(mirai)
library(tibble)
library(dplyr)
library(purrr)

cores <- parallel::detectCores()
if (is.null(cores) || is.na(cores) || cores < 2L) {
  cores <- 2L
}

candidate_workers <- unique(c(0L, max(1L, cores - 1L)))
benchmark_runs <- 3L

spec <- gbm_spec(
  initial_value = 100,
  drift = 0.05,
  volatility = 0.2
)

n_paths <- 5000
n_steps <- 126

benchmark_once <- function(workers, run_id) {
  mirai::daemons(0)
  if (workers > 0) {
    mirai::daemons(workers)
  }
  timing <- system.time({
    compare_with_finite_diff(
      process_spec = spec,
      strikes = seq(80, 120, by = 10),
      maturity = 1,
      risk_free_rate = 0.03,
      volatility = 0.2,
      option_type = "call",
      bump_spot = 1,
      bump_vol = 0.01,
      n_paths = n_paths,
      n_steps = n_steps,
      seed = 123
    )
  })[["elapsed"]]
  mirai::daemons(0)
  tibble(
    workers = workers,
    run = run_id,
    elapsed = timing
  )
}

run_results <- purrr::map(candidate_workers, function(workers) {
  purrr::map(seq_len(benchmark_runs), function(run_id) {
    benchmark_once(workers, run_id)
  }) |> purrr::list_rbind()
}) |> purrr::list_rbind()

summary_results <- run_results |>
  dplyr::group_by(workers) |>
  dplyr::summarise(
    mean_elapsed = mean(elapsed),
    median_elapsed = median(elapsed),
    sd_elapsed = sd(elapsed),
    .groups = "drop"
  )

print(run_results)
print(summary_results)
