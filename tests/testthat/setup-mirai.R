if (requireNamespace("mirai", quietly = TRUE)) {
  cores <- parallel::detectCores()
  if (is.null(cores) || is.na(cores)) {
    cores <- 1L
  }
  workers <- max(1L, as.integer(cores) - 1L)
  if (workers > 0) {
    mirai::daemons(workers)
  }
}
