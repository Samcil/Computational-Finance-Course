#' Plot Simulated Paths
#'
#' Creates professional visualizations of simulated stochastic process paths
#' using ggplot2. Automatically detects the process type and applies appropriate
#' styling.
#'
#' @param paths_data Tibble. Output from \code{simulate_paths()}.
#' @param n_paths_plot Integer. Maximum number of paths to plot. Default is NULL (all paths).
#' @param alpha Numeric. Transparency level for path lines (0 to 1). Default is 0.3.
#' @param theme Character. ggplot2 theme to use. One of "minimal", "classic", "bw". Default is "minimal".
#'
#' @return A ggplot object.
#'
#' @details
#' The function automatically detects the process type from the paths_data attributes
#' and applies appropriate labels and styling. For GBM processes, it includes the
#' SDE equation in the subtitle.
#'
#' @examples
#' # Simulate and plot GBM
#' gbm_spec(100, 0.05, 0.2) |>
#'   simulate_paths(n_paths = 50, n_steps = 252, maturity = 1) |>
#'   plot_paths(n_paths_plot = 10)
#'
#' @export
plot_paths <- function(paths_data,
                       n_paths_plot = NULL,
                       alpha = 0.3,
                       theme = c("minimal", "classic", "bw")) {
  
  # Validation
  checkmate::assert_data_frame(paths_data)
  checkmate::assert_int(n_paths_plot, lower = 1, null.ok = TRUE)
  checkmate::assert_number(alpha, lower = 0, upper = 1)
  theme <- match.arg(theme)
  
  # Get process type from attributes
  process_type <- attr(paths_data, "process_type")
  spec <- attr(paths_data, "spec")
  
  # Filter paths if requested
  if (!is.null(n_paths_plot)) {
    max_path_id <- max(paths_data$path_id)
    paths_to_plot <- min(n_paths_plot, max_path_id)
    paths_data <- dplyr::filter(paths_data, path_id <= paths_to_plot)
  }
  
  # Determine value column name
  value_col <- if ("stock_price" %in% names(paths_data)) {
    "stock_price"
  } else {
    "value"
  }
  
  # Create appropriate labels based on process type
  plot_config <- get_plot_config(process_type, spec)
  
  # Create base plot
  p <- ggplot2::ggplot(
    paths_data,
    ggplot2::aes(x = time, y = .data[[value_col]], group = path_id)
  ) +
    ggplot2::geom_line(
      alpha = alpha,
      linewidth = 0.5,
      color = plot_config$color
    ) +
    ggplot2::labs(
      title = plot_config$title,
      subtitle = plot_config$subtitle,
      x = "Time (years)",
      y = plot_config$y_label
    )
  
  # Apply theme
  p <- p + switch(
    theme,
    minimal = ggplot2::theme_minimal(),
    classic = ggplot2::theme_classic(),
    bw = ggplot2::theme_bw()
  )
  
  # Additional styling
  p <- p +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 14),
      plot.subtitle = ggplot2::element_text(size = 10, color = "gray30"),
      axis.title = ggplot2::element_text(face = "bold"),
      panel.grid.minor = ggplot2::element_blank()
    )
  
  p
}


#' Get Plot Configuration for Process Type
#'
#' Returns appropriate plot labels and colors based on process type.
#'
#' @param process_type Character. Type of stochastic process.
#' @param spec Process specification object.
#'
#' @return List with title, subtitle, y_label, and color.
#'
#' @keywords internal
get_plot_config <- function(process_type, spec) {
  
  if (is.null(process_type) || process_type == "gbm") {
    list(
      title = "Geometric Brownian Motion Paths",
      subtitle = sprintf(
        "dS(t) = %.4f·S(t)·dt + %.4f·S(t)·dW(t) | S(0) = %.2f",
        spec$drift, spec$volatility, spec$initial_value
      ),
      y_label = "Stock Price",
      color = "steelblue"
    )
  } else if (process_type == "abm") {
    list(
      title = "Arithmetic Brownian Motion Paths",
      subtitle = sprintf(
        "dX(t) = %.4f·dt + %.4f·dW(t) | X(0) = %.2f",
        spec$drift, spec$volatility, spec$initial_value
      ),
      y_label = "Process Value",
      color = "darkred"
    )
  } else {
    # Default
    list(
      title = "Simulated Paths",
      subtitle = "",
      y_label = "Value",
      color = "black"
    )
  }
}


#' Quick Demonstration of Path Simulation
#'
#' Convenience function to quickly simulate and visualize paths with sensible defaults.
#'
#' @param process_spec Process specification object.
#' @param n_paths Integer. Number of paths. Default is 25.
#' @param ... Additional arguments passed to \code{simulate_paths()}.
#'
#' @return List with paths (tibble) and plot (ggplot object).
#'
#' @examples
#' # Quick GBM demo
#' demo <- demo_paths(gbm_spec(100, 0.05, 0.3), n_paths = 50)
#' print(demo$plot)
#'
#' @export
demo_paths <- function(process_spec, n_paths = 25, ...) {
  
  # Simulate with defaults
  paths <- simulate_paths(
    process_spec = process_spec,
    n_paths = n_paths,
    n_steps = 252,
    maturity = 1.0,
    ...
  )
  
  # Create plot
  p <- plot_paths(paths, n_paths_plot = min(25, n_paths))
  
  # Return both
  list(
    paths = paths,
    plot = p
  )
}
