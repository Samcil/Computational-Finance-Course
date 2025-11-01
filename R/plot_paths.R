#' Plot Simulated Paths
#'
#' Creates professional visualizations of simulated stochastic process paths
#' using ggplot2. Automatically detects the process type and applies appropriate
#' styling.
#'
#' @param paths_data Tibble. Output from \code{simulate_paths()}.
#' @param n_paths_plot Integer. Maximum number of paths to plot. Default is NULL (all paths).
#' @param plot_type Character. For Poisson: "standard" or "compensated". Ignored for other processes.
#' @param alpha Numeric. Transparency level for path lines (0 to 1). Default is 0.3.
#' @param theme Character. ggplot2 theme to use. One of "minimal", "classic", "bw". Default is "minimal".
#'
#' @return A ggplot object.
#'
#' @details
#' The function automatically detects the process type from the paths_data attributes
#' and applies appropriate labels and styling. For GBM processes, it includes the
#' SDE equation in the subtitle. For Poisson processes, use `plot_type` to choose
#' between standard or compensated processes.
#'
#' @examples
#' # Simulate and plot GBM
#' gbm_spec(100, 0.05, 0.2) |>
#'   simulate_paths(n_paths = 50, n_steps = 252, maturity = 1) |>
#'   plot_paths(n_paths_plot = 10)
#'
#' # Simulate and plot Poisson
#' poisson_spec(1.0) |>
#'   simulate_paths(25, 500, 30) |>
#'   plot_paths(plot_type = "compensated")
#'
#' @export
plot_paths <- function(paths_data,
                       n_paths_plot = NULL,
                       plot_type = c("standard", "compensated"),
                       alpha = 0.3,
                       theme = c("minimal", "classic", "bw")) {
  # Validation
  # Use assert_integerish to allow numeric values coercible to integers (e.g., 100 instead of 100L)
  checkmate::assert_data_frame(paths_data)
  checkmate::assert_integerish(n_paths_plot, lower = 1, len = 1, null.ok = TRUE, any.missing = FALSE)
  checkmate::assert_number(alpha, lower = 0, upper = 1)
  plot_type <- match.arg(plot_type)
  theme <- match.arg(theme)

  # Coerce to integer for internal use
  if (!is.null(n_paths_plot)) {
    n_paths_plot <- as.integer(n_paths_plot)
  }

  # Get process type from attributes
  process_type <- attr(paths_data, "process_type")
  intensity <- attr(paths_data, "intensity")

  # Filter paths if requested
  if (!is.null(n_paths_plot)) {
    max_path_id <- max(paths_data$path_id)
    paths_to_plot <- min(n_paths_plot, max_path_id)
    paths_data <- dplyr::filter(paths_data, path_id <= paths_to_plot)
  }

  # Determine value column name based on process type and plot_type
  if (process_type == "poisson") {
    value_col <- if (plot_type == "compensated") "compensated_count" else "count"
  } else if ("stock_price" %in% names(paths_data)) {
    value_col <- "stock_price"
  } else {
    value_col <- "value"
  }

  # Create appropriate labels based on process type
  plot_config <- get_plot_config(process_type, intensity, plot_type)
  has_component <- "component" %in% names(paths_data)
  if (has_component) {
    base_aes <- ggplot2::aes(
      x = time,
      y = .data[[value_col]],
      group = interaction(path_id, component),
      color = component
    )
  } else {
    base_aes <- ggplot2::aes(
      x = time,
      y = .data[[value_col]],
      group = path_id
    )
  }

  # Create base plot
  if (has_component) {
    line_layer <- ggplot2::geom_line(
      alpha = alpha,
      linewidth = 0.5
    )
  } else if (is.null(plot_config$color)) {
    line_layer <- ggplot2::geom_line(
      alpha = alpha,
      linewidth = 0.5
    )
  } else {
    line_layer <- ggplot2::geom_line(
      alpha = alpha,
      linewidth = 0.5,
      color = plot_config$color
    )
  }

  p <- ggplot2::ggplot(
    paths_data,
    base_aes
  ) +
    line_layer +
    ggplot2::labs(
      title = plot_config$title,
      subtitle = plot_config$subtitle,
      x = "Time (years)",
      y = plot_config$y_label
    )

  # Add expected value line for Poisson processes
  if (process_type == "poisson") {
    if (plot_type == "standard") {
      # E[N(t)] = λt
      p <- p + ggplot2::geom_line(
        data = data.frame(time = unique(paths_data$time)),
        ggplot2::aes(x = time, y = intensity * time),
        color = "red",
        linetype = "dashed",
        linewidth = 1,
        inherit.aes = FALSE
      )
    } else {
      # E[Ñ(t)] = 0
      p <- p + ggplot2::geom_hline(
        yintercept = 0,
        color = "red",
        linetype = "dashed",
        linewidth = 1
      )
    }
  } else if (process_type == "cir") {
    long_term_mean <- attr(paths_data, "long_term_mean")
    if (!is.null(long_term_mean)) {
      p <- p + ggplot2::geom_hline(
        yintercept = long_term_mean,
        color = "red",
        linetype = "dashed",
        linewidth = 1
      )
    }
  } else if (process_type == "correlated_bm" && has_component) {
    p <- p + ggplot2::guides(color = ggplot2::guide_legend(title = "Component"))
  }

  # Apply theme
  p <- p + switch(theme,
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
#' @param intensity Numeric. Intensity parameter for Poisson process.
#' @param plot_type Character. Type of plot (for Poisson: "standard" or "compensated").
#'
#' @return List with title, subtitle, y_label, and color.
#'
#' @keywords internal
get_plot_config <- function(process_type, intensity = NULL, plot_type = "standard") {
  if (is.null(process_type)) {
    process_type <- "gbm"
  }

  config_generators <- list(
    gbm = function(intensity = NULL, plot_type = "standard") {
      list(
        title = "Geometric Brownian Motion Paths",
        subtitle = "dS(t) = \u03bc S(t) dt + \u03c3 S(t) dW(t)",
        y_label = "Stock Price S(t)",
        color = "steelblue"
      )
    },
    abm = function(intensity = NULL, plot_type = "standard") {
      list(
        title = "Arithmetic Brownian Motion Paths",
        subtitle = "dX(t) = \u03bc dt + \u03c3 dW(t)",
        y_label = "Process Value X(t)",
        color = "darkred"
      )
    },
    poisson = function(intensity = NULL, plot_type = "standard") {
      if (plot_type == "compensated") {
        list(
          title = "Compensated Poisson Process Paths",
          subtitle = sprintf(
            "\u00d1(t) = N(t) - \u03bbt  |  \u03bb = %.2f (martingale)",
            intensity
          ),
          y_label = "Compensated Count \u00d1(t)",
          color = "darkgreen"
        )
      } else {
        list(
          title = "Poisson Process Paths",
          subtitle = sprintf(
            "N(t) ~ Poisson(\u03bbt)  |  \u03bb = %.2f",
            intensity
          ),
          y_label = "Count N(t)",
          color = "darkorange"
        )
      }
    },
    cir = function(intensity = NULL, plot_type = "standard") {
      list(
        title = "CIR Variance Paths",
        subtitle = "dv(t) = \u03ba(\u03b8 - v(t))dt + \u03c3\u221av(t)dW(t)",
        y_label = "Variance v(t)",
        color = "purple"
      )
    },
    correlated_bm = function(intensity = NULL, plot_type = "standard") {
      list(
        title = "Correlated Brownian Motion Paths",
        subtitle = "dX(t) = \u03bc dt + L dW(t)",
        y_label = "Component value",
        color = NULL
      )
    },
    heston = function(intensity = NULL, plot_type = "standard") {
      list(
        title = "Heston Stochastic Volatility Paths",
        subtitle = "Joint stock price and variance under correlated Brownian motions",
        y_label = "Stock price S(t)",
        color = "navy"
      )
    },
    merton = function(intensity = NULL, plot_type = "standard") {
      list(
        title = "Merton Jump-Diffusion Paths",
        subtitle = "dS/S = (r - \u03bb \u03ba) dt + \u03c3 dW + J dN",
        y_label = "Stock price S(t)",
        color = "firebrick"
      )
    },
    default = function(intensity = NULL, plot_type = "standard") {
      list(
        title = "Simulated Paths",
        subtitle = "",
        y_label = "Value",
        color = "black"
      )
    }
  )

  generator <- config_generators[[process_type]]
  if (is.null(generator)) {
    generator <- config_generators$default
  }

  generator(intensity = intensity, plot_type = plot_type)
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
