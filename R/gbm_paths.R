#' Generate Geometric Brownian Motion and Arithmetic Brownian Motion Paths
#'
#' Simulates stock price paths under the risk-neutral measure using the
#' Euler-Maruyama discretization scheme. The ABM process X(t) follows:
#'   dX(t) = (r - 0.5*sigma^2)*dt + sigma*dW(t)
#' And the GBM process S(t) is obtained as:
#'   S(t) = exp(X(t))
#'
#' @param n_paths Integer. Number of Monte Carlo paths to simulate.
#' @param n_steps Integer. Number of time steps for discretization.
#' @param maturity Numeric. Time to maturity in years.
#' @param interest_rate Numeric. Risk-free interest rate (annualized).
#' @param volatility Numeric. Volatility parameter (annualized).
#' @param initial_price Numeric. Initial stock price at time t=0.
#' @param seed Integer. Random seed for reproducibility. Default is 123.
#' @param return_format Character. Format for returned data: "tidy" (default) returns
#'   a tibble in long format; "matrix" returns a list with matrices.
#'
#' @return If return_format="tidy", a tibble with columns:
#'   \describe{
#'     \item{path_id}{Path identifier (1 to n_paths)}
#'     \item{time}{Time point}
#'     \item{log_price}{Log of stock price (ABM process X(t))}
#'     \item{stock_price}{Simulated stock price (GBM process S(t))}
#'   }
#'   If return_format="matrix", a list with:
#'   \describe{
#'     \item{time}{Time grid vector of length (n_steps + 1)}
#'     \item{log_price}{Matrix of ABM paths (n_paths x n_steps+1)}
#'     \item{stock_price}{Matrix of GBM paths (n_paths x n_steps+1)}
#'   }
#'
#' @details
#' The function implements the Euler-Maruyama discretization scheme for the ABM:
#'   X[i+1] = X[i] + (r - 0.5*sigma^2)*dt + sigma*sqrt(dt)*Z[i]
#' where Z[i] ~ N(0,1) are standard normal random variables.
#'
#' For multiple paths, the normal samples are standardized at each time step
#' to ensure mean 0 and variance 1.
#'
#' Uses purrr for functional programming and checkmate for robust input validation.
#'
#' @examples
#' # Simulate 100 paths over 1 year with 252 time steps
#' paths_tidy <- generate_gbm_abm_paths(
#'   n_paths = 100,
#'   n_steps = 252,
#'   maturity = 1.0,
#'   interest_rate = 0.05,
#'   volatility = 0.2,
#'   initial_price = 100
#' )
#'
#' # Return in matrix format for faster computation
#' paths_matrix <- generate_gbm_abm_paths(
#'   n_paths = 1000,
#'   n_steps = 500,
#'   maturity = 1.0,
#'   interest_rate = 0.05,
#'   volatility = 0.4,
#'   initial_price = 100,
#'   return_format = "matrix"
#' )
#'
#' @export
generate_gbm_abm_paths <- function(n_paths,
                                    n_steps,
                                    maturity,
                                    interest_rate,
                                    volatility,
                                    initial_price,
                                    seed = 123,
                                    return_format = c("tidy", "matrix")) {
  
  # Input validation using checkmate
  # Use assert_integerish to allow numeric values coercible to integers (e.g., 100 instead of 100L)
  checkmate::assert_integerish(n_paths, lower = 1, len = 1, any.missing = FALSE)
  checkmate::assert_integerish(n_steps, lower = 1, len = 1, any.missing = FALSE)
  checkmate::assert_number(maturity, lower = 0, finite = TRUE)
  checkmate::assert_number(interest_rate, finite = TRUE)
  checkmate::assert_number(volatility, lower = 0, finite = TRUE)
  checkmate::assert_number(initial_price, lower = 0, finite = TRUE)
  checkmate::assert_integerish(seed, len = 1, any.missing = FALSE)
  return_format <- match.arg(return_format)
  
  # Coerce to integer for internal use
  n_paths <- as.integer(n_paths)
  n_steps <- as.integer(n_steps)
  seed <- as.integer(seed)
  
  # Set random seed for reproducibility
  set.seed(seed)
  
  # Time step size
  dt <- maturity / n_steps
  
  # Generate random normal samples
  Z <- matrix(stats::rnorm(n_paths * n_steps, mean = 0, sd = 1),
              nrow = n_paths, ncol = n_steps)
  
  # Standardize each column using purrr::map_dfc
  if (n_paths > 1) {
    Z <- purrr::map_dfc(
      seq_len(n_steps),
      ~ {
        z_col <- Z[, .x]
        (z_col - mean(z_col)) / stats::sd(z_col)
      }
    ) |>
      as.matrix()
  }
  
  # Initialize with initial condition
  X_initial <- rep(log(initial_price), n_paths)
  
  # Compute increments using purrr::accumulate
  # Each step adds drift + diffusion
  X_list <- purrr::accumulate(
    seq_len(n_steps),
    function(x_prev, step_idx) {
      x_prev + 
        (interest_rate - 0.5 * volatility^2) * dt + 
        volatility * sqrt(dt) * Z[, step_idx]
    },
    .init = X_initial
  )
  
  # Convert list to matrix (each element is a vector of all paths at time t)
  X <- do.call(cbind, X_list)
  
  # Convert ABM to GBM: S(t) = exp(X(t))
  S <- exp(X)
  
  # Time grid using purrr::map_dbl
  time_grid <- purrr::map_dbl(0:n_steps, ~ .x * dt)
  
  # Return in requested format
  if (return_format == "matrix") {
    return(list(
      time = time_grid,
      log_price = X,
      stock_price = S
    ))
  } else {
    # Convert to tidy format using purrr::map_dfr
    paths_tidy <- purrr::map_dfr(
      seq_len(n_paths),
      ~ tibble::tibble(
        path_id = .x,
        time = time_grid,
        log_price = X[.x, ],
        stock_price = S[.x, ]
      )
    )
    
    return(paths_tidy)
  }
}


#' Plot GBM and ABM Paths
#'
#' Create visualizations of simulated GBM and ABM paths using ggplot2.
#'
#' @param paths_data Either a tibble from generate_gbm_abm_paths() with return_format="tidy",
#'   or a list with return_format="matrix". If list, it will be converted to tidy format.
#' @param plot_type Character. Type of plot: "both" (default), "gbm", or "abm".
#' @param max_paths Integer. Maximum number of paths to plot. Default is 25.
#'   If more paths exist, a random sample will be plotted.
#' @param theme Character. ggplot2 theme to use: "minimal" (default), "classic", or "bw".
#'
#' @return A ggplot2 object or list of ggplot2 objects.
#'
#' @examples
#' paths <- generate_gbm_abm_paths(
#'   n_paths = 100,
#'   n_steps = 252,
#'   maturity = 1.0,
#'   interest_rate = 0.05,
#'   volatility = 0.2,
#'   initial_price = 100
#' )
#'
#' # Plot both GBM and ABM
#' plot_gbm_abm_paths(paths, plot_type = "both")
#'
#' # Plot only GBM paths
#' plot_gbm_abm_paths(paths, plot_type = "gbm", max_paths = 10)
#'
#' @export
plot_gbm_abm_paths <- function(paths_data,
                                plot_type = c("both", "gbm", "abm"),
                                max_paths = 25,
                                theme = c("minimal", "classic", "bw")) {
  
  # Input validation using checkmate
  # Use assert_integerish to allow numeric values coercible to integers (e.g., 100 instead of 100L)
  checkmate::assert_integerish(max_paths, lower = 1, len = 1, any.missing = FALSE)
  plot_type <- match.arg(plot_type)
  theme <- match.arg(theme)
  
  # Coerce to integer for internal use
  max_paths <- as.integer(max_paths)
  
  # Convert to tidy format if needed
  if (is.list(paths_data) && !tibble::is_tibble(paths_data)) {
    # Convert matrix format to tidy using purrr::map_dfr
    n_paths <- nrow(paths_data$log_price)
    
    paths_data <- purrr::map_dfr(
      seq_len(n_paths),
      ~ tibble::tibble(
        path_id = .x,
        time = paths_data$time,
        log_price = paths_data$log_price[.x, ],
        stock_price = paths_data$stock_price[.x, ]
      )
    )
  }
  
  # Sample paths if too many
  unique_paths <- unique(paths_data$path_id)
  if (length(unique_paths) > max_paths) {
    sampled_paths <- sample(unique_paths, max_paths)
    paths_data <- paths_data |>
      dplyr::filter(path_id %in% sampled_paths)
  }
  
  # Select theme
  theme_func <- switch(theme,
                       "minimal" = ggplot2::theme_minimal(),
                       "classic" = ggplot2::theme_classic(),
                       "bw" = ggplot2::theme_bw())
  
  # Create plots using purrr::map
  plot_specs <- list(
    abm = list(
      y_var = "log_price",
      y_lab = "Log Price X(t)",
      color = "steelblue",
      title = "Arithmetic Brownian Motion (ABM) Paths",
      subtitle = "Simulated paths: X(t) with dX(t) = (r - 0.5σ²)dt + σdW(t)"
    ),
    gbm = list(
      y_var = "stock_price",
      y_lab = "Stock Price S(t)",
      color = "darkred",
      title = "Geometric Brownian Motion (GBM) Paths",
      subtitle = "Stock price: S(t) = exp(X(t))"
    )
  )
  
  # Filter plot specs based on plot_type
  if (plot_type == "abm") {
    plot_specs <- plot_specs["abm"]
  } else if (plot_type == "gbm") {
    plot_specs <- plot_specs["gbm"]
  }
  
  # Generate plots using purrr::map
  plots <- purrr::map(
    plot_specs,
    ~ ggplot2::ggplot(
        paths_data, 
        ggplot2::aes(x = time, y = .data[[.x$y_var]], group = path_id)
      ) +
      ggplot2::geom_line(alpha = 0.5, color = .x$color) +
      ggplot2::labs(
        title = .x$title,
        subtitle = .x$subtitle,
        x = "Time (years)",
        y = .x$y_lab
      ) +
      theme_func +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 14),
        plot.subtitle = ggplot2::element_text(size = 10)
      )
  )
  
  # Return single plot or list
  if (length(plots) == 1) {
    return(plots[[1]])
  } else {
    return(plots)
  }
}


#' Generate Example GBM/ABM Paths with Visualization
#'
#' Convenience function that generates and plots GBM/ABM paths with default parameters.
#' Useful for quick demonstrations and examples.
#'
#' @param n_paths Integer. Number of paths. Default is 25.
#' @param n_steps Integer. Number of steps. Default is 500.
#' @param maturity Numeric. Time to maturity. Default is 1.0.
#' @param interest_rate Numeric. Interest rate. Default is 0.05.
#' @param volatility Numeric. Volatility. Default is 0.4.
#' @param initial_price Numeric. Initial price. Default is 100.
#' @param plot Logical. Whether to create plots. Default is TRUE.
#'
#' @return A list containing:
#'   \describe{
#'     \item{paths}{The simulated paths (tibble)}
#'     \item{plots}{List of ggplot2 objects (if plot=TRUE)}
#'   }
#'
#' @examples
#' # Generate and plot with defaults
#' result <- demo_gbm_abm_paths()
#'
#' # Access paths
#' head(result$paths)
#'
#' # Display plots
#' print(result$plots$gbm)
#' print(result$plots$abm)
#'
#' @export
demo_gbm_abm_paths <- function(n_paths = 25,
                                n_steps = 500,
                                maturity = 1.0,
                                interest_rate = 0.05,
                                volatility = 0.4,
                                initial_price = 100,
                                plot = TRUE) {
  
  # Input validation using checkmate
  # Use assert_integerish to allow numeric values coercible to integers (e.g., 100 instead of 100L)
  checkmate::assert_integerish(n_paths, lower = 1, len = 1, any.missing = FALSE)
  checkmate::assert_integerish(n_steps, lower = 1, len = 1, any.missing = FALSE)
  checkmate::assert_number(maturity, lower = 0, finite = TRUE)
  checkmate::assert_number(interest_rate, finite = TRUE)
  checkmate::assert_number(volatility, lower = 0, finite = TRUE)
  checkmate::assert_number(initial_price, lower = 0, finite = TRUE)
  checkmate::assert_logical(plot, len = 1)
  
  # Coerce to integer for internal use
  n_paths <- as.integer(n_paths)
  n_steps <- as.integer(n_steps)
  
  # Generate paths
  paths <- generate_gbm_abm_paths(
    n_paths = n_paths,
    n_steps = n_steps,
    maturity = maturity,
    interest_rate = interest_rate,
    volatility = volatility,
    initial_price = initial_price,
    return_format = "tidy"
  )
  
  # Create result list
  result <- list(paths = paths)
  
  # Add plots if requested
  if (plot) {
    result$plots <- plot_gbm_abm_paths(paths, plot_type = "both")
  }
  
  return(result)
}
