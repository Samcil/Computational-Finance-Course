# Mathematical Notation and Rendering

## Overview

The R package uses professional mathematical notation throughout, ensuring that formulas are properly rendered in documentation, help files, and console output.

## LaTeX Notation in roxygen2 Documentation

### Inline Equations

Use `\eqn{formula}` for inline mathematical expressions:

```r
#' @param drift Numeric. Drift parameter \eqn{\mu}.
#' @param volatility Numeric. Volatility parameter \eqn{\sigma}.
```

**Renders as:** "Drift parameter μ" and "Volatility parameter σ"

### Display Equations

Use `\deqn{formula}` for centered, display-style equations:

```r
#' @details
#' The process follows:
#' \deqn{dS(t) = \mu S(t) dt + \sigma S(t) dW(t)}
```

**Renders as:** A centered, larger equation in help files

### Common Mathematical Notation

```r
# Greek letters
\eqn{\mu}      # μ (mu)
\eqn{\sigma}   # σ (sigma)
\eqn{\lambda}  # λ (lambda)
\eqn{\theta}   # θ (theta)

# Subscripts and superscripts
\eqn{S_0}      # S₀
\eqn{X_t}      # Xₜ
\eqn{\sigma^2} # σ²

# Fractions
\eqn{\frac{\sigma^2}{2}}  # σ²/2

# Exponentials and functions
\eqn{\exp(x)}
\eqn{\log(x)}
\eqn{\sqrt{x}}
```

## Unicode Symbols in cli Output

### Greek Letters (Direct Unicode)

```r
cli::cli_text("Drift (\u03bc): {value}")        # μ
cli::cli_text("Volatility (\u03c3): {value}")   # σ
cli::cli_text("Lambda (\u03bb): {value}")       # λ
```

### Subscripts (Unicode)

```r
cli::cli_text("Initial value (S\u2080): {value}")  # S₀
cli::cli_text("Initial value (X\u2080): {value}")  # X₀
```

### Common Unicode Math Symbols

| Symbol | Unicode | Description |
|--------|---------|-------------|
| μ | `\u03bc` | Greek mu |
| σ | `\u03c3` | Greek sigma |
| λ | `\u03bb` | Greek lambda |
| θ | `\u03b8` | Greek theta |
| ρ | `\u03c1` | Greek rho |
| ν | `\u03bd` | Greek nu |
| ₀ | `\u2080` | Subscript 0 |
| ₁ | `\u2081` | Subscript 1 |
| ² | `\u00b2` | Superscript 2 |
| ≈ | `\u2248` | Approximately equal |
| ≤ | `\u2264` | Less than or equal |
| ≥ | `\u2265` | Greater than or equal |
| ∞ | `\u221e` | Infinity |

## cli Package Features

### Semantic Formatting

```r
# Headers
cli::cli_h1("Main Title")
cli::cli_h2("Section Title")
cli::cli_h3("Subsection Title")

# Lists
cli::cli_dl(c(
  "Parameter" = "Value",
  "Another" = "Value"
))

# Alerts
cli::cli_alert_info("Information message")
cli::cli_alert_success("Success message")
cli::cli_alert_warning("Warning message")
cli::cli_alert_danger("Error message")
```

### Inline Styling

```r
# Function references (clickable in RStudio)
cli::cli_text("Use {.fn function_name}")

# Variables
cli::cli_text("The {.var variable_name} is important")

# Field names
cli::cli_text("The {.field field_name} contains data")

# Emphasis
cli::cli_text("This is {.emph emphasized}")

# Strong emphasis
cli::cli_text("This is {.strong important}")
```

### Colors

```r
# Predefined colors
cli::col_cyan("Cyan text")
cli::col_green("Green text")
cli::col_blue("Blue text")
cli::col_red("Red text")
cli::col_yellow("Yellow text")
cli::col_magenta("Magenta text")

# Usage in definition lists
cli::cli_dl(c(
  "Initial value" = cli::col_cyan("{value}"),
  "Drift" = cli::col_green("{value}"),
  "Volatility" = cli::col_blue("{value}")
))
```

## ggplot2 Plot Annotations

### Mathematical Expressions in Plot Text

For plot titles and labels, use Unicode directly or `expression()`:

```r
# Using Unicode (simpler, works everywhere)
subtitle <- sprintf(
  "dS(t) = \u03bc S(t) dt + \u03c3 S(t) dW(t)  |  \u03bc = %.4f",
  drift_value
)

# Using expression() for more complex math
p + labs(
  subtitle = expression(
    dS(t) == mu * S(t) * dt + sigma * S(t) * dW(t)
  )
)
```

## Documentation Rendering Platforms

### R Help System (`?function`)

- **LaTeX equations**: Fully rendered
- **Inline math**: `\eqn{}` rendered as unicode/text
- **Display math**: `\deqn{}` rendered centered
- **Platform**: Works in RStudio help pane, terminal, HTML help

### pkgdown Websites

- **LaTeX equations**: Rendered with MathJax (beautiful HTML)
- **Display equations**: Centered, professional typesetting
- **Inline equations**: Properly styled inline
- **Live examples**: Can show interactive demonstrations

### PDF Documentation

- **LaTeX equations**: Professional typeset (native LaTeX)
- **Best quality**: Full LaTeX rendering engine
- **Print-ready**: Publication quality output

## Best Practices

### 1. Consistency

Use the same notation throughout:
- `\mu` for drift (not `r`, `m`, or `drift`)
- `\sigma` for volatility (not `vol` or `sd`)
- `S_0` or `S(0)` for initial value (consistent choice)

### 2. Documentation

Always include:
- Parameter descriptions with notation: `@param drift Numeric. Drift \eqn{\mu}.`
- Full SDE in `@details` section
- Analytical solutions when available

### 3. Print Methods

Use cli for all console output:
```r
print.my_spec <- function(x, ...) {
  cli::cli_h2("My Specification")
  cli::cli_text("Process: {.emph [SDE with unicode]}")
  cli::cli_dl(c(
    "Parameter (\u03bc)" = cli::col_cyan("{x$param}")
  ))
  cli::cli_alert_info("Use {.fn next_function} to continue")
  invisible(x)
}
```

### 4. Testing

Verify rendering in:
- RStudio help pane (`?function`)
- Terminal R session
- pkgdown website (if deployed)
- PDF manual (R CMD check)

## Examples

### Complete Function Documentation

```r
#' Generate Geometric Brownian Motion Paths
#'
#' Simulates sample paths under the risk-neutral measure.
#'
#' @param initial_value Numeric. Initial stock price \eqn{S_0}.
#' @param drift Numeric. Drift parameter \eqn{\mu} (risk-free rate).
#' @param volatility Numeric. Volatility \eqn{\sigma} (annualized).
#'
#' @return A tibble with columns `path_id`, `time`, and `value`.
#'
#' @details
#' ## Stochastic Differential Equation
#'
#' The GBM process follows:
#' \deqn{dS(t) = \mu S(t) dt + \sigma S(t) dW(t)}
#'
#' Where:
#' \itemize{
#'   \item \eqn{S(t)} = Stock price at time \eqn{t}
#'   \item \eqn{\mu} = Drift (risk-free rate)
#'   \item \eqn{\sigma} = Volatility
#'   \item \eqn{W(t)} = Standard Brownian motion
#' }
#'
#' ## Analytical Solution
#'
#' \deqn{S(t) = S_0 \exp\left[(\mu - \frac{\sigma^2}{2})t + \sigma W(t)\right]}
#'
#' @examples
#' paths <- generate_gbm_paths(
#'   initial_value = 100,
#'   drift = 0.05,
#'   volatility = 0.2
#' )
#'
#' @export
```

### Complete Print Method

```r
#' @export
print.gbm_spec <- function(x, ...) {
  cli::cli_h2("Geometric Brownian Motion Specification")
  cli::cli_text("Process: {.emph dS(t) = \u03bc S(t) dt + \u03c3 S(t) dW(t)}")
  cli::cli_text("")
  cli::cli_dl(c(
    "Initial value (S\u2080)" = cli::col_cyan("{x$initial_value}"),
    "Drift (\u03bc)" = cli::col_green("{x$drift}"),
    "Volatility (\u03c3)" = cli::col_blue("{x$volatility}")
  ))
  cli::cli_text("")
  cli::cli_alert_info("Use {.fn simulate_paths} to generate sample paths")
  invisible(x)
}
```

## References

- [roxygen2 documentation](https://roxygen2.r-lib.org/)
- [cli package](https://cli.r-lib.org/)
- [LaTeX Mathematics](https://en.wikibooks.org/wiki/LaTeX/Mathematics)
- [Unicode Mathematical Symbols](https://en.wikipedia.org/wiki/Mathematical_operators_and_symbols_in_Unicode)
