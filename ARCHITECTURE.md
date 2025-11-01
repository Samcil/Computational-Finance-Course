# R Package Architecture: tidymodels/hardhat Design Principles

## Overview

The CompFinanceR package implements a modern, modular architecture inspired by tidymodels and hardhat design principles. This document explains the architectural choices and how they benefit both users and developers.

## Design Philosophy

### 1. **Separation of Interface and Implementation (hardhat principle)**

Following hardhat's philosophy, we separate:
- **User Interface**: Simple, consistent functions that users interact with
- **Engine Implementation**: Complex computational logic hidden from users
- **Bridge Layer**: Methods that connect interface to engines

### 2. **Spec-Fit-Predict Pattern (tidymodels principle)**

All stochastic processes follow a consistent three-step workflow:

1. **Specify**: Create a process specification with parameters
   ```r
   spec <- gbm_spec(initial_value = 100, drift = 0.05, volatility = 0.2)
   ```

2. **Fit/Simulate**: Generate sample paths using the specification
   ```r
   paths <- simulate_paths(spec, n_paths = 100, n_steps = 252, maturity = 1.0)
   ```

3. **Predict/Analyze**: Visualize or extract results
   ```r
   plot_paths(paths)
   ```

### 3. **Functional Programming with purrr**

- **Zero for loops or apply functions** throughout the codebase
- Use modern purrr functions:
  - `purrr::accumulate()` - For cumulative operations (replaces for loops)
  - `purrr::map()` - For transformations
  - `purrr::list_rbind()` - For combining data frames (modern replacement for `map_dfr`)
  - Anonymous function syntax `\(x)` (R >= 4.1.0)

### 4. **Robust Input Validation with checkmate**

- Every user-facing function validates inputs with checkmate
- Provides clear, actionable error messages
- Fast C-based validation
- Self-documenting function contracts

## Architecture Layers

### Layer 1: User Interface (Public API)

**Files**: `R/simulate_paths.R`, `R/plot_paths.R`

**Purpose**: Provide simple, consistent interface for users

**Key Functions**:
- `gbm_spec()`, `abm_spec()`, etc. - Create process specifications
- `simulate_paths()` - Generic function for all processes
- `plot_paths()` - Universal plotting function
- `demo_paths()` - Quick demonstrations

**Design Principles**:
- Consistent argument ordering (process_spec first, then simulation params)
- Pipeable with `|>`
- Clear, descriptive names (snake_case)
- Comprehensive roxygen2 documentation
- Input validation at the entry point

**Example**:
```r
# Simple, consistent interface regardless of process complexity
gbm_spec(100, 0.05, 0.2) |>
  simulate_paths(n_paths = 50, n_steps = 252, maturity = 1.0) |>
  plot_paths(n_paths_plot = 10)
```

### Layer 2: Engine Implementation (Internal Logic)

**Files**: `R/simulate_paths_engines.R`

**Purpose**: Handle computational complexity

**Key Components**:
- `simulate_paths.gbm_spec()` - GBM engine
- `simulate_paths.abm_spec()` - ABM engine
- `generate_standardized_normals()` - Shared utility
- `compute_cumulative_paths()` - Functional cumulative computation
- `matrix_to_tidy()` - Format conversion

**Design Principles**:
- S3 method dispatch based on spec class
- Pure functional programming (no side effects except RNG)
- Modular helper functions (DRY principle)
- Hidden from users (not exported)
- Optimized for performance

**Example** (internal):
```r
# Engine method called via dispatch
simulate_paths.gbm_spec <- function(process_spec, n_paths, n_steps, maturity, seed) {
  # Complex implementation hidden from users
  # Uses functional programming throughout
  # Returns standardized tidy tibble
}
```

### Layer 3: Utilities and Helpers

**Purpose**: Reusable components shared across engines

**Key Functions**:
- `generate_standardized_normals()` - Random number generation
- `compute_cumulative_paths()` - Functional path computation
- `matrix_to_tidy()` - Format conversion
- `get_plot_config()` - Visualization configuration

**Design Principles**:
- Pure functions (no side effects)
- Highly reusable
- Well-tested
- Documented internally

## Modern R Practices

### 1. No Deprecated Functions

**Avoided**: `purrr::map_dfr()`, `purrr::map_dfc()` (superseded)

**Used Instead**:
- `purrr::list_rbind()` - Modern row-binding
- `do.call(cbind, list)` - For column-binding with clear intent
- `purrr::map()` + explicit binding

### 2. Anonymous Function Syntax

Using R >= 4.1.0 syntax for clarity:

```r
# Old style
purrr::map(x, function(i) i^2)

# New style (clearer, more concise)
purrr::map(x, \(i) i^2)
```

### 3. Type-Stable Operations

Using purrr's typed variants:
- `map_dbl()` - Returns numeric vector
- `map_int()` - Returns integer vector
- `map_chr()` - Returns character vector

Ensures outputs are always the expected type, preventing surprises.

### 4. Tidy Evaluation

Using `rlang` for dynamic column names:
```r
tibble::tibble(
  path_id = i,
  time = time_grid,
  !!value_name := value_matrix[i, ]  # Dynamic column naming
)
```

## Consistency Across Processes

All stochastic processes share the same interface:

```r
# GBM
gbm_spec(initial_value, drift, volatility) |> simulate_paths(...)

# ABM
abm_spec(initial_value, drift, volatility) |> simulate_paths(...)

# Future: CIR
cir_spec(initial_value, mean_reversion, long_term_mean, volatility) |> simulate_paths(...)

# Future: Heston
heston_spec(initial_price, initial_variance, ...) |> simulate_paths(...)
```

**Benefits**:
- Easy to learn (learn once, apply everywhere)
- Consistent return format (always tidy tibbles)
- Same plotting function works for all processes
- Easy to swap processes in analysis pipelines

## DRY Principle Implementation

### Avoiding Repetition

1. **Shared validation logic**: checkmate assertions in user-facing functions
2. **Shared random number generation**: `generate_standardized_normals()`
3. **Shared path computation**: `compute_cumulative_paths()`
4. **Shared format conversion**: `matrix_to_tidy()`
5. **Shared plotting**: Single `plot_paths()` for all processes

### Modularity Benefits

- **Easier testing**: Test each component in isolation
- **Easier maintenance**: Fix bugs in one place
- **Easier extension**: Add new processes by implementing one method
- **Better performance**: Optimize shared utilities once

## Extension Pattern

Adding a new process requires:

1. **Create specification function**:
   ```r
   new_process_spec <- function(param1, param2, ...) {
     # Validation
     # Return spec object with appropriate class
   }
   ```

2. **Implement engine method**:
   ```r
   simulate_paths.new_process_spec <- function(process_spec, ...) {
     # Use shared utilities
     # Return tidy tibble
   }
   ```

3. **Add plot configuration**:
   ```r
   # Update get_plot_config() in plot_paths.R
   ```

That's it! The rest works automatically through dispatch and shared utilities.

## Testing Strategy

### Unit Tests

- Test each specification function independently
- Test each engine method independently  
- Test shared utilities independently
- Test error handling and validation

### Integration Tests

- Test full pipelines (spec → simulate → plot)
- Test reproducibility (same seed → same results)
- Test pipeable workflows

### Property-Based Tests

- Test mathematical properties (initial conditions, dimensions, etc.)
- Test tidy data invariants

## Performance Considerations

1. **Vectorization**: Operations are vectorized where possible
2. **Functional programming**: No loop overhead, compiler-friendly
3. **Memory efficiency**: Avoid unnecessary copies
4. **Lazy evaluation**: Compute only what's needed

## Documentation Standards

### roxygen2 Tags

Every exported function includes:
- `@param` - All parameters with types and constraints
- `@return` - Return value structure
- `@examples` - Working examples
- `@details` - Implementation details and equations
- `@export` - Only for public API

### Mathematical Notation

Include SDEs and equations in documentation:
```r
#' dS(t) = r*S(t)*dt + sigma*S(t)*dW(t)
```

### Usage Examples

Show both basic and pipeline usage:
```r
#' @examples
#' # Basic usage
#' spec <- gbm_spec(100, 0.05, 0.2)
#' paths <- simulate_paths(spec, 100, 252, 1.0)
#'
#' # Pipeline usage
#' gbm_spec(100, 0.05, 0.2) |>
#'   simulate_paths(100, 252, 1.0) |>
#'   plot_paths()
```

## Benefits of This Architecture

### For Users

- **Simple, consistent interface** across all processes
- **Pipeline-friendly** design works naturally with tidyverse
- **Clear error messages** from checkmate validation
- **Comprehensive documentation** with examples
- **Type-safe operations** prevent surprises

### For Developers

- **Modular design** makes code easy to understand and modify
- **DRY principle** means fixing bugs in one place
- **Pure functions** are easy to test
- **Clear separation** between interface and implementation
- **Extensible** - adding new processes is straightforward

### For Maintainers

- **Consistent patterns** across codebase
- **Well-tested** at multiple levels
- **Modern R practices** ensure longevity
- **Following community standards** (tidyverse, tidymodels)

## References

- [hardhat: Construct Modeling Packages](https://hardhat.tidymodels.org/)
- [tidymodels principles](https://www.tidymodels.org/learn/)
- [Tidyverse style guide](https://style.tidyverse.org/)
- [purrr documentation](https://purrr.tidyverse.org/)
- [checkmate documentation](https://mllg.github.io/checkmate/)
