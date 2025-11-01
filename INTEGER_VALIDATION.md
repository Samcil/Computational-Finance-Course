# Integer Validation with checkmate

## Overview

The package uses `checkmate::assert_integerish()` instead of `checkmate::assert_int()` for integer parameter validation. This provides a more user-friendly experience while maintaining type safety.

## Rationale

### User-Friendly Behavior

R users commonly write numeric literals without the `L` suffix (e.g., `100` instead of `100L`). Requiring strict integer types would force users to remember to add `L`:

```r
# Without relaxed validation - requires L suffix
simulate_paths(spec, n_paths = 100L, n_steps = 252L, maturity = 1.0)  # Annoying!

# With relaxed validation - natural numeric literals work
simulate_paths(spec, n_paths = 100, n_steps = 252, maturity = 1.0)   # Natural!
```

### Type Safety Maintained

`assert_integerish()` still validates that values can be safely coerced to integers:

```r
# These work (can be coerced without loss)
n_paths = 100    # numeric -> integer (OK)
n_paths = 100.0  # numeric -> integer (OK)
n_paths = 100L   # already integer (OK)

# These fail (cannot be coerced safely)
n_paths = 100.5  # fractional (ERROR)
n_paths = NA     # missing (ERROR)
n_paths = "100"  # character (ERROR)
```

## Implementation Pattern

All functions that accept integer parameters follow this pattern:

```r
my_function <- function(n_paths, n_steps, seed = 123) {
  # Validation: Allow coercible numeric values
  checkmate::assert_integerish(n_paths, lower = 1, len = 1, any.missing = FALSE)
  checkmate::assert_integerish(n_steps, lower = 1, len = 1, any.missing = FALSE)
  checkmate::assert_integerish(seed, len = 1, any.missing = FALSE)
  
  # Coercion: Convert to integer for internal use
  n_paths <- as.integer(n_paths)
  n_steps <- as.integer(n_steps)
  seed <- as.integer(seed)
  
  # Function body uses integer values
  # ...
}
```

## Functions Updated

All functions with integer parameters have been updated:

- `simulate_paths()` - n_paths, n_steps, seed
- `gbm_spec()`, `abm_spec()` - no integer parameters
- `plot_paths()` - n_paths_plot
- `generate_gbm_abm_paths()` - n_paths, n_steps, seed
- `plot_gbm_abm_paths()` - max_paths
- `demo_gbm_abm_paths()` - n_paths, n_steps

## Benefits

1. **Natural R syntax**: Users can write numbers naturally without `L` suffix
2. **Type safety**: Invalid values (fractional, NA, non-numeric) still caught
3. **Clear errors**: checkmate provides informative error messages
4. **Consistent behavior**: All functions use same validation pattern
5. **No silent errors**: Coercion is explicit, not implicit

## Example Error Messages

```r
# Fractional value
simulate_paths(spec, n_paths = 100.5, n_steps = 252, maturity = 1.0)
# Error: Assertion on 'n_paths' failed: Must be of type 'integerish', but has no integer values.

# Character value
simulate_paths(spec, n_paths = "100", n_steps = 252, maturity = 1.0)
# Error: Assertion on 'n_paths' failed: Must be of type 'integerish', not 'character'.

# Missing value
simulate_paths(spec, n_paths = NA, n_steps = 252, maturity = 1.0)
# Error: Assertion on 'n_paths' failed: Contains missing values.
```

## References

- checkmate documentation: `?checkmate::assert_integerish`
- tidymodels hardhat vignette on input validation
- R Packages book on function arguments
