#' Normalize Unicode Strings to ASCII
#'
#' Converts Unicode-heavy lecture annotations (for example Greek letters or
#' superscripts) into ASCII equivalents so that documentation and tests remain
#' portable across locales. The helper is primarily used when importing text
#' snippets from the course's Python notebooks.
#'
#' @param x Character vector to normalise.
#'
#' @return Character vector with ASCII-only content.
#'
#' @examples
#' normalize_to_ascii("X ~ N(\u03bc, \u03c3\u00b2)")
#'
#' normalize_to_ascii(c(
#'   "\u03ba\u0304(T\u2081) = 4\u03ba v\u0304 / \u03b3\u00b2",
#'   "\u03bb\u03bc = \u03bb(e^(\u03bc_J + 0.5\u03c3_J\u00b2) - 1)"
#' ))
#' @export
normalize_to_ascii <- function(x) {
  checkmate::assert_character(x, any.missing = FALSE)

  normalized <- x

  multi_char_replacements <- c(
    "\u03ba\u0304" = "kappa_bar",
    "\u0076\u0304" = "v_bar",
    "\u0056\u0304" = "V_bar",
    "\u03b3\u00b2" = "gamma^2",
    "\u03c3\u00b2" = "sigma^2",
    "\u03bc\u2093" = "mu_x"
  )

  normalized <- purrr::reduce2(
    .x = names(multi_char_replacements),
    .y = multi_char_replacements,
    .init = normalized,
    .f = \(state, pattern, replacement) gsub(pattern, replacement, state, fixed = TRUE)
  )

  single_char_replacements <- c(
    "\u03bc" = "mu",
    "\u03c3" = "sigma",
    "\u03bb" = "lambda",
    "\u03c7" = "chi",
    "\u03ba" = "kappa",
    "\u03b3" = "gamma",
    "\u03c4" = "tau",
    "\u03c6" = "phi",
    "\u03c9" = "omega",
    "\u03b1" = "alpha",
    "\u03b2" = "beta",
    "\u03b8" = "theta",
    "\u03c1" = "rho",
    "\u0394" = "Delta",
    "\u00b2" = "^2",
    "\u00b3" = "^3",
    "\u2074" = "^4",
    "\u2080" = "_0",
    "\u2081" = "_1",
    "\u2082" = "_2",
    "\u2083" = "_3",
    "\u00b7" = "*",
    "\u2192" = "->",
    "\u2212" = "-",
    "\u00d7" = "*"
  )

  normalized <- purrr::reduce2(
    .x = names(single_char_replacements),
    .y = single_char_replacements,
    .init = normalized,
    .f = \(state, pattern, replacement) gsub(pattern, replacement, state, fixed = TRUE)
  )

  normalized <- iconv(normalized, from = "UTF-8", to = "ASCII//TRANSLIT", sub = "")
  normalized[is.na(normalized)] <- ""
  normalized <- gsub(" {2,}", " ", normalized)
  trimws(normalized)
}
