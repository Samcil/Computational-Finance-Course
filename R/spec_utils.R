#' Specification Utilities
#'
#' Internal constructors and engine helpers used to build model specifications
#' that mimic tidymodels-style objects.
#'
#' @keywords internal
new_model_spec <- function(class, args, mode) {
  checkmate::assert_list(args, names = "unique")
  checkmate::assert_string(mode)

  structure(
    list(
      args = args,
      eng_args = list(),
      mode = mode,
      method = list(engine = NULL),
      metadata = list()
    ),
    class = c(class, "model_spec")
  )
}


set_engine_base <- function(object, engine, eng_args) {
  object$method$engine <- engine
  object$eng_args <- eng_args
  object
}


#' Create a Process Specification with Hierarchical Classes
#'
#' Internal helper building on the lightweight model spec pattern but exposing
#' arguments at the top level so simulation engines can access parameters
#' directly. Additional inheritance classes ensure that specs share behaviour
#' and documentation via their families (e.g. diffusion, jump_diffusion).
#'
#' @param class Character scalar for the concrete spec class.
#' @param args Named list of validated arguments stored directly on the spec.
#' @param process_type Character scalar identifying the process for downstream
#'   printing and plotting.
#' @param inheritance Character vector of parent classes to prepend (e.g.
#'   "abm_spec" so GBM inherits ABM behaviour).
#' @param family Character scalar naming the process family (e.g.
#'   "diffusion_spec"). When supplied the corresponding family class is added
#'   to the class vector.
#' @param mode Character scalar describing the modelling mode. Defaults to
#'   "stochastic_simulation" for process specs.
#'
#' @keywords internal
new_process_spec <- function(class,
                             args,
                             process_type,
                             inheritance = character(),
                             family = NULL,
                             mode = "stochastic_simulation") {
  checkmate::assert_string(class)
  checkmate::assert_list(args, names = "unique")
  checkmate::assert_string(process_type)
  checkmate::assert_character(inheritance, null.ok = TRUE, any.missing = FALSE)
  if (!is.null(family)) {
    checkmate::assert_string(family)
  }
  checkmate::assert_string(mode)

  spec <- args
  spec$process_type <- process_type
  spec$mode <- mode
  spec$eng_args <- list()
  spec$method <- list(engine = NULL)
  spec$metadata <- list()
  spec$args <- args

  classes <- c(class, inheritance)
  if (!is.null(family)) {
    family_class <- if (grepl("_spec$", family)) family else paste0(family, "_spec")
    classes <- c(classes, family_class)
  }

  classes <- unique(c(classes, "process_spec", "model_spec"))
  class(spec) <- classes
  spec
}


#' Create a Diffusion Process Specification
#'
#' Convenience wrapper above `new_process_spec()` for continuous diffusions.
#'
#' @keywords internal
new_diffusion_spec <- function(class,
                               args,
                               process_type,
                               inheritance = character()) {
  new_process_spec(
    class = class,
    args = args,
    process_type = process_type,
    inheritance = inheritance,
    family = "diffusion"
  )
}


#' Create a Jump Process Specification
#'
#' @keywords internal
new_jump_spec <- function(class,
                          args,
                          process_type,
                          inheritance = character()) {
  new_process_spec(
    class = class,
    args = args,
    process_type = process_type,
    inheritance = inheritance,
    family = "jump"
  )
}


#' Create a Jump-Diffusion Specification
#'
#' Ensures jump diffusions inherit diffusion behaviour while recording jump
#' metadata.
#'
#' @keywords internal
new_jump_diffusion_spec <- function(class,
                                    args,
                                    process_type,
                                    inheritance = character()) {
  new_process_spec(
    class = class,
    args = args,
    process_type = process_type,
    inheritance = c(inheritance, "diffusion_spec"),
    family = "jump_diffusion"
  )
}


#' Create a Stochastic Volatility Specification
#'
#' @keywords internal
new_stochastic_vol_spec <- function(class,
                                    args,
                                    process_type,
                                    inheritance = character()) {
  new_process_spec(
    class = class,
    args = args,
    process_type = process_type,
    inheritance = c(inheritance, "diffusion_spec"),
    family = "stochastic_vol"
  )
}


#' Create a Hybrid Process Specification
#'
#' Convenience wrapper for hybrid equity/interest-rate specifications that
#' combine diffusion components. Ensures hybrid specs share the `hybrid_spec`
#' family class for downstream dispatching and documentation.
#'
#' @keywords internal
new_hybrid_spec <- function(class,
                            args,
                            process_type,
                            inheritance = character()) {
  new_process_spec(
    class = class,
    args = args,
    process_type = process_type,
    inheritance = inheritance,
    family = "hybrid"
  )
}


#' Attach Metadata to a Process Specification
#'
#' @keywords internal
set_process_metadata <- function(spec, ...) {
  checkmate::assert_true(inherits(spec, "process_spec"))
  metadata <- rlang::list2(...)
  spec$metadata <- utils::modifyList(spec$metadata, metadata)
  spec
}


#' Attach Metadata to a Model Specification
#'
#' Generic helper for specs that do not inherit from `process_spec` but still
#' need to publish engine information or auxiliary descriptors.
#'
#' @keywords internal
set_spec_metadata <- function(spec, ...) {
  checkmate::assert_true(inherits(spec, "model_spec"))
  metadata <- rlang::list2(...)
  current <- spec$metadata
  if (is.null(current)) {
    current <- list()
  }
  spec$metadata <- utils::modifyList(current, metadata)
  spec
}


#' Process Specification Lineage
#'
#' Returns the ordered class lineage for a process specification, including the
#' concrete spec, parent specs, and shared family classes such as
#' `diffusion_spec` or `jump_diffusion_spec`.
#'
#' @param spec A process specification inheriting from `process_spec`.
#' @return Character vector of class names in precedence order.
#' @export
spec_lineage <- function(spec) {
  checkmate::assert_true(inherits(spec, "process_spec"))
  lineage <- class(spec)
  lineage[!lineage %in% c("list", "model_spec")]
}


format_spec_lineage <- function(spec, collapse = " -> ") {
  paste(spec_lineage(spec), collapse = collapse)
}


#' Set the Engine for a Model Specification
#'
#' Lightweight analogue of tidymodels' `set_engine()` tailored for
#' CompFinanceR specs.
#'
#' @param object A specification inheriting from `model_spec`.
#' @param engine Character string identifying the engine.
#' @param ... Engine-specific options stored with the specification.
#'
#' @return Updated specification with engine metadata recorded.
#'
#' @export
set_engine <- function(object, engine, ...) {
  UseMethod("set_engine")
}


#' @export
set_engine.model_spec <- function(object, engine, ...) {
  checkmate::assert_string(engine)
  eng_args <- rlang::list2(...)
  set_engine_base(object, engine = engine, eng_args = eng_args)
}


#' @export
set_engine.default <- function(object, engine, ...) {
  rlang::abort(
    message = "No set_engine() method exists for this specification",
    class = "set_engine_not_implemented"
  )
}
