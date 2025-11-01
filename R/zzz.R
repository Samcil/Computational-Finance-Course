register_print_method <- function(pkgname, class) {
  method_name <- paste0("print.", class)
  method <- get(method_name, envir = asNamespace(pkgname))
  base::registerS3method("print", class, method, envir = asNamespace(pkgname))
}

.onLoad <- function(libname, pkgname) {
  register_print_method(pkgname, "gbm_spec")
  register_print_method(pkgname, "abm_spec")
  register_print_method(pkgname, "poisson_spec")
  register_print_method(pkgname, "heston_spec")
  register_print_method(pkgname, "merton_spec")
}
