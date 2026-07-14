#' Effective sample size
#'
#' @description
#' `ess()` is the generic for the effective sample size, a summary of how much
#' information a set of weights retains relative to an unweighted sample. When
#' weights vary substantially, the effective sample size can be much smaller
#' than the number of observations.
#'
#' @details
#' This package defines only the generic. Methods live in the packages that own
#' the relevant objects. A method for a weight vector typically returns a single
#' number, for example the Kish effective sample size
#' \eqn{(\sum w)^2 / \sum w^2}, while a method for a fitted model may return a
#' tibble with one row per group. The generic is intentionally minimal:
#' method-specific arguments, such as whether to remove missing values, are
#' passed through `...`.
#'
#' @param x An object to compute the effective sample size for, such as a
#'   weight vector or a fitted model. `ess()` dispatches on this argument.
#' @param ... Arguments passed to methods.
#'
#' @return A method-defined value. Weight-vector methods generally return a
#'   single number; fitted-model methods may return a tibble.
#'
#' @seealso The `halfmoon` and `balancing` packages for methods.
#'
#' @export
ess <- function(x, ...) {
  UseMethod("ess")
}

#' @export
ess.default <- function(x, ...) {
  stop_no_method("ess", x)
}
