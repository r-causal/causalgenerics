#' Causal weight accessors
#'
#' @description
#' Generics for inspecting and setting the metadata carried by causal weight
#' vectors.
#'
#' * `is_causal_wt()` tests whether an object is a causal weight vector.
#' * `estimand()` returns the causal estimand the weights target, such as
#'   `"ate"` or `"att"`.
#' * `estimand<-()` sets the causal estimand.
#'
#' @details
#' This package defines only the generics. The weight classes and their methods
#' live in the packages that own them, such as `propensity`. `is_causal_wt()`
#' defaults to `FALSE` so that any object that is not a recognized causal weight
#' vector reports as much; `estimand()` and `estimand<-()` have no meaningful
#' default and signal an error when no method is registered for the object.
#'
#' @param x An object to inspect or modify. These generics dispatch on this
#'   argument.
#' @param value The causal estimand to assign.
#' @param ... Arguments passed to methods.
#'
#' @return
#' `is_causal_wt()` returns a single logical value. `estimand()` returns the
#' estimand, typically a character string, or `NULL` when none is recorded.
#' `estimand<-()` returns `x` with the estimand updated.
#'
#' @seealso The `propensity` package for methods.
#'
#' @name causal-weights
NULL

#' @rdname causal-weights
#' @export
is_causal_wt <- function(x, ...) {
  UseMethod("is_causal_wt")
}

#' @export
is_causal_wt.default <- function(x, ...) {
  FALSE
}

#' @rdname causal-weights
#' @export
estimand <- function(x, ...) {
  UseMethod("estimand")
}

#' @export
estimand.default <- function(x, ...) {
  stop_no_method("estimand", x)
}

#' @rdname causal-weights
#' @export
`estimand<-` <- function(x, ..., value) {
  UseMethod("estimand<-")
}

#' @export
`estimand<-.default` <- function(x, ..., value) {
  stop_no_method("estimand<-", x)
}
