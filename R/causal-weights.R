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
#' This package defines the generics and the abstract weight class they are
#' written against. It supplies a method for each generic on `causal_wts`, so a
#' concrete class built with [new_causal_wts()] inherits all three. The concrete
#' classes themselves live in the packages that own them, such as `propensity`,
#' as do methods for any object that is not a `causal_wts` vector.
#' `is_causal_wt()` defaults to `FALSE` so that any object that is not a
#' recognized causal weight vector reports as much; `estimand()` and
#' `estimand<-()` have no meaningful default and signal an error when no method
#' is registered for the object.
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
#' @seealso [new_causal_wts()] for the class these generics have methods for,
#'   and the `propensity` package for concrete weight classes.
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
