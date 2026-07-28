#' Inverse probability weighted estimation
#'
#' @description
#' `ipw()` is the generic for bring-your-own-model inverse probability weighted
#' estimation of causal effects. A method takes a fitted weighting or propensity
#' score model together with a fitted weighted outcome model and returns causal
#' effect estimates with standard errors that account for the two-step
#' estimation process.
#'
#' @details
#' This package defines the generic and the shared result class its methods
#' return. The methods themselves live in the packages that own the relevant
#' model classes, such as `propensity` for propensity score models and
#' `balancing` for balancing weight fits. Dispatch is on `wt_mod`, the weighting
#' object; each method documents the outcome model classes and further arguments
#' it accepts.
#'
#' A method builds its return value with [new_ipw()] rather than a result object
#' of its own. The field names and their order are a cross-package contract, so
#' that an IPW estimate reads the same way whichever package produced it, and
#' constructing through [new_ipw()] is also what gives a method the shared
#' `print()` and `as.data.frame()` methods.
#'
#' @param wt_mod The weighting object that produced the weights, for example a
#'   fitted propensity score model. `ipw()` dispatches on this argument.
#' @param outcome_mod A fitted weighted outcome model.
#' @param ... Arguments passed to methods.
#'
#' @return An object of class `ipw`, holding the causal effect estimates and
#'   their standard errors alongside the models they came from. See [new_ipw()]
#'   for the components and their order.
#'
#' @seealso [new_ipw()], the constructor every method returns through, and the
#'   `propensity` and `balancing` packages for methods.
#'
#' @export
ipw <- function(wt_mod, outcome_mod, ...) {
  UseMethod("ipw")
}
