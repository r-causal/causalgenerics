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
#' This package defines only the generic. Methods live in the packages that own
#' the relevant model classes, such as `propensity` for propensity score models
#' and `balancing` for balancing weight fits. Dispatch is on `ps_mod`, the
#' weighting object; each method documents the outcome model classes and further
#' arguments it accepts.
#'
#' @param ps_mod The weighting object that produced the weights, for example a
#'   fitted propensity score model. `ipw()` dispatches on this argument.
#' @param outcome_mod A fitted weighted outcome model.
#' @param ... Arguments passed to methods.
#'
#' @return A method-defined object holding causal effect estimates and their
#'   standard errors.
#'
#' @seealso The `propensity` and `balancing` packages for methods.
#'
#' @export
ipw <- function(ps_mod, outcome_mod, ...) {
  UseMethod("ipw")
}
