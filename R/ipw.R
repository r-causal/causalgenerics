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
#' The weighting argument was named `ps_mod` in propensity 0.1.0, where `ipw()`
#' was a plain function rather than a generic. That spelling is deprecated in
#' favor of `wt_mod`, which reads for a weighting object of any kind rather than
#' only a propensity score model. `ps_mod` is not an argument of the generic,
#' but a call that still names the argument that way warns once per session and
#' then behaves in every respect as though the object had been passed as
#' `wt_mod`, dispatch included. New code should use `wt_mod`.
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
  # The deprecated `ps_mod` spelling has to be repaired here, before dispatch,
  # and repairing it means calling `ipw()` again with the object genuinely
  # supplied as `wt_mod`. Assigning the object to `wt_mod` in this body and
  # falling through to `UseMethod()` does not work: when the argument dispatch
  # happens on is missing from the call, `UseMethod()` dispatches on the first
  # argument of the call whatever that argument is named, and never consults the
  # value assigned here. A caller who wrote `ipw(outcome_mod = o, ps_mod = w)`
  # would then silently get the method for the outcome model's class. A method
  # of its own cannot fix that either, since the re-entered call is what decides
  # which method runs.
  if (missing(wt_mod) && "ps_mod" %in% ...names()) {
    warn_ipw_ps_mod()

    dots <- list(...)
    # A logical index rather than `dots$ps_mod`, because `ps_mod = NULL` is an
    # argument the caller supplied and `$` cannot tell it from one they left
    # out. Indexing also removes every element the name matches, so nothing of
    # the old spelling reaches the method.
    from_ps_mod <- names(dots) == "ps_mod"
    args <- c(
      list(wt_mod = dots[from_ps_mod][[1L]]),
      # Left out rather than passed along when the caller omitted it, so that
      # the method still sees it as missing.
      if (!missing(outcome_mod)) list(outcome_mod = outcome_mod),
      dots[!from_ps_mod]
    )

    return(do.call("ipw", check_ipw_reentry(args)))
  }

  UseMethod("ipw")
}

# Check the argument list assembled for the re-entered `ipw()` call, and return
# it unchanged.
#
# The re-entry is single level by construction: the new call supplies `wt_mod`
# and carries no `ps_mod`, so `missing(wt_mod)` is `FALSE` the second time
# through and the shim cannot fire again. That invariant is checked rather than
# assumed, because breaking it would recurse without bound rather than give a
# wrong answer.
check_ipw_reentry <- function(args, call = sys.call(-1)) {
  arg_names <- names(args)

  if (!("wt_mod" %in% arg_names) || "ps_mod" %in% arg_names) {
    stop(errorCondition(
      "The re-entered `ipw()` call must supply `wt_mod` and drop `ps_mod`.",
      class = c(
        "causalgenerics_ipw_reentry",
        "causalgenerics_internal_error"
      ),
      call = call
    ))
  }

  args
}

# Package-local state. The one entry is the flag below.
deprecation_state <- new.env(parent = emptyenv())

# Warn that `ipw()` was called with the deprecated `ps_mod` spelling, once per
# session rather than once per call. Repeating it on every call would bury the
# output of a script that loops over models, which is the reason
# `lifecycle::deprecate_warn()` throttles the same way; a flag is all that
# behavior needs here, and this package takes no dependency to get it.
#
# The classes follow the shape the error helpers in `R/utils.R` use: one keyed
# to this deprecation, and one general class for a handler that wants any
# deprecation this package signals.
warn_ipw_ps_mod <- function(call = sys.call(-1)) {
  if (isTRUE(deprecation_state$ipw_ps_mod)) {
    return(invisible(NULL))
  }

  deprecation_state$ipw_ps_mod <- TRUE

  warning(warningCondition(
    "The `ps_mod` argument of `ipw()` is deprecated; name it `wt_mod` instead.",
    class = c(
      "causalgenerics_deprecated_ipw_ps_mod",
      "causalgenerics_deprecated"
    ),
    call = call
  ))

  invisible(NULL)
}
