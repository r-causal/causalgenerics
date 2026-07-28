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
#' Keeping that promise takes more than reading the object out from under the
#' old name, because the name a call gives one argument decides where its other
#' arguments land. A call naming `ps_mod` is therefore matched against the
#' 0.1.0 formals before it is routed into the current ones, so every argument
#' around it stays where 0.1.0 put it. In `ipw(ps_mod = model, fit)`, `fit` is
#' the outcome model, as it was in 0.1.0, rather than the weighting object the
#' current formals would read it as.
#'
#' Naming the weighting argument both ways in one call is an error rather than a
#' deprecation: the two names are one argument, and a call that spells it twice
#' has not said which object the method should use. Naming `ps_mod` twice is an
#' error for the same reason, and was one in 0.1.0 as well.
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
  # `ps_mod` is not an argument of the generic, so a call that names it lands in
  # `...`. Reading the names rather than forcing `list(...)` keeps a call with
  # nothing to repair from paying to evaluate its dots.
  from_ps_mod <- ...names() == "ps_mod"

  if (any(from_ps_mod)) {
    # Two objects under the one name, and nothing in the call to say which of
    # them the method should weight by. 0.1.0 rejected this too, where the
    # second object matched a formal the first had already taken.
    if (sum(from_ps_mod) > 1L) {
      stop_invalid_argument("ps_mod", "be supplied at most once")
    }

    # The deprecated spelling has to be repaired before dispatch, and repairing
    # it means re-entering `ipw()` with the object genuinely supplied as
    # `wt_mod`. Assigning the object to `wt_mod` in this body and falling
    # through to `UseMethod()` does not work: when the argument dispatch happens
    # on is missing from the call, `UseMethod()` dispatches on the first
    # argument of the call whatever that argument is named, and never consults
    # the value assigned here. A caller who wrote `ipw(outcome_mod = o,
    # ps_mod = w)` would then silently get the method for the outcome model's
    # class. A method of its own cannot fix that either, since the re-entered
    # call is what decides which method runs.
    #
    # Re-entry goes through a function carrying the 0.1.0 formals, and it is the
    # call itself that is handed to it rather than an argument list rebuilt from
    # this frame. Swapping only the function part leaves every argument with the
    # position, the name, and the promise it already has, so matching it against
    # those formals puts each one where 0.1.0 put it. That is what makes the
    # promise above true for a call such as `ipw(ps_mod = w, o)`, whose unnamed
    # argument is the outcome model and not, as the current formals would read
    # it, the weighting object.
    #
    # Rewriting the argument names in the call instead does not work. A call
    # that forwards its own dots arrives here as `ipw(...)`, whose one argument
    # name is `...`: there is nothing to rewrite, and re-evaluating it recurses
    # without bound. Splicing the dots in with `match.call()` does not work
    # either, since that substitutes the argument expressions and loses the
    # environments they have to be evaluated in.
    call <- sys.call()
    call[[1L]] <- ipw_ps_mod_signature

    return(eval(call, parent.frame()))
  }

  UseMethod("ipw")
}

# The 0.1.0 signature of `ipw()`, kept so that a call using the deprecated
# `ps_mod` spelling can be matched against the formals it was written for and
# then routed into the current ones. Nothing calls this directly; `ipw()`
# reaches it by swapping it into the function part of the call it was given.
ipw_ps_mod_signature <- function(ps_mod, outcome_mod, ...) {
  # The call as the caller wrote it, with the function part put back. The swap
  # above is this package's doing, and no condition should report it as though
  # the caller had made it.
  call <- sys.call()
  call[[1L]] <- quote(ipw)

  # `wt_mod` is not a formal here, so a call that also spells the weighting
  # argument that way leaves it in `...`, which is what makes it visible.
  # `pmatch()` rather than an equality test, because R matches a named argument
  # partially against the formals that precede `...`: `ipw(wt = w, ps_mod = m)`
  # supplies `wt_mod` as surely as writing the name out does, and the call is
  # ambiguous either way.
  if (!all(is.na(pmatch(...names(), "wt_mod", duplicates.ok = TRUE)))) {
    stop_invalid_argument(
      "ps_mod",
      "not be supplied together with `wt_mod`, which names the same argument",
      call = call
    )
  }

  warn_ipw_ps_mod(call)

  args <- c(
    list(wt_mod = ps_mod),
    # Left out rather than passed along when the caller omitted it, so that the
    # method still sees it as missing.
    if (!missing(outcome_mod)) list(outcome_mod = outcome_mod),
    list(...)
  )

  do.call("ipw", check_ipw_reentry(args, call))
}

# Check the argument list assembled for the re-entered `ipw()` call, and return
# it unchanged.
#
# The re-entry is single level by construction: the list always names `wt_mod`,
# and it can never carry `ps_mod`, which matched a formal of
# `ipw_ps_mod_signature()` and so never reached its dots. `missing(wt_mod)` is
# therefore `FALSE` the second time through and the shim cannot fire again. That
# invariant is checked rather than assumed, because breaking it would recurse
# without bound rather than give a wrong answer.
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
#
# `call` is passed rather than taken from `sys.call(-1)`, because the only
# caller runs under the 0.1.0 formals that `ipw()` swapped in, and its own call
# carries that substituted function part.
warn_ipw_ps_mod <- function(call) {
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
