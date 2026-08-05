#' The presentation mode of an inverse probability weighted result
#'
#' @description
#' The two readings of a result, defined for the class [new_ipw()] constructs.
#'
#' * `as_marginal()` returns the result reporting the causal contrast estimates.
#' * `as_conditional()` returns the result presenting the outcome model's
#'   coefficient surface.
#'
#' @details
#' Both surfaces exist on every result, so these generics record which one the
#' result presents rather than computing anything. They set the `effects` field
#' of the [new_ipw()] contract and read nothing else, which makes them the
#' supported way to move a result between the two readings: a caller writes
#' `as_conditional(res)` rather than assigning to the field.
#'
#' The methods on `ipw` are total. Every result has one of the two modes, so
#' asking for either is always answerable: they never error, asking twice says
#' what asking once said, and a result that goes out to the other reading and
#' back is the result that went in. A result built before the field existed
#' carries six fields rather than seven and reads as marginal, which is the mode
#' every method produced then.
#'
#' The generics live here for the reason `print()` does. Two packages each
#' registering `as_conditional.ipw()` would collide in the shared S3 method
#' table, and a caller writing against a result would then get whichever package
#' was installed last rather than the contract. There is no marginal or
#' conditional reading of an object that is not an IPW result, so the default
#' method signals an error rather than inventing one.
#'
#' @param x An `ipw` object. These generics dispatch on this argument.
#' @param ... Arguments passed to methods.
#'
#' @return `x` with its presentation mode set to the one asked for. The methods
#'   on `ipw` change the `effects` field and nothing else, so every other field
#'   comes back as it went in, the covariance attached to `estimates` included.
#'
#' @seealso [new_ipw()] for the result class and the field these generics set.
#'
#' @examples
#' dat <- data.frame(
#'   x = rep(c(-1.5, -0.5, 0.5, 1.5), each = 5),
#'   z = rep(c(0, 1), 10),
#'   y = rep(c(0, 1, 1, 0, 1), 4)
#' )
#'
#' # Written out literally, in the shape the `ipw()` return contract documents.
#' estimates <- data.frame(
#'   effect = "rd",
#'   estimate = 0.199882,
#'   std.err = 0.092425,
#'   z = 2.1626,
#'   ci.lower = 0.018732,
#'   ci.upper = 0.381032,
#'   conf.level = 0.95,
#'   p.value = 0.030570
#' )
#'
#' res <- new_ipw(
#'   estimand = "ate",
#'   wt_mod = glm(z ~ x, family = binomial(), data = dat),
#'   outcome_mod = glm(y ~ z, family = quasibinomial(), data = dat),
#'   estimates = estimates,
#'   se_method = "linearization",
#'   fit = NULL
#' )
#'
#' # A method that names no mode reports marginal effects.
#' res$effects
#'
#' as_conditional(res)$effects
#'
#' # Asking for the reading a result already reports gives the result back, and
#' # the two readings round-trip.
#' identical(as_marginal(res), res)
#' identical(as_marginal(as_conditional(res)), res)
#'
#' # Neither reading exists for an object that is not an IPW result.
#' try(as_marginal(1:3))
#'
#' @name ipw-modes
NULL

#' @rdname ipw-modes
#' @export
as_marginal <- function(x, ...) {
  UseMethod("as_marginal")
}

#' @export
as_marginal.ipw <- function(x, ...) {
  x$effects <- "marginal"
  x
}

#' @export
as_marginal.default <- function(x, ...) {
  stop_no_method("as_marginal", x)
}

#' @rdname ipw-modes
#' @export
as_conditional <- function(x, ...) {
  UseMethod("as_conditional")
}

#' @export
as_conditional.ipw <- function(x, ...) {
  x$effects <- "conditional"
  x
}

#' @export
as_conditional.default <- function(x, ...) {
  stop_no_method("as_conditional", x)
}

#' The presentation mode a result records
#'
#' The mode is the `effects` field of the [new_ipw()] contract, read through one
#' helper so that a result which does not carry it is read the same way
#' everywhere. Two shapes do not carry it: a result stored before the field
#' existed has six fields, and a list built with `effects = NULL` has seven, one
#' of them empty. Neither records a mode, and marginal is the reading every
#' method produced when neither did, so both are read as marginal rather than
#' refused.
#'
#' A field holding anything else was assigned to directly, which is outside the
#' contract the constructor and the mode generics keep. It is refused here
#' rather than passed on, so that a caller reading the mode gets one of the two
#' readings or an error and never a third thing to branch on.
#'
#' @param object An `ipw` object.
#' @param call The call to report an invalid stored mode against, which is the
#'   accessor's rather than this helper's.
#'
#' @return A single string, either `"marginal"` or `"conditional"`.
#'
#' @noRd
ipw_effects <- function(object, call = sys.call(-1)) {
  effects <- object$effects
  if (is.null(effects)) {
    return("marginal")
  }
  check_ipw_effects(effects, call = call)
  effects
}

#' The mode an accessor's `effects` argument asks a result for
#'
#' An accessor that takes an `effects` argument reports the mode the caller
#' names, and the mode the result records when the caller names none. `NULL` is
#' therefore how a caller declines to override, and any other value has to meet
#' the contract the constructor holds its own argument to.
#'
#' A named mode says which surface to report, so the stored field is not read at
#' all in that case. A result whose field was assigned to directly is therefore
#' still readable by naming a mode, and reading it without naming one is refused
#' where the field is read.
#'
#' @param object An `ipw` object.
#' @param effects The `effects` argument as the caller supplied it, or `NULL`.
#' @param call The call to report an invalid value against, which is the
#'   accessor's rather than this helper's.
#'
#' @return A single string, either `"marginal"` or `"conditional"`.
#'
#' @noRd
resolve_ipw_effects <- function(object, effects, call = sys.call(-1)) {
  if (is.null(effects)) {
    return(ipw_effects(object, call = call))
  }
  check_ipw_effects(effects, call = call)
  effects
}

#' Refuse an `effects` value that is not a mode
#'
#' The field has two readings and no third, and it holds the string itself, so
#' anything that is not one of the two exactly is refused. `NULL` is refused
#' along with the rest: it is not a single string, and a result that stored it
#' would record no mode at all. The callers that do treat `NULL` as an answer,
#' such as `resolve_ipw_effects()`, deal with it before they get here.
#'
#' @param effects The `effects` argument as the caller supplied it.
#' @param call The call to report the error against, which is the constructor's
#'   or the accessor's rather than this helper's.
#'
#' @return `effects`, invisibly, when it is a mode.
#'
#' @noRd
check_ipw_effects <- function(effects, call = sys.call(-1)) {
  # `%in%` rather than `==` so that a missing string is refused along with the
  # rest rather than making the test itself missing.
  is_mode <- is.character(effects) &&
    length(effects) == 1L &&
    effects %in% c("marginal", "conditional")

  if (!is_mode) {
    stop_invalid_argument(
      "effects",
      'be a single string, either "marginal" or "conditional"',
      call = call
    )
  }

  invisible(effects)
}
