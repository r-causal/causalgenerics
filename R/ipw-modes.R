#' The presentation mode of an inverse probability weighted result
#'
#' @description
#' The two readings of a result, defined for the class [new_ipw()] constructs
#' and for the pooled results [pool_ipw()] returns.
#'
#' * `as_marginal()` returns the result reporting the causal contrast estimates.
#' * `as_conditional()` returns the result presenting the outcome model's
#'   coefficient surface.
#'
#' @details
#' Both surfaces exist on an `ipw` result that supports both readings, so on one
#' of those these generics record which one the result presents rather than
#' computing anything. They set the `effects` field of the [new_ipw()] contract
#' and read the `readings` field beside it, which makes them the supported way
#' to move a result between the two readings: a caller writes
#' `as_conditional(res)` rather than assigning to the field, and a result that
#' does not support the reading asked for is refused rather than left recording
#' one it has no surface for.
#'
#' The methods on `ipw` hold for a result that supports both readings, which is
#' what a result records unless the package that built it said otherwise. Asking
#' such a result for either reading is always answerable: the methods never
#' error on one, asking twice says what asking once said, and a result that goes
#' out to the other reading and back is the result that went in.
#'
#' A result that supports one reading is refused the other rather than moved to
#' it. The `readings` field says which ones a result can present, and a fitting
#' package records one of them when the other has no meaning for the analysis it
#' ran. Asking for a reading outside that set raises an error of class
#' `causalgenerics_unsupported_reading_marginal` or
#' `causalgenerics_unsupported_reading_conditional`, and of the general class
#' `causalgenerics_unsupported_reading`, which carries the reading asked for
#' under `effects` and the set the result records under `readings`. The reading
#' such a result does support is the no-op it is on a result carrying both.
#'
#' A result built before the fields existed carries six fields rather than
#' eight. It reads as marginal, which is the mode every method produced then,
#' and as supporting both readings, which is what every result was assumed to
#' support when none of them recorded otherwise.
#'
#' The generics live here for the reason `print()` does. Two packages each
#' registering `as_conditional.ipw()` would collide in the shared S3 method
#' table, and a caller writing against a result would then get whichever package
#' was installed last rather than the contract. There is no marginal or
#' conditional reading of an object that is not an IPW result, so the default
#' method signals an error rather than inventing one.
#'
#' # The readings a pooled result carries
#'
#' [pool_ipw()] pools both readings of one set of results whenever it can
#' compute both, and stores the one the call did not name whole, under the
#' `alternate` component. The methods on `ipw_pooled` therefore swap which
#' reading the result presents rather than setting a field: the pooled
#' estimates, the pooling diagnostics, and the recorded mode move together, and
#' the components that describe the pooled analyses rather than a reading of
#' them stay where they are.
#'
#' Where both readings were pooled, the properties above hold. Asking a pooled
#' result for the reading it already presents gives that result back, and a
#' result taken out to the other reading and back is the result that went in.
#' Totality is what these methods cannot keep. A reading the pooling could not
#' compute is recorded as unavailable rather than computed, and asking for it
#' raises an error of class `causalgenerics_pool_missing_surface_marginal` or
#' `causalgenerics_pool_missing_surface_conditional`, and of the general class
#' `causalgenerics_pool_missing_surface`. That condition carries the reading
#' under `effects` and, under `reason`, the refusal that reading raised when it
#' was pooled, which is the same wording the caller would have seen from
#' [pool_ipw()] had they asked for it there.
#'
#' A result pooled before both readings were kept records no alternate at all.
#' Asking such a result for the other reading is refused with the same two
#' classes and a `reason` of `NULL`, and the message says that pooling the
#' results again gives a result carrying both.
#'
#' @param x An `ipw` or `ipw_pooled` object. These generics dispatch on this
#'   argument.
#' @param ... Arguments passed to methods.
#'
#' @return `x` presenting the reading asked for. The methods on `ipw` change the
#'   `effects` field and nothing else, so every other field comes back as it
#'   went in, the covariance attached to `estimates` and the set of readings the
#'   result supports included. A result that does not support the reading asked
#'   for raises an error rather than returning one. The methods on
#'   `ipw_pooled` exchange `estimates`, `pooling`, and the recorded mode with the
#'   reading stored under `alternate`, so the reading that was presented is what
#'   the returned result stores there, and the components shared by both readings
#'   come back as they went in. A pooled result that does not carry the reading
#'   asked for raises an error rather than returning one.
#'
#' @seealso [new_ipw()] for the result class and the field these generics set,
#'   and [pool_ipw()] for the pooled result and the reading it stores beside the
#'   one it presents.
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
#' # A result that supports one reading is refused the other rather than moved
#' # to it.
#' marginal_only <- new_ipw(
#'   estimand = "ate",
#'   wt_mod = glm(z ~ x, family = binomial(), data = dat),
#'   outcome_mod = glm(y ~ z, family = quasibinomial(), data = dat),
#'   estimates = estimates,
#'   se_method = "linearization",
#'   fit = NULL,
#'   readings = "marginal"
#' )
#'
#' try(as_conditional(marginal_only))
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
  # The set the result records is read before the field is written, so that a
  # result which cannot present this reading is refused rather than left
  # recording a mode it has no surface for.
  check_ipw_reading(x, "marginal")
  x$effects <- "marginal"
  x
}

#' @export
as_marginal.ipw_pooled <- function(x, ...) {
  flip_ipw_pooled(x, "marginal")
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
  check_ipw_reading(x, "conditional")
  x$effects <- "conditional"
  x
}

#' @export
as_conditional.ipw_pooled <- function(x, ...) {
  flip_ipw_pooled(x, "conditional")
}

#' @export
as_conditional.default <- function(x, ...) {
  stop_no_method("as_conditional", x)
}

#' Move a pooled result between the two readings it carries
#'
#' [pool_ipw()] pools both readings of one set of results and stores the one the
#' call did not name under `alternate`, so a pooled result is moved between them
#' by swapping the three fields that differ rather than by setting one. The
#' shared fields are properties of the analyses rather than of a reading of
#' them, and they stay where they are.
#'
#' The swap is by assignment to fields the result already has, which is what
#' makes the round trip the object that went in: the field order is the
#' contract, and the reading that was active goes back under `alternate` in the
#' shape it was recorded in. Asking for the reading the result already presents
#' is answered with the result itself, so the generics are idempotent on a
#' pooled result and on one stored before both readings were kept alike.
#'
#' Unlike the methods on `ipw`, these cannot be total. A reading that could not
#' be pooled is recorded as unavailable, and a result stored before both
#' readings were kept records nothing at all; both are refused here rather than
#' answered with a surface the result does not hold. The recorded reason is the
#' refusal the same request raised at pooling time, so it is passed on as it
#' stands.
#'
#' @param x An `ipw_pooled` object.
#' @param effects The reading to present.
#' @param call The call to report a missing reading against, which is the
#'   generic's or the accessor's rather than this helper's.
#'
#' @return `x` presenting the reading asked for.
#'
#' @noRd
flip_ipw_pooled <- function(x, effects, call = sys.call(-1)) {
  current <- ipw_effects(x, call = call)
  if (current == effects) {
    return(x)
  }

  alternate <- x$alternate
  if (is.null(alternate)) {
    stop_pool_missing_surface(effects, reason = NULL, call = call)
  }
  if (is.null(alternate$estimates)) {
    stop_pool_missing_surface(effects, alternate$reason, call = call)
  }

  # Read out before either side is written, since the two halves of a swap
  # cannot both be assigned first.
  stashed <- list(
    effects = current,
    estimates = x$estimates,
    pooling = x$pooling
  )

  x$effects <- alternate$effects
  x$estimates <- alternate$estimates
  x$pooling <- alternate$pooling
  x$alternate <- stashed
  x
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
#' A named mode is checked against the readings the result supports, and one
#' outside them is refused rather than reported from a surface the result has no
#' reading of. A mode read off the result needs no such check: the constructor
#' refuses a result whose mode is not one of its readings, so the stored mode is
#' one of them for every result a constructor built, and a result whose fields
#' were assigned to directly is outside the contract either way.
#'
#' That last part is about resolving the mode, which is all this helper does. It
#' holds for a caller that goes on to report the resolved mode directly, as the
#' accessors on `ipw` do. A pooled result is reported by moving it to the reading
#' asked for, and the move has to read the stored mode to know whether that
#' reading is the one the result already presents, so naming a mode there does
#' not keep the stored field from being read.
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
  check_ipw_reading(object, effects, call = call)
  effects
}

#' The readings a result supports
#'
#' The set is the `readings` field of the [new_ipw()] contract, read through one
#' helper so that a result which does not carry it is read the same way
#' everywhere. Three shapes do not carry it: a result stored before either field
#' existed has six fields, one stored between them has seven, and a list built
#' with `readings = NULL` has eight, one of them empty. None of them records a
#' set, and both readings were assumed to exist on every result when none of
#' them recorded otherwise, so all three are read as supporting both rather than
#' refused.
#'
#' A field holding anything else was assigned to directly, which is outside the
#' contract the constructor keeps, and it is refused here rather than passed on.
#' The mode the field has to include is not checked here, since the field being
#' read says nothing about which reading the result presents and the constructor
#' has already refused a result whose mode is outside its readings.
#'
#' @param object An `ipw` object.
#' @param call The call to report an invalid stored set against, which is the
#'   accessor's rather than this helper's.
#'
#' @return A character vector of one or both of `"marginal"` and
#'   `"conditional"`.
#'
#' @noRd
ipw_readings <- function(object, call = sys.call(-1)) {
  readings <- object$readings
  if (is.null(readings)) {
    return(c("marginal", "conditional"))
  }
  check_ipw_readings(readings, call = call)
  readings
}

#' Refuse a reading the result does not support
#'
#' Three surfaces ask the same question of a result: the mode generics, which
#' would otherwise record a reading the result cannot present, and the
#' accessors and `as.data.frame()`, which would otherwise report one. The
#' question is asked here once so that all of them refuse in the same words and
#' with the same classes.
#'
#' @param object An `ipw` object.
#' @param effects The reading that was asked for.
#' @param call The call to report the refusal against, which is the generic's or
#'   the method's rather than this helper's.
#'
#' @return `effects`, invisibly, when the result supports it.
#'
#' @noRd
check_ipw_reading <- function(object, effects, call = sys.call(-1)) {
  readings <- ipw_readings(object, call = call)
  if (!effects %in% readings) {
    stop_unsupported_reading(effects, readings, call = call)
  }

  invisible(effects)
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

#' Refuse a set of readings that is not one
#'
#' There are two readings and no third, so a set is one or both of them. The
#' empty set is refused along with the rest, since a result supporting no
#' reading reports nothing at all, and so is `NULL`: that is how a result stored
#' before the field existed records no set, and it is read as both wherever it
#' is read, so naming it here would ask the constructor to write that shape
#' deliberately. A reading named twice says nothing a set naming it once does
#' not, and it would leave the field with more entries than there are readings.
#'
#' The set is checked against the mode as well when there is a mode to check it
#' against. A result cannot record a reading it does not support: the mode says
#' which surface the result presents, and a result presenting one it has no
#' reading of is a state every method downstream would have to have an answer
#' for. `effects` is `NULL` where a stored set is read back rather than
#' supplied, which is the one case with no mode in hand to compare it to.
#'
#' @param readings The `readings` argument as the caller supplied it, or the
#'   field as a result stores it.
#' @param effects The mode the result records, or `NULL` when there is none to
#'   check the set against.
#' @param call The call to report the error against, which is the constructor's
#'   or the accessor's rather than this helper's.
#'
#' @return `readings`, invisibly, when it is a set of readings.
#'
#' @noRd
check_ipw_readings <- function(readings, effects = NULL, call = sys.call(-1)) {
  # `%in%` rather than `==` so that a missing string is refused along with the
  # rest rather than making the test itself missing.
  # Membership and distinctness bound the length between them, so a set of the
  # right kind cannot be longer than there are readings.
  is_set <- is.character(readings) &&
    length(readings) >= 1L &&
    all(readings %in% c("marginal", "conditional")) &&
    anyDuplicated(readings) == 0L

  if (!is_set) {
    stop_invalid_argument(
      "readings",
      paste0(
        'be a character vector of one or both of "marginal" and ',
        '"conditional", each named at most once'
      ),
      call = call
    )
  }

  if (!is.null(effects) && !effects %in% readings) {
    stop_invalid_argument(
      "readings",
      paste0(
        "include the reading the result records, since a result cannot record ",
        "a reading it does not support, and this one records the ",
        encodeString(effects, quote = '"'),
        " reading"
      ),
      call = call
    )
  }

  invisible(readings)
}
