# Read-only methods for the abstract weight class.
#
# Every method here unwraps to the underlying double and hands the work to the
# base or vctrs implementation. None of them names a metadata field, so they are
# correct for any concrete class built on `causal_wts` regardless of what that
# class records. Operations that give back something that is still a weight
# vector delegate rather than reconstruct, which leaves the concrete class's own
# `vec_restore()` in charge of the metadata.

#' @importFrom vctrs vec_math vec_math_base vec_data vec_restore
#' @export
vec_math.causal_wts <- function(.fn, .x, ...) {
  # Some functions like cumsum/cumprod should preserve the weight class
  if (.fn %in% c("cumsum", "cumprod", "cummin", "cummax")) {
    result <- vctrs::vec_math_base(.fn, vctrs::vec_data(.x), ...)
    return(vctrs::vec_restore(result, .x))
  }
  # Other functions like log, sqrt return numeric
  vctrs::vec_math_base(.fn, vctrs::vec_data(.x), ...)
}

#' @export
Summary.causal_wts <- function(..., na.rm = FALSE) {
  .fn <- .Generic
  args <- list(...)
  numeric_args <- lapply(args, vctrs::vec_data)
  do.call(.fn, c(numeric_args, list(na.rm = na.rm)))
}

#' @export
min.causal_wts <- function(..., na.rm = FALSE) {
  args <- list(...)
  numeric_args <- lapply(args, vctrs::vec_data)
  do.call("min", c(numeric_args, list(na.rm = na.rm)))
}

#' @export
max.causal_wts <- function(..., na.rm = FALSE) {
  args <- list(...)
  numeric_args <- lapply(args, vctrs::vec_data)
  do.call("max", c(numeric_args, list(na.rm = na.rm)))
}

#' @export
range.causal_wts <- function(..., na.rm = FALSE) {
  args <- list(...)
  numeric_args <- lapply(args, vctrs::vec_data)
  do.call("range", c(numeric_args, list(na.rm = na.rm)))
}

# median() and quantile() are not members of the Summary group generic, so the
# abstract class needs its own methods to reach them. Both summarize the
# underlying double and return a plain numeric.
#' @importFrom stats median
#' @export
median.causal_wts <- function(x, na.rm = FALSE, ...) {
  stats::median(vctrs::vec_data(x), na.rm = na.rm, ...)
}

#' @importFrom stats quantile
#' @export
quantile.causal_wts <- function(
  x,
  probs = seq(0, 1, 0.25),
  na.rm = FALSE,
  ...
) {
  stats::quantile(vctrs::vec_data(x), probs = probs, na.rm = na.rm, ...)
}

#' @export
summary.causal_wts <- function(object, ...) {
  # `vec_data()` rather than `as.numeric()`, which routes through `vec_cast()`
  # and so errors for a concrete class that supplies no cast to double.
  summary(vctrs::vec_data(object), ...)
}

#' @export
anyDuplicated.causal_wts <- function(x, incomparables = FALSE, ...) {
  # `anyDuplicated.vctrs_vctr()` answers `TRUE` or `FALSE`; callers expect the
  # index of the first duplicate that the base default gives for a double.
  anyDuplicated(vctrs::vec_data(x), incomparables = incomparables, ...)
}

#' @export
diff.causal_wts <- function(x, lag = 1L, differences = 1L, ...) {
  diff(vctrs::vec_data(x), lag = lag, differences = differences, ...)
}

#' @export
`[.causal_wts` <- function(x, i, ...) {
  # Bare `x[]` must return the whole weight vector unchanged; reading `i` here
  # would force the missing argument and error.
  if (missing(i)) {
    return(NextMethod())
  }
  if (is.matrix(i) || is.array(i)) {
    return(vctrs::vec_data(x)[i, ...])
  }
  NextMethod()
}

#' Comparison operators on `causal_wts` short-circuit `vec_equal()`/`vec_compare()`.
#'
#' Without these methods, `==.vctrs_vctr` and friends route through
#' `vec_equal()` -> `vec_cast_common()` -> the concrete class's
#' `vec_ptype2()` method, which fires a class downgrade warning once per call.
#' `glm.fit()` evaluates `weights == 0` and `weights > 0` repeatedly inside
#' `profile.glm()`, so a single `tidy(glm, conf.int = TRUE)` call can emit 100+
#' identical warnings. Comparing weights returns a logical vector, so no class is
#' actually being downgraded from the user's perspective; the warning is spurious
#' here.
#'
#' Strict vctrs size semantics are preserved: `vec_recycle_common()` enforces the
#' same N-or-1 size rule that `vec_equal()` would, so length-mismatched
#' comparisons error rather than silently recycling per base R rules. Combine
#' paths (`vec_c()`, `vec_cast()`) still go through the warning-emitting
#' `vec_ptype2` methods, so the user is still informed when the weight class
#' really is dropped.
#'
#' `call` is what the size error blames. Left to its own default,
#' `vec_recycle_common()` blames its caller, which is this helper, and the user
#' is told an unexported name they never wrote. `sys.call(-1)` records the
#' operator method instead, matching how `stop_no_method()` and
#' `stop_invalid_argument()` take a call in `R/utils.R`. An environment would be
#' the more usual thing to hand rlang, but not here: rlang treats a dispatch
#' frame as uninformative and walks past it, so the caller's environment names
#' whatever called the comparison rather than the comparison.
#' @importFrom vctrs vec_recycle_common
#' @noRd
causal_wts_compare <- function(e1, e2, call = sys.call(-1)) {
  args <- vctrs::vec_recycle_common(e1, e2, .call = call)
  if (inherits(args[[1]], "causal_wts")) {
    args[[1]] <- vctrs::vec_data(args[[1]])
  }
  if (inherits(args[[2]], "causal_wts")) {
    args[[2]] <- vctrs::vec_data(args[[2]])
  }
  list(e1 = args[[1]], e2 = args[[2]])
}

#' @export
`==.causal_wts` <- function(e1, e2) {
  args <- causal_wts_compare(e1, e2)
  args$e1 == args$e2
}

#' @export
`!=.causal_wts` <- function(e1, e2) {
  args <- causal_wts_compare(e1, e2)
  args$e1 != args$e2
}

#' @export
`<.causal_wts` <- function(e1, e2) {
  args <- causal_wts_compare(e1, e2)
  args$e1 < args$e2
}

#' @export
`>.causal_wts` <- function(e1, e2) {
  args <- causal_wts_compare(e1, e2)
  args$e1 > args$e2
}

#' @export
`<=.causal_wts` <- function(e1, e2) {
  args <- causal_wts_compare(e1, e2)
  args$e1 <= args$e2
}

#' @export
`>=.causal_wts` <- function(e1, e2) {
  args <- causal_wts_compare(e1, e2)
  args$e1 >= args$e2
}
