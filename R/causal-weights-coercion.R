#' Common type of two causal weight vectors of the same class
#'
#' @description
#' `causal_wts_ptype2()` supplies the one coercion rule that concrete weight
#' classes share: two causal weight vectors combine only when their estimands
#' are identical, and otherwise the common type is a plain double. Concrete
#' classes call it from inside their own `vctrs::vec_ptype2()` method for the
#' pairing of two vectors of that class.
#'
#' @details
#' This is a helper rather than a method, and it cannot become one. vctrs
#' resolves `vec_ptype2()` with an exact lookup on `class(x)[[1]]` and never
#' walks the rest of the class vector, so a method registered on `causal_wts` is
#' never found for a concrete subclass. Each package therefore keeps its own
#' `vec_ptype2.<cls>.<cls>` method and calls this from inside it:
#'
#' ```r
#' vec_ptype2.bw.bw <- function(x, y, ...) {
#'   causal_wts_ptype2(
#'     x,
#'     y,
#'     new_bw(estimand = estimand(x)),
#'     on_downgrade = function() warn_bw_downgrade("bw")
#'   )
#' }
#' ```
#'
#' The estimands are read through [estimand()], so a class that computes its
#' estimand rather than storing it is compared on the value its method returns.
#' Two vectors that record no estimand are compatible, which is the ordinary
#' case for continuous-exposure weights rather than an edge case.
#'
#' `ptype` is evaluated only on the compatible path. A caller's prototype is
#' usually a live constructor call over metadata the two vectors do not agree
#' on, so the downgrade path never forces it.
#'
#' The condition belongs to the caller. This function signals nothing of its
#' own: `on_downgrade()` runs exactly once when the estimands differ, and
#' whatever it signals reaches the user unchanged. The concrete classes
#' downstream disagree on the class of that condition and on whether the
#' downgrade warns at all, so choosing one here would break them.
#'
#' The scope is one pairing and one rule, which is as much as the concrete
#' classes share. A class that checks nothing beyond the estimand hands over its
#' whole method, as `balancing` does. A class that checks further metadata can
#' use this for the estimand alone and keep its own control flow: `propensity`
#' compares five more fields after the estimand, each with its own message and
#' its own bail-out. Such a caller can pass a sentinel as `ptype` and carry on
#' with its remaining checks when that sentinel comes back.
#'
#' @param x,y The two weight vectors being combined, as passed to
#'   `vec_ptype2()`.
#' @param ptype The common type to return when the estimands agree. Returned
#'   exactly as supplied and never inspected, and evaluated only on that path.
#' @param on_downgrade A function of no arguments, called once when the
#'   estimands differ. Its value is discarded; it exists so the caller can
#'   signal a condition of its own choosing.
#'
#' @return `ptype` when `estimand(x)` and `estimand(y)` are identical, and
#'   `double()` otherwise.
#'
#' @seealso [new_causal_wts()] for the class this rule is written against, and
#'   [estimand()] for the accessor it compares.
#'
#' @export
#'
#' @examples
#' ate <- new_causal_wts(c(1.2, 0.8), subclass = "my_wts", estimand = "ate")
#' att <- new_causal_wts(c(1.0, 1.5), subclass = "my_wts", estimand = "att")
#'
#' # Matching estimands give back the prototype the caller built.
#' causal_wts_ptype2(
#'   ate,
#'   ate,
#'   new_causal_wts(subclass = "my_wts", estimand = estimand(ate)),
#'   on_downgrade = function() message("the estimands differ")
#' )
#'
#' # Differing estimands run the callback and drop to a bare double.
#' causal_wts_ptype2(
#'   ate,
#'   att,
#'   new_causal_wts(subclass = "my_wts", estimand = estimand(ate)),
#'   on_downgrade = function() message("the estimands differ")
#' )
causal_wts_ptype2 <- function(x, y, ptype, on_downgrade) {
  # The generic rather than `attr(x, "estimand")`, because a concrete class is
  # free to compute its estimand instead of storing one.
  if (!identical(estimand(x), estimand(y))) {
    on_downgrade()
    return(double())
  }

  # `ptype` is still a promise here, and everything above is written so that it
  # stays one on the downgrade path. Forcing it, even to name it in a message,
  # would evaluate a constructor call over metadata the two vectors just failed
  # to agree on.
  ptype
}
