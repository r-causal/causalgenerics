#' Construct a causal weight vector
#'
#' @description
#' `new_causal_wts()` is the low-level constructor for the abstract weight class
#' that the r-causal ecosystem shares. It is intended for package developers
#' building a concrete weight class, not for end users, and it assumes `x` is
#' already a double vector.
#'
#' @details
#' Weights built here are vctrs vectors over a double, so the class vector is
#' `c(subclass, "causal_wts", "vctrs_vctr", "double")`. Callers may rely on that
#' shape and on its order: `subclass` first, so that methods for the concrete
#' class take precedence, then `causal_wts`, so that the concrete class inherits
#' every method registered against the abstract one. `subclass` must therefore be
#' a single non-empty string. Anything else would change the length of the class
#' vector or name a class that no method can be registered against.
#'
#' The dots carry whatever metadata the concrete class records, such as the
#' estimand or the exposure groups. This constructor neither names nor types
#' those fields. A field passed as `NULL` leaves its attribute absent rather than
#' recording a `NULL`, which is what lets a concrete constructor default its
#' metadata to `NULL` and pass it straight through.
#'
#' The abstract layer reaches well past construction. A concrete class inherits
#' the metadata accessors [is_causal_wt()], [estimand()], and [estimand<-()],
#' and it inherits every read-only operation that names no metadata field and so
#' gives the answer the underlying double would: `vec_math()`, the `Summary`
#' group generic, `min()`, `max()`, `range()`, `median()`, `quantile()`,
#' `summary()`, `anyDuplicated()`, `diff()`, `[`, and the six comparison
#' operators. [ess()] arrives as well, through a default that computes the Kish
#' effective sample size for any numeric vector.
#'
#' Do not write any of those again for a concrete class. A method registered on
#' the subclass takes precedence over the shared one, and the comparison
#' operators in particular do more than delegate: they short-circuit
#' `vec_equal()` and `vec_compare()`, which would otherwise route an expression
#' such as `weights > 0` through the concrete class's own `vec_ptype2()` method
#' and signal whatever that method signals on a downgrade, once per call.
#' `glm.fit()` evaluates comparisons of that shape repeatedly inside
#' `profile.glm()`, so a single profiled confidence interval on a weighted fit
#' would emit the same warning a hundred times over.
#'
#' What a concrete class does still own is its `vec_ptype2()`, `vec_cast()`,
#' `vec_ptype_abbr()`, `vec_ptype_full()`, `vec_restore()`, and `vec_arith()`
#' methods, for three distinct reasons.
#'
#' `vec_ptype2()` and `vec_cast()` cannot be inherited at all. vctrs resolves
#' them through `s3_method_specific()`, which keys on `class(x)[[1]]` and never
#' walks the rest of the class vector, so a `causal_wts` method is never found.
#'
#' `vec_ptype_abbr()` and `vec_ptype_full()` are likewise not reached from the
#' abstract class. A concrete class that omits either one falls back to the vctrs
#' default and prints a bare class name where the estimand should appear.
#'
#' `vec_restore()` and `vec_arith()` do dispatch through `causal_wts`, but they
#' still have to be written for the concrete class, because they reconcile
#' metadata that this package knows nothing about. An abstract `vec_restore()`
#' that copies attributes wholesale errors on named weights and re-attaches
#' index-typed metadata at the wrong length after slicing, and an abstract
#' `vec_arith()` would silently legalize arithmetic between two different
#' concrete weight classes.
#'
#' @param x A double vector of weights.
#' @param subclass A single non-empty string naming the concrete weight class.
#' @param ... Metadata to store as attributes on the result. A `NULL` value
#'   leaves the attribute absent.
#'
#' @return A vector of class `c(subclass, "causal_wts", "vctrs_vctr", "double")`.
#'
#' @seealso [estimand()] and [is_causal_wt()], the accessors this class supplies
#'   methods for, and [causal_wts_ptype2()], the coercion rule a concrete class
#'   calls from its own `vec_ptype2()` method.
#'
#' @importFrom vctrs new_vctr vec_assert
#' @export
#'
#' @examples
#' wts <- new_causal_wts(c(1.2, 0.8), subclass = "my_wts", estimand = "ate")
#' class(wts)
#' estimand(wts)
new_causal_wts <- function(x = double(), subclass, ...) {
  # Checked before `x` so that the cheap scalar check runs first. A type and
  # length check on its own would admit `NA_character_` and `""`, neither of
  # which names a class anything can dispatch on.
  valid_subclass <- !missing(subclass) &&
    is.character(subclass) &&
    length(subclass) == 1L &&
    !is.na(subclass) &&
    nzchar(subclass)

  if (!valid_subclass) {
    stop_invalid_argument("subclass", "be a single non-empty string")
  }

  vctrs::vec_assert(x, ptype = double())

  vctrs::new_vctr(
    x,
    ...,
    class = c(subclass, "causal_wts"),
    inherit_base_type = TRUE
  )
}

#' @export
is_causal_wt.causal_wts <- function(x, ...) {
  TRUE
}

#' @export
estimand.causal_wts <- function(x, ...) {
  attr(x, "estimand")
}

#' @export
`estimand<-.causal_wts` <- function(x, ..., value) {
  attr(x, "estimand") <- value
  x
}
