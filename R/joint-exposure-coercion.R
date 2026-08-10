# Coercion rules for the joint exposure class.
#
# vctrs resolves `vec_ptype2()` and `vec_cast()` with an exact lookup on the
# first element of each class vector, so every pairing needs its own method and
# nothing here is reached by inheritance from `"factor"`. The rule they spell out
# is one rule: a combined vector carries a crossing only when both sides declare
# the same one, and otherwise the declaration is given up out loud. Casting into
# the class is refused rather than degraded, because a crossing cannot be
# recovered from labels; a factor whose labels merely look like one would record
# a declaration nobody made.

#' Report whether two joint exposures declare the same crossing
#'
#' Three ways they can disagree, one per clause: the components are named
#' differently, their levels are different values, or their levels are the same
#' values in a different order. Each one means the two vectors were drawn from
#' different declarations, and a combine that kept the class would hand
#' downstream code one declaration for data from two. The cells are compared as
#' well as the components, so that a vector built by hand rather than through
#' [joint_exposure()] cannot agree on the components and disagree on the labels.
#'
#' @param x,y Two joint exposures.
#'
#' @return `TRUE` when the two declarations agree, `FALSE` otherwise.
#'
#' @noRd
joint_declarations_agree <- function(x, y) {
  identical(joint_components(x), joint_components(y)) &&
    identical(levels(x), levels(y))
}

#' The factor the cells of two vectors degrade into
#'
#' Every cell either side declares stays a level, in the order the left side
#' declares them and then whatever the right side adds. A degraded vector reports
#' the same observations under the same labels; what it no longer says is which
#' two treatments those labels cross.
#'
#' @param x,y The two vectors being combined.
#'
#' @return A zero-length factor over the union of the two level sets.
#'
#' @noRd
joint_degraded_ptype <- function(x, y) {
  factor(levels = union(levels(x), levels(y)))
}

#' @importFrom vctrs vec_ptype2 vec_cast vec_ptype_abbr vec_ptype_full
#' @export
vec_ptype2.joint_exposure.joint_exposure <- function(x, y, ...) {
  if (joint_declarations_agree(x, y)) {
    return(new_joint_exposure(integer(), levels(x), joint_components(x)))
  }
  warn_joint_exposure_incompatible_metadata()
  joint_degraded_ptype(x, y)
}

#' @export
vec_ptype2.joint_exposure.factor <- function(x, y, ...) {
  warn_joint_exposure_foreign_type("factor")
  joint_degraded_ptype(x, y)
}

#' @export
vec_ptype2.factor.joint_exposure <- function(x, y, ...) {
  warn_joint_exposure_foreign_type("factor")
  joint_degraded_ptype(x, y)
}

#' @export
vec_ptype2.joint_exposure.character <- function(x, y, ...) {
  warn_joint_exposure_foreign_type("character")
  character()
}

#' @export
vec_ptype2.character.joint_exposure <- function(x, y, ...) {
  warn_joint_exposure_foreign_type("character")
  character()
}

#' @export
vec_cast.joint_exposure.joint_exposure <- function(
  x,
  to,
  ...,
  x_arg = "",
  to_arg = ""
) {
  # Reached from a combine only after `vec_ptype2()` has agreed the two
  # declarations match, so the codes already index the target's cells. A caller
  # who reaches it directly with two declarations that disagree is asking for a
  # recoding this class has no rule for.
  if (!joint_declarations_agree(x, to)) {
    vctrs::stop_incompatible_cast(x, to, x_arg = x_arg, to_arg = to_arg)
  }
  x
}

#' @export
vec_cast.factor.joint_exposure <- function(x, to, ...) {
  # A prototype with no levels of its own names no cells, and the ones the
  # crossing declares are the right answer: falling back to the observed values
  # would drop a cell that this vector happens not to contain.
  if (nlevels(to) == 0L) {
    return(joint_bare_factor(x))
  }
  vctrs::vec_cast(as.character(x), to)
}

#' @export
vec_cast.character.joint_exposure <- function(x, to, ...) {
  as.character(x)
}

# `table()` reaches the codes through `as.integer()`, which for a vctrs vector is
# this cast. Without it, counting the cells of a joint exposure errors.
#' @export
vec_cast.integer.joint_exposure <- function(x, to, ...) {
  vctrs::vec_data(x)
}

#' @export
vec_ptype_full.joint_exposure <- function(x, ...) {
  paste0("joint_exposure<", toString(names(joint_components(x))), ">")
}

#' @export
vec_ptype_abbr.joint_exposure <- function(x, ...) {
  "jnt_exp"
}
