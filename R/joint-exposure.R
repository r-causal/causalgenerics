#' Declare an exposure that crosses two treatments
#'
#' @description
#' `joint_exposure()` crosses two discrete treatments into one categorical
#' exposure and records the declaration on the result: which two treatments were
#' crossed, in what order, which levels each of them takes, and which cell the
#' effects are reported against. The package that builds the weights, the
#' package that checks balance over the cells, and the package that reports the
#' effects all read that one declaration, so they agree on which cells exist
#' without any of them owning it.
#'
#' @details
#' Conditioning on a variable and jointly intervening on two treatments are
#' different causal questions, and this class serves the second. Effect
#' modification asks how the effect of one treatment differs across the levels
#' of something else, which need not be a treatment and need not be anything an
#' analysis could set; what comes back is one effect of that one treatment per
#' subgroup, and a result reporting those carries the `group` column
#' [new_ipw()] documents. Interaction asks what would happen under an
#' intervention that sets both treatments at once; what comes back is one
#' exposure with a cell per combination, reported against the reference cell.
#' The packages that implement [ipw()] methods develop the distinction, in the
#' weights they fit and in the effects they report. What is settled here is the
#' declaration all of them read.
#'
#' The result is a factor. Its class vector is
#' `c("joint_exposure", "factor", "vctrs_vctr", "integer")`, with `"factor"`
#' ahead of `"vctrs_vctr"` so that the formula machinery treats a joint exposure
#' natively: `model.matrix()` gives exactly what a plain factor over the same
#' cells gives, `assign` and `contrasts` attributes included.
#'
#' The cells are the crossing of the two components' levels, labelled
#' `"name1 = level1, name2 = level2"` and ordered with the first component
#' varying fastest, which is the order [interaction()] produces. The reference
#' cell crosses the two components' own reference levels, and that ordering is
#' what puts it first, where the contrasts a model builds are read against it. A
#' factor component's reference level is the first of its `levels()`; any other
#' component's is the smallest of its values, which makes `0` the reference of a
#' 0/1 treatment and `FALSE` the reference of a logical one.
#'
#' Every cell of the crossing must be populated. A cell nobody falls in is a
#' positivity violation rather than a small sample: the joint effect of the two
#' treatments is not identified in those data, and no weight, contrast, or
#' balance check over the crossing can be computed. Construction refuses it, as
#' it refuses a component that never varies and a component with a missing
#' value.
#'
#' The two components must also be named distinctly. The names are what tell the
#' treatments apart in the cells, so one name used twice would label every cell
#' with the same variable on either side of it and leave nothing downstream able
#' to tell which component a cell varies. Construction refuses that as well.
#'
#' The declaration survives every operation that keeps the full set of cells,
#' including `[`, `vctrs::vec_slice()`, `sort()`, and a round trip through a
#' data frame, and it survives combining two joint exposures that declare the
#' same crossing. Operations that narrow or rewrite the level set give it up and
#' say so: `droplevels()`, `x[i, drop = TRUE]`, and `levels<-` each warn and
#' return a plain factor. So does combining a joint exposure with a factor, with
#' a character vector, or with a joint exposure that declares a different
#' crossing. Casting the other way is refused outright, because a crossing
#' cannot be recovered from labels that merely look like one.
#'
#' Two ways of combining escape all of that and give a bare factor without a
#' word, because neither reaches a method this package can register.
#' `unlist()` combines the underlying codes in base C code without dispatching
#' at all, and `c()` dispatches on its first argument, so a combine that begins
#' with a plain factor never reaches this class. Reach for `vctrs::vec_c()`
#' where either could apply.
#'
#' @param ... Exactly two named vectors, one per treatment, recycled to a common
#'   length. Each may be a factor, a character vector, or a numeric or logical
#'   vector. The names are the treatment names the cell labels are written from,
#'   and the two must be distinct.
#'
#' @return A vector of class
#'   `c("joint_exposure", "factor", "vctrs_vctr", "integer")` with one level per
#'   cell of the crossing, carrying the component levels and the cell labels the
#'   accessors read back.
#'
#' @seealso [is_joint_exposure()], [joint_components()], and
#'   [joint_reference()], which read the declaration back off the vector.
#'
#' @importFrom vctrs new_vctr vec_data vec_recycle_common
#' @export
#'
#' @examples
#' x <- joint_exposure(
#'   qsmk = c(0, 1, 0, 1),
#'   exercise = c("no", "no", "yes", "yes")
#' )
#' x
#'
#' levels(x)
#' joint_components(x)
#' joint_reference(x)
#'
#' # The formula machinery sees the factor it is, and reports every cell
#' # against the reference one.
#' colnames(model.matrix(~x, data.frame(x = x)))
joint_exposure <- function(...) {
  components <- list(...)

  if (length(components) != 2L) {
    stop_joint_exposure_two_components(length(components))
  }

  component_names <- names(components)
  no_name <- is.null(component_names) ||
    anyNA(component_names) ||
    !all(nzchar(component_names))
  if (no_name) {
    stop_joint_exposure_unnamed_component()
  }

  # Settled before the crossing is built, so that a caller who passed one column
  # twice is sent to name the second treatment rather than told that the cells
  # the duplicate leaves empty violate positivity.
  if (identical(component_names[[1L]], component_names[[2L]])) {
    stop_joint_exposure_shared_name(component_names[[1L]])
  }

  # Recycling before the per-component checks so that the checks and the codes
  # below see the same lengths. Two exposures measured on the same units already
  # agree, so in practice this is where a caller learns that theirs do not.
  components <- do.call(vctrs::vec_recycle_common, components)

  # Indexed by position, which pairs each component with the name it arrived
  # under, so that a fault is reported against the component carrying it
  # whichever of the two positions that component is in.
  for (i in seq_along(components)) {
    component <- components[[i]]
    name <- component_names[[i]]
    # The levels as well as the values. A factor built by `addNA()` or by
    # `factor(exclude = NULL)` declares `NA` as a level, which gives its missing
    # observations an ordinary code and hides them from `anyNA()`. Those rows
    # still have no known exposure, and the cell the crossing would write for
    # them names an absence rather than a treatment; when nobody falls in that
    # cell, the check below would refuse the crossing for a positivity failure
    # and send the caller to fix a component that is not the problem.
    declares_na <- is.factor(component) && anyNA(levels(component))
    if (anyNA(component) || declares_na) {
      stop_joint_exposure_missing_value(name)
    }
    # Observed values rather than declared levels. A factor that declares two
    # levels and shows one is the same fault as a constant numeric, and a check
    # written against `levels()` alone would let it through and then refuse the
    # cells it leaves empty, which names the wrong thing.
    if (length(unique(component)) < 2L) {
      stop_joint_exposure_constant_component(name)
    }
  }

  crossed <- lapply(components, joint_component_factor)
  component_levels <- lapply(crossed, levels)
  cells <- joint_cell_labels(component_names, component_levels)
  codes <- joint_cell_codes(crossed)

  empty <- cells[tabulate(codes, nbins = length(cells)) == 0L]
  if (length(empty) > 0L) {
    stop_joint_exposure_empty_cell(empty)
  }

  new_joint_exposure(codes, cells, component_levels)
}

#' Read the declaration a joint exposure carries
#'
#' @description
#' These three functions read back what [joint_exposure()] recorded, so that a
#' package holding a column of exposures can act on the crossing without
#' rebuilding it.
#'
#' `is_joint_exposure()` reports whether an object declares a crossing.
#' `joint_components()` reports which two treatments were crossed and the levels
#' each of them takes. `joint_reference()` reports the cell the effects are
#' reported against.
#'
#' @details
#' They are plain functions rather than generics. The declaration lives in the
#' attributes of one class, and a package that meets a joint exposure needs the
#' same answer this package would give rather than a place to register a
#' different one.
#'
#' `joint_components()` gives the component names in the order the crossing was
#' declared in, and each component's levels in the order the crossing varies
#' them, which is the order that settles the cell labels. `joint_reference()`
#' gives the first level of the vector, which is the cell crossing the two
#' components' reference levels.
#'
#' @param x A joint exposure, as built by [joint_exposure()]. `is_joint_exposure()`
#'   accepts any object.
#'
#' @return `is_joint_exposure()` returns a single `TRUE` or `FALSE`.
#'   `joint_components()` returns a named list of two character vectors, one per
#'   component, named for the treatments. `joint_reference()` returns the
#'   reference cell's label as a single string.
#'
#' @seealso [joint_exposure()], which records what these read back.
#'
#' @name joint-exposure-accessors
#'
#' @examples
#' x <- joint_exposure(
#'   qsmk = c(0, 1, 0, 1),
#'   exercise = c("no", "no", "yes", "yes")
#' )
#'
#' is_joint_exposure(x)
#' joint_components(x)
#' joint_reference(x)
NULL

#' @rdname joint-exposure-accessors
#' @export
is_joint_exposure <- function(x) {
  inherits(x, "joint_exposure")
}

#' @rdname joint-exposure-accessors
#' @export
joint_components <- function(x) {
  check_joint_exposure(x)
  attr(x, "components")
}

#' @rdname joint-exposure-accessors
#' @export
joint_reference <- function(x) {
  check_joint_exposure(x)
  # Read off the levels rather than stored separately. The reference cell is the
  # first cell by construction, and a second copy of it would be a second thing
  # to keep true through every operation that rebuilds the vector.
  levels(x)[[1L]]
}

#' Construct a joint exposure from validated parts
#'
#' Assumes the crossing has already been checked, which is why the subsetting
#' and coercion methods can call it: they rebuild a vector whose declaration was
#' validated when it was first declared, and revalidating there would refuse a
#' slice that happens to leave a cell empty, which is exactly the case the class
#' is meant to carry through.
#'
#' @param codes An integer vector of cell indices into `cells`.
#' @param cells A character vector of cell labels, in the declared order.
#' @param components A named list of the components' levels, in the declared
#'   order.
#'
#' @return A vector of class
#'   `c("joint_exposure", "factor", "vctrs_vctr", "integer")`.
#'
#' @noRd
new_joint_exposure <- function(codes, cells, components) {
  vctrs::new_vctr(
    codes,
    levels = cells,
    components = components,
    class = c("joint_exposure", "factor"),
    inherit_base_type = TRUE
  )
}

#' Refuse an object that carries no crossing
#'
#' The accessors read attributes that only this class carries, so an object
#' without them would answer `NULL` and let a caller carry on as though it had
#' read a declaration.
#'
#' @param x The object an accessor was called on.
#' @param call The call to blame, which is the accessor's own.
#'
#' @return `NULL`, invisibly, when `x` is a joint exposure.
#'
#' @noRd
check_joint_exposure <- function(x, call = sys.call(-1)) {
  if (!is_joint_exposure(x)) {
    stop_invalid_argument("x", "be a joint exposure", call = call)
  }
  invisible(NULL)
}

#' Read a component as the factor whose levels the crossing varies
#'
#' A factor states its own reference level and its own level order, and a
#' crossing that re-sorted them would silently change which comparison the
#' analysis reports. Anything else is read the way `factor()` reads it, which
#' sorts the observed values and so makes `0` the reference of a 0/1 treatment.
#'
#' @param x One component, already checked for missing values.
#'
#' @return A factor.
#'
#' @noRd
joint_component_factor <- function(x) {
  if (is.factor(x)) {
    return(x)
  }
  factor(x)
}

#' Write the cell labels of a crossing
#'
#' The first component varies fastest, which is the order `interaction()`
#' produces with its default `lex.order = FALSE` and the order that puts the
#' cell crossing both reference levels first.
#'
#' @param names The two component names, in the declared order.
#' @param levels The two components' levels, in the declared order.
#'
#' @return A character vector with one label per cell.
#'
#' @noRd
joint_cell_labels <- function(names, levels) {
  first <- rep(levels[[1L]], times = length(levels[[2L]]))
  second <- rep(levels[[2L]], each = length(levels[[1L]]))
  paste0(names[[1L]], " = ", first, ", ", names[[2L]], " = ", second)
}

#' Index each observation into the cells of a crossing
#'
#' @param crossed The two components as factors, in the declared order.
#'
#' @return An integer vector of cell indices, one per observation.
#'
#' @noRd
joint_cell_codes <- function(crossed) {
  first <- as.integer(crossed[[1L]])
  second <- as.integer(crossed[[2L]])
  first + (second - 1L) * nlevels(crossed[[1L]])
}

#' Give up the declaration and keep the cells
#'
#' The vector a degrading operation hands back. It holds the same observations
#' labelled the same way and says nothing about the crossing, which is what
#' makes it a plain factor rather than a joint exposure with a field removed.
#'
#' @param x A joint exposure.
#'
#' @return A factor of class `"factor"` with the joint exposure's cells as its
#'   levels.
#'
#' @noRd
joint_bare_factor <- function(x) {
  factor(as.character(x), levels = levels(x))
}

# `levels.vctrs_vctr()` answers `NULL`, and `"vctrs_vctr"` is in the class vector
# ahead of anything base would reach for a factor. Everything that reads the
# cells goes through here, `print.factor()` and `table()` included, so without
# this method the class reports no levels at all.
#' @export
levels.joint_exposure <- function(x) {
  attr(x, "levels")
}

#' @export
`levels<-.joint_exposure` <- function(x, value) {
  warn_joint_exposure_replaced_levels()
  bare <- joint_bare_factor(x)
  levels(bare) <- value
  bare
}

# `"factor"` precedes `"vctrs_vctr"`, so `[.factor` would run first: it
# re-attaches the levels and the class and leaves `drop = TRUE` to `factor()`,
# which errors on a vctrs vector. The slice is built here instead, from the codes
# and the declaration the vector already carries.
#' @export
`[.joint_exposure` <- function(x, i, ..., drop = FALSE) {
  if (isTRUE(drop)) {
    warn_joint_exposure_dropped_levels()
    bare <- joint_bare_factor(x)
    return(droplevels(if (missing(i)) bare else bare[i]))
  }
  if (missing(i)) {
    return(x)
  }
  new_joint_exposure(vctrs::vec_data(x)[i], levels(x), joint_components(x))
}

#' @export
`[[.joint_exposure` <- function(x, i, ...) {
  new_joint_exposure(vctrs::vec_data(x)[[i]], levels(x), joint_components(x))
}

# `"factor"` precedes `"vctrs_vctr"`, so `c.factor()` would run first and combine
# the codes into a bare factor without a word, which is the one degradation that
# would happen silently. Routing through `vec_c()` puts `c()` on the same
# coercion rules every other combine takes: two vectors declaring the same
# crossing keep it, and anything else warns on its way down.
#' @export
c.joint_exposure <- function(...) {
  vctrs::vec_c(...)
}

#' @export
droplevels.joint_exposure <- function(x, ...) {
  warn_joint_exposure_dropped_levels()
  droplevels(joint_bare_factor(x), ...)
}

#' @export
format.joint_exposure <- function(x, ...) {
  format(as.character(x), ...)
}

# The header is what tells a reader that these are cells of a declared crossing
# rather than labels somebody pasted together, and it names the two treatments so
# that the crossing can be read without reading every level. The reference cell
# closes the display because it is the comparison every reported effect is
# against.
#' @export
print.joint_exposure <- function(x, ...) {
  cat(
    "<joint_exposure[",
    length(x),
    "]: ",
    toString(names(joint_components(x))),
    ">\n",
    sep = ""
  )
  if (length(x) > 0L) {
    print(format(x), quote = FALSE)
  }
  cat("Reference: ", joint_reference(x), "\n", sep = "")
  invisible(x)
}
