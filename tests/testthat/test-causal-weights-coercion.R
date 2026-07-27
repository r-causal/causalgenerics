# `causal_wts_ptype2()` is not a method and never will be. vctrs resolves
# `vec_ptype2()` through `s3_method_specific()`, which keys on `class(x)[[1]]`
# alone, so `vec_ptype2.causal_wts.causal_wts` is never found for a concrete
# subclass. Each package therefore keeps its own
# `vec_ptype2.<cls>.<cls>` method and calls this helper from inside it. Every
# assertion below calls the helper directly for that reason; there is no
# dispatch to isolate.
#
# The helper covers the pairing of two vectors of the same concrete class, which
# is the only pairing whose rule is shared. `vec_ptype2.<cls>.double` and the
# character and integer pairings differ between the two packages down to whether
# they warn at all, and none of them has an estimand on both sides to compare.
#
# What is shared is one rule: two causal weight vectors combine only when their
# estimands are identical, and otherwise the common type is a plain double. The
# rest of the behavior belongs to the caller, and most of this file exists to pin
# that boundary. The caller supplies the prototype for the compatible path, the
# caller decides whether the downgrade warns, and the caller owns the class of
# whatever condition it signals. propensity warns
# `propensity_coercion_warning` with a per-field message and balancing warns
# `balancing_class_downgrade_warning` with one generic message; both classes are
# pinned by tests in those packages. A helper that signalled anything of its own,
# or that unified those classes, would break them.

cg_coercion_wts <- function(estimand = NULL, x = double()) {
  new_causal_wts(x, subclass = "cg_coercion_wts", estimand = estimand)
}

# The prototype the caller hands to the compatible path. It carries metadata
# beyond the estimand so that returning it and rebuilding an equivalent vector
# around the same estimand are distinguishable: `expect_identical()` compares the
# whole attribute set, and a helper that reconstructed the prototype would lose
# the fields it cannot name.
cg_coercion_ptype <- function(estimand = NULL) {
  new_causal_wts(
    subclass = "cg_coercion_wts",
    estimand = estimand,
    stabilized = TRUE,
    groups = list(exposed = 1L, unexposed = 2L)
  )
}

# A downgrade callback that records how many times it ran.
#
# `expect_warning()` proves a warning was signalled at least once, which is not
# the claim. A caller sizes its own warning output on the promise that one
# combine produces one call, and only a count can hold the helper to that.
#
# The callback also returns a value. An implementation written as
# `return(on_downgrade())` produces the right warning and the right call count
# and still hands the caller the wrong answer, so every assertion that uses this
# checks the result as well.
new_call_counter <- function(action = function() NULL) {
  state <- new.env(parent = emptyenv())
  state$n <- 0L
  state$callback <- function() {
    state$n <- state$n + 1L
    action()
    "the callback's return value"
  }
  state
}

# Signal a warning of `class` and nothing else, standing in for a downstream
# package's own classed warning helper.
warn_probe <- function(class) {
  function() {
    warning(warningCondition("the estimands differ", class = class))
  }
}

# Every condition `expr` signals, in order, along with its value.
#
# `expect_warning(class = )` proves that a condition of that class was signalled
# but says nothing about what else was, and what else was is the claim this file
# is built around: the helper must add nothing on top of what the caller's
# callback signals. Warnings and messages are muffled so that a stray one fails
# an assertion here rather than leaking into the test reporter.
capture_conditions <- function(expr) {
  conditions <- list()
  result <- withCallingHandlers(
    expr,
    condition = function(cnd) {
      conditions[[length(conditions) + 1L]] <<- cnd
      if (inherits(cnd, "warning")) {
        invokeRestart("muffleWarning")
      }
      if (inherits(cnd, "message")) {
        invokeRestart("muffleMessage")
      }
    }
  )
  list(result = result, conditions = conditions)
}

test_that("identical estimands return the prototype the caller supplied", {
  counter <- new_call_counter()
  ptype <- cg_coercion_ptype("ate")

  expect_identical(
    causal_wts_ptype2(
      cg_coercion_wts("ate"),
      cg_coercion_wts("ate"),
      ptype,
      on_downgrade = counter$callback
    ),
    ptype
  )
  expect_identical(counter$n, 0L)
})

test_that("two unset estimands are compatible", {
  # `estimand()` returns `NULL` when nothing was recorded, and
  # `identical(NULL, NULL)` is `TRUE`, so two vectors that never named an
  # estimand combine and keep their class. That is the ordinary case for
  # continuous-exposure weights rather than an edge case, and a rule written as
  # "both estimands must be set and equal" would send it down the downgrade path
  # instead.
  counter <- new_call_counter()
  ptype <- cg_coercion_ptype()

  expect_identical(
    causal_wts_ptype2(
      cg_coercion_wts(),
      cg_coercion_wts(),
      ptype,
      on_downgrade = counter$callback
    ),
    ptype
  )
  expect_identical(counter$n, 0L)
})

test_that("differing estimands return a bare double", {
  # `expect_identical()` against `double()` pins the type, the zero length, and
  # the absence of any attribute in one assertion. A downgrade that handed back
  # a zero-length weight vector would still combine without error and would
  # quietly keep the class the downgrade exists to drop.
  expect_identical(
    causal_wts_ptype2(
      cg_coercion_wts("ate"),
      cg_coercion_wts("att"),
      cg_coercion_ptype("ate"),
      on_downgrade = function() NULL
    ),
    double()
  )
})

test_that("an estimand set on one side only is a downgrade", {
  # `NULL` against a string is the asymmetric case, and it has to fail in both
  # orders. `identical()` already answers this, but a rule written with
  # `is.null()` branches or with `==`, which returns `logical(0)` against `NULL`
  # and would never be `FALSE`, does not.
  expect_identical(
    causal_wts_ptype2(
      cg_coercion_wts(),
      cg_coercion_wts("ate"),
      cg_coercion_ptype("ate"),
      on_downgrade = function() NULL
    ),
    double()
  )
  expect_identical(
    causal_wts_ptype2(
      cg_coercion_wts("ate"),
      cg_coercion_wts(),
      cg_coercion_ptype("ate"),
      on_downgrade = function() NULL
    ),
    double()
  )
})

test_that("the downgrade callback runs exactly once", {
  counter <- new_call_counter()

  expect_identical(
    causal_wts_ptype2(
      cg_coercion_wts("ate"),
      cg_coercion_wts("att"),
      cg_coercion_ptype("ate"),
      on_downgrade = counter$callback
    ),
    double()
  )
  expect_identical(counter$n, 1L)
})

test_that("the prototype is never evaluated on the downgrade path", {
  # The prototype a caller supplies is a live constructor call, such as
  # `new_bw(estimand = estimand(x))`, and the downgrade path exists precisely
  # because that vector cannot be built. Lazy evaluation gives the helper this
  # for free, and any `force(ptype)` or use of the prototype in a message takes
  # it away again, so it is worth an assertion that cannot pass by accident: a
  # prototype that errors when touched.
  counter <- new_call_counter()

  expect_identical(
    causal_wts_ptype2(
      cg_coercion_wts("ate"),
      cg_coercion_wts("att"),
      stop("the prototype must not be evaluated"),
      on_downgrade = counter$callback
    ),
    double()
  )
  expect_identical(counter$n, 1L)
})

test_that("the helper signals nothing of its own on either path", {
  # A callback that does nothing leaves the helper as the only thing that could
  # signal, on the path where a downgrade is happening and a warning would be
  # most tempting to add.
  silent <- function() NULL

  expect_no_condition(
    causal_wts_ptype2(
      cg_coercion_wts("ate"),
      cg_coercion_wts("att"),
      cg_coercion_ptype("ate"),
      on_downgrade = silent
    )
  )
  expect_no_condition(
    causal_wts_ptype2(
      cg_coercion_wts("ate"),
      cg_coercion_wts("ate"),
      cg_coercion_ptype("ate"),
      on_downgrade = silent
    )
  )
})

test_that("the callback's condition reaches the caller unchanged", {
  # Two callbacks differing only in the class they signal, standing in for
  # propensity's `warn_incompatible_metadata()` and balancing's
  # `warn_bw_downgrade()`. The helper has to be agnostic to both: exactly one
  # condition comes out, it is the one the callback signalled, and the helper
  # neither wraps it, re-classes it, nor adds a second condition alongside it.
  first <- capture_conditions(causal_wts_ptype2(
    cg_coercion_wts("ate"),
    cg_coercion_wts("att"),
    cg_coercion_ptype("ate"),
    on_downgrade = warn_probe("cg_probe_downgrade_one")
  ))

  expect_identical(first$result, double())
  expect_length(first$conditions, 1L)
  expect_s3_class(
    first$conditions[[1]],
    c("cg_probe_downgrade_one", "warning", "condition"),
    exact = TRUE
  )
  expect_identical(
    conditionMessage(first$conditions[[1]]),
    "the estimands differ"
  )

  second <- capture_conditions(causal_wts_ptype2(
    cg_coercion_wts("ate"),
    cg_coercion_wts("att"),
    cg_coercion_ptype("ate"),
    on_downgrade = warn_probe("cg_probe_downgrade_two")
  ))

  expect_identical(second$result, double())
  expect_length(second$conditions, 1L)
  expect_s3_class(
    second$conditions[[1]],
    c("cg_probe_downgrade_two", "warning", "condition"),
    exact = TRUE
  )
})

test_that("the estimand is read through estimand() rather than the attribute", {
  # A concrete class is free to compute its estimand rather than store it, and
  # both classes downstream reach it through the generic. `attr(x, "estimand")`
  # agrees with the generic for every other fixture in this file, so a
  # registered method that disagrees with the stored attribute is the only thing
  # that separates the two implementations.
  local_s3_method("estimand", "cg_coercion_wts", function(x, ...) "ate")

  counter <- new_call_counter()
  ptype <- cg_coercion_ptype("ate")

  expect_identical(
    causal_wts_ptype2(
      cg_coercion_wts("att"),
      cg_coercion_wts("atu"),
      ptype,
      on_downgrade = counter$callback
    ),
    ptype
  )
  expect_identical(counter$n, 0L)
})
