test_that("ipw() dispatches on the weighting object", {
  local_s3_method("ipw", "cg_weight", function(wt_mod, outcome_mod, ...) {
    "dispatched"
  })

  wt_mod <- structure(list(), class = "cg_weight")
  expect_identical(
    dispatch_from_baseenv(ipw, wt_mod, outcome_mod = NULL),
    "dispatched"
  )
})

test_that("ipw() errors for an object with no registered method", {
  # `ipw()` has no default method, so the failure comes from base R's dispatch.
  # Pin that message: a bare `expect_error()` here would also be satisfied by a
  # call that failed before it ever reached dispatch.
  x <- structure(list(), class = "cg_unregistered")
  expect_error(
    dispatch_from_baseenv(ipw, x, outcome_mod = NULL),
    "no applicable method for 'ipw'"
  )
})

test_that("ipw() reaches the same method through all four call forms", {
  local_s3_method("ipw", "cg_weight", function(wt_mod, outcome_mod, ...) {
    list(
      weighting = class(wt_mod),
      outcome = class(outcome_mod),
      dots = list(...)
    )
  })
  local_ps_mod_warning_reset()

  wt_mod <- structure(list(), class = "cg_weight")
  outcome_mod <- structure(list(), class = "cg_outcome")

  # The empty `dots` is the other half of what the deprecated forms have to
  # show: the old name is consumed by the repair rather than passed along to
  # the method.
  expected <- list(
    weighting = "cg_weight",
    outcome = "cg_outcome",
    dots = list()
  )

  expect_identical(
    dispatch_from_baseenv(ipw, wt_mod = wt_mod, outcome_mod = outcome_mod),
    expected
  )
  expect_identical(
    dispatch_from_baseenv(ipw, wt_mod, outcome_mod),
    expected
  )
  # `suppressWarnings()` rather than `expect_warning()` here, because
  # `expect_warning()` returns the condition rather than the value under test.
  # The warning has tests of its own below.
  expect_identical(
    suppressWarnings(dispatch_from_baseenv(
      ipw,
      ps_mod = wt_mod,
      outcome_mod = outcome_mod
    )),
    expected
  )
  expect_identical(
    suppressWarnings(dispatch_from_baseenv(
      ipw,
      outcome_mod = outcome_mod,
      ps_mod = wt_mod
    )),
    expected
  )
})

test_that("ipw() dispatches on `ps_mod` even when `outcome_mod` is named first", {
  # The regression the shim exists for, and the reason it cannot be written as a
  # method. `UseMethod()` falls back to the first argument of the call when the
  # argument it dispatches on is missing from that call, so without the
  # re-entry this call runs the method registered for the outcome model's class
  # and returns its answer with no error at all.
  local_s3_method("ipw", "cg_weight", function(wt_mod, outcome_mod, ...) {
    "weighting method"
  })
  local_s3_method("ipw", "cg_outcome", function(wt_mod, outcome_mod, ...) {
    "outcome method"
  })
  local_ps_mod_warning_reset()

  wt_mod <- structure(list(), class = "cg_weight")
  outcome_mod <- structure(list(), class = "cg_outcome")

  expect_identical(
    suppressWarnings(dispatch_from_baseenv(
      ipw,
      outcome_mod = outcome_mod,
      ps_mod = wt_mod
    )),
    "weighting method"
  )
})

test_that("ipw() errors when the weighting object is named both ways", {
  # The two spellings are one argument, so a call that uses both has said the
  # same thing twice and there is no reading of it that is safe to guess at.
  # Before this error the extra `ps_mod` rode through `...` to the method, which
  # ignored it, and nothing anywhere said the argument had been spelled twice.
  local_s3_method("ipw", "cg_weight", function(wt_mod, outcome_mod, ...) {
    "weighting method"
  })
  local_ps_mod_warning_reset()

  wt_mod <- structure(list(), class = "cg_weight")
  other <- structure(list(), class = "cg_weight")

  expect_error(
    dispatch_from_baseenv(ipw, wt_mod = wt_mod, ps_mod = other),
    class = "causalgenerics_invalid_argument_ps_mod"
  )
  # The general class as well, since the specific one alone would not show that
  # a handler for any argument this package rejects catches it.
  expect_error(
    dispatch_from_baseenv(ipw, wt_mod = wt_mod, ps_mod = other),
    class = "causalgenerics_invalid_argument"
  )
  # Both spellings are named, since the whole content of the error is that these
  # two names are the same argument. Snapshotted for the wording and asserted
  # here as well, since snapshots do not run on CRAN.
  expect_error(
    dispatch_from_baseenv(ipw, wt_mod = wt_mod, ps_mod = other),
    "`ps_mod`.+`wt_mod`"
  )
  expect_snapshot(error = TRUE, ipw(wt_mod = wt_mod, ps_mod = other))
})

test_that("ipw() matches a `ps_mod` call against the 0.1.0 formals", {
  # The forms that make the deprecation's promise true or false. Under the 0.1.0
  # formals `(ps_mod, outcome_mod, ...)`, an unnamed argument beside a named
  # `ps_mod` was the outcome model, and that call form is on CRAN. Read against
  # the current formals it fills `wt_mod` instead, where it takes over dispatch;
  # since outcome models are `glm`s and propensity registers its weighting
  # method for `glm`, that dispatch succeeds and answers with the wrong numbers.
  #
  # Methods are registered for the outcome and stray classes as well, so a call
  # that dispatched on either would return that method's answer rather than fail
  # for want of one.
  local_s3_method("ipw", "cg_weight", function(wt_mod, outcome_mod, ...) {
    list(
      weighting = class(wt_mod),
      outcome = if (missing(outcome_mod)) "<missing>" else class(outcome_mod),
      dots = list(...)
    )
  })
  local_s3_method("ipw", "cg_outcome", function(wt_mod, outcome_mod, ...) {
    "outcome method"
  })
  local_s3_method("ipw", "cg_stray", function(wt_mod, outcome_mod, ...) {
    "stray method"
  })
  local_ps_mod_warning_reset()

  wt_mod <- structure(list(), class = "cg_weight")
  outcome_mod <- structure(list(), class = "cg_outcome")
  stray <- structure(list(), class = "cg_stray")

  # `suppressWarnings()` rather than `expect_warning()`, which returns the
  # condition rather than the value under test. The warning has tests of its own
  # below.
  matched <- function(...) suppressWarnings(dispatch_from_baseenv(ipw, ...))

  # The unnamed argument is the outcome model, whichever side of `ps_mod` it
  # falls on.
  expect_identical(
    matched(ps_mod = wt_mod, outcome_mod),
    list(weighting = "cg_weight", outcome = "cg_outcome", dots = list())
  )
  expect_identical(
    matched(outcome_mod, ps_mod = wt_mod),
    list(weighting = "cg_weight", outcome = "cg_outcome", dots = list())
  )

  # Only the first unnamed argument is the outcome model. Once both formals are
  # taken, by name or by position, the rest go to `...` as they did in 0.1.0.
  expect_identical(
    matched(ps_mod = wt_mod, outcome_mod, stray),
    list(weighting = "cg_weight", outcome = "cg_outcome", dots = list(stray))
  )
  expect_identical(
    matched(ps_mod = wt_mod, outcome_mod = outcome_mod, stray),
    list(weighting = "cg_weight", outcome = "cg_outcome", dots = list(stray))
  )

  # And with nothing to fill it, `outcome_mod` is still missing rather than
  # filled from `...`.
  expect_identical(
    matched(ps_mod = wt_mod, conf_level = 0.9),
    list(
      weighting = "cg_weight",
      outcome = "<missing>",
      dots = list(conf_level = 0.9)
    )
  )
})

test_that("ipw() repairs a `ps_mod` call forwarded through `...`", {
  # A wrapper that forwards its dots reaches the generic as `ipw(...)`. Its one
  # argument name is `...`, so a repair that rewrote argument names in the call
  # would find nothing to rewrite and re-evaluate the same call forever. A
  # repair that spliced the forwarded arguments in with `match.call()` would
  # substitute their expressions and then evaluate them in the wrapper's frame,
  # where the names below do not exist. Re-entering through the 0.1.0 formals
  # leaves the promises where they are, which is what this asserts: the objects
  # arrive, and they arrive in the slots 0.1.0 put them in.
  local_s3_method("ipw", "cg_weight", function(wt_mod, outcome_mod, ...) {
    list(weighting = class(wt_mod), outcome = class(outcome_mod))
  })
  local_ps_mod_warning_reset()

  forward <- function(...) ipw(...)
  caller <- function() {
    caller_weighting <- structure(list(), class = "cg_weight")
    caller_outcome <- structure(list(), class = "cg_outcome")
    forward(ps_mod = caller_weighting, caller_outcome)
  }

  expect_identical(
    suppressWarnings(caller()),
    list(weighting = "cg_weight", outcome = "cg_outcome")
  )
})

test_that("ipw() errors when `ps_mod` is supplied twice", {
  # Two objects under the one name, and nothing in the call to say which of them
  # to weight by. 0.1.0 rejected it outright, since the second object matched a
  # formal the first had taken. Before this check the shim read the first,
  # dropped the second, and warned as though the call had been unambiguous.
  local_s3_method("ipw", "cg_weight", function(wt_mod, outcome_mod, ...) {
    "weighting method"
  })
  local_s3_method("ipw", "cg_other", function(wt_mod, outcome_mod, ...) {
    "other method"
  })
  local_ps_mod_warning_reset()

  wt_mod <- structure(list(), class = "cg_weight")
  other <- structure(list(), class = "cg_other")
  outcome_mod <- structure(list(), class = "cg_outcome")

  supplied_twice <- function() {
    # jarl-ignore duplicated_arguments: the duplicated name is the call under test
    dispatch_from_baseenv(
      ipw,
      ps_mod = wt_mod,
      ps_mod = other,
      outcome_mod = outcome_mod
    )
  }

  expect_error(
    supplied_twice(),
    class = "causalgenerics_invalid_argument_ps_mod"
  )
  # The general class as well, since the specific one alone would not show that
  # a handler for any argument this package rejects catches it.
  expect_error(supplied_twice(), class = "causalgenerics_invalid_argument")
  # Asserted here as well as snapshotted, since snapshots do not run on CRAN.
  expect_error(supplied_twice(), "`ps_mod` must be supplied at most once")
  expect_snapshot(
    error = TRUE,
    # jarl-ignore duplicated_arguments: the duplicated name is the call under test
    ipw(ps_mod = wt_mod, ps_mod = other, outcome_mod = outcome_mod)
  )
})

test_that("ipw() errors when `wt_mod` is spelled partially beside `ps_mod`", {
  # R matches a named argument partially against the formals that precede `...`,
  # so `wt` supplies `wt_mod` exactly as the full spelling does and the call is
  # the same ambiguity. A check that compared names for equality would miss it
  # and then hand two `wt_mod` arguments to the re-entered call.
  local_s3_method("ipw", "cg_weight", function(wt_mod, outcome_mod, ...) {
    "weighting method"
  })
  local_ps_mod_warning_reset()

  wt_mod <- structure(list(), class = "cg_weight")
  other <- structure(list(), class = "cg_weight")

  expect_error(
    dispatch_from_baseenv(ipw, wt = wt_mod, ps_mod = other),
    class = "causalgenerics_invalid_argument_ps_mod"
  )
})

test_that("ipw() dispatches once when it repairs the deprecated name", {
  # The re-entered call supplies `wt_mod`, so the shim cannot fire a second
  # time. Counting the method's runs is what shows that from the outside.
  runs <- 0L
  local_s3_method("ipw", "cg_weight", function(wt_mod, outcome_mod, ...) {
    runs <<- runs + 1L
    "dispatched"
  })
  local_ps_mod_warning_reset()

  wt_mod <- structure(list(), class = "cg_weight")

  expect_identical(
    suppressWarnings(dispatch_from_baseenv(
      ipw,
      ps_mod = wt_mod,
      outcome_mod = NULL
    )),
    "dispatched"
  )
  expect_identical(runs, 1L)
})

test_that("check_ipw_reentry() refuses an argument list that would recurse", {
  # The generic builds neither of these, which is the point: the guard is what
  # turns a mistake in that construction into an error rather than unbounded
  # recursion, so it is asserted on its own.
  expect_error(
    check_ipw_reentry(list(outcome_mod = NULL)),
    class = "causalgenerics_ipw_reentry"
  )
  expect_error(
    check_ipw_reentry(list(wt_mod = NULL, ps_mod = NULL)),
    class = "causalgenerics_internal_error"
  )
  expect_snapshot(error = TRUE, check_ipw_reentry(list(outcome_mod = NULL)))

  expect_identical(check_ipw_reentry(list(wt_mod = 1)), list(wt_mod = 1))
})

test_that("ipw() leaves a missing `outcome_mod` missing", {
  local_s3_method("ipw", "cg_weight", function(wt_mod, outcome_mod, ...) {
    missing(outcome_mod)
  })
  local_ps_mod_warning_reset()

  wt_mod <- structure(list(), class = "cg_weight")

  expect_true(suppressWarnings(dispatch_from_baseenv(ipw, ps_mod = wt_mod)))
  expect_false(suppressWarnings(dispatch_from_baseenv(
    ipw,
    ps_mod = wt_mod,
    outcome_mod = NULL
  )))
})

test_that("ipw() passes further arguments on unchanged", {
  local_s3_method("ipw", "cg_weight", function(wt_mod, outcome_mod, ...) {
    list(...)
  })
  local_ps_mod_warning_reset()

  wt_mod <- structure(list(), class = "cg_weight")

  expect_identical(
    suppressWarnings(dispatch_from_baseenv(
      ipw,
      ps_mod = wt_mod,
      outcome_mod = NULL,
      conf_level = 0.9,
      se_method = "linearization"
    )),
    list(conf_level = 0.9, se_method = "linearization")
  )
})

test_that("the deprecated `ps_mod` spelling warns", {
  local_s3_method("ipw", "cg_weight", function(wt_mod, outcome_mod, ...) {
    "dispatched"
  })
  wt_mod <- structure(list(), class = "cg_weight")

  local_ps_mod_warning_reset()
  expect_warning(
    ipw(ps_mod = wt_mod, outcome_mod = NULL),
    class = "causalgenerics_deprecated_ipw_ps_mod"
  )

  # The general class, on a fresh flag, since the specific one alone would not
  # show that a handler for any deprecation this package signals catches it.
  local_ps_mod_warning_reset()
  expect_warning(
    ipw(ps_mod = wt_mod, outcome_mod = NULL),
    class = "causalgenerics_deprecated"
  )

  # The message names the argument to use instead. Snapshotted for the wording
  # and asserted here as well, since snapshots do not run on CRAN and this is
  # the half that has to hold there.
  local_ps_mod_warning_reset()
  expect_warning(ipw(ps_mod = wt_mod, outcome_mod = NULL), "`wt_mod` instead")
})

test_that("the `ps_mod` deprecation warning reads as intended", {
  local_s3_method("ipw", "cg_weight", function(wt_mod, outcome_mod, ...) {
    "dispatched"
  })
  local_ps_mod_warning_reset()

  wt_mod <- structure(list(), class = "cg_weight")

  # `expect_snapshot_warning()` rather than `expect_snapshot()`, which lets the
  # warning reach the default handler and so fails outright under the
  # `warn = 2` the test gate sets. The wording this records is asserted for real
  # in the test above, since snapshots do not run on CRAN.
  expect_snapshot_warning(ipw(ps_mod = wt_mod, outcome_mod = NULL))
})

test_that("the `ps_mod` deprecation warns once per session", {
  local_s3_method("ipw", "cg_weight", function(wt_mod, outcome_mod, ...) {
    "dispatched"
  })
  local_ps_mod_warning_reset()

  wt_mod <- structure(list(), class = "cg_weight")

  expect_warning(
    ipw(ps_mod = wt_mod, outcome_mod = NULL),
    class = "causalgenerics_deprecated"
  )
  # Still repaired, just no longer announced.
  expect_no_warning(
    expect_identical(ipw(ps_mod = wt_mod, outcome_mod = NULL), "dispatched")
  )

  # And the flag is what holds it back, rather than anything about the call.
  local_ps_mod_warning_reset()
  expect_warning(
    ipw(ps_mod = wt_mod, outcome_mod = NULL),
    class = "causalgenerics_deprecated"
  )
})

test_that("ipw() does not warn when the weighting object is passed as `wt_mod`", {
  local_s3_method("ipw", "cg_weight", function(wt_mod, outcome_mod, ...) {
    "dispatched"
  })
  local_ps_mod_warning_reset()

  wt_mod <- structure(list(), class = "cg_weight")

  expect_no_warning(dispatch_from_baseenv(
    ipw,
    wt_mod = wt_mod,
    outcome_mod = NULL
  ))
  expect_no_warning(dispatch_from_baseenv(ipw, wt_mod, NULL))
})
