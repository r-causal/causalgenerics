# `ess()` gains a real default here: the Kish effective sample size,
# `sum(w)^2 / sum(w^2)`. That formula is correct for any numeric weight vector,
# so the default is the whole implementation and no concrete weight class needs
# a method of its own. That is the point of putting it here. A `psw`
# and a `bw` both walk past their own class to `ess.default()` and compute
# correctly, which is what lets halfmoon stop masking the generic.
#
# Two facts about dispatch shape what follows. A weight vector's class is
# `c("psw", "causal_wts", "vctrs_vctr", "double")`, so `UseMethod()` walks that
# vector and then `default`; it never appends an implicit `"numeric"`, and an
# `ess.numeric()` method would therefore never fire for one. And `is.numeric()`
# is `TRUE` for a weight vector, so the input check has to be `!is.numeric(x)`
# rather than anything keyed on the class.

# A weight vector with nothing registered for its own class, and in particular
# no `vec_arith()` method, which the abstract layer deliberately does not
# supply. That is what a concrete class is allowed to look like partway through
# being built, and it is also the instrument several of these tests rely on:
# `wts^2` raises `vctrs_error_incompatible_op` on this fixture, so a value
# coming back at all proves the default reached the underlying double before
# squaring rather than squaring the weight vector itself.
cg_probe_wts <- function(x, ...) {
  new_causal_wts(x, subclass = "cg_probe_wts", ...)
}

test_that("ess() computes the Kish effective sample size for a bare double", {
  # Equal weights carry the same information as an unweighted sample, so the
  # effective size is the actual size; unequal weights make it smaller. Both
  # fixtures are dyadic, which is what lets `expect_identical()` compare against
  # a literal: `rep(1.2, 5)` gives `5 + 8.9e-16` rather than `5`.
  expect_identical(dispatch_from_baseenv(ess, rep(2, 5)), 5)
  expect_identical(dispatch_from_baseenv(ess, c(1, 1, 1, 1, 2)), 4.5)
})

test_that("ess() agrees across a double, an integer, and a weight vector", {
  # The three-way agreement is the contract that lets every concrete weight
  # class skip writing a method. It also pins that the input check is
  # `!is.numeric(x)` rather than a class test: `is.numeric()` is `TRUE` for a
  # weight vector, asserted below, so a single predicate admits all three of
  # these, while a check keyed on `causal_wts` would reject the first two and a
  # check keyed on `double` would reject the third.
  #
  # The weights are deliberately unequal, so the answer is not the length of the
  # input and an implementation that returned `length(x)` cannot pass.
  wts <- cg_probe_wts(c(1, 1, 1, 1, 2))
  expect_true(is.numeric(wts))

  expect_identical(dispatch_from_baseenv(ess, c(1, 1, 1, 1, 2)), 4.5)
  expect_identical(dispatch_from_baseenv(ess, c(1L, 1L, 1L, 1L, 2L)), 4.5)
  expect_identical(dispatch_from_baseenv(ess, wts), 4.5)
})

test_that("ess() returns a bare double for a weight vector", {
  # `sum()` on a `causal_wts` vector unwraps through the inherited `Summary`
  # method, so the quotient is an ordinary number. A result that still carried
  # the weight class, or the metadata that came in with it, would be a weight
  # vector of length one, which an effective sample size is not.
  result <- dispatch_from_baseenv(
    ess,
    cg_probe_wts(rep(2, 5), estimand = "ate")
  )

  expect_type(result, "double")
  expect_null(attributes(result))
})

test_that("ess() drops missing values when na.rm passes through the dots", {
  # The generic is `ess(x, ...)`, so `na.rm` reaches the default through the
  # dots rather than through a formal of the generic.
  expect_identical(dispatch_from_baseenv(ess, c(2, 2, NA), na.rm = TRUE), 2)
  expect_identical(
    dispatch_from_baseenv(ess, cg_probe_wts(c(2, 2, NA)), na.rm = TRUE),
    2
  )
})

test_that("ess() takes na.rm as an argument of its own, not out of the dots", {
  # Passing by name works under either signature, so nothing above separates
  # `ess.default(x, na.rm = FALSE, ...)` from a version that reaches into the
  # dots for `na.rm`. A positional call does separate them: under the second
  # shape `TRUE` stays in the dots, is never read, and the caller gets `NA`
  # back with no indication that the argument was ignored. The version being
  # moved up out of halfmoon is `ess(wts, na.rm = FALSE)`, so its callers can
  # already write this, and the move has to keep working for them.
  expect_identical(dispatch_from_baseenv(ess, c(2, 2, NA), TRUE), 2)
})

test_that("ess() keeps missing values by default", {
  # Both sums are missing, so the quotient is missing too. Dropping the missing
  # weights unasked would report an effective sample size computed from the
  # weights that happen to be present and say nothing about the ones that are
  # not, so the default has to be the conservative answer.
  expect_identical(dispatch_from_baseenv(ess, c(2, 2, NA)), NA_real_)
  expect_identical(
    dispatch_from_baseenv(ess, cg_probe_wts(c(2, 2, NA))),
    NA_real_
  )
})

test_that("a registered method takes precedence over the default", {
  # The fixture is a classed double, so `is.numeric()` is `TRUE` for it and the
  # default would return a Kish value rather than refuse it. The method returns
  # a sentinel instead, which is what separates the two: an assertion written
  # against the Kish value would be satisfied by whichever of them ran. The
  # second half computes what the default would have given, so the value the
  # first half must not be is named rather than implied.
  local_s3_method("ess", "cg_weight", function(x, ...) {
    "from the method"
  })

  wts <- structure(rep(2, 5), class = "cg_weight")

  expect_identical(dispatch_from_baseenv(ess, wts), "from the method")
  expect_identical(dispatch_from_baseenv(ess, unclass(wts)), 5)
})

test_that("ess() passes further arguments through to a registered method", {
  # The default is not the only thing the dots have to reach. Downstream methods
  # take arguments of their own through the same route, so this asserts the
  # generic's signature against a method rather than against the default.
  local_s3_method("ess", "cg_weight", function(x, na.rm = FALSE, ...) {
    if (na.rm) "dropped" else "kept"
  })

  wts <- structure(c(1, 1, NA), class = "cg_weight")

  expect_identical(dispatch_from_baseenv(ess, wts, na.rm = TRUE), "dropped")
  expect_identical(dispatch_from_baseenv(ess, wts), "kept")
})

test_that("ess() rejects input that is not numeric", {
  expect_error(
    dispatch_from_baseenv(ess, c("a", "b")),
    class = "causalgenerics_invalid_argument_x"
  )
  expect_error(
    dispatch_from_baseenv(ess, list(1, 2)),
    class = "causalgenerics_invalid_argument_x"
  )
  expect_error(
    dispatch_from_baseenv(ess, data.frame(wts = c(1, 2))),
    class = "causalgenerics_invalid_argument_x"
  )
})

test_that("ess() rejects input that is numeric underneath but not numeric", {
  # `is.numeric()` is the whole check, and it is `FALSE` for both of these even
  # though a factor stores integer level codes and a logical coerces to zero and
  # one. An implementation that unwrapped the payload before checking it would
  # return a Kish value computed from level codes, which is not an effective
  # sample size of anything.
  expect_error(
    dispatch_from_baseenv(ess, factor(c("a", "b"))),
    class = "causalgenerics_invalid_argument_x"
  )
  expect_error(
    dispatch_from_baseenv(ess, c(TRUE, TRUE)),
    class = "causalgenerics_invalid_argument_x"
  )
})

test_that("ess() rejects NULL because NULL is not numeric", {
  # A deliberate change from the version being moved up out of halfmoon, which
  # unwraps `NULL` to `NULL` and then evaluates `0 / 0`. A caller who passed
  # nothing at all gets an error naming the argument instead of a `NaN` that
  # travels on into whatever it is compared against.
  #
  # This is the type check firing, and nothing more. It is not a policy against
  # `NaN`: `NULL` is refused for the same reason a character vector is, and
  # numeric input that produces `NaN` honestly still gets `NaN` back, which the
  # next block pins.
  expect_error(
    dispatch_from_baseenv(ess, NULL),
    class = "causalgenerics_invalid_argument_x"
  )
})

test_that("ess() returns NaN for numeric input whose quotient is 0 / 0", {
  # The two numeric inputs the Kish formula has no answer for: no weights at
  # all, and weights that are all zero. Each makes both the sum and the sum of
  # squares zero, and each returns `NaN` from the version being moved up out of
  # halfmoon, whose own suite asserts the all-zero case directly.
  #
  # Refusing them the way `NULL` is refused would read as a natural extension
  # of the check above, and it would break that suite with nothing here to
  # catch it, so the answer is pinned rather than left to follow from the
  # `NULL` rule. It also settles what an empty weight vector returns.
  expect_identical(dispatch_from_baseenv(ess, rep(0, 5)), NaN)
  expect_identical(dispatch_from_baseenv(ess, numeric(0)), NaN)
})

test_that("the non-numeric error carries the general condition class", {
  # The two-class shape that `stop_no_method()` establishes and
  # `new_causal_wts()` reuses: a specific class for handlers that want this one
  # error, and a general one for callers that want any argument this package
  # rejected.
  expect_error(
    dispatch_from_baseenv(ess, c("a", "b")),
    class = "causalgenerics_invalid_argument"
  )
})

test_that("the non-numeric error states the contract", {
  # Two inputs rather than one, and they are expected to record the same
  # message. That sameness is the claim: `NULL` earns no wording of its own,
  # because it reaches the error through the same type check as any other
  # non-numeric input.
  expect_snapshot(error = TRUE, ess(c("a", "b")))
  expect_snapshot(error = TRUE, ess(NULL))
})
