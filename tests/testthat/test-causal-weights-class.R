# The class vector `new_causal_wts()` builds is a cross-package contract, not an
# implementation detail. Both packages that build on it pin their own vector
# with `exact = TRUE`, so any change to the order or the length of what this
# constructor produces breaks them at once.
#
# The accessor assertions call the generics through `dispatch_from_baseenv()` so
# that they exercise the `S3method()` registration that downstream packages
# reach these methods through, rather than anything visible from the test frame.

test_that("new_causal_wts() builds the documented class vector", {
  wts <- new_causal_wts(c(0.5, 1, 1.5), subclass = "cg_weight")

  expect_s3_class(
    wts,
    c("cg_weight", "causal_wts", "vctrs_vctr", "double"),
    exact = TRUE
  )

  # `as.vector()` drops the class and the metadata and leaves the double
  # payload. `as.double()` would route through `vec_cast()`, which needs a
  # method for the concrete class that the abstract layer deliberately does not
  # supply, so it errors here rather than returning the payload.
  expect_identical(as.vector(wts), c(0.5, 1, 1.5))
})

test_that("new_causal_wts() defaults to a zero-length double", {
  wts <- new_causal_wts(subclass = "cg_weight")

  expect_s3_class(
    wts,
    c("cg_weight", "causal_wts", "vctrs_vctr", "double"),
    exact = TRUE
  )
  expect_length(wts, 0L)
  expect_identical(as.vector(wts), double())
})

test_that("new_causal_wts() stores propensity-shaped metadata from the dots", {
  # The dots carry whatever the concrete class needs. This is the full set of
  # fields a propensity score weight vector records.
  wts <- new_causal_wts(
    c(1.2, 0.8),
    subclass = "cg_weight",
    estimand = "ate",
    stabilized = TRUE,
    trimmed = FALSE,
    truncated = FALSE,
    calibrated = FALSE,
    stabilization_score = 0.94
  )

  expect_identical(attr(wts, "estimand"), "ate")
  expect_identical(attr(wts, "stabilized"), TRUE)
  expect_identical(attr(wts, "trimmed"), FALSE)
  expect_identical(attr(wts, "truncated"), FALSE)
  expect_identical(attr(wts, "calibrated"), FALSE)
  expect_identical(attr(wts, "stabilization_score"), 0.94)
  expect_s3_class(
    wts,
    c("cg_weight", "causal_wts", "vctrs_vctr", "double"),
    exact = TRUE
  )
})

test_that("new_causal_wts() stores balancing-shaped metadata from the dots", {
  # A balancing weight vector records the exposure groups as a named list of
  # index vectors alongside the estimand. Nothing about the metadata is typed or
  # named by this constructor, so an attribute holding indices has to survive
  # unchanged.
  groups <- list(untreated = c(1L, 3L), treated = c(2L, 4L))
  wts <- new_causal_wts(
    c(1, 1, 1, 1),
    subclass = "cg_weight",
    estimand = "att",
    groups = groups
  )

  expect_identical(attr(wts, "estimand"), "att")
  expect_identical(attr(wts, "groups"), groups)
  expect_s3_class(
    wts,
    c("cg_weight", "causal_wts", "vctrs_vctr", "double"),
    exact = TRUE
  )
})

test_that("estimand() reads the estimand a causal_wts vector records", {
  wts <- new_causal_wts(c(1.2, 0.8), subclass = "cg_weight", estimand = "ate")

  expect_identical(dispatch_from_baseenv(estimand, wts), "ate")
})

test_that("estimand() returns NULL when no estimand was recorded", {
  # Downstream `vec_ptype_abbr()` methods branch on the absent case, so it has
  # to be `NULL` rather than an error or an empty string.
  wts <- new_causal_wts(c(1.2, 0.8), subclass = "cg_weight")

  expect_null(dispatch_from_baseenv(estimand, wts))
})

test_that("estimand<-() writes the estimand attribute", {
  # The replacement form `estimand(wts) <- "att"` looks `estimand<-()` up where
  # it is written, which `baseenv()` cannot offer. Calling the replacement
  # generic by value keeps the assertion isolated from the test frame; the
  # surface syntax it expands to is base R's, not this package's.
  wts <- new_causal_wts(c(1.2, 0.8), subclass = "cg_weight", estimand = "ate")
  wts <- dispatch_from_baseenv(`estimand<-`, wts, value = "att")

  expect_identical(attr(wts, "estimand"), "att")
  expect_s3_class(
    wts,
    c("cg_weight", "causal_wts", "vctrs_vctr", "double"),
    exact = TRUE
  )
})

test_that("estimand<-() records an estimand on a vector that had none", {
  wts <- new_causal_wts(c(1.2, 0.8), subclass = "cg_weight")
  wts <- dispatch_from_baseenv(`estimand<-`, wts, value = "ato")

  expect_identical(dispatch_from_baseenv(estimand, wts), "ato")
})

test_that("is_causal_wt() is TRUE for a causal_wts vector", {
  wts <- new_causal_wts(c(1.2, 0.8), subclass = "cg_weight")

  expect_true(dispatch_from_baseenv(is_causal_wt, wts))
})

test_that("is_causal_wt() is FALSE for objects outside the class", {
  expect_false(dispatch_from_baseenv(is_causal_wt, c(1.2, 0.8)))
  expect_false(dispatch_from_baseenv(is_causal_wt, list(1.2, 0.8)))
})

test_that("new_causal_wts() rejects an x that is not a double", {
  # The message belongs to vctrs and can change with it, so only the condition
  # class is asserted and no snapshot is taken.
  expect_error(
    new_causal_wts(1:3, subclass = "cg_weight"),
    class = "vctrs_error_assert_ptype"
  )
  expect_error(
    new_causal_wts(c("a", "b"), subclass = "cg_weight"),
    class = "vctrs_error_assert_ptype"
  )
})

test_that("new_causal_wts() rejects a subclass that is not length one", {
  # A `subclass` of length two builds
  # `c("a", "b", "causal_wts", "vctrs_vctr", "double")`, a five-element class
  # vector, and both packages downstream pin the four-element shape with
  # `expect_s3_class(..., exact = TRUE)`. Catching that at construction is the
  # reason the low-level constructor takes the subclass at all.
  expect_error(
    new_causal_wts(1, subclass = c("a", "b")),
    class = "causalgenerics_invalid_argument_subclass"
  )
  expect_error(
    new_causal_wts(1, subclass = character()),
    class = "causalgenerics_invalid_argument_subclass"
  )
})

test_that("new_causal_wts() rejects a subclass that is not a character vector", {
  expect_error(
    new_causal_wts(1, subclass = 1),
    class = "causalgenerics_invalid_argument_subclass"
  )
  expect_error(
    new_causal_wts(1, subclass = TRUE),
    class = "causalgenerics_invalid_argument_subclass"
  )
  expect_error(
    new_causal_wts(1, subclass = NULL),
    class = "causalgenerics_invalid_argument_subclass"
  )
})

test_that("new_causal_wts() rejects a missing or empty subclass", {
  # Both of these are single character strings, so a type and length check on
  # its own lets them through. Either one names a class that no method can be
  # registered against, and `""` also makes `class(x)[[1]]` empty, which is what
  # vctrs keys `vec_ptype2()` and `vec_cast()` dispatch on.
  expect_error(
    new_causal_wts(1, subclass = NA_character_),
    class = "causalgenerics_invalid_argument_subclass"
  )
  expect_error(
    new_causal_wts(1, subclass = ""),
    class = "causalgenerics_invalid_argument_subclass"
  )
})

test_that("the subclass error carries the general condition class", {
  # The two-class shape is what `stop_no_method()` already establishes: a
  # specific class for tests and downstream handlers that want this one error,
  # and a general class for callers that want any argument this package
  # rejected.
  expect_error(
    new_causal_wts(1, subclass = c("a", "b")),
    class = "causalgenerics_invalid_argument"
  )
})
