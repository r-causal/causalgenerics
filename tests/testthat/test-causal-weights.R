test_that("is_causal_wt() defaults to FALSE", {
  expect_false(dispatch_from_baseenv(is_causal_wt, 1:10))
  expect_false(
    dispatch_from_baseenv(
      is_causal_wt,
      structure(list(), class = "cg_unregistered")
    )
  )
})

test_that("is_causal_wt() dispatches to a registered method", {
  local_s3_method("is_causal_wt", "cg_weight", function(x, ...) {
    TRUE
  })

  expect_true(
    dispatch_from_baseenv(is_causal_wt, structure(1, class = "cg_weight"))
  )
})

test_that("estimand() dispatches to a registered method", {
  local_s3_method("estimand", "cg_weight", function(x, ...) {
    attr(x, "estimand")
  })

  wts <- structure(1, class = "cg_weight", estimand = "ate")
  expect_identical(dispatch_from_baseenv(estimand, wts), "ate")
})

test_that("estimand<-() dispatches to a registered method", {
  local_s3_method("estimand<-", "cg_weight", function(x, ..., value) {
    attr(x, "estimand") <- value
    x
  })

  # The replacement form `estimand(wts) <- "att"` needs `estimand<-()` to be
  # visible where it is written, which `baseenv()` cannot offer. Calling the
  # replacement generic by value keeps the assertion isolated from the test
  # frame; the surface syntax it expands to is base R's, not this package's.
  wts <- structure(1, class = "cg_weight")
  wts <- dispatch_from_baseenv(`estimand<-`, wts, value = "att")
  expect_identical(attr(wts, "estimand"), "att")
})

test_that("estimand() and estimand<-() error for an unregistered class", {
  x <- structure(list(), class = "cg_unregistered")
  expect_error(
    dispatch_from_baseenv(estimand, x),
    class = "causalgenerics_no_method"
  )
  expect_error(
    dispatch_from_baseenv(`estimand<-`, x, value = "ate"),
    class = "causalgenerics_no_method"
  )
  expect_snapshot(error = TRUE, estimand(x))
})
