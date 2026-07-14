test_that("is_causal_wt() defaults to FALSE and dispatches when registered", {
  expect_false(is_causal_wt(1:10))
  expect_false(is_causal_wt(structure(list(), class = "cg_unregistered")))

  is_causal_wt.cg_weight <- function(x, ...) {
    TRUE
  }
  expect_true(is_causal_wt(structure(1, class = "cg_weight")))
})

test_that("estimand() dispatches to a registered method", {
  estimand.cg_weight <- function(x, ...) {
    attr(x, "estimand")
  }

  wts <- structure(1, class = "cg_weight", estimand = "ate")
  expect_identical(estimand(wts), "ate")
})

test_that("estimand<-() dispatches to a registered method", {
  `estimand<-.cg_weight` <- function(x, ..., value) {
    attr(x, "estimand") <- value
    x
  }

  wts <- structure(1, class = "cg_weight")
  estimand(wts) <- "att"
  expect_identical(attr(wts, "estimand"), "att")
})

test_that("estimand() and estimand<-() error for an unregistered class", {
  x <- structure(list(), class = "cg_unregistered")
  expect_error(estimand(x), class = "causalgenerics_no_method")
  expect_error(estimand(x) <- "ate", class = "causalgenerics_no_method")
  expect_snapshot(error = TRUE, estimand(x))
})
