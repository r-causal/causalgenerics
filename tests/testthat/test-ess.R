test_that("ess() dispatches to a registered method", {
  ess.cg_weight <- function(x, ...) {
    sum(unclass(x))^2 / sum(unclass(x)^2)
  }

  wts <- structure(rep(2, 5), class = "cg_weight")
  expect_identical(ess(wts), 5)
})

test_that("ess() passes further arguments through the dots", {
  ess.cg_weight <- function(x, na.rm = FALSE, ...) {
    w <- unclass(x)
    sum(w, na.rm = na.rm)^2 / sum(w^2, na.rm = na.rm)
  }

  wts <- structure(c(1, 1, NA), class = "cg_weight")
  expect_identical(ess(wts, na.rm = TRUE), 2)
})

test_that("ess() errors informatively for an unregistered class", {
  x <- structure(list(), class = "cg_unregistered")
  expect_error(ess(x), class = "causalgenerics_no_method")
  expect_snapshot(error = TRUE, ess(x))
})
