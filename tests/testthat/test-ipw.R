test_that("ipw() dispatches on the weighting object", {
  ipw.cg_weight <- function(ps_mod, outcome_mod, ...) {
    "dispatched"
  }

  ps_mod <- structure(list(), class = "cg_weight")
  expect_identical(ipw(ps_mod, outcome_mod = NULL), "dispatched")
})

test_that("ipw() errors for an object with no registered method", {
  x <- structure(list(), class = "cg_unregistered")
  expect_error(ipw(x, outcome_mod = NULL))
})
