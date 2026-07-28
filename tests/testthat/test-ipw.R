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
