test_that("local_s3_method() registers a method in the namespace method table", {
  local_s3_method("ess", "cg_probe", function(x, ...) {
    "registered"
  })

  expect_true(
    exists("ess.cg_probe", envir = s3_methods_table(), inherits = FALSE)
  )
  expect_identical(
    dispatch_from_baseenv(ess, structure(1, class = "cg_probe")),
    "registered"
  )
})

test_that("local_s3_method() removes the method when the test that registered it exits", {
  # This test runs after the one above, which registered `ess.cg_probe`. A
  # failure in the first half means the cleanup in `local_s3_method()` left the
  # entry behind, and the second half means a later registration for the same
  # class would be reading a stale method.
  x <- structure(1, class = "cg_probe")
  expect_false(
    exists("ess.cg_probe", envir = s3_methods_table(), inherits = FALSE)
  )
  expect_error(
    dispatch_from_baseenv(ess, x),
    class = "causalgenerics_no_method"
  )

  local_s3_method("ess", "cg_probe", function(x, ...) {
    "second"
  })
  expect_identical(dispatch_from_baseenv(ess, x), "second")
})

test_that("dispatch_from_baseenv() cannot see a method defined in the test frame", {
  # The guard on the guard. `UseMethod()` searches the frame the generic was
  # called from, so the direct call below succeeds on nothing but a local
  # function. If the isolated call ever succeeds too, every dispatch assertion
  # in this suite has stopped proving that registration happened.
  ess.cg_probe <- function(x, ...) {
    "frame local"
  }
  x <- structure(1, class = "cg_probe")

  expect_identical(ess(x), "frame local")
  expect_error(
    dispatch_from_baseenv(ess, x),
    class = "causalgenerics_no_method"
  )
})
