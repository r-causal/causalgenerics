test_that("local_s3_method() registers a method in the namespace method table", {
  local_s3_method("ess", "cg_probe", function(x, ...) {
    "registered"
  })

  expect_identical(
    s3_methods_table("ess"),
    asNamespace("causalgenerics")$.__S3MethodsTable__.
  )
  expect_true(
    exists("ess.cg_probe", envir = s3_methods_table("ess"), inherits = FALSE)
  )
  expect_identical(
    dispatch_from_baseenv(ess, structure(1, class = "cg_probe")),
    "registered"
  )
})

test_that("local_s3_method() removes the method when the scope that registered it exits", {
  x <- structure(1, class = "cg_probe")

  # `local_s3_method()` defers its cleanup to the frame it was called from, so
  # the registration happens inside a function this test can watch return. That
  # also keeps the test from depending on a sibling test having registered the
  # class first, which would leave it vacuous whenever the file is run filtered
  # to this test alone.
  register_first <- function() {
    local_s3_method("ess", "cg_probe", function(x, ...) {
      "first"
    })
    expect_identical(dispatch_from_baseenv(ess, x), "first")
  }
  expect_no_warning(register_first())

  expect_false(
    exists("ess.cg_probe", envir = s3_methods_table("ess"), inherits = FALSE)
  )
  expect_error(
    dispatch_from_baseenv(ess, x),
    class = "causalgenerics_no_method"
  )

  # `registerS3method()` overwrites, so this half cannot fail by reading a stale
  # method. What it establishes is that registering the same class again after a
  # cleanup works, so one test can never block a later one.
  local_s3_method("ess", "cg_probe", function(x, ...) {
    "second"
  })
  expect_identical(dispatch_from_baseenv(ess, x), "second")
})

test_that("local_s3_method() uses the table belonging to a foreign generic", {
  # `registerS3method()` stores a method in the table of the environment the
  # generic lives in, which for `median()` is stats' namespace rather than this
  # package's. A lookup pinned to this package would leave the entry in stats'
  # table permanently and warn while trying to remove something that was never
  # there. Later work here registers methods for several generics this package
  # does not own.
  x <- structure(1, class = "cg_probe")
  stats_table <- asNamespace("stats")$.__S3MethodsTable__.

  register_median <- function() {
    local_s3_method("median", "cg_probe", function(x, ...) {
      "foreign"
    })
    expect_identical(s3_methods_table("median"), stats_table)
    expect_true(
      exists("median.cg_probe", envir = stats_table, inherits = FALSE)
    )
    expect_false(
      exists(
        "median.cg_probe",
        envir = s3_methods_table("ess"),
        inherits = FALSE
      )
    )
    expect_identical(dispatch_from_baseenv(median, x), "foreign")
  }
  expect_no_warning(register_median())

  expect_false(exists("median.cg_probe", envir = stats_table, inherits = FALSE))
})

test_that("local_s3_method() uses base's table for primitives and group generics", {
  # `[` is a primitive rather than a closure, and `Summary` is a group generic
  # with no function of its own to take an environment from. `registerS3method()`
  # files both under base, and the helper has to agree.
  base_table <- asNamespace("base")$.__S3MethodsTable__.

  register_base_generics <- function() {
    local_s3_method("[", "cg_probe", function(x, ...) {
      "subset"
    })
    local_s3_method("Summary", "cg_probe", function(..., na.rm = FALSE) {
      "summarized"
    })
    expect_identical(s3_methods_table("["), base_table)
    expect_identical(s3_methods_table("Summary"), base_table)
    expect_true(exists("[.cg_probe", envir = base_table, inherits = FALSE))
    expect_true(exists(
      "Summary.cg_probe",
      envir = base_table,
      inherits = FALSE
    ))
  }
  expect_no_warning(register_base_generics())

  expect_false(exists("[.cg_probe", envir = base_table, inherits = FALSE))
  expect_false(exists("Summary.cg_probe", envir = base_table, inherits = FALSE))
})

test_that("local_s3_method() restores a method it replaced", {
  # Registering over a method the NAMESPACE already supplies must put the
  # original back. Removing the entry instead would take `ess.default()` out of
  # the method table for the rest of the session, after which an unregistered
  # class would raise base R's dispatch failure rather than this package's
  # condition, and every test of the default behavior would be asserting the
  # wrong thing.
  x <- structure(list(), class = "cg_unregistered")

  replace_default <- function() {
    local_s3_method("ess", "default", function(x, ...) {
      "replacement default"
    })
    expect_identical(dispatch_from_baseenv(ess, x), "replacement default")
  }
  expect_no_warning(replace_default())

  expect_identical(
    get("ess.default", envir = s3_methods_table("ess"), inherits = FALSE),
    ess.default
  )
  expect_error(
    dispatch_from_baseenv(ess, x),
    class = "causalgenerics_no_method"
  )
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
