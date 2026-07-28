# The condition contract is pinned here, at the two helpers that build it, and
# nowhere else. Every other assertion in the suite matches these conditions by
# class membership, which is the right altitude for a call site: it says the
# error a function signals is the one the caller is told to handle, and it stays
# true when the helper's message or fields change. What membership cannot say is
# what the helpers themselves promise, and that is what the two tests below are
# for.
#
# Class order is part of that promise. `class(cnd)` is the S3 dispatch order for
# the condition object, so `print()`, `format()`, `conditionMessage()`, and
# rlang's `cnd_header()`, `cnd_body()`, and `cnd_footer()` all consult it in
# sequence: with the specific class first, a method registered for one argument
# or one generic wins over a method registered for the general class, and with
# the order reversed the specific method would never run. This is not about
# `tryCatch()`, which picks a handler by the order the handlers were supplied
# rather than by the order of the condition's classes.

test_that("stop_invalid_argument() builds the documented condition", {
  # A named wrapper rather than a bare call, so that the call the helper records
  # is one this test can write down. `sys.call(-1)` records the caller's call,
  # and for a helper invoked straight from `tryCatch()` that is one of
  # `tryCatch()`'s own frames.
  reject_arg <- function(value) {
    stop_invalid_argument("value", "be a numeric vector")
  }

  cnd <- tryCatch(reject_arg(1), error = identity)

  expect_identical(
    class(cnd),
    c(
      "causalgenerics_invalid_argument_value",
      "causalgenerics_invalid_argument",
      "error",
      "condition"
    )
  )
  # The field, which is what a handler reads to report which argument was
  # rejected without parsing the message for it.
  expect_identical(cnd$arg, "value")
  # `must` completes the sentence, so the helper owns the backticks, the leading
  # "must", and the period.
  expect_identical(
    conditionMessage(cnd),
    "`value` must be a numeric vector."
  )
  expect_identical(conditionCall(cnd), quote(reject_arg(1)))
})

test_that("stop_no_method() builds the documented condition", {
  refuse_method <- function(x) stop_no_method("estimand", x)

  # Two classes, so that the message pins which one it names rather than
  # matching whatever a single-class object would give either way.
  x <- structure(list(), class = c("cg_probe", "cg_probe_base"))
  cnd <- tryCatch(refuse_method(x), error = identity)

  expect_identical(
    class(cnd),
    c(
      "causalgenerics_no_method_estimand",
      "causalgenerics_no_method",
      "error",
      "condition"
    )
  )
  expect_identical(cnd$generic, "estimand")
  expect_identical(
    conditionMessage(cnd),
    "No `estimand()` method for an object of class <cg_probe>."
  )
  expect_identical(conditionCall(cnd), quote(refuse_method(x)))
})
