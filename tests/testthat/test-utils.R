# The condition contract is pinned here, at the helpers that build it, and
# nowhere else. Every other assertion in the suite matches these conditions by
# class membership, which is the right altitude for a call site: it says the
# error a function signals is the one the caller is told to handle, and it stays
# true when the helper's message or fields change. What membership cannot say is
# what the helpers themselves promise, and that is what the tests below are for,
# one per helper.
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

test_that("stop_no_vcov() builds the documented condition", {
  # A named wrapper rather than a bare call, for the reason the test above uses
  # one: `sys.call(-1)` records the caller's call, and this is the caller.
  refuse_vcov <- function(result) stop_no_vcov(result)

  cnd <- tryCatch(refuse_vcov("ipw"), error = identity)

  expect_identical(
    class(cnd),
    c(
      "causalgenerics_no_vcov_ipw",
      "causalgenerics_no_vcov",
      "error",
      "condition"
    )
  )
  # The field, which is what a handler reads to report which kind of object
  # carried no covariance without parsing the message for it.
  expect_identical(cnd$result, "ipw")
  # The message names the object in backticks, joins the two clauses with a
  # semicolon, and writes the attribute the contract is keyed to as code. It
  # says "object" rather than "result" because the helper answers for a result
  # and for a component model of one alike, and the sentence has to read
  # correctly whichever of the two is named.
  expect_identical(
    conditionMessage(cnd),
    paste0(
      "This `ipw` object records no covariance to report; the package that ",
      "produced it attaches one when it supports the `ipw_vcov` contract."
    )
  )
  expect_identical(conditionCall(cnd), quote(refuse_vcov("ipw")))

  # The name is interpolated rather than written into the sentence, which is
  # what lets the component model name itself here too.
  model_cnd <- tryCatch(refuse_vcov("ipw_model"), error = identity)
  expect_match(conditionMessage(model_cnd), "ipw_model", fixed = TRUE)
})

test_that("stop_no_conditional_vcov() builds the documented condition", {
  # A named wrapper rather than a bare call, for the reason the test above uses
  # one: `sys.call(-1)` records the caller's call, and this is the caller.
  refuse_conditional <- function() stop_no_conditional_vcov()

  cnd <- tryCatch(refuse_conditional(), error = identity)

  # The specific class first and the general one behind it. A handler written
  # for any covariance this package cannot report catches this one alongside
  # `stop_no_vcov()`, and a handler written for the conditional reading tells
  # the two apart: the first says the result records no covariance of its
  # effects, the second that the outcome model carries no corrected one.
  expect_identical(
    class(cnd),
    c(
      "causalgenerics_no_conditional_vcov",
      "causalgenerics_no_vcov",
      "error",
      "condition"
    )
  )

  # The shape of the message rather than its wording, which the snapshots at the
  # accessors record: a single string naming the reading that cannot be
  # answered, ending in a period the way every condition here does.
  message <- conditionMessage(cnd)

  expect_type(message, "character")
  expect_length(message, 1L)
  expect_match(message, "conditional", fixed = TRUE)
  expect_match(message, "\\.$")
  expect_identical(conditionCall(cnd), quote(refuse_conditional()))
})

test_that("stop_conditional_vcov_mismatch() builds the documented condition", {
  # A named wrapper rather than a bare call, for the reason the tests above use
  # one: `sys.call(-1)` records the caller's call, and this is the caller.
  refuse_pairing <- function() {
    stop_conditional_vcov_mismatch(
      c("theta1", "theta2"),
      c("(Intercept)", "z")
    )
  }

  cnd <- tryCatch(refuse_pairing(), error = identity)

  # The specific class first and the general one behind it, which is the shape
  # every condition here has. A handler written for any covariance this package
  # cannot report catches this one alongside `stop_no_vcov()` and
  # `stop_no_conditional_vcov()`.
  expect_identical(
    class(cnd),
    c(
      "causalgenerics_conditional_vcov_mismatch",
      "causalgenerics_no_vcov",
      "error",
      "condition"
    )
  )

  # Not the reading's own refusal, and the class vector above is not the whole
  # of that claim: `inherits()` is what a handler asks, and `print.ipw()` asks
  # it of every condition the covariance lookup raises. Answering yes would put
  # the note that says to wrap the model over a block that is already there and
  # cannot be read against these coefficients.
  expect_false(inherits(cnd, "causalgenerics_no_conditional_vcov"))

  # The fields, which are what a handler reads to report the two label sets
  # without parsing the message for them.
  expect_identical(cnd$block_labels, c("theta1", "theta2"))
  expect_identical(cnd$coef_labels, c("(Intercept)", "z"))

  # The message names both sets, since the pairing is what failed and neither
  # set on its own says why. The lists are quoted and comma separated the way
  # `select_effects()` writes the labels of a surface, the two clauses are
  # joined with a semicolon, and the remedy names the constructor as code.
  expect_identical(
    conditionMessage(cnd),
    paste0(
      "The conditional covariance is labelled \"theta1\", \"theta2\" and the ",
      "outcome model reports coefficients named \"(Intercept)\", \"z\"; the ",
      "package that produced the result attaches the block labelled by ",
      "coefficient name with `new_ipw_model()`."
    )
  )
  expect_identical(conditionCall(cnd), quote(refuse_pairing()))
})

test_that("stop_conditional_vcov_mismatch() names coefficients that have none", {
  # A model whose coefficients carry no names is one of the ways the pairing
  # cannot be made, and the clause that would list them has nothing to list.
  # Naming an empty set would read as a model reporting no coefficients at all,
  # which is a different fact about a different model.
  refuse_unnamed <- function() {
    stop_conditional_vcov_mismatch(c("theta1", "theta2"), NULL)
  }

  cnd <- tryCatch(refuse_unnamed(), error = identity)

  expect_identical(
    class(cnd),
    c(
      "causalgenerics_conditional_vcov_mismatch",
      "causalgenerics_no_vcov",
      "error",
      "condition"
    )
  )
  expect_null(cnd$coef_labels)
  expect_identical(
    conditionMessage(cnd),
    paste0(
      "The conditional covariance is labelled \"theta1\", \"theta2\" and the ",
      "outcome model reports unnamed coefficients; the package that produced ",
      "the result attaches the block labelled by coefficient name with ",
      "`new_ipw_model()`."
    )
  )
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
