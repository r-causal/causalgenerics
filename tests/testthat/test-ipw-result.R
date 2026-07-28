# `new_ipw()`, `print.ipw()`, `format_model_call()`, and `as.data.frame.ipw()`
# are the shared result layer for `ipw()`. Every package that supplies an
# `ipw()` method builds its return with this constructor and inherits these
# methods, so the field names, their order, and the printed and data frame forms
# are a cross-package contract rather than a presentation detail. Two packages
# each defining their own `print.ipw()` would collide in the shared S3 method
# table, which is the situation this package exists to prevent.
#
# The estimates frames below are written out literally, in the shape the `ipw()`
# return contract documents, rather than produced by fitting a model. Nothing
# here needs a real estimator, and this package holds no model-fitting
# dependency to build one with.

# A binary-exposure estimates frame: one row per effect measure, no `comparison`
# column. The risk difference is on the raw scale and the two ratios are on the
# log scale, which is what makes `exponentiate` mean anything.
binary_estimates <- function() {
  data.frame(
    effect = c("rd", "log(rr)", "log(or)"),
    estimate = c(0.199882, 0.560414, 0.878313),
    std.err = c(0.092425, 0.273519, 0.418661),
    z = c(2.1626, 2.0489, 2.0979),
    ci.lower = c(0.018732, 0.024326, 0.057753),
    ci.upper = c(0.381032, 1.096502, 1.698873),
    conf.level = 0.95,
    p.value = c(0.030570, 0.040470, 0.035910)
  )
}

# A categorical-exposure estimates frame. The effect labels repeat across
# contrasts, so a `comparison` column sits immediately after `effect` and names
# the non-reference and reference level of each one.
categorical_estimates <- function() {
  data.frame(
    effect = rep(c("rd", "log(rr)", "log(or)"), times = 2),
    comparison = rep(c("b vs a", "c vs a"), each = 3),
    estimate = c(0.081945, 0.168870, 0.328762, 0.166939, 0.318293, 0.676435),
    std.err = c(0.050387, 0.104633, 0.203058, 0.045182, 0.091898, 0.185786),
    z = c(1.6263, 1.6139, 1.6191, 3.6948, 3.4635, 3.6409),
    ci.lower = c(-0.016811, -0.036207, -0.069225, 0.078384, 0.138176, 0.312300),
    ci.upper = c(0.180701, 0.373947, 0.726749, 0.255494, 0.498410, 1.040570),
    conf.level = 0.95,
    p.value = c(
      0.1038822,
      0.1065433,
      0.1054363,
      0.0002200,
      0.0005331,
      0.0002717
    )
  )
}

# A continuous-outcome estimates frame. There is only a difference in means, so
# neither ratio row is present and `exponentiate = TRUE` has nothing to act on.
continuous_estimates <- function() {
  data.frame(
    effect = "diff",
    estimate = 2.25255,
    std.err = 0.17524,
    z = 12.854,
    ci.lower = 1.909083,
    ci.upper = 2.596017,
    conf.level = 0.95,
    p.value = 4.5e-38
  )
}

# Fixed rather than simulated data, so that the deparsed calls the snapshots pin
# never move between runs.
ipw_result_data <- function() {
  data.frame(
    x = rep(c(-1.5, -0.5, 0.5, 1.5), each = 5),
    z = rep(c(0, 1), 10),
    y = rep(c(0, 1, 1, 0, 1), 4)
  )
}

# The model slots hold real fits, so that the `Call:` lines the snapshots pin are
# the ones a caller actually sees rather than a stand-in.
ipw_result_models <- function() {
  dat <- ipw_result_data()
  list(
    ps_mod = glm(z ~ x, family = binomial(), data = dat),
    outcome_mod = glm(y ~ z, family = quasibinomial(), data = dat)
  )
}

# Assemble a result the way a method would, so that each test names only the
# part it is about.
ipw_result <- function(estimates) {
  mods <- ipw_result_models()
  new_ipw(
    estimand = "ate",
    ps_mod = mods$ps_mod,
    outcome_mod = mods$outcome_mod,
    estimates = estimates,
    se_method = "mestimation",
    fit = NULL
  )
}

# ---- new_ipw() ---------------------------------------------------------------

test_that("new_ipw() builds the documented six-field list", {
  estimates <- binary_estimates()
  fit <- structure(list(theta = c(1, 2)), class = "cg_mestimator")

  res <- new_ipw(
    estimand = "ate",
    ps_mod = "a propensity score model",
    outcome_mod = "an outcome model",
    estimates = estimates,
    se_method = "mestimation",
    fit = fit
  )

  expect_type(res, "list")
  expect_s3_class(res, "ipw", exact = TRUE)
  expect_identical(
    names(res),
    c("estimand", "ps_mod", "outcome_mod", "estimates", "se_method", "fit")
  )

  expect_identical(res$estimand, "ate")
  expect_identical(res$ps_mod, "a propensity score model")
  expect_identical(res$outcome_mod, "an outcome model")
  expect_identical(res$estimates, estimates)
  expect_identical(res$se_method, "mestimation")
  expect_identical(res$fit, fit)
})

test_that("new_ipw() keeps a NULL fit as a named field", {
  # The linearization path has no M-estimator object to report. `fit` still has
  # to be present and `NULL`, because callers read the field by name and a list
  # that dropped it would have five elements rather than six.
  res <- new_ipw(
    estimand = "att",
    ps_mod = NULL,
    outcome_mod = NULL,
    estimates = binary_estimates(),
    se_method = "linearization",
    fit = NULL
  )

  expect_length(res, 6L)
  expect_identical(
    names(res),
    c("estimand", "ps_mod", "outcome_mod", "estimates", "se_method", "fit")
  )
  expect_null(res$fit)
  expect_identical(res$se_method, "linearization")
})

test_that("new_ipw() takes its arguments in the documented order", {
  # Both packages downstream call this positionally in places, so the signature
  # order is part of the contract and not just the field order.
  estimates <- binary_estimates()

  positional <- new_ipw("ate", "ps", "outcome", estimates, "mestimation", "fit")
  named <- new_ipw(
    estimand = "ate",
    ps_mod = "ps",
    outcome_mod = "outcome",
    estimates = estimates,
    se_method = "mestimation",
    fit = "fit"
  )

  expect_identical(positional, named)
})

# ---- print.ipw() -------------------------------------------------------------

# Each print test snapshots the whole output and then asserts the few lines that
# carry the contract. The snapshot alone would not hold it, for two reasons. On
# the run that first writes `_snaps/ipw-result.md` the implementation supplies
# its own baseline, so a printed form that was wrong from the outset is recorded
# as correct rather than rejected. And testthat skips snapshot comparison on
# CRAN, so on the checks that matter most the snapshot verifies nothing at all.
# The named assertions below are what make these tests fail rather than record,
# and they are the only print coverage that runs on CRAN. This pairing mirrors
# the one the error tests use, where the condition class is asserted alongside
# the snapshotted message.

# Where a model's heading sits in the captured output. A model's `Call:` line is
# the line directly below its heading, so reading the call through the heading
# index is what ties it to the model it reports and fixes the order of the two
# blocks. Asserting that both headings and both calls merely appear somewhere
# holds just as well when the two models are printed under each other's heading,
# which misstates the caller's own fits rather than formatting them differently.
heading_index <- function(out, heading) {
  index <- which(out == heading)
  expect_length(index, 1L)
  index
}

test_that("print() summarizes a binary-exposure result", {
  res <- ipw_result(binary_estimates())

  expect_snapshot(print(res))

  out <- capture.output(print(res))

  expect_match(out, "^Inverse Probability Weight Estimator$", all = FALSE)
  # The estimand is reported in upper case, so `ate` reads as ATE.
  expect_match(out, "^Estimand: ATE\\s*$", all = FALSE)
  # The propensity score model is reported first, and each `Call:` line belongs
  # to the heading above it. The two fixtures differ in both formula and family,
  # which is what makes the two lines tell each other apart.
  ps_heading <- heading_index(out, "Propensity Score Model:")
  outcome_heading <- heading_index(out, "Outcome Model:")

  expect_lt(ps_heading, outcome_heading)
  expect_match(out[ps_heading + 1L], "glm(formula = z ~ x", fixed = TRUE)
  expect_match(out[outcome_heading + 1L], "glm(formula = y ~ z", fixed = TRUE)

  expect_match(out, "^Estimates:$", all = FALSE)

  # Rows are labelled by effect, and the character `effect` column is gone
  # rather than formatted as a number.
  expect_match(out, "^rd +0\\.199882 ", all = FALSE)
  expect_match(out, "^log\\(rr\\) +0\\.560414 ", all = FALSE)
  expect_match(out, "^log\\(or\\) +0\\.878313 ", all = FALSE)
  expect_false(any(grepl("effect", out, fixed = TRUE)))

  # Significance stars are left on, which is what `has.Pvalue` is for.
  expect_match(out, "^rd .*\\*$", all = FALSE)
  expect_match(out, "^Signif\\. codes:", all = FALSE)
})

test_that("print() keys rows by effect and comparison for a categorical result", {
  # The categorical shape is the only place the two print paths differ. Row
  # labels become `effect` and `comparison` together, and both character columns
  # have to be gone before `printCoefmat()` sees the frame.
  res <- ipw_result(categorical_estimates())

  expect_snapshot(print(res))

  out <- capture.output(print(res))

  expect_match(out, "^rd b vs a +0\\.081945 ", all = FALSE)
  expect_match(out, "^log\\(rr\\) b vs a +0\\.168870 ", all = FALSE)
  expect_match(out, "^log\\(or\\) c vs a +0\\.676435 ", all = FALSE)

  # Dropping `effect` alone is the failure worth guarding against, because it
  # does not error: `printCoefmat()` sends the character `comparison` column
  # through `data.matrix()`, which factor-codes it into a column reading
  # `1.000000` and `2.000000` alongside the real estimates.
  expect_false(any(grepl("comparison", out, fixed = TRUE)))
  expect_false(any(grepl("effect", out, fixed = TRUE)))

  # With six rows the stars do not fit beside the estimates and
  # `printCoefmat()` wraps them into a block of their own, still keyed by the
  # same row labels.
  expect_match(out, "^rd c vs a +\\*\\*\\*$", all = FALSE)
  expect_match(out, "^Signif\\. codes:", all = FALSE)
})

test_that("print() summarizes a continuous-outcome result", {
  # A continuous outcome has one row and no ratio rows, so it takes the same
  # path as the binary shape with a single effect label.
  res <- ipw_result(continuous_estimates())

  expect_snapshot(print(res))

  out <- capture.output(print(res))

  expect_match(out, "^diff +2\\.25255 ", all = FALSE)
  expect_match(out, "^Signif\\. codes:", all = FALSE)
  expect_false(any(grepl("effect", out, fixed = TRUE)))
})

test_that("print() returns its input invisibly", {
  res <- ipw_result(binary_estimates())

  # The printed form is snapshotted above, so here the output only needs to go
  # somewhere other than the test report.
  withr::local_output_sink(withr::local_tempfile())
  printed <- withVisible(print(res))

  expect_false(printed$visible)
  expect_identical(printed$value, res)
})

# ---- format_model_call() -----------------------------------------------------

test_that("format_model_call() deparses a call the object exposes", {
  dat <- ipw_result_data()
  mod <- lm(y ~ x, data = dat)

  expect_identical(format_model_call(mod), "lm(formula = y ~ x, data = dat)")
})

test_that("format_model_call() joins a call that deparses to several lines", {
  # `deparse()` breaks a long call into several elements. The result has to be a
  # single string, because `print()` writes it as one field of one line.
  dat <- ipw_result_data()
  mod <- lm(y ~ x + I(x^2) + I(x^3) + I(x^4) + I(x^5) + I(x^6), data = dat)

  formatted <- format_model_call(mod)

  expect_length(formatted, 1L)
  expect_match(formatted, "\n", fixed = TRUE)
  expect_match(formatted, "lm(formula = y ~ x + I(x^2)", fixed = TRUE)
  expect_match(formatted, "data = dat)", fixed = TRUE)
})

test_that("format_model_call() falls back to a class label with no call", {
  # `getCall()` returns `NULL` for an object that is subsettable but records no
  # call, which is the ordinary case for a weighting object that was not fitted
  # by a modelling function.
  expect_identical(
    format_model_call(structure(list(), class = "cg_no_call")),
    "<cg_no_call>"
  )
  expect_identical(
    format_model_call(structure(list(), class = c("cg_fit", "cg_base"))),
    "<cg_fit/cg_base>"
  )
})

test_that("format_model_call() falls back for an object that cannot be subset", {
  # `getCall()` reaches the call through `getElement()`, which subsets. An
  # object that cannot be subset therefore raises a condition instead of
  # returning `NULL`, and the fallback has to catch that as well. This is the
  # path a fit built on an object system with no `$` method takes, and it is why
  # `print()` works for a weighting object of any shape.
  mod <- structure(1:3, class = "cg_unsubsettable")

  expect_error(stats::getCall(mod))
  expect_identical(format_model_call(mod), "<cg_unsubsettable>")
})

test_that("print() shows the class label for a model with no accessible call", {
  res <- new_ipw(
    estimand = "ate",
    ps_mod = structure(1:3, class = "cg_unsubsettable"),
    outcome_mod = structure(list(), class = "cg_no_call"),
    estimates = binary_estimates(),
    se_method = "mestimation",
    fit = NULL
  )

  out <- capture.output(print(res))

  ps_heading <- heading_index(out, "Propensity Score Model:")
  outcome_heading <- heading_index(out, "Outcome Model:")

  expect_match(out[ps_heading + 1L], "Call: <cg_unsubsettable>", fixed = TRUE)
  expect_match(out[outcome_heading + 1L], "Call: <cg_no_call>", fixed = TRUE)
})

# ---- as.data.frame.ipw() -----------------------------------------------------

test_that("as.data.frame() returns the estimates unchanged by default", {
  estimates <- binary_estimates()
  res <- ipw_result(estimates)

  expect_identical(as.data.frame(res), estimates)
})

test_that("as.data.frame() returns a categorical estimates frame unchanged", {
  estimates <- categorical_estimates()
  res <- ipw_result(estimates)

  expect_identical(as.data.frame(res), estimates)
})

test_that("as.data.frame() passes row.names through", {
  res <- ipw_result(binary_estimates())

  df <- as.data.frame(res, row.names = c("first", "second", "third"))

  expect_identical(rownames(df), c("first", "second", "third"))
})

test_that("as.data.frame(exponentiate = TRUE) moves only the ratio rows", {
  estimates <- binary_estimates()
  res <- ipw_result(estimates)

  df <- as.data.frame(res, exponentiate = TRUE)

  # The two log rows are relabelled and their point estimate and interval move
  # to the natural scale. The risk difference was never on the log scale, so it
  # is untouched.
  expect_identical(df$effect, c("rd", "rr", "or"))
  expect_identical(
    df$estimate,
    c(estimates$estimate[1], exp(estimates$estimate[2:3]))
  )
  expect_identical(
    df$ci.lower,
    c(estimates$ci.lower[1], exp(estimates$ci.lower[2:3]))
  )
  expect_identical(
    df$ci.upper,
    c(estimates$ci.upper[1], exp(estimates$ci.upper[2:3]))
  )

  # Inference stays on the log scale, where it is done. Exponentiating a
  # standard error would state a quantity nobody asked for.
  expect_identical(df$std.err, estimates$std.err)
  expect_identical(df$z, estimates$z)
  expect_identical(df$p.value, estimates$p.value)
  expect_identical(df$conf.level, estimates$conf.level)

  expect_identical(names(df), names(estimates))
})

test_that("as.data.frame(exponentiate = TRUE) relabels every comparison", {
  estimates <- categorical_estimates()
  res <- ipw_result(estimates)

  df <- as.data.frame(res, exponentiate = TRUE)

  expect_identical(df$effect, rep(c("rd", "rr", "or"), times = 2))
  # The comparison labels identify the contrast, not the scale, so they survive
  # untouched and stay in place after `effect`.
  expect_identical(df$comparison, estimates$comparison)
  expect_identical(names(df), names(estimates))

  ratios <- df$effect %in% c("rr", "or")
  expect_identical(df$estimate[ratios], exp(estimates$estimate[ratios]))
  expect_identical(df$ci.lower[ratios], exp(estimates$ci.lower[ratios]))
  expect_identical(df$ci.upper[ratios], exp(estimates$ci.upper[ratios]))
  expect_identical(df$std.err, estimates$std.err)
})

test_that("as.data.frame(exponentiate = TRUE) matches the log labels exactly", {
  # The rows to move are the ones labelled `log(rr)` and `log(or)`, not every
  # row whose label happens to contain `rr` or `or`. An already exponentiated
  # frame is the fixture that separates the two: its labels are `rr` and `or`,
  # so a substring match would exponentiate it a second time while an equality
  # match leaves it alone.
  once <- as.data.frame(ipw_result(binary_estimates()), exponentiate = TRUE)
  twice <- as.data.frame(ipw_result(once), exponentiate = TRUE)

  expect_identical(twice, once)
})

test_that("as.data.frame(exponentiate = TRUE) is a no-op with no ratio rows", {
  # A continuous outcome reports a difference in means and nothing else, so
  # asking for the natural scale has to leave the frame alone rather than error
  # on the absent rows.
  estimates <- continuous_estimates()
  res <- ipw_result(estimates)

  expect_identical(as.data.frame(res, exponentiate = TRUE), estimates)
})

# ---- registration and export -------------------------------------------------

test_that("the print and data frame methods are registered against ipw", {
  # Every assertion above reaches these methods from the test frame, and
  # `UseMethod()` searches the calling frame before it consults the method
  # table. All of them therefore pass whether or not the NAMESPACE carries an
  # `S3method()` directive. Downstream packages call `print()` and
  # `as.data.frame()` from their own namespaces, where the table is the only
  # route, so table membership is the claim that matters and nothing else here
  # makes it.
  #
  # Both generics are base closures, so both entries belong to base's table
  # rather than to this package's.
  generics <- c("print", "as.data.frame")

  registered <- vapply(
    generics,
    function(generic) {
      table <- s3_methods_table(generic)
      !is.null(table) &&
        exists(paste0(generic, ".ipw"), envir = table, inherits = FALSE)
    },
    logical(1)
  )

  expect_identical(generics[!registered], character())
})

test_that("new_ipw() is exported and format_model_call() is not", {
  # `new_ipw()` is the reason the result layer moved here: a package supplying
  # an `ipw()` method has to be able to build the object. `format_model_call()`
  # is a helper of `print()` and stays internal.
  #
  # `pkgload::load_all()` attaches every object in the package to the search
  # path, so an assertion made against `package:causalgenerics` would hold for
  # an unexported function too. The namespace's own export registry is not
  # widened that way, and the second assertion is what keeps that true: it names
  # a function that is present in the namespace and must be absent from the
  # registry, so if the registry ever did report everything this test would fail
  # rather than pass for the wrong reason.
  #
  # A namespace's parent chain runs out through its imports environment and
  # `namespace:base` to the global environment and the search path, so the
  # default `inherits = TRUE` would report a name found anywhere on that chain,
  # including the `package:causalgenerics` copy that `load_all()` attaches.
  # `inherits = FALSE` is what makes the assertion about the namespace itself.
  exports <- getNamespaceExports("causalgenerics")

  expect_true("new_ipw" %in% exports)

  expect_true(
    exists(
      "format_model_call",
      envir = asNamespace("causalgenerics"),
      inherits = FALSE
    )
  )
  expect_false("format_model_call" %in% exports)
})
