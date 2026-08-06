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

# The categorical frame with the column naming its contrasts called `contrast`,
# which is the name the result layer reports it under. It is derived from the
# frame above rather than written out again, so that the two differ in the column
# name and in nothing else. The section at the foot of this file says why the
# name matters.
contrast_estimates <- function() {
  estimates <- categorical_estimates()
  names(estimates)[names(estimates) == "comparison"] <- "contrast"
  estimates
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
    wt_mod = glm(z ~ x, family = binomial(), data = dat),
    outcome_mod = glm(y ~ z, family = quasibinomial(), data = dat)
  )
}

# Assemble a result the way a method would, so that each test names only the
# part it is about.
ipw_result <- function(estimates) {
  mods <- ipw_result_models()
  new_ipw(
    estimand = "ate",
    wt_mod = mods$wt_mod,
    outcome_mod = mods$outcome_mod,
    estimates = estimates,
    se_method = "mestimation",
    fit = NULL
  )
}

# A result built with a given presentation mode, so that a test naming a mode
# does not restate the other six fields. The two model slots hold strings rather
# than fits, since nothing about the mode reads them.
ipw_with_effects <- function(effects) {
  new_ipw(
    estimand = "ate",
    wt_mod = "a propensity score model",
    outcome_mod = "an outcome model",
    estimates = binary_estimates(),
    se_method = "mestimation",
    fit = NULL,
    effects = effects
  )
}

# A result of the shape the constructor built before the mode existed: six
# fields and no `effects`. Results stored from an earlier version of a fitting
# package have this shape, as does one built by hand against the earlier
# contract. It is derived from the seven-field result rather than written out
# again, so that the two differ in the field and in nothing else.
legacy_result <- function(estimates = binary_estimates()) {
  fields <- unclass(ipw_result(estimates))
  fields$effects <- NULL
  structure(fields, class = "ipw")
}

# Data for the conditional fixtures, where the outcome moves with the exposure:
# three of the ten unexposed observations and seven of the ten exposed ones have
# `y == 1`. The frame the other fixtures use is balanced against both `x` and
# `z`, so its outcome model has a slope of zero to machine noise, and a number
# that small is not the same in its last digits from one platform to the next.
# The conditional table prints the coefficients themselves, so it needs
# coefficients worth printing.
conditional_data <- function() {
  data.frame(
    x = rep(c(-1.5, -0.5, 0.5, 1.5), each = 5),
    z = rep(c(0, 1), 10),
    y = c(rep(c(0, 1), 7), rep(c(1, 0), 3))
  )
}

# The models a conditional result is built around, fitted the way the marginal
# fixtures' are so that the metadata block above the table reads the same in
# both modes.
conditional_models <- function() {
  dat <- conditional_data()
  list(
    wt_mod = glm(z ~ x, family = binomial(), data = dat),
    outcome_mod = glm(y ~ z, family = quasibinomial(), data = dat)
  )
}

# The corrected covariance of that model's coefficients, in the shape a fitting
# package hands to `new_ipw_model()`: the outcome block of the stacked sandwich,
# labelled with the coefficient names. It is built from the model's own
# covariance so that the dimnames agree without being restated, and scaled so
# that it differs from it. Accounting for the weights having been estimated
# ordinarily widens the block, and the scaling is what separates "the printed
# standard errors are the corrected ones" from "the printed standard errors are
# whatever the outcome model computed for itself".
corrected_outcome_vcov <- function(mod) {
  stats::vcov(mod) * 1.4
}

# The corrected block labelled in the reverse of the coefficient order, which is
# the shape a fitting package produces when its stacked system holds the outcome
# parameters the other way round. The two diagonal entries are far apart, so a
# table that took the standard errors off the diagonal by position rather than
# by label would print the intercept's against the slope and the slope's against
# the intercept, and it would look exactly like a correct table while doing it.
reversed_outcome_vcov <- function() {
  matrix(
    c(0.0625, 0.03125, 0.03125, 0.25),
    nrow = 2,
    dimnames = list(c("z", "(Intercept)"), c("z", "(Intercept)"))
  )
}

# A block labelled with the parameter names of a stacked system rather than with
# the model's coefficient names. `new_ipw_model()` takes it, since only the
# fitting package knows which block of the sandwich belongs to which model, and
# there is no reading of it against these coefficients.
theta_outcome_vcov <- function() {
  matrix(
    c(0.0625, 0.03125, 0.03125, 0.25),
    nrow = 2,
    dimnames = list(c("theta1", "theta2"), c("theta1", "theta2"))
  )
}

# A result that presents the conditional reading. `wrap = FALSE` leaves the
# outcome model as it was fitted, carrying no covariance from the joint
# estimation, which is the shape `print()` has to report rather than refuse.
# `vcov` names the block to wrap it with, for the tests that are about which
# block the table is built from; the default is the corrected one.
conditional_result <- function(wrap = TRUE, vcov = NULL) {
  mods <- conditional_models()
  outcome_mod <- mods$outcome_mod
  if (wrap) {
    if (is.null(vcov)) {
      vcov <- corrected_outcome_vcov(outcome_mod)
    }
    outcome_mod <- new_ipw_model(outcome_mod, vcov)
  }

  new_ipw(
    estimand = "ate",
    wt_mod = mods$wt_mod,
    outcome_mod = outcome_mod,
    estimates = binary_estimates(),
    se_method = "mestimation",
    fit = NULL,
    effects = "conditional"
  )
}

# ---- new_ipw() ---------------------------------------------------------------

test_that("new_ipw() builds the documented seven-field list", {
  estimates <- binary_estimates()
  fit <- structure(list(theta = c(1, 2)), class = "cg_mestimator")

  res <- new_ipw(
    estimand = "ate",
    wt_mod = "a propensity score model",
    outcome_mod = "an outcome model",
    estimates = estimates,
    se_method = "mestimation",
    fit = fit,
    effects = "conditional"
  )

  expect_type(res, "list")
  expect_s3_class(res, "ipw", exact = TRUE)
  expect_identical(
    names(res),
    c(
      "estimand",
      "wt_mod",
      "outcome_mod",
      "estimates",
      "se_method",
      "fit",
      "effects"
    )
  )

  expect_identical(res$estimand, "ate")
  expect_identical(res$wt_mod, "a propensity score model")
  expect_identical(res$outcome_mod, "an outcome model")
  expect_identical(res$estimates, estimates)
  expect_identical(res$se_method, "mestimation")
  expect_identical(res$fit, fit)
  expect_identical(res$effects, "conditional")
})

test_that("new_ipw() keeps a NULL fit as a named field", {
  # The linearization path has no M-estimator object to report. `fit` still has
  # to be present and `NULL`, because callers read the field by name and a list
  # that dropped it would have six elements rather than seven.
  res <- new_ipw(
    estimand = "att",
    wt_mod = NULL,
    outcome_mod = NULL,
    estimates = binary_estimates(),
    se_method = "linearization",
    fit = NULL
  )

  expect_length(res, 7L)
  expect_identical(
    names(res),
    c(
      "estimand",
      "wt_mod",
      "outcome_mod",
      "estimates",
      "se_method",
      "fit",
      "effects"
    )
  )
  expect_null(res$fit)
  expect_identical(res$se_method, "linearization")
})

test_that("new_ipw() takes its arguments in the documented order", {
  # Both packages downstream call this positionally in places, so the signature
  # order is part of the contract and not just the field order.
  estimates <- binary_estimates()

  positional <- new_ipw(
    "ate",
    "wt",
    "outcome",
    estimates,
    "mestimation",
    "fit",
    "conditional"
  )
  named <- new_ipw(
    estimand = "ate",
    wt_mod = "wt",
    outcome_mod = "outcome",
    estimates = estimates,
    se_method = "mestimation",
    fit = "fit",
    effects = "conditional"
  )

  expect_identical(positional, named)

  # `effects` is the last argument as well as the last field, so a call written
  # against the six-argument signature is still a call this one answers, and it
  # still means the marginal reading.
  earlier <- new_ipw("ate", "wt", "outcome", estimates, "mestimation", "fit")

  expect_identical(names(earlier), names(positional))
  expect_identical(earlier$effects, "marginal")
})

# ---- the effects mode --------------------------------------------------------

test_that("new_ipw() defaults the effects mode to marginal", {
  # A method that names no mode reports marginal effects, which is what every
  # method downstream computes today. The default is stored rather than left
  # absent, so that a caller reading `res$effects` gets the mode whichever way
  # the result was built.
  res <- new_ipw(
    estimand = "ate",
    wt_mod = "a propensity score model",
    outcome_mod = "an outcome model",
    estimates = binary_estimates(),
    se_method = "mestimation",
    fit = NULL
  )

  expect_identical(res$effects, "marginal")
  expect_identical(names(res)[[7]], "effects")
  expect_identical(res[[7]], "marginal")
})

test_that("new_ipw() stores either mode as given", {
  expect_identical(ipw_with_effects("marginal")$effects, "marginal")
  expect_identical(ipw_with_effects("conditional")$effects, "conditional")
})

test_that("new_ipw() rejects an effects value that is not a mode", {
  # The field has two readings and no third. A misspelling stored unchecked
  # would sit in the result until something downstream branched on the mode and
  # took the branch neither reading names.
  expect_error(
    ipw_with_effects("conditonal"),
    class = "causalgenerics_invalid_argument_effects"
  )
  expect_error(
    ipw_with_effects("conditonal"),
    class = "causalgenerics_invalid_argument"
  )
  expect_error(
    ipw_with_effects("everything"),
    class = "causalgenerics_invalid_argument_effects"
  )
  # The value is the string, so case is part of it.
  expect_error(
    ipw_with_effects("Marginal"),
    class = "causalgenerics_invalid_argument_effects"
  )
})

test_that("new_ipw() rejects an effects value that is not a single string", {
  # A factor prints as its label and a length-two vector holds a mode in its
  # first element, so both would pass a check that read the first element of
  # whatever it was given, and neither is a mode.
  expect_error(
    ipw_with_effects(factor("marginal")),
    class = "causalgenerics_invalid_argument_effects"
  )
  expect_error(
    ipw_with_effects(1),
    class = "causalgenerics_invalid_argument_effects"
  )
  expect_error(
    ipw_with_effects(TRUE),
    class = "causalgenerics_invalid_argument_effects"
  )
  expect_error(
    ipw_with_effects(c("marginal", "conditional")),
    class = "causalgenerics_invalid_argument_effects"
  )
  expect_error(
    ipw_with_effects(c("marginal", "marginal")),
    class = "causalgenerics_invalid_argument_effects"
  )
  # `NULL` is not a single string either, and a result that stored it would
  # record no mode at all.
  expect_error(
    ipw_with_effects(NULL),
    class = "causalgenerics_invalid_argument_effects"
  )
})

test_that("the effects error states the contract", {
  expect_snapshot(error = TRUE, ipw_with_effects("everything"))
  expect_snapshot(error = TRUE, ipw_with_effects(c("marginal", "conditional")))
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

# Whether `print()` writes a row labelled exactly `label`.
#
# The estimate columns are numeric, so a row label is a line whose remainder
# after the label is spaces and then a number. Testing only that a line starts
# with the label would hold for a prefix of the real one, and `printCoefmat()`
# also wraps significance stars onto lines of their own under the same labels
# when the frame is wide. Requiring a number keeps both out.
labels_a_printed_row <- function(out, label) {
  candidates <- out[startsWith(out, label)]
  remainder <- substring(candidates, nchar(label) + 1L)
  any(grepl("^ +-?[0-9]", remainder))
}

# The numbers `print()` writes on the row labelled `label`.
#
# Which columns a coefficient table carries is the implementation's to choose,
# and the estimate and its standard error are the first two whichever others are
# there, so reading the row by position is what lets an assertion name those two
# quantities without fixing the rest of the table. Tokens that are not numbers,
# such as significance stars and the `<` of a small p-value, are dropped rather
# than parsed. A row that is not there gives back nothing, so that the caller's
# assertions fail on an empty vector rather than on a subscript error.
printed_numbers <- function(out, label) {
  rows <- out[startsWith(out, label)]
  expect_length(rows, 1L)
  if (length(rows) != 1L) {
    return(numeric())
  }

  tokens <- strsplit(trimws(substring(rows, nchar(label) + 1L)), " +")[[1]]
  as.numeric(tokens[grepl("^-?[0-9]", tokens)])
}

test_that("print() summarizes a binary-exposure result", {
  res <- ipw_result(binary_estimates())

  expect_snapshot(print(res))

  out <- capture.output(print(res))

  expect_match(out, "^Inverse Probability Weight Estimator$", all = FALSE)
  # The estimand is reported in upper case, so `ate` reads as ATE.
  expect_match(out, "^Estimand: ATE\\s*$", all = FALSE)
  # The weight estimator is reported first, and each `Call:` line belongs to the
  # heading above it. The two fixtures differ in both formula and family, which
  # is what makes the two lines tell each other apart.
  wt_heading <- heading_index(out, "Weight Estimator:")
  outcome_heading <- heading_index(out, "Outcome Model:")

  expect_lt(wt_heading, outcome_heading)
  expect_match(out[wt_heading + 1L], "glm(formula = z ~ x", fixed = TRUE)
  expect_match(out[outcome_heading + 1L], "glm(formula = y ~ z", fixed = TRUE)

  expect_match(out, "^Marginal estimates:$", all = FALSE)

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

# ---- print.ipw() and the presentation mode -----------------------------------

# The two readings are different tables of different numbers, and a printed
# result that did not say which one it was showing would be read as the other.
# So the mode is reported twice over: once in the metadata block, where a reader
# looking for what the result is finds it beside the estimand, and once in the
# heading of the table itself, where a reader looking at the numbers finds it
# without scrolling back. The two are asserted separately below for that reason.

test_that("print() names the marginal reading in its metadata and its table", {
  res <- ipw_result(binary_estimates())

  out <- capture.output(print(res))

  # The mode is a fact about the result rather than about the table, so it is
  # reported with the estimand, and it is reported in both readings rather than
  # only in the one a reader might not expect.
  expect_match(
    out,
    "^Effects: marginal \\(population-averaged\\)\\s*$",
    all = FALSE
  )

  effects <- grep("^Effects:", out)
  expect_length(effects, 1L)
  expect_true(all(effects < heading_index(out, "Weight Estimator:")))

  # The heading names the reading the table reports. The unqualified heading is
  # gone rather than kept alongside: a result printed under it says nothing
  # about which of the two surfaces its numbers came from.
  expect_match(out, "^Marginal estimates:$", all = FALSE)
  expect_false(any(grepl("^Estimates:$", out)))

  # The table under that heading is the effects table it always was, keyed by
  # effect label and not by the outcome model's coefficient names.
  expect_true(labels_a_printed_row(out, "rd"))
  expect_true(labels_a_printed_row(out, "log(rr)"))
  expect_false(labels_a_printed_row(out, "(Intercept)"))
})

test_that("print() names the conditional reading in its metadata and its table", {
  res <- conditional_result()

  expect_snapshot(print(res))

  out <- capture.output(print(res))

  # The parenthetical is the implementation's to word, but it has to name the
  # outcome model: "conditional" on its own says the estimates are conditional
  # on something without saying on what.
  expect_match(
    out,
    "^Effects: conditional \\(.*outcome model.*\\)\\s*$",
    all = FALSE
  )

  effects <- grep("^Effects:", out)
  expect_length(effects, 1L)
  expect_true(all(effects < heading_index(out, "Weight Estimator:")))

  expect_match(out, "^Conditional estimates \\(outcome model\\):$", all = FALSE)
  expect_false(any(grepl("^Estimates:$", out)))
  expect_false(any(grepl("^Marginal estimates:$", out)))

  # Everything above the table is the same block in either reading: the
  # estimand, and the two component models with each call under its own
  # heading. The mode changes which surface is tabulated and nothing else about
  # what the result is.
  expect_match(out, "^Inverse Probability Weight Estimator$", all = FALSE)
  expect_match(out, "^Estimand: ATE\\s*$", all = FALSE)

  wt_heading <- heading_index(out, "Weight Estimator:")
  outcome_heading <- heading_index(out, "Outcome Model:")

  expect_lt(wt_heading, outcome_heading)
  expect_match(out[wt_heading + 1L], "glm(formula = z ~ x", fixed = TRUE)
  expect_match(out[outcome_heading + 1L], "glm(formula = y ~ z", fixed = TRUE)

  # The table is the outcome model's coefficient surface. The effects the other
  # reading reports are still on the object, and they are not what this reading
  # prints.
  expect_true(labels_a_printed_row(out, "(Intercept)"))
  expect_true(labels_a_printed_row(out, "z"))
  expect_false(labels_a_printed_row(out, "rd"))
  expect_false(labels_a_printed_row(out, "log(rr)"))
  expect_false(labels_a_printed_row(out, "log(or)"))
})

test_that("print() builds the conditional table from the corrected covariance", {
  # The estimate column is the outcome model's coefficients and the column
  # beside it is the standard error implied by the corrected block, which is
  # the covariance the joint estimation of the weights and the outcome gives.
  # The covariance the model computed for itself is the wrong answer rather
  # than an approximate one, and it is right there to be printed instead: it
  # treats the estimated weights as fixed and reports an uncertainty the
  # coefficients do not have.
  mod <- conditional_models()$outcome_mod
  res <- conditional_result()

  out <- capture.output(print(res))

  coefficients <- coef(mod)
  corrected_se <- sqrt(diag(corrected_outcome_vcov(mod)))
  own_se <- sqrt(diag(stats::vcov(mod)))

  intercept <- printed_numbers(out, "(Intercept)")
  slope <- printed_numbers(out, "z")

  expect_equal(
    intercept[1],
    unname(coefficients["(Intercept)"]),
    tolerance = 1e-4
  )
  expect_equal(slope[1], unname(coefficients["z"]), tolerance = 1e-4)
  expect_equal(
    intercept[2],
    unname(corrected_se["(Intercept)"]),
    tolerance = 1e-4
  )
  expect_equal(slope[2], unname(corrected_se["z"]), tolerance = 1e-4)

  # The fixture is discriminating: the two sets of standard errors are far
  # enough apart that the assertions above tell them apart at the printed
  # precision.
  expect_gt(abs(corrected_se[["z"]] - own_se[["z"]]), 1e-3)
})

test_that("print() pairs each coefficient with its own standard error", {
  # The block a fitting package attaches is labelled, and the labels are what
  # say which coefficient each variance belongs to. A table that took the
  # standard errors off the diagonal by position instead prints one
  # coefficient's uncertainty against another's and looks exactly like a
  # correct table: the column is where it belongs, the numbers are plausible,
  # and nothing in the output says they were crossed.
  res <- conditional_result(vcov = reversed_outcome_vcov())

  expect_snapshot(print(res))

  out <- capture.output(print(res))

  # The block's diagonal is 0.25 for the intercept and 0.0625 for the slope, so
  # the standard errors are 0.5 and 0.25. Crossed, they would be 0.25 and 0.5.
  expect_equal(printed_numbers(out, "(Intercept)")[2], 0.5, tolerance = 1e-4)
  expect_equal(printed_numbers(out, "z")[2], 0.25, tolerance = 1e-4)

  # The z statistic is the estimate over that standard error, so the rest of
  # the row moves with the pairing and reports a significance that belongs to
  # neither coefficient.
  coefficients <- coef(conditional_models()$outcome_mod)
  expect_equal(
    printed_numbers(out, "z")[3],
    unname(coefficients["z"]) / 0.25,
    tolerance = 1e-4
  )
  expect_equal(
    printed_numbers(out, "(Intercept)")[3],
    unname(coefficients["(Intercept)"]) / 0.5,
    tolerance = 1e-4
  )
})

test_that("print() reports a conditional result with no corrected covariance", {
  # An outcome model that was not wrapped carries no covariance from the joint
  # estimation, and `vcov()` refuses the reading for that reason. `print()` is
  # the summary of whatever the caller is holding, so refusing here would leave
  # a result that cannot be looked at at all. It reports the coefficients and
  # says what is missing instead.
  res <- conditional_result(wrap = FALSE)

  expect_error(vcov(res), class = "causalgenerics_no_conditional_vcov")
  expect_no_error(capture.output(print(res)))

  expect_snapshot(print(res))

  out <- capture.output(print(res))

  expect_match(out, "^Conditional estimates \\(outcome model\\):$", all = FALSE)
  expect_true(labels_a_printed_row(out, "(Intercept)"))
  expect_true(labels_a_printed_row(out, "z"))

  # The coefficients are the ones the model reports, whatever else the table
  # can carry with no covariance to compute it from.
  expect_equal(
    printed_numbers(out, "z")[1],
    unname(coef(conditional_models()$outcome_mod)["z"]),
    tolerance = 1e-4
  )

  # The standard errors the model computed for itself are not among the numbers
  # on the row. A column of them would satisfy everything above while reporting
  # the uncertainty the coefficients would have had if the weights had been
  # fixed rather than estimated, which is the reading this table exists to
  # correct, and a reader has no way to tell such a column from the corrected
  # one.
  own_se <- sqrt(diag(stats::vcov(conditional_models()$outcome_mod)))
  expect_false(any(abs(printed_numbers(out, "z") - own_se[["z"]]) < 1e-4))

  # The note states the fact. Its wording is the implementation's, but a reader
  # has to be told that the covariance is the missing part rather than left to
  # notice that a column is absent. It is part of the printed form rather than
  # a condition, so it travels with the output a caller captures.
  expect_true(any(grepl("covariance", out, ignore.case = TRUE)))
})

test_that("print() refuses a conditional result whose wrapper carries nothing", {
  # The test above is the fitting package that attached no block, which is a
  # result a caller can still look at: the note says what is missing and the
  # coefficients are printed without it. This is a different object. The outcome
  # model claims the wrapper class and carries no covariance behind it, which
  # `new_ipw_model()` cannot produce, so the note would tell a caller to wrap a
  # model that is already wrapped. It is reported as the error it is instead.
  # The handler `print()` puts around the covariance names the reading's
  # condition, and this one is the model's, so it travels past.
  #
  # `capture.output()` because the error is raised partway down the summary, and
  # the lines written before it belong in neither the test report nor the
  # assertion.
  res <- conditional_result()
  attr(res$outcome_mod, "ipw_vcov") <- NULL

  expect_error(
    capture.output(print(res)),
    class = "causalgenerics_no_vcov_ipw_model"
  )
  expect_error(capture.output(print(res)), class = "causalgenerics_no_vcov")
})

test_that("print() refuses a conditional block labelled for other parameters", {
  # A model that was never wrapped is reported rather than refused: the note
  # says what is missing and the coefficients are still worth looking at. This
  # is not that case. The block is there and cannot be read against these
  # coefficients, so there is no table to write and the note has nothing to
  # tell the caller, whose model is wrapped already. The condition travels past
  # the handler and out of `print()` the way a wrapper carrying nothing does.
  #
  # `capture.output()` because the error is raised partway down the summary,
  # and the lines written before it belong in neither the test report nor the
  # assertion.
  res <- conditional_result(vcov = theta_outcome_vcov())

  expect_error(
    capture.output(print(res)),
    class = "causalgenerics_conditional_vcov_mismatch"
  )
  expect_error(capture.output(print(res)), class = "causalgenerics_no_vcov")

  # The frame the condition names is the method the caller invoked, which is
  # what the accessors' conditions name too. The covariance is asked for inside
  # a `tryCatch()` here, so a call resolved when the helper's own default was
  # forced would name `doTryCatch()`: an internal frame of base R, in neither
  # this package nor the code the caller wrote.
  cnd <- tryCatch(capture.output(print(res)), error = identity)
  expect_identical(conditionCall(cnd), quote(print.ipw(res)))
})

test_that("print() says a conditional result has no coefficients to tabulate", {
  # An outcome model that reports no coefficients has no rows to tabulate under
  # this heading. Everything above the table is still the summary of a result
  # the caller is holding, so it is written, the heading names the reading, and
  # the reason the table is absent takes its place.
  res <- new_ipw(
    estimand = "ate",
    wt_mod = conditional_models()$wt_mod,
    outcome_mod = structure(list(), class = "cg_no_coef"),
    estimates = binary_estimates(),
    se_method = "mestimation",
    fit = NULL,
    effects = "conditional"
  )

  expect_no_error(capture.output(print(res)))

  expect_snapshot(print(res))

  out <- capture.output(print(res))

  expect_match(out, "^Inverse Probability Weight Estimator$", all = FALSE)
  expect_match(out, "^Conditional estimates \\(outcome model\\):$", all = FALSE)

  # The sentence itself, read out of the output with its line breaks collapsed.
  # Where the lines are broken is the implementation's to choose; that the
  # reader is told why the table is missing, in those words, is not.
  expect_match(
    paste(trimws(out), collapse = " "),
    paste0(
      "The outcome model reports no coefficients, so there is no ",
      "conditional table to print."
    ),
    fixed = TRUE
  )

  # The metadata block is the one both readings write, so a result with no
  # table is still a result a caller can look at.
  expect_match(out[heading_index(out, "Outcome Model:") + 1L], "cg_no_coef")
})

test_that("print() reports no coefficients without asking for a covariance", {
  # The guard is on the coefficients, and it runs before the covariance is
  # looked up. A wrapper carrying nothing raises an error of its own from that
  # lookup, and a result with no rows to tabulate has no use for a covariance
  # either way, so reaching for one would refuse a result over a part of it
  # that nothing was going to read.
  stripped <- new_ipw(
    estimand = "ate",
    wt_mod = conditional_models()$wt_mod,
    outcome_mod = structure(list(), class = c("ipw_model", "cg_no_coef")),
    estimates = binary_estimates(),
    se_method = "mestimation",
    fit = NULL,
    effects = "conditional"
  )

  expect_no_error(capture.output(print(stripped)))

  # `coef()` answers `NULL` for the fixture above, through the default method
  # reading a list element that is not there. A model with a method of its own
  # reports the same fact as an empty numeric vector, and both are a model with
  # no coefficients, so the guard is on the length rather than on the value.
  local_s3_method("coef", "cg_empty_coef", function(object, ...) numeric())
  empty <- new_ipw(
    estimand = "ate",
    wt_mod = conditional_models()$wt_mod,
    outcome_mod = structure(list(), class = "cg_empty_coef"),
    estimates = binary_estimates(),
    se_method = "mestimation",
    fit = NULL,
    effects = "conditional"
  )

  expect_no_error(capture.output(print(empty)))
  expect_match(
    paste(trimws(capture.output(print(empty))), collapse = " "),
    paste0(
      "The outcome model reports no coefficients, so there is no ",
      "conditional table to print."
    ),
    fixed = TRUE
  )
})

test_that("print() reads a result with no mode as marginal", {
  # A result stored before the field existed carries six fields. Marginal is
  # the reading every method produced then, so it is the one printed, and what
  # is printed is the seven-field marginal result's form exactly rather than a
  # third thing that says the mode some other way.
  legacy <- legacy_result()

  expect_length(legacy, 6L)
  expect_null(legacy$effects)

  out <- capture.output(print(legacy))

  expect_match(
    out,
    "^Effects: marginal \\(population-averaged\\)\\s*$",
    all = FALSE
  )
  expect_match(out, "^Marginal estimates:$", all = FALSE)
  expect_identical(out, capture.output(print(ipw_result(binary_estimates()))))
})

test_that("print() refuses a stored mode that is not a reading", {
  # The field is set by the constructor, `as_marginal()`, and
  # `as_conditional()`, none of which can put anything but a reading there. A
  # field that says something else was assigned to directly, and there is no
  # reading to fall back on: printing the marginal table would present a result
  # that names neither reading as though it named that one. The accessors
  # refuse the same value, so a caller handles one class wherever it came from.
  res <- ipw_result(binary_estimates())
  res$effects <- "banana"

  expect_error(print(res), class = "causalgenerics_invalid_argument_effects")
  expect_error(print(res), class = "causalgenerics_invalid_argument")

  # The supported route back is the one that put a reading there in the first
  # place, and taking it leaves nothing behind to print around.
  expect_no_error(capture.output(print(as_marginal(res))))

  expect_snapshot(error = TRUE, print(res))
})

test_that("print() returns its input invisibly in the conditional mode", {
  res <- conditional_result()

  # The conditional table is what ran, so the invisible return asserted below is
  # that branch's rather than one taken by a `print()` that never read the mode.
  expect_match(
    capture.output(print(res)),
    "^Conditional estimates \\(outcome model\\):$",
    all = FALSE
  )

  # The printed form is snapshotted above, so from here the output only needs to
  # go somewhere other than the test report.
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
    wt_mod = structure(1:3, class = "cg_unsubsettable"),
    outcome_mod = structure(list(), class = "cg_no_call"),
    estimates = binary_estimates(),
    se_method = "mestimation",
    fit = NULL
  )

  out <- capture.output(print(res))

  wt_heading <- heading_index(out, "Weight Estimator:")
  outcome_heading <- heading_index(out, "Outcome Model:")

  expect_match(out[wt_heading + 1L], "Call: <cg_unsubsettable>", fixed = TRUE)
  expect_match(out[outcome_heading + 1L], "Call: <cg_no_call>", fixed = TRUE)
})

# ---- as.data.frame.ipw() -----------------------------------------------------

# `as.data.frame()` is the result's tidier-shaped table rather than a copy of the
# `estimates` frame. The storage names of that frame are the constructor's
# contract with the packages that build a result; the names below are the
# contract with the packages that report one, and they are the column names a
# tidier uses. A fitting package's `tidy()` method is `as_tibble()` over this
# table and nothing else, so what that method returns is settled here rather than
# restated in every package that has one.
#
# The interval is the one part of the table this method computes rather than
# reads. At the level the frame stores, the stored bounds come back as they
# stand: they need not be the normal pair at all, a bootstrap or profile interval
# being asymmetric about the estimate, and even a normal one rounded on its way
# into the frame is not the number recomputing gives. At any other level the
# normal approximation is built from the estimate and its standard error.

# The bounds the method builds at a level the frame does not store: one half
# width from the normal approximation, added and subtracted. `qnorm()` is not
# exactly antisymmetric, so taking the lower bound from the lower tail instead
# would put the two a little apart from each other.
normal_bounds <- function(estimates, level) {
  half_width <- stats::qnorm(1 - (1 - level) / 2) * estimates$std.err
  list(
    lower = estimates$estimate - half_width,
    upper = estimates$estimate + half_width
  )
}

# A binary-exposure frame whose rows do not share one confidence level. A frame
# like this is what tells a rule reading the stored level row by row from one
# reading it as a property of the frame, and the two disagree about the rows
# whose stored level is the one asked for.
mixed_level_estimates <- function() {
  estimates <- binary_estimates()
  estimates$conf.level <- c(0.95, 0.90, 0.95)
  estimates
}

# A binary-exposure frame whose ratios are already on the natural scale, which is
# what separates matching the `log(rr)` and `log(or)` labels exactly from
# matching every label holding `rr` or `or` somewhere. Only the labels are
# changed: what the numbers under them mean is beside the point for a test about
# which rows are picked out.
exponentiated_estimates <- function() {
  estimates <- binary_estimates()
  estimates$effect <- c("rd", "rr", "or")
  estimates
}

# The covariance of the binary fixture's three effects, in the shape a method
# attaches it: a square matrix labelled on both margins by effect label, with
# off-diagonal entries, since effects computed from the same weighted means are
# correlated.
binary_vcov <- function() {
  labels <- binary_estimates()$effect
  matrix(
    c(
      0.008542,
      0.022753,
      0.034812,
      0.022753,
      0.074813,
      0.113106,
      0.034812,
      0.113106,
      0.175277
    ),
    nrow = 3,
    dimnames = list(labels, labels)
  )
}

# An estimates frame carrying the covariance a method attaches to it.
with_vcov <- function(estimates, covariance = binary_vcov()) {
  attr(estimates, "ipw_vcov") <- covariance
  estimates
}

test_that("as.data.frame() renames the estimates into the tidier columns", {
  estimates <- binary_estimates()
  res <- ipw_result(estimates)

  df <- as.data.frame(res)

  expect_s3_class(df, "data.frame", exact = TRUE)
  expect_identical(
    names(df),
    c("term", "estimate", "std.error", "statistic", "p.value")
  )

  # The label of a row is its term, and the label is what the frame's `effect`
  # column holds. The comparison a categorical result reports is a column of its
  # own rather than part of this one.
  expect_type(df$term, "character")
  expect_identical(df$term, c("rd", "log(rr)", "log(or)"))

  # Renamed and selected, not recomputed. The statistic is the frame's `z` and
  # the p-value is the frame's own, so a result reporting a statistic that is
  # not the estimate over its standard error still reports that statistic here.
  expect_identical(df$estimate, estimates$estimate)
  expect_identical(df$std.error, estimates$std.err)
  expect_identical(df$statistic, estimates$z)
  expect_identical(df$p.value, estimates$p.value)

  # The storage names are gone rather than kept alongside. A table carrying both
  # spellings of one quantity would let a caller read the same column twice
  # under two names and would not stack with a tidier's output.
  expect_false(
    any(
      c("effect", "std.err", "z", "ci.lower", "ci.upper", "conf.level") %in%
        names(df)
    )
  )
})

test_that("as.data.frame() puts comparison after term for a categorical result", {
  # A categorical exposure reports each effect measure once per contrast, so the
  # term alone names several rows the same thing. The column naming the contrast
  # follows the term it qualifies rather than sitting among the numbers.
  estimates <- categorical_estimates()
  res <- ipw_result(estimates)

  df <- as.data.frame(res)

  expect_identical(
    names(df),
    c("term", "comparison", "estimate", "std.error", "statistic", "p.value")
  )
  expect_identical(df$term, rep(c("rd", "log(rr)", "log(or)"), times = 2))
  expect_identical(df$comparison, estimates$comparison)
  expect_identical(df$estimate, estimates$estimate)
  expect_identical(df$std.error, estimates$std.err)
})

test_that("as.data.frame() omits comparison when the result reports none", {
  # A binary or continuous exposure has one contrast, so there is nothing for a
  # comparison column to say. It is absent rather than present and constant: a
  # column of one repeated value would stack with a categorical result's table
  # and read as a contrast that was named.
  df <- as.data.frame(ipw_result(continuous_estimates()))

  expect_identical(
    names(df),
    c("term", "estimate", "std.error", "statistic", "p.value")
  )
  expect_identical(nrow(df), 1L)
  expect_identical(df$term, "diff")
  expect_identical(df$estimate, continuous_estimates()$estimate)
})

test_that("as.data.frame() reports no interval and no level by default", {
  # The bounds are opt-in, and the level is an argument rather than a column.
  # A `conf.level` column would repeat one number down every row and would be
  # read as a column of the table rather than as the level the two bounds beside
  # it were built at.
  df <- as.data.frame(ipw_result(binary_estimates()))

  expect_false(any(c("conf.low", "conf.high", "conf.level") %in% names(df)))
})

test_that("as.data.frame(conf.int = TRUE) appends the bounds last", {
  estimates <- binary_estimates()
  res <- ipw_result(estimates)

  df <- as.data.frame(res, conf.int = TRUE)

  expect_identical(
    names(df),
    c(
      "term",
      "estimate",
      "std.error",
      "statistic",
      "p.value",
      "conf.low",
      "conf.high"
    )
  )

  # The level the bounds were built at is not a column of the table.
  expect_false("conf.level" %in% names(df))
})

test_that("as.data.frame(conf.int = TRUE) keeps comparison after term", {
  # The bounds go on the end, so the categorical table adds them after the same
  # five columns the binary one has and leaves the contrast where it was.
  df <- as.data.frame(ipw_result(categorical_estimates()), conf.int = TRUE)

  expect_identical(
    names(df),
    c(
      "term",
      "comparison",
      "estimate",
      "std.error",
      "statistic",
      "p.value",
      "conf.low",
      "conf.high"
    )
  )
  expect_identical(df$comparison, categorical_estimates()$comparison)
})

test_that("as.data.frame(conf.int = TRUE) returns the stored bounds", {
  # At the level the frame stores, the frame's own bounds come back. They need
  # not be the normal pair, and the fixture is discriminating on that point: its
  # bounds are rounded to six places, so recomputing them from the estimate and
  # the standard error gives numbers that agree to the eye and differ in the
  # digits an assertion reads.
  estimates <- binary_estimates()
  res <- ipw_result(estimates)

  df <- as.data.frame(res, conf.int = TRUE, conf.level = 0.95)

  expect_identical(df$conf.low, estimates$ci.lower)
  expect_identical(df$conf.high, estimates$ci.upper)

  recomputed <- normal_bounds(estimates, 0.95)
  expect_false(any(df$conf.low == recomputed$lower))
  expect_false(any(df$conf.high == recomputed$upper))
})

test_that("as.data.frame(conf.int = TRUE) recomputes at another level", {
  # A level the frame does not store has no bounds to return, so the normal
  # approximation is built from the estimate and its standard error.
  estimates <- binary_estimates()
  res <- ipw_result(estimates)

  df <- as.data.frame(res, conf.int = TRUE, conf.level = 0.9)

  expected <- normal_bounds(estimates, 0.9)
  expect_identical(df$conf.low, expected$lower)
  expect_identical(df$conf.high, expected$upper)

  # The stored bounds describe the level the frame stores, and reporting them
  # under another one would label a 95% interval as a 90% one.
  expect_false(any(df$conf.low == estimates$ci.lower))
  expect_false(any(df$conf.high == estimates$ci.upper))

  # The rest of the table is what it is without an interval, so asking for one
  # changes nothing but the two columns it adds.
  expect_identical(
    df[names(as.data.frame(res))],
    as.data.frame(res)
  )
})

test_that("as.data.frame(conf.int = TRUE) recomputes every row when the stored levels differ", {
  # The stored level is a property of the frame rather than of a row: the bounds
  # come back as they stand only when every row was reported at the level asked
  # for. A frame whose rows disagree has no one level to return bounds for, and
  # a table mixing stored bounds with recomputed ones would report two intervals
  # under one column heading.
  estimates <- mixed_level_estimates()
  res <- ipw_result(estimates)

  df <- as.data.frame(res, conf.int = TRUE, conf.level = 0.95)

  expected <- normal_bounds(estimates, 0.95)
  expect_identical(df$conf.low, expected$lower)
  expect_identical(df$conf.high, expected$upper)

  # The two rows whose stored level is the one asked for are recomputed with the
  # rest, which is what makes this a claim about the frame and not about a row.
  expect_false(any(df$conf.low[c(1, 3)] == estimates$ci.lower[c(1, 3)]))
})

test_that("as.data.frame(conf.int = TRUE) recomputes when the frame stores no level", {
  # A frame with no `conf.level` column says nothing about what its bounds
  # describe, so there is no level for the stored pair to be returned at and the
  # normal approximation is built for whatever was asked for.
  estimates <- binary_estimates()
  estimates$conf.level <- NULL
  res <- ipw_result(estimates)

  df <- as.data.frame(res, conf.int = TRUE, conf.level = 0.95)

  expected <- normal_bounds(estimates, 0.95)
  expect_identical(df$conf.low, expected$lower)
  expect_identical(df$conf.high, expected$upper)
})

test_that("as.data.frame(exponentiate = TRUE) moves only the ratio rows", {
  estimates <- binary_estimates()
  res <- ipw_result(estimates)

  df <- as.data.frame(res, exponentiate = TRUE)

  # The two log rows move to the natural scale and their terms move with them,
  # since the label names the scale. The risk difference was never on the log
  # scale, so it is untouched.
  expect_identical(df$term, c("rd", "rr", "or"))
  expect_identical(
    df$estimate,
    c(estimates$estimate[1], exp(estimates$estimate[2:3]))
  )

  # Inference stays on the log scale, where it is done. Exponentiating a
  # standard error would state a quantity nobody asked for, and the statistic
  # and the p-value describe the log scale estimate.
  expect_identical(df$std.error, estimates$std.err)
  expect_identical(df$statistic, estimates$z)
  expect_identical(df$p.value, estimates$p.value)

  # The columns are the ones the table always has: exponentiating changes what
  # the numbers mean and not which columns report them.
  expect_identical(names(df), names(as.data.frame(res)))
})

test_that("as.data.frame(exponentiate = TRUE) moves the bounds it reports", {
  estimates <- binary_estimates()
  res <- ipw_result(estimates)

  stored <- as.data.frame(res, conf.int = TRUE, exponentiate = TRUE)

  expect_identical(
    stored$conf.low,
    c(estimates$ci.lower[1], exp(estimates$ci.lower[2:3]))
  )
  expect_identical(
    stored$conf.high,
    c(estimates$ci.upper[1], exp(estimates$ci.upper[2:3]))
  )

  # The interval is settled before the scale is: the bounds a level the frame
  # does not store gets are built on the log scale, where the standard error
  # lives, and exponentiated afterwards. Building them from the exponentiated
  # estimate instead would add a half width on the log scale to a number on the
  # natural one.
  recomputed <- as.data.frame(
    res,
    conf.int = TRUE,
    conf.level = 0.9,
    exponentiate = TRUE
  )
  expected <- normal_bounds(estimates, 0.9)

  expect_identical(
    recomputed$conf.low,
    c(expected$lower[1], exp(expected$lower[2:3]))
  )
  expect_identical(
    recomputed$conf.high,
    c(expected$upper[1], exp(expected$upper[2:3]))
  )
})

test_that("as.data.frame(exponentiate = TRUE) relabels every comparison", {
  estimates <- categorical_estimates()
  res <- ipw_result(estimates)

  df <- as.data.frame(res, conf.int = TRUE, exponentiate = TRUE)

  expect_identical(df$term, rep(c("rd", "rr", "or"), times = 2))

  # The comparison labels identify the contrast, not the scale, so they survive
  # untouched and stay in place after the term.
  expect_identical(df$comparison, estimates$comparison)
  expect_identical(names(df)[1:2], c("term", "comparison"))

  ratios <- df$term %in% c("rr", "or")
  expect_identical(df$estimate[ratios], exp(estimates$estimate[ratios]))
  expect_identical(df$conf.low[ratios], exp(estimates$ci.lower[ratios]))
  expect_identical(df$conf.high[ratios], exp(estimates$ci.upper[ratios]))
  expect_identical(df$std.error, estimates$std.err)
})

test_that("as.data.frame(exponentiate = TRUE) matches the log labels exactly", {
  # The rows to move are the ones labelled `log(rr)` and `log(or)`, not every
  # row whose label happens to contain `rr` or `or`. A frame whose ratios are
  # already on the natural scale is the fixture that separates the two: its
  # labels are `rr` and `or`, so a substring match would exponentiate it a
  # second time while an equality match leaves it alone.
  res <- ipw_result(exponentiated_estimates())

  df <- as.data.frame(res, conf.int = TRUE, exponentiate = TRUE)

  expect_identical(df$term, c("rd", "rr", "or"))
  expect_identical(df, as.data.frame(res, conf.int = TRUE))
})

test_that("as.data.frame(exponentiate = TRUE) is a no-op with no ratio rows", {
  # A continuous outcome reports a difference in means and nothing else, so
  # asking for the natural scale has to leave the table alone rather than error
  # on the absent rows.
  res <- ipw_result(continuous_estimates())

  df <- as.data.frame(res, conf.int = TRUE, exponentiate = TRUE)

  expect_identical(df$term, "diff")
  expect_identical(df, as.data.frame(res, conf.int = TRUE))
})

test_that("as.data.frame() carries the covariance the result stores", {
  # The covariance belongs to the estimates rather than to a column of them, so
  # it travels on the table as the attribute the `new_ipw()` contract names it
  # under. A caller who reads the table and then asks it for the covariance of
  # the rows it holds gets the matrix the fitting package attached.
  covariance <- binary_vcov()
  estimates <- with_vcov(binary_estimates(), covariance)
  res <- ipw_result(estimates)

  df <- as.data.frame(res)

  expect_identical(
    names(df),
    c("term", "estimate", "std.error", "statistic", "p.value")
  )
  expect_identical(attr(df, "ipw_vcov", exact = TRUE), covariance)

  # Asking for the interval adds two columns and says nothing about the
  # covariance of the rows.
  expect_identical(
    attr(as.data.frame(res, conf.int = TRUE), "ipw_vcov", exact = TRUE),
    covariance
  )

  # Including at a level the frame does not store, where the bounds are built
  # here rather than read off. Recomputing an interval says nothing about the
  # covariance of the rows either: the estimates are the ones the result holds
  # whichever level the bounds beside them report.
  expect_identical(
    attr(
      as.data.frame(res, conf.int = TRUE, conf.level = 0.9),
      "ipw_vcov",
      exact = TRUE
    ),
    covariance
  )

  # A result whose method attached none reports none rather than a matrix built
  # from the standard errors, which would report the effects as uncorrelated.
  expect_null(attr(as.data.frame(ipw_result(binary_estimates())), "ipw_vcov"))
})

test_that("as.data.frame(exponentiate = TRUE) drops the covariance", {
  # The covariance describes the estimates on the scale they were estimated on.
  # Exponentiating two of the rows leaves it describing neither the table it is
  # attached to nor anything else, and a caller has no way to tell such a matrix
  # from one that still matches. There is no correct matrix to put in its place
  # either, since the delta method answer is not what the fitting package
  # computed, so it is dropped.
  covariance <- binary_vcov()
  res <- ipw_result(with_vcov(binary_estimates(), covariance))

  expect_null(attr(as.data.frame(res, exponentiate = TRUE), "ipw_vcov"))
  expect_null(
    attr(
      as.data.frame(res, conf.int = TRUE, exponentiate = TRUE),
      "ipw_vcov"
    )
  )

  # The result is the caller's own and reading a table off it does not change
  # it, whichever scale the table was asked for.
  expect_identical(attr(res$estimates, "ipw_vcov", exact = TRUE), covariance)
})

test_that("as.data.frame() returns the same table in either mode", {
  # The mode says which surface `print()` and the accessors report. The table is
  # built from the estimates the result stores, which both readings hold, so
  # flipping a result to read its coefficients does not change what
  # `as.data.frame()` means. A caller who does not want the effects table asks
  # the outcome model for its own.
  marginal <- ipw_result(binary_estimates())
  conditional <- as_conditional(marginal)

  expect_identical(conditional$effects, "conditional")
  expect_identical(as.data.frame(conditional), as.data.frame(marginal))
  expect_identical(
    as.data.frame(conditional, conf.int = TRUE, conf.level = 0.9),
    as.data.frame(marginal, conf.int = TRUE, conf.level = 0.9)
  )
  expect_identical(
    as.data.frame(conditional, exponentiate = TRUE),
    as.data.frame(marginal, exponentiate = TRUE)
  )
})

test_that("as.data.frame() returns the same table for a result with no mode", {
  # A result stored before the field existed carries six fields. The table is
  # the estimates it holds either way, so it is the seven-field result's table
  # exactly rather than a third thing.
  legacy <- legacy_result()

  expect_length(legacy, 6L)
  expect_null(legacy$effects)

  marginal <- ipw_result(binary_estimates())

  expect_identical(as.data.frame(legacy), as.data.frame(marginal))
  expect_identical(
    as.data.frame(legacy, conf.int = TRUE, conf.level = 0.9),
    as.data.frame(marginal, conf.int = TRUE, conf.level = 0.9)
  )
  expect_identical(
    as.data.frame(legacy, conf.int = TRUE, exponentiate = TRUE),
    as.data.frame(marginal, conf.int = TRUE, exponentiate = TRUE)
  )
})

test_that("as.data.frame() sets the row names from row.names", {
  res <- ipw_result(binary_estimates())

  df <- as.data.frame(res, row.names = c("first", "second", "third"))

  expect_identical(rownames(df), c("first", "second", "third"))
  expect_identical(
    df[names(as.data.frame(res))],
    `rownames<-`(as.data.frame(res), c("first", "second", "third"))
  )
})

test_that("as.data.frame() takes row.names and optional positionally", {
  # The generic in \pkg{base} takes `row.names` and `optional` in that order
  # before its dots, and this method matches it, so a positional call written
  # against the generic means here what it means there. The three arguments this
  # method adds sit after the dots and are named at every call site, which is
  # what keeps a second positional argument from being read as one of them.
  res <- ipw_result(binary_estimates())

  expect_identical(
    as.data.frame(res, c("first", "second", "third")),
    as.data.frame(res, row.names = c("first", "second", "third"))
  )

  # Every column this method builds is named, so there is nothing for `optional`
  # to make optional. It is accepted for the generic's sake and changes nothing.
  expect_identical(as.data.frame(res, NULL, TRUE), as.data.frame(res))
  expect_identical(
    as.data.frame(res, optional = TRUE, conf.int = TRUE),
    as.data.frame(res, conf.int = TRUE)
  )

  # Partial matching reaches only the arguments before the dots, so an
  # abbreviation of one that sits after them is swallowed by the dots rather
  # than matched. `expon` is therefore not `exponentiate`, and a table asked for
  # that way is the table asked for with no arguments at all.
  expect_identical(as.data.frame(res, expon = TRUE), as.data.frame(res))
})

test_that("as.data.frame() ignores further arguments", {
  # The dots are the generic's, and this method reads nothing out of them. An
  # argument it does not know is neither an error nor a column of the table.
  res <- ipw_result(binary_estimates())

  expect_identical(
    as.data.frame(res, cg_unused = "ignored"),
    as.data.frame(res)
  )
})

test_that("as.data.frame() refuses a conf.level that is not a probability", {
  res <- ipw_result(binary_estimates())

  expect_error(
    as.data.frame(res, conf.int = TRUE, conf.level = 1),
    class = "causalgenerics_invalid_argument_conf.level"
  )
  expect_error(
    as.data.frame(res, conf.int = TRUE, conf.level = 1),
    class = "causalgenerics_invalid_argument"
  )

  # The bounds of the open interval are refused along with everything outside
  # it: an interval at level 0 or 1 is not one the normal approximation gives a
  # finite pair of numbers for.
  expect_error(
    as.data.frame(res, conf.int = TRUE, conf.level = 0),
    class = "causalgenerics_invalid_argument_conf.level"
  )
  expect_error(
    as.data.frame(res, conf.int = TRUE, conf.level = -0.1),
    class = "causalgenerics_invalid_argument_conf.level"
  )
  expect_error(
    as.data.frame(res, conf.int = TRUE, conf.level = 95),
    class = "causalgenerics_invalid_argument_conf.level"
  )

  # A missing level names no interval, and a vector of levels names several
  # where the table reports one.
  expect_error(
    as.data.frame(res, conf.int = TRUE, conf.level = NA_real_),
    class = "causalgenerics_invalid_argument_conf.level"
  )
  expect_error(
    as.data.frame(res, conf.int = TRUE, conf.level = c(0.9, 0.95)),
    class = "causalgenerics_invalid_argument_conf.level"
  )
  expect_error(
    as.data.frame(res, conf.int = TRUE, conf.level = numeric()),
    class = "causalgenerics_invalid_argument_conf.level"
  )
  expect_error(
    as.data.frame(res, conf.int = TRUE, conf.level = NULL),
    class = "causalgenerics_invalid_argument_conf.level"
  )

  # A level is a number. A string that reads as one is not, and neither is the
  # logical a comparison would produce.
  expect_error(
    as.data.frame(res, conf.int = TRUE, conf.level = "0.95"),
    class = "causalgenerics_invalid_argument_conf.level"
  )
  expect_error(
    as.data.frame(res, conf.int = TRUE, conf.level = TRUE),
    class = "causalgenerics_invalid_argument_conf.level"
  )
})

test_that("as.data.frame() checks conf.level with no interval asked for", {
  # The argument means the same thing whether or not the bounds are reported,
  # and a call that names a level the method cannot build an interval at has
  # said something wrong. Checking only on the branch that reads it would accept
  # the call now and refuse the same value the moment `conf.int` was turned on.
  res <- ipw_result(binary_estimates())

  expect_error(
    as.data.frame(res, conf.level = 95),
    class = "causalgenerics_invalid_argument_conf.level"
  )
  expect_error(
    as.data.frame(res, conf.int = FALSE, conf.level = NA_real_),
    class = "causalgenerics_invalid_argument_conf.level"
  )
})

test_that("as.data.frame() refuses a conf.int that is not a flag", {
  res <- ipw_result(binary_estimates())

  expect_error(
    as.data.frame(res, conf.int = NA),
    class = "causalgenerics_invalid_argument_conf.int"
  )
  expect_error(
    as.data.frame(res, conf.int = NA),
    class = "causalgenerics_invalid_argument"
  )
  expect_error(
    as.data.frame(res, conf.int = c(TRUE, FALSE)),
    class = "causalgenerics_invalid_argument_conf.int"
  )
  expect_error(
    as.data.frame(res, conf.int = logical()),
    class = "causalgenerics_invalid_argument_conf.int"
  )
  expect_error(
    as.data.frame(res, conf.int = NULL),
    class = "causalgenerics_invalid_argument_conf.int"
  )

  # A number and a string are not flags. Both would take the true branch of a
  # bare `if`, so a method that did not check them would report an interval for
  # `conf.int = "no"`.
  expect_error(
    as.data.frame(res, conf.int = 1),
    class = "causalgenerics_invalid_argument_conf.int"
  )
  expect_error(
    as.data.frame(res, conf.int = "TRUE"),
    class = "causalgenerics_invalid_argument_conf.int"
  )
})

test_that("as.data.frame() refuses an exponentiate that is not a flag", {
  res <- ipw_result(binary_estimates())

  expect_error(
    as.data.frame(res, exponentiate = NA),
    class = "causalgenerics_invalid_argument_exponentiate"
  )
  expect_error(
    as.data.frame(res, exponentiate = NA),
    class = "causalgenerics_invalid_argument"
  )
  expect_error(
    as.data.frame(res, exponentiate = c(TRUE, TRUE)),
    class = "causalgenerics_invalid_argument_exponentiate"
  )
  expect_error(
    as.data.frame(res, exponentiate = NULL),
    class = "causalgenerics_invalid_argument_exponentiate"
  )
  expect_error(
    as.data.frame(res, exponentiate = 1),
    class = "causalgenerics_invalid_argument_exponentiate"
  )
  expect_error(
    as.data.frame(res, exponentiate = "yes"),
    class = "causalgenerics_invalid_argument_exponentiate"
  )
})

test_that("the as.data.frame() argument errors state the contract", {
  # One value per argument. The message names the argument and the contract it
  # broke rather than the value that broke it, so a second value of the same
  # argument would record the same sentence again.
  #
  # Each value below is one the method reports on rather than one base R trips
  # over on its way past. A logical argument that reaches an `if` unchecked
  # errors from the `if` itself, and a snapshot of that message would record the
  # absence of the check as though it were the check.
  res <- ipw_result(binary_estimates())

  expect_snapshot(error = TRUE, as.data.frame(res, conf.level = 95))
  expect_snapshot(error = TRUE, as.data.frame(res, conf.int = NA))
  expect_snapshot(error = TRUE, as.data.frame(res, exponentiate = 1))
})

# ---- the contrast column -----------------------------------------------------

# The column a categorical result names its contrasts with is `contrast`, in the
# `estimates` frame the constructor takes and in the tidier-shaped table
# `as.data.frame()` reports. Two things downstream read that table by column name
# and both spell it that way. `mice::pool()` groups the rows it pools by a closed
# whitelist of column names that holds `contrast` and not `comparison`, so a
# table under the other name loses the contrast it was keyed by. And
# \pkg{marginaleffects} names the same column `contrast`, so a table under the
# other name does not stack with one of theirs.
#
# The labels the rows carry are not the column. A label is `effect` and the
# contrast pasted together, such as `"rd b vs a"`, and those strings are what
# `print()` writes down the side of its table and what `coef()`, `vcov()`, and
# `confint()` name their results with. The column the labels are built from is
# renamed and the labels themselves are not, so every one of those surfaces reads
# the same before and after. That is the second half of what this section pins,
# and it is the half a rename is most likely to break quietly.

# The labels the categorical fixture's rows carry, written out rather than pasted
# together. The label rule is what these assertions are for, so restating it with
# the `paste()` the implementation uses would assert nothing.
contrast_labels <- function() {
  c(
    "rd b vs a",
    "log(rr) b vs a",
    "log(or) b vs a",
    "rd c vs a",
    "log(rr) c vs a",
    "log(or) c vs a"
  )
}

# The covariance of those six effects, in the shape a method attaches it: a
# square matrix labelled on both margins by effect label, with off-diagonal
# entries, since effects computed from the same weighted means are correlated.
contrast_vcov <- function() {
  std_err <- contrast_estimates()$std.err
  index <- seq_along(std_err)
  covariance <- 0.9^abs(outer(index, index, "-")) * outer(std_err, std_err)
  dimnames(covariance) <- list(contrast_labels(), contrast_labels())
  covariance
}

test_that("as.data.frame() puts contrast after term for a categorical result", {
  # A categorical exposure reports each effect measure once per contrast, so the
  # term alone names several rows the same thing. The column naming the contrast
  # follows the term it qualifies rather than sitting among the numbers, and it
  # is the only column that names it: a table carrying two spellings of one
  # quantity would let a caller read the same column twice under two names.
  estimates <- contrast_estimates()
  res <- ipw_result(estimates)

  df <- as.data.frame(res)

  expect_identical(
    names(df),
    c("term", "contrast", "estimate", "std.error", "statistic", "p.value")
  )
  expect_identical(names(df)[1:2], c("term", "contrast"))
  expect_false("comparison" %in% names(df))

  expect_identical(df$contrast, estimates$contrast)
  expect_identical(df$term, rep(c("rd", "log(rr)", "log(or)"), times = 2))
  expect_identical(df$estimate, estimates$estimate)
  expect_identical(df$std.error, estimates$std.err)
})

test_that("as.data.frame(conf.int = TRUE) keeps contrast after term", {
  # The bounds go on the end, so the categorical table adds them after the same
  # six columns and leaves the contrast where it was.
  estimates <- contrast_estimates()

  df <- as.data.frame(ipw_result(estimates), conf.int = TRUE)

  expect_identical(
    names(df),
    c(
      "term",
      "contrast",
      "estimate",
      "std.error",
      "statistic",
      "p.value",
      "conf.low",
      "conf.high"
    )
  )
  expect_identical(df$contrast, estimates$contrast)
  expect_false("comparison" %in% names(df))
})

test_that("as.data.frame(exponentiate = TRUE) keeps every contrast label", {
  # The contrast labels identify the contrast, not the scale, so they survive
  # untouched and stay in place after the term while the two ratio rows move.
  estimates <- contrast_estimates()

  df <- as.data.frame(
    ipw_result(estimates),
    conf.int = TRUE,
    exponentiate = TRUE
  )

  expect_identical(df$term, rep(c("rd", "rr", "or"), times = 2))
  expect_identical(df$contrast, estimates$contrast)
  expect_identical(names(df)[1:2], c("term", "contrast"))
})

test_that("print() keys rows by effect and contrast", {
  # Row labels are `effect` and `contrast` together, and both character columns
  # have to be gone before `printCoefmat()` sees the frame. Leaving the contrast
  # column in is the failure worth guarding against, because it does not error:
  # `printCoefmat()` sends a character column through `data.matrix()`, which
  # factor-codes it into a column reading `1.000000` and `2.000000` alongside the
  # real estimates, under its own heading.
  res <- ipw_result(contrast_estimates())

  out <- capture.output(print(res))

  expect_match(out, "^rd b vs a +0\\.081945 ", all = FALSE)
  expect_match(out, "^log\\(rr\\) b vs a +0\\.168870 ", all = FALSE)
  expect_match(out, "^log\\(or\\) c vs a +0\\.676435 ", all = FALSE)

  expect_false(any(grepl("contrast", out, fixed = TRUE)))
  expect_false(any(grepl("effect", out, fixed = TRUE)))
})

test_that("the contrast column leaves the effect labels as they were", {
  # Renaming the column renames no row. `"rd b vs a"` is the line `print()`
  # writes, the name `coef()` gives, the dimname `vcov()` carries, and the row
  # `confint()` both labels and matches a character `parm` against, and it was
  # all four of those before the column was renamed as well.
  covariance <- contrast_vcov()
  res <- ipw_result(with_vcov(contrast_estimates(), covariance))

  expect_identical(names(coef(res)), contrast_labels())
  expect_identical(rownames(confint(res)), contrast_labels())
  expect_identical(
    dimnames(vcov(res)),
    list(contrast_labels(), contrast_labels())
  )

  # Asserted against each other as well as against the literal, so that a change
  # reaching one surface and not the rest is caught here rather than by a caller
  # who reads a covariance out by the name `coef()` gave.
  expect_identical(rownames(vcov(res)), names(coef(res)))
  expect_identical(colnames(vcov(res)), names(coef(res)))
  expect_identical(rownames(confint(res)), names(coef(res)))

  # The label selects a row as well as naming one, which is the surface a caller
  # addresses `confint()` through.
  expect_identical(
    rownames(confint(res, parm = "log(or) c vs a")),
    "log(or) c vs a"
  )
})

# A frame stored by an earlier version of a fitting package names the same column
# `comparison`, and results built against that contract are still results a
# caller is holding. The frame is read the way `legacy_result()` above is read:
# the older spelling is an alias for the canonical one rather than a shape the
# result layer has no reading for.
#
# What the alias buys is only on the way in. Every surface on the way out is the
# canonical one whichever name the stored frame used, so a caller cannot tell
# from a table, a label, or a printed row which vintage of frame produced it.
# That is the whole of the guarantee, and each half of it needs its own
# assertion: a frame whose column the reader does not recognize reaches the
# output surfaces as no contrast at all, and the damage is silent rather than
# loud. `as.data.frame()` drops the column, leaving six rows whose terms repeat
# in pairs, which is exactly the shape `mice::pool()` groups across; and the
# labels collapse to three names used twice, which `coef()` and `confint()`
# report as duplicates and which `print()` cannot use as row names at all.
#
# The three tests after the first therefore hold today. They are here for the
# implementation they are written against rather than the one they are written
# on, since the reading that satisfies the canonical pins above satisfies them
# too only if it keeps the alias.

test_that("as.data.frame() reports a stored comparison column as contrast", {
  # The column name is the caller's contract, and it is the canonical one for a
  # frame of either vintage. A table that reported the stored spelling would
  # hand `mice::pool()` a column outside the whitelist it groups by and would
  # not stack with a table built from a frame that used the other name.
  estimates <- categorical_estimates()
  res <- ipw_result(estimates)

  df <- as.data.frame(res)

  expect_identical(
    names(df),
    c("term", "contrast", "estimate", "std.error", "statistic", "p.value")
  )
  expect_identical(names(df)[1:2], c("term", "contrast"))
  expect_false("comparison" %in% names(df))

  expect_identical(df$contrast, estimates$comparison)
  expect_identical(df$term, rep(c("rd", "log(rr)", "log(or)"), times = 2))

  # The table is the one the canonical frame gives, down to the covariance it
  # carries. Which column name the result was built with is not a fact about the
  # table, so there is nothing in it for a caller to branch on.
  expect_identical(df, as.data.frame(ipw_result(contrast_estimates())))
})

test_that("a stored comparison column keeps the effect labels as they were", {
  # The labels are built from the contrast whatever the column holding it is
  # called, so they are the six full ones rather than three names used twice.
  # Duplicated labels are the failure this guards against: `coef()` would return
  # two elements answering to `"rd"`, and `confint(parm = "rd b vs a")` would
  # refuse a label the result reports.
  covariance <- contrast_vcov()
  res <- ipw_result(with_vcov(categorical_estimates(), covariance))

  expect_identical(names(coef(res)), contrast_labels())
  expect_identical(rownames(confint(res)), contrast_labels())
  expect_identical(
    dimnames(vcov(res)),
    list(contrast_labels(), contrast_labels())
  )

  expect_identical(anyDuplicated(names(coef(res))), 0L)
  expect_identical(rownames(vcov(res)), names(coef(res)))
  expect_identical(rownames(confint(res)), names(coef(res)))
  expect_identical(
    rownames(confint(res, parm = "log(or) c vs a")),
    "log(or) c vs a"
  )
})

test_that("print() keys rows by effect and a stored comparison column", {
  # Row labels are the same six, and both character columns are gone before
  # `printCoefmat()` sees the frame. A reader that did not know the stored
  # spelling would leave that column in the numeric matrix and would have three
  # row names to give six rows, which `printCoefmat()` cannot be handed at all.
  res <- ipw_result(categorical_estimates())

  out <- capture.output(print(res))

  expect_match(out, "^rd b vs a +0\\.081945 ", all = FALSE)
  expect_match(out, "^log\\(rr\\) b vs a +0\\.168870 ", all = FALSE)
  expect_match(out, "^log\\(or\\) c vs a +0\\.676435 ", all = FALSE)

  expect_false(any(grepl("comparison", out, fixed = TRUE)))
  expect_false(any(grepl("contrast", out, fixed = TRUE)))
  expect_false(any(grepl("effect", out, fixed = TRUE)))
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
