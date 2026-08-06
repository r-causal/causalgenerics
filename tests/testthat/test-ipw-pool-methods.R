# The methods an `ipw_pooled` result carries: `print()`, `coef()`, `vcov()`,
# `confint()`, `nobs()`, and `as.data.frame()`. `pool_ipw()` itself, and the
# construction of the class, are in `test-ipw-pool.R`, the way `new_ipw()` and
# its methods are split across `test-ipw-result.R` and `test-ipw-accessors.R`.
# The test files in this suite mirror the files under `R/` one for one, and
# these methods are a surface of their own.
#
# They live in this package for the reason the `ipw` methods do. Two packages
# each registering `print.ipw_pooled()` would collide in the shared S3 method
# table, and a caller writing against a pooled result would then get whichever
# package was installed last rather than the contract.
#
# Everything the methods read is a field `pool_ipw()` writes: the `estimates`
# frame, the `ipw_vcov` attribute attached to it, `estimand`, `effects`, `m`,
# `dfcom`, `nobs`, and `outcome_link`. They recompute nothing, so a pooled
# result read back out of a saved file answers exactly as one just built does.
#
# The pooled surfaces differ from the unpooled ones in one way that runs through
# every method here: the inference is referred to t on each effect's own pooled
# degrees of freedom rather than to the normal. A method that reused the `ipw`
# spelling would agree on the estimates and be wrong about every interval and
# every p-value, by more the fewer imputations there are.
#
# The fixtures are rebuilt here rather than shared with `test-ipw-pool.R`, since
# testthat sources each file in an environment of its own. They hold the same
# numbers, so the pooled quantities that file works out by hand are the ones
# asserted against here: the `rd` row of the binary fixture pools to an estimate
# of 0.30 and a standard error of `sqrt(0.06)`, on 9639 / 1048 degrees of
# freedom when the complete-data count is 17.

# ---- fixtures ----------------------------------------------------------------

pool_data <- function(rows = 20) {
  dat <- data.frame(
    x = rep(c(-1.5, -0.5, 0.5, 1.5), each = 5),
    z = rep(c(0, 1), 10),
    y = rep(c(0, 1, 1, 0, 1), 4)
  )
  dat[seq_len(rows), , drop = FALSE]
}

pool_imputed_data <- function(i) {
  dat <- pool_data()
  dat$y[i] <- 1 - dat$y[i]
  dat
}

pool_outcome_model <- function(rows = 20) {
  glm(y ~ z, family = quasibinomial(), data = pool_data(rows))
}

pool_wt_model <- function(rows = 20) {
  glm(z ~ x, family = binomial(), data = pool_data(rows))
}

pool_binary_estimates <- function(estimate, std.err, conf.level = 0.95) {
  half_width <- stats::qnorm(1 - (1 - conf.level) / 2) * std.err
  data.frame(
    effect = c("rd", "log(rr)", "log(or)"),
    estimate = estimate,
    std.err = std.err,
    z = estimate / std.err,
    ci.lower = estimate - half_width,
    ci.upper = estimate + half_width,
    conf.level = conf.level,
    p.value = 2 * stats::pnorm(-abs(estimate / std.err))
  )
}

# The same three imputations `test-ipw-pool.R` uses. The `rd` column is what
# carries the literals: estimates of 0.20, 0.30, and 0.40 against standard
# errors of 0.10, 0.20, and 0.30.
pool_binary_estimate <- function() {
  rbind(
    c(0.20, 0.55, 0.85),
    c(0.30, 0.60, 0.95),
    c(0.40, 0.65, 1.05)
  )
}

pool_binary_std_err <- function() {
  rbind(
    c(0.10, 0.25, 0.40),
    c(0.20, 0.30, 0.45),
    c(0.30, 0.35, 0.50)
  )
}

pool_binary_labels <- function() {
  c("rd", "log(rr)", "log(or)")
}

pool_categorical_labels <- function() {
  c(
    "rd b vs a",
    "log(rr) b vs a",
    "log(or) b vs a",
    "rd c vs a",
    "log(rr) c vs a",
    "log(or) c vs a"
  )
}

pool_categorical_estimates <- function(estimate, std.err) {
  half_width <- stats::qnorm(0.975) * std.err
  data.frame(
    effect = rep(c("rd", "log(rr)", "log(or)"), times = 2),
    contrast = rep(c("b vs a", "c vs a"), each = 3),
    estimate = estimate,
    std.err = std.err,
    z = estimate / std.err,
    ci.lower = estimate - half_width,
    ci.upper = estimate + half_width,
    conf.level = 0.95,
    p.value = 2 * stats::pnorm(-abs(estimate / std.err))
  )
}

pool_categorical_estimate <- function() {
  rbind(
    c(0.08, 0.17, 0.33, 0.16, 0.31, 0.67),
    c(0.10, 0.20, 0.38, 0.20, 0.35, 0.75),
    c(0.12, 0.23, 0.43, 0.24, 0.39, 0.83)
  )
}

pool_categorical_std_err <- function() {
  rbind(
    c(0.05, 0.10, 0.20, 0.04, 0.09, 0.18),
    c(0.06, 0.11, 0.21, 0.05, 0.10, 0.19),
    c(0.07, 0.12, 0.22, 0.06, 0.11, 0.20)
  )
}

pool_effects_vcov <- function(std_err, labels, rho = 0.9) {
  index <- seq_along(std_err)
  covariance <- rho^abs(outer(index, index, "-")) * outer(std_err, std_err)
  dimnames(covariance) <- list(labels, labels)
  covariance
}

pool_fit <- function(
  estimates,
  vcov = NULL,
  effects = "marginal",
  outcome_mod = pool_outcome_model()
) {
  if (!is.null(vcov)) {
    attr(estimates, "ipw_vcov") <- vcov
  }
  new_ipw(
    estimand = "ate",
    wt_mod = pool_wt_model(),
    outcome_mod = outcome_mod,
    estimates = estimates,
    se_method = "linearization",
    fit = NULL,
    effects = effects
  )
}

pool_binary_fits <- function(vcov = TRUE) {
  estimate <- pool_binary_estimate()
  std_err <- pool_binary_std_err()
  lapply(1:3, function(i) {
    pool_fit(
      pool_binary_estimates(estimate[i, ], std_err[i, ]),
      vcov = if (vcov) {
        pool_effects_vcov(std_err[i, ], pool_binary_labels())
      } else {
        NULL
      }
    )
  })
}

pool_categorical_fits <- function(vcov = TRUE) {
  estimate <- pool_categorical_estimate()
  std_err <- pool_categorical_std_err()
  lapply(1:3, function(i) {
    pool_fit(
      pool_categorical_estimates(estimate[i, ], std_err[i, ]),
      vcov = if (vcov) {
        pool_effects_vcov(std_err[i, ], pool_categorical_labels())
      } else {
        NULL
      }
    )
  })
}

# The outcome models of three imputations, wrapped with the corrected covariance
# the joint estimation implies. The conditional reading pools that surface, so
# every pooled conditional result carries a covariance and the labels are
# coefficient names rather than effect labels.
pool_conditional_models <- function() {
  lapply(1:3, function(i) {
    mod <- glm(y ~ z, family = quasibinomial(), data = pool_imputed_data(i))
    new_ipw_model(mod, stats::vcov(mod) * (1 + i / 10))
  })
}

# The same, on an identity link. A coefficient there is not on a scale an
# exponential undoes, which is what separates the two conditional
# `exponentiate` paths below.
pool_identity_models <- function() {
  lapply(1:3, function(i) {
    mod <- lm(y ~ z, data = pool_imputed_data(i))
    new_ipw_model(mod, stats::vcov(mod) * (1 + i / 10))
  })
}

pool_conditional_fits <- function(models = pool_conditional_models()) {
  estimate <- pool_binary_estimate()
  std_err <- pool_binary_std_err()
  Map(
    function(i, mod) {
      pool_fit(
        pool_binary_estimates(estimate[i, ], std_err[i, ]),
        outcome_mod = mod,
        effects = "conditional"
      )
    },
    1:3,
    models
  )
}

# The pooled results every assertion below is written against. The complete-data
# count is named rather than read off the fits, so the degrees of freedom are
# the ones the literals were worked out for.
pooled_binary <- function(vcov = TRUE) {
  pool_ipw(pool_binary_fits(vcov = vcov), dfcom = 17)
}

pooled_categorical <- function() {
  pool_ipw(pool_categorical_fits(), dfcom = 17)
}

pooled_conditional <- function(models = pool_conditional_models()) {
  pool_ipw(pool_conditional_fits(models), dfcom = 18)
}

# A pooled result whose stored bounds are not the ones a recomputation gives:
# asymmetric about the estimate, the way a bootstrap interval is, and recorded
# at 0.9 rather than the 0.95 a caller who names no level asks for. Both facts
# are what make "the stored bounds come back at the stored level" observable.
pooled_stored_interval <- function() {
  res <- pooled_binary()
  res$estimates$ci.lower <- c(0.061304, 0.128951, 0.207738)
  res$estimates$ci.upper <- c(0.372015, 1.043772, 1.622904)
  res$estimates$conf.level <- 0.9
  res
}

# The pooled degrees of freedom of the binary fixture at a complete-data count
# of 17, worked out by hand in `test-ipw-pool.R`. Written here rather than read
# off the result, since the bounds and the p-values below are the claim and
# reading their own `df` column back would assert nothing about it.
pooled_rd_df <- function() {
  9639 / 1048
}

# Whether `print()` writes a row labelled exactly `label`. The estimate columns
# are numeric, so a row label is a line whose remainder after the label is
# spaces and then a number. Testing only that a line starts with the label would
# hold for a prefix of the real one, and `printCoefmat()` wraps significance
# stars onto lines of their own under the same labels when the frame is wide.
labels_a_printed_row <- function(out, label) {
  candidates <- out[startsWith(out, label)]
  remainder <- substring(candidates, nchar(label) + 1L)
  any(grepl("^ +-?[0-9]", remainder))
}

# The numbers `print()` writes on the row labelled `label`. Which columns the
# table carries is the implementation's to choose, and the estimate and its
# standard error are the first two whichever others are there, so reading by
# position lets an assertion name those two without fixing the rest.
printed_numbers <- function(out, label) {
  rows <- out[startsWith(out, label)]
  rows <- rows[grepl("[0-9]", substring(rows, nchar(label) + 1L))]
  if (length(rows) != 1L) {
    return(numeric())
  }
  tokens <- strsplit(trimws(substring(rows, nchar(label) + 1L)), " +")[[1]]
  as.numeric(tokens[grepl("^-?[0-9]", tokens)])
}

# ---- print.ipw_pooled() ------------------------------------------------------

# Each print test asserts the lines that carry the contract. The green round
# adds the snapshot of the whole output beside them, the way every other print
# test in this suite pairs the two: the snapshot alone would not hold the
# contract, since the first run records whatever the implementation produced and
# testthat skips snapshot comparison on CRAN, so the named assertions are the
# only print coverage that runs there.

test_that("print() summarizes a pooled binary result", {
  res <- pooled_binary()

  out <- capture.output(print(res))

  expect_match(
    out,
    "^Pooled Inverse Probability Weight Estimator$",
    all = FALSE
  )
  # The estimand is reported in upper case, the way an unpooled result reports
  # it, so `ate` reads as ATE.
  expect_match(out, "^Estimand: ATE\\s*$", all = FALSE)
  # The reading is named the way `print.ipw()` names it, since a pooled result
  # presents one of the same two surfaces.
  expect_match(out, "^Effects: marginal", all = FALSE)
  # The two facts an unpooled result has no counterpart for. How many analyses
  # went in decides how much of the interval is between-imputation variance, and
  # the complete-data count is what the small-sample adjustment spent, so a
  # reader cannot judge the degrees of freedom in the table without both.
  expect_match(out, "Imputations: 3", fixed = TRUE, all = FALSE)
  expect_match(out, "Complete-data df: 17", fixed = TRUE, all = FALSE)

  expect_match(out, "^Pooled marginal estimates:$", all = FALSE)

  # Rows are keyed by effect label, and the character `effect` column is gone
  # rather than formatted as a number.
  for (label in pool_binary_labels()) {
    expect_true(labels_a_printed_row(out, label))
  }
  expect_false(any(grepl("effect", out, fixed = TRUE)))

  # The estimate and the standard error are the first two numbers on the row.
  expect_equal(printed_numbers(out, "rd")[1], 0.3, tolerance = 1e-4)
  expect_equal(printed_numbers(out, "rd")[2], sqrt(0.06), tolerance = 1e-4)
})

test_that("print() keys pooled rows by effect and contrast", {
  # The categorical shape is the only place the two marginal print paths differ.
  # Row labels become `effect` and `contrast` together, and both character
  # columns have to be gone before `printCoefmat()` sees the frame: dropping
  # `effect` alone does not error, it factor-codes `contrast` into a column
  # reading `1.000000` and `2.000000` beside the real estimates.
  res <- pooled_categorical()

  out <- capture.output(print(res))

  for (label in pool_categorical_labels()) {
    expect_true(labels_a_printed_row(out, label))
  }
  expect_false(any(grepl("contrast", out, fixed = TRUE)))
  expect_false(any(grepl("effect", out, fixed = TRUE)))
})

test_that("print() names the conditional reading of a pooled result", {
  # The conditional reading tabulates the pooled outcome-model coefficients, so
  # the rows are coefficient names and the heading says which surface is being
  # shown. A reader handed the table alone would otherwise take a pooled
  # `(Intercept)` for a causal effect.
  res <- pooled_conditional()

  out <- capture.output(print(res))

  expect_match(out, "^Effects: conditional", all = FALSE)
  expect_match(out, "^Pooled conditional estimates", all = FALSE)

  for (label in c("(Intercept)", "z")) {
    expect_true(labels_a_printed_row(out, label))
  }
  # Not the marginal surface, which these same results also carry.
  expect_false(labels_a_printed_row(out, "log(or)"))
})

test_that("print() leaves the pooling diagnostics out of the table", {
  # `riv`, `lambda`, and `fmi` describe how much of the uncertainty came from
  # the imputation rather than from the data. They are a diagnostic rather than
  # a result, they are in the `pooling` field for a caller who wants them, and
  # five more columns beside the estimates would push the table past the width
  # where `printCoefmat()` keeps a row on one line.
  out <- capture.output(print(pooled_binary()))

  expect_false(any(grepl("riv", out, fixed = TRUE)))
  expect_false(any(grepl("lambda", out, fixed = TRUE)))
  expect_false(any(grepl("fmi", out, fixed = TRUE)))
  expect_false(any(grepl("ubar", out, fixed = TRUE)))
  # The per-effect degrees of freedom do appear: they are what the interval and
  # the p-value beside them were computed on.
  expect_true(any(grepl("df", out, fixed = TRUE)))
})

test_that("print() returns its input invisibly", {
  res <- pooled_binary()

  # The pooled table is what ran, so the invisible return asserted below is that
  # method's rather than one taken by `print.default()`, which returns its input
  # invisibly too and would satisfy the rest of this test on any list.
  expect_match(
    capture.output(print(res)),
    "^Pooled marginal estimates:$",
    all = FALSE
  )

  # The printed form is asserted above and snapshotted in the green round, so
  # from here the output only needs to go somewhere other than the test report.
  withr::local_output_sink(withr::local_tempfile())
  printed <- withVisible(print(res))

  expect_false(printed$visible)
  expect_identical(printed$value, res)
})

# ---- coef.ipw_pooled() -------------------------------------------------------

test_that("coef() returns the pooled estimates named by effect", {
  res <- pooled_binary()

  expect_identical(
    coef(res),
    setNames(
      res$estimates$estimate,
      c(
        "rd",
        "log(rr)",
        "log(or)"
      )
    )
  )
  expect_equal(unname(coef(res)[["rd"]]), 0.3)
  expect_length(coef(res), 3L)
})

test_that("coef() keys a pooled categorical result by effect and contrast", {
  # The effect labels repeat across contrasts, so `effect` alone would name
  # three of the six rows twice over.
  res <- pooled_categorical()

  expect_identical(names(coef(res)), pool_categorical_labels())
  expect_identical(unname(coef(res)), res$estimates$estimate)
})

test_that("coef() reports the coefficient names in the conditional reading", {
  # The pooled conditional surface is the outcome models' coefficients, so the
  # names are theirs. There is no `effects` argument here: a pooled result holds
  # one surface, the one `pool_ipw()` pooled, rather than both.
  res <- pooled_conditional()

  expect_identical(names(coef(res)), c("(Intercept)", "z"))
  expect_identical(unname(coef(res)), res$estimates$estimate)
  expect_false(any(names(coef(res)) %in% pool_binary_labels()))
})

# ---- vcov.ipw_pooled() -------------------------------------------------------

test_that("vcov() returns the pooled covariance", {
  res <- pooled_binary()
  covariance <- attr(res$estimates, "ipw_vcov", exact = TRUE)

  expect_identical(vcov(res), covariance)
  expect_identical(
    dimnames(vcov(res)),
    list(pool_binary_labels(), pool_binary_labels())
  )
  expect_identical(vcov(res), t(vcov(res)))
  # The diagonal is the square of the pooled standard errors, so a caller
  # reading a variance out by the name `coef()` gave gets that effect's.
  expect_equal(diag(vcov(res)), res$estimates$std.err^2, ignore_attr = TRUE)
  expect_identical(rownames(vcov(res)), names(coef(res)))
})

test_that("vcov() returns the pooled covariance of a categorical result", {
  res <- pooled_categorical()

  expect_identical(dim(vcov(res)), c(6L, 6L))
  expect_identical(rownames(vcov(res)), pool_categorical_labels())
})

test_that("vcov() names the conditional block by coefficient", {
  res <- pooled_conditional()

  expect_identical(dim(vcov(res)), c(2L, 2L))
  expect_identical(
    dimnames(vcov(res)),
    list(c("(Intercept)", "z"), c("(Intercept)", "z"))
  )
  expect_identical(rownames(vcov(res)), names(coef(res)))
})

test_that("vcov() refuses a pooled result that carries no covariance", {
  # `pool_ipw()` attaches no covariance when any of the results it pooled
  # carried none, since the average of the within-imputation covariances needs
  # one from every imputation. There is no honest fallback here for the reason
  # there is none on an unpooled result: the standard errors in the frame give
  # the diagonal and say nothing about the off-diagonal entries, which are far
  # from zero for effects estimated from the same weighted means.
  res <- pooled_binary(vcov = FALSE)

  expect_null(attr(res$estimates, "ipw_vcov", exact = TRUE))
  expect_error(vcov(res), class = "causalgenerics_no_vcov_ipw_pooled")
  expect_error(vcov(res), class = "causalgenerics_no_vcov")

  # The message names the object it is about, the way it does for an `ipw` and
  # for an `ipw_model`, so a caller holding several kinds of result knows which
  # one refused.
  cnd <- tryCatch(vcov(res), error = identity)
  expect_match(conditionMessage(cnd), "ipw_pooled", fixed = TRUE)
})

# ---- confint.ipw_pooled() ----------------------------------------------------

test_that("confint() bounds the pooled estimates on their own df", {
  # The half width is `qt(1 - (1 - level) / 2, df)` times the standard error,
  # taken at each row's own pooled degrees of freedom. Those are in single
  # figures here, where t and the normal are far apart, so a method that reused
  # `confint.ipw()`'s `qnorm()` would be visibly too narrow.
  res <- pooled_binary()
  estimates <- res$estimates

  ci <- confint(res, level = 0.8)
  half_width <- stats::qt(0.9, estimates$df) * estimates$std.err

  expect_true(is.matrix(ci))
  expect_identical(dim(ci), c(3L, 2L))
  expect_identical(dimnames(ci), list(pool_binary_labels(), c("10 %", "90 %")))
  expect_equal(
    ci[, 1],
    setNames(
      estimates$estimate - half_width,
      pool_binary_labels()
    )
  )
  expect_equal(
    ci[, 2],
    setNames(
      estimates$estimate + half_width,
      pool_binary_labels()
    )
  )

  # Wider than the normal interval at the same level, which is what says the
  # degrees of freedom were used at all.
  normal <- stats::qnorm(0.9) * estimates$std.err
  expect_true(all(half_width > normal))
})

test_that("confint() bounds the rd row at a number that can be written down", {
  # The `rd` row pools to an estimate of 0.30 and a standard error of
  # `sqrt(0.06)` on 9639 / 1048 degrees of freedom, which `test-ipw-pool.R`
  # works out by hand. The bounds follow from those three numbers and the t
  # quantile, with nothing read back off the result.
  res <- pooled_binary()

  ci <- confint(res, parm = "rd", level = 0.95)
  half_width <- stats::qt(0.975, pooled_rd_df()) * sqrt(0.06)

  expect_identical(rownames(ci), "rd")
  expect_equal(ci[1, 1], 0.3 - half_width)
  expect_equal(ci[1, 2], 0.3 + half_width)
})

test_that("confint() returns the stored bounds at the stored level", {
  # The bounds in the frame are the ones `pool_ipw()` reported, and at the level
  # it reported them for they are returned rather than rebuilt. The fixture's
  # stored bounds are asymmetric about the estimate, so no recomputation
  # produces them, and the level they are stored at is not the one a caller who
  # names none asks for.
  res <- pooled_stored_interval()
  estimates <- res$estimates

  stored <- confint(res, level = 0.9)

  expect_identical(
    dimnames(stored),
    list(pool_binary_labels(), c("5 %", "95 %"))
  )
  expect_identical(
    stored[, 1],
    setNames(estimates$ci.lower, pool_binary_labels())
  )
  expect_identical(
    stored[, 2],
    setNames(estimates$ci.upper, pool_binary_labels())
  )

  # At any other level, including the 0.95 a caller who names none gets, the
  # bounds are rebuilt from the estimate, the standard error, and the row's df.
  default <- confint(res)
  half_width <- stats::qt(0.975, estimates$df) * estimates$std.err

  expect_identical(
    dimnames(default),
    list(pool_binary_labels(), c("2.5 %", "97.5 %"))
  )
  expect_equal(
    default[, 1],
    setNames(estimates$estimate - half_width, pool_binary_labels())
  )
  expect_false(isTRUE(all.equal(
    unname(default[, 1]),
    estimates$ci.lower
  )))
})

test_that("confint() labels its columns the way stats does", {
  res <- pooled_binary()

  expect_identical(colnames(confint(res, level = 0.99)), c("0.5 %", "99.5 %"))
  expect_identical(colnames(confint(res, level = 0.8)), c("10 %", "90 %"))
})

test_that("confint() selects pooled rows by label and by position", {
  # `parm` means what it means everywhere else in this package: a character
  # vector names the rows and a numeric one indexes them, and either way the
  # rows come back in the order they were asked for.
  res <- pooled_binary()

  expect_identical(rownames(confint(res, parm = "log(rr)")), "log(rr)")
  expect_identical(dim(confint(res, parm = "log(rr)")), c(1L, 2L))
  expect_identical(
    rownames(confint(res, parm = c("log(or)", "rd"))),
    c("log(or)", "rd")
  )
  expect_identical(rownames(confint(res, parm = 2)), "log(rr)")
  expect_identical(rownames(confint(res, parm = c(3, 1))), c("log(or)", "rd"))
})

test_that("confint() selects a pooled categorical row by its full label", {
  # The label is effect and contrast together, so `parm = "rd"` names nothing in
  # a categorical result even though `rd` is a value of the `effect` column.
  res <- pooled_categorical()

  expect_identical(rownames(confint(res, parm = "rd c vs a")), "rd c vs a")
  expect_error(
    confint(res, parm = "rd"),
    class = "causalgenerics_invalid_argument_parm"
  )
})

test_that("confint() refuses a parm the pooled result does not report", {
  # The same guard `confint.ipw()` keeps, reached through the same helper.
  # Dropping the unmatched labels and returning the rest would answer a question
  # that was not asked, with a matrix of the wrong number of rows and nothing to
  # say why.
  res <- pooled_binary()

  expect_error(
    confint(res, parm = "rr"),
    class = "causalgenerics_invalid_argument_parm"
  )
  expect_error(
    confint(res, parm = "rr"),
    class = "causalgenerics_invalid_argument"
  )
  expect_error(
    confint(res, parm = c("rd", "rr")),
    class = "causalgenerics_invalid_argument_parm"
  )
  expect_error(
    confint(res, parm = 99),
    class = "causalgenerics_invalid_argument_parm"
  )
  # `NA` is the subscript that would go furthest unnoticed: a logical `NA`
  # recycles to the length of the surface, so indexing with it gives one `NA`
  # row per effect, a matrix of exactly the shape a caller expects holding
  # nothing.
  expect_error(
    confint(res, parm = NA),
    class = "causalgenerics_invalid_argument_parm"
  )
  expect_error(
    confint(res, parm = NA_integer_),
    class = "causalgenerics_invalid_argument_parm"
  )
})

test_that("confint() gives a numeric parm its ordinary subscript meaning", {
  # Zero selects nothing and a negative position drops the row it names, the way
  # a subscript reads everywhere else. An empty selection is still a labelled
  # matrix, so a caller who built one from a filter that matched nothing can
  # read its columns or bind it to another without treating the case apart.
  res <- pooled_binary()

  none <- confint(res, parm = 0)

  expect_true(is.matrix(none))
  expect_identical(dim(none), c(0L, 2L))
  expect_identical(colnames(none), c("2.5 %", "97.5 %"))

  dropped <- confint(res, parm = -1)

  expect_identical(dim(dropped), c(2L, 2L))
  expect_identical(rownames(dropped), c("log(rr)", "log(or)"))
})

test_that("confint() bounds the conditional reading the same way", {
  res <- pooled_conditional()
  estimates <- res$estimates

  ci <- confint(res)
  half_width <- stats::qt(0.975, estimates$df) * estimates$std.err

  expect_identical(rownames(ci), c("(Intercept)", "z"))
  expect_equal(
    ci[, 1],
    setNames(estimates$estimate - half_width, c("(Intercept)", "z"))
  )
})

# ---- nobs.ipw_pooled() -------------------------------------------------------

test_that("nobs() reports the smallest observation count, as an integer", {
  # The count `pool_ipw()` stored, which is the smallest any pooled analysis was
  # estimated from.
  #
  # This is the one method here whose answer is already right without it:
  # `nobs.default()` reads `object[["nobs"]]` from any list, and the pooled
  # result has a field of that name. A registered method is what makes the
  # answer a contract rather than a coincidence of two names agreeing, so the
  # registration is asserted alongside the value, and the value is read through
  # a call that cannot see this frame.
  res <- pooled_binary()

  expect_identical(nobs(res), 20L)
  expect_type(nobs(res), "integer")
  expect_identical(nobs(res), res$nobs)

  expect_true(exists(
    "nobs.ipw_pooled",
    envir = s3_methods_table("nobs"),
    inherits = FALSE
  ))
  expect_identical(dispatch_from_baseenv(nobs, res), 20L)
})

# ---- df.residual() -------------------------------------------------------

test_that("df.residual() has no method for a pooled result", {
  # Deliberately absent. Residual degrees of freedom are a property of one fit,
  # and a pooled result has as many as it had imputations. What it does have is
  # a pooled degrees of freedom per effect, which is not a residual count and
  # differs from row to row, so it lives in the `df` column of the `estimates`
  # frame where it sits beside the interval and the p-value it was used for.
  # Answering with any single number here would be answering a question the
  # object cannot be asked.
  #
  # `df.residual.default()` reads a `df.residual` element, and a pooled result
  # has none, so the generic gives `NULL`. This test passes today and has to
  # keep passing: it pins the absence rather than a behavior to be built.
  res <- pooled_binary()

  expect_null(stats::df.residual(res))
  expect_false(exists(
    "df.residual.ipw_pooled",
    envir = s3_methods_table("df.residual"),
    inherits = FALSE
  ))
  expect_true("df" %in% names(res$estimates))
  expect_length(res$estimates$df, 3L)
})

# ---- as.data.frame.ipw_pooled() ----------------------------------------------

test_that("as.data.frame() reports the tidier-shaped pooled table", {
  # The tidier column names, with `df` after the statistic: the pooled inference
  # is referred to t on it, so a tidier surface that dropped it would leave a
  # statistic no reader could refer to anything.
  res <- pooled_binary()
  df <- as.data.frame(res)

  expect_s3_class(df, "data.frame")
  expect_identical(
    names(df),
    c("term", "estimate", "std.error", "statistic", "df", "p.value")
  )
  expect_identical(df$term, pool_binary_labels())
  expect_identical(df$estimate, res$estimates$estimate)
  expect_identical(df$std.error, res$estimates$std.err)
  expect_identical(df$statistic, res$estimates$t)
  expect_identical(df$df, res$estimates$df)
  expect_identical(df$p.value, res$estimates$p.value)
  # The storage names are gone, so code reading the table by them is reading a
  # column that is not there rather than the wrong one.
  expect_false(any(c("effect", "std.err", "z", "conf.level") %in% names(df)))
})

test_that("as.data.frame() names the contrast column of a pooled table", {
  res <- pooled_categorical()
  df <- as.data.frame(res)

  expect_identical(
    names(df),
    c(
      "term",
      "contrast",
      "estimate",
      "std.error",
      "statistic",
      "df",
      "p.value"
    )
  )
  expect_identical(df$contrast, rep(c("b vs a", "c vs a"), each = 3))
  expect_identical(nrow(df), 6L)
})

test_that("as.data.frame() heads the conditional table by coefficient", {
  # A pooled conditional table names one contrast, so it carries no `contrast`
  # column and its terms are coefficient names.
  res <- pooled_conditional()
  df <- as.data.frame(res)

  expect_identical(df$term, c("(Intercept)", "z"))
  expect_false("contrast" %in% names(df))
})

test_that("as.data.frame() appends the bounds only when they are asked for", {
  res <- pooled_binary()

  expect_false(any(c("conf.low", "conf.high") %in% names(as.data.frame(res))))

  df <- as.data.frame(res, conf.int = TRUE)

  # Last, after the columns the table always carries, so asking for an interval
  # adds to the table rather than rearranging it.
  expect_identical(
    names(df),
    c(
      "term",
      "estimate",
      "std.error",
      "statistic",
      "df",
      "p.value",
      "conf.low",
      "conf.high"
    )
  )
  # At the stored level the stored bounds come back.
  expect_identical(df$conf.low, res$estimates$ci.lower)
  expect_identical(df$conf.high, res$estimates$ci.upper)
})

test_that("as.data.frame() takes the stored level when none is named", {
  # `conf.level = NULL` is the default here rather than `0.95`, since a pooled
  # result records the level its bounds were built at and a fixed default would
  # rebuild them for a result stored at any other one. The fixture stores 0.9,
  # so a hardcoded default would fail this.
  res <- pooled_stored_interval()

  df <- as.data.frame(res, conf.int = TRUE)

  expect_identical(df$conf.low, res$estimates$ci.lower)
  expect_identical(df$conf.high, res$estimates$ci.upper)
})

test_that("as.data.frame() recomputes the bounds at another level", {
  # The frame-level all-or-nothing rule `as.data.frame.ipw()` keeps: the stored
  # bounds come back only when every row records the level asked for, and
  # otherwise the pair is rebuilt, here from t on each row's own df.
  res <- pooled_binary()
  estimates <- res$estimates

  df <- as.data.frame(res, conf.int = TRUE, conf.level = 0.8)
  half_width <- stats::qt(0.9, estimates$df) * estimates$std.err

  expect_equal(df$conf.low, estimates$estimate - half_width)
  expect_equal(df$conf.high, estimates$estimate + half_width)
  expect_false(isTRUE(all.equal(df$conf.low, estimates$ci.lower)))

  # A frame whose rows disagree about the level records none, so both bounds are
  # rebuilt even at a level some of the rows name. The fixture is the one whose
  # stored bounds are asymmetric about the estimate: the bounds a rebuild gives
  # are not the stored ones for any row, which is what makes the frame-level
  # rule separable from a per-row one here. Built on `pooled_binary()` the two
  # rules agree to the bit, since its stored bounds are the recomputation.
  mixed <- pooled_stored_interval()
  mixed$estimates$conf.level <- c(0.9, 0.95, 0.9)

  rebuilt <- as.data.frame(mixed, conf.int = TRUE, conf.level = 0.9)
  mixed_half <- stats::qt(0.95, mixed$estimates$df) * mixed$estimates$std.err

  expect_equal(rebuilt$conf.low, mixed$estimates$estimate - mixed_half)
  expect_equal(rebuilt$conf.high, mixed$estimates$estimate + mixed_half)
  # Not the stored pair for the two rows that do record 0.9, which is what a
  # per-row rule would return for them.
  expect_false(isTRUE(all.equal(rebuilt$conf.low, mixed$estimates$ci.lower)))
  expect_false(isTRUE(all.equal(
    rebuilt$conf.low[c(1, 3)],
    mixed$estimates$ci.lower[c(1, 3)]
  )))
})

test_that("as.data.frame() carries the pooled covariance on the plain table", {
  # The covariance belongs to the estimates rather than to a column of them, so
  # it travels on the returned table under the same attribute, exactly as
  # `as.data.frame.ipw()` carries an unpooled one. A result whose fits did not
  # all carry one carries none: assigning `NULL` sets no attribute.
  res <- pooled_binary()
  df <- as.data.frame(res)

  expect_identical(
    attr(df, "ipw_vcov", exact = TRUE),
    attr(res$estimates, "ipw_vcov", exact = TRUE)
  )
  expect_identical(
    dimnames(attr(df, "ipw_vcov", exact = TRUE)),
    list(pool_binary_labels(), pool_binary_labels())
  )

  none <- as.data.frame(pooled_binary(vcov = FALSE))
  expect_null(attr(none, "ipw_vcov", exact = TRUE))
})

test_that("as.data.frame() refuses arguments that are not what they must be", {
  # The same three checks `as.data.frame.ipw()` makes, with the same classes. A
  # number and a string both take the true branch of a bare `if`, so an
  # unchecked `conf.int = "no"` would report an interval.
  res <- pooled_binary()

  expect_error(
    as.data.frame(res, conf.int = "no"),
    class = "causalgenerics_invalid_argument_conf.int"
  )
  expect_error(
    as.data.frame(res, conf.int = NA),
    class = "causalgenerics_invalid_argument_conf.int"
  )
  expect_error(
    as.data.frame(res, exponentiate = 1),
    class = "causalgenerics_invalid_argument_exponentiate"
  )
  expect_error(
    as.data.frame(res, conf.level = 2),
    class = "causalgenerics_invalid_argument_conf.level"
  )
  expect_error(
    as.data.frame(res, conf.level = c(0.9, 0.95)),
    class = "causalgenerics_invalid_argument_conf.level"
  )
  expect_error(
    as.data.frame(res, conf.level = 2),
    class = "causalgenerics_invalid_argument"
  )
  # `conf.level` is checked whether or not the bounds were asked for, since a
  # call naming a level no interval can be built at has said something wrong
  # either way.
  expect_error(
    as.data.frame(res, conf.int = FALSE, conf.level = 0),
    class = "causalgenerics_invalid_argument_conf.level"
  )
})

test_that("as.data.frame() exponentiates the pooled ratio rows", {
  # Exactly `as.data.frame.ipw()`'s semantics on a marginal table. The rows to
  # move are matched on the labels `log(rr)` and `log(or)` exactly, so a table
  # whose ratios are already on the natural scale is left alone rather than
  # exponentiated twice, and the two terms are relabelled with the scale they
  # now carry.
  res <- pooled_binary()
  plain <- as.data.frame(res)

  df <- as.data.frame(res, exponentiate = TRUE)

  expect_identical(df$term, c("rd", "rr", "or"))
  expect_identical(df$estimate[[1]], plain$estimate[[1]])
  expect_equal(df$estimate[[2]], exp(plain$estimate[[2]]))
  expect_equal(df$estimate[[3]], exp(plain$estimate[[3]]))
  # The inference is done on the log scale and stays there.
  expect_identical(df$std.error, plain$std.error)
  expect_identical(df$statistic, plain$statistic)
  expect_identical(df$p.value, plain$p.value)
  expect_identical(df$df, plain$df)
})

test_that("as.data.frame() settles the interval before the scale", {
  # The bounds are built on the log scale and exponentiated afterwards, so a
  # recomputed bound is a t half width on the log scale added to an estimate on
  # the log scale. Exponentiating first and then adding a half width would give
  # an interval that is not the image of the log-scale one and can run below
  # zero for a ratio.
  res <- pooled_binary()
  estimates <- res$estimates

  df <- as.data.frame(
    res,
    conf.int = TRUE,
    conf.level = 0.8,
    exponentiate = TRUE
  )
  half_width <- stats::qt(0.9, estimates$df) * estimates$std.err
  lower <- estimates$estimate - half_width
  upper <- estimates$estimate + half_width

  expect_equal(df$conf.low, c(lower[1], exp(lower[2]), exp(lower[3])))
  expect_equal(df$conf.high, c(upper[1], exp(upper[2]), exp(upper[3])))
  expect_true(all(df$conf.low[2:3] > 0))
})

test_that("as.data.frame() drops the covariance when it exponentiates", {
  # The matrix describes the estimates on the scale they were estimated on, and
  # once the ratio rows have moved it describes neither the table it sits on nor
  # anything else. There is no correct matrix to put in its place, since the
  # delta method answer is not what the pooling computed.
  res <- pooled_binary()

  expect_null(attr(
    as.data.frame(res, exponentiate = TRUE),
    "ipw_vcov",
    exact = TRUE
  ))
  expect_null(attr(
    as.data.frame(res, conf.int = TRUE, exponentiate = TRUE),
    "ipw_vcov",
    exact = TRUE
  ))
})

test_that("as.data.frame() exponentiates a conditional table on a log link", {
  # A conditional table has no rows labelled as ratios to pick out, so the
  # stored link settles whether there is anything for an exponential to undo. A
  # logit link puts every coefficient on the log odds scale, which is such a
  # scale, so every estimate moves and no term is relabelled: the terms are
  # coefficient names, and a coefficient does not change its name with the scale
  # its estimate is reported on.
  res <- pooled_conditional()
  plain <- as.data.frame(res)

  expect_identical(res$outcome_link, "logit")

  df <- as.data.frame(res, conf.int = TRUE, exponentiate = TRUE)

  expect_identical(df$term, c("(Intercept)", "z"))
  expect_equal(df$estimate, exp(plain$estimate))
  expect_equal(df$conf.low, exp(as.data.frame(res, conf.int = TRUE)$conf.low))
  expect_identical(df$std.error, plain$std.error)
  expect_identical(df$statistic, plain$statistic)
  expect_identical(df$p.value, plain$p.value)
})

test_that("as.data.frame() refuses to exponentiate an identity-link table", {
  # A coefficient on the identity scale exponentiates to a number describing
  # nothing, and there is no subset of rows to move instead, so the call is
  # refused rather than answered in part. The refusal names the link, since that
  # is the fact about the result the caller has to act on.
  res <- pooled_conditional(pool_identity_models())

  expect_identical(res$outcome_link, "identity")

  expect_error(
    as.data.frame(res, exponentiate = TRUE),
    class = "causalgenerics_invalid_argument_exponentiate"
  )
  expect_error(
    as.data.frame(res, exponentiate = TRUE),
    class = "causalgenerics_invalid_argument"
  )

  cnd <- tryCatch(
    as.data.frame(res, exponentiate = TRUE),
    error = identity
  )
  expect_match(conditionMessage(cnd), "identity", fixed = TRUE)

  # The table itself is still readable, and the refusal is about the one
  # argument rather than about the result.
  expect_no_error(as.data.frame(res))
  expect_no_error(as.data.frame(res, conf.int = TRUE))
  # A marginal pooled table is unaffected by the link, since it picks its rows
  # out by label.
  expect_no_error(as.data.frame(pooled_binary(), exponentiate = TRUE))
})

# ---- registration ------------------------------------------------------------

test_that("the pooled result's methods are registered for dispatch", {
  # Downstream packages reach these through the S3 method table a NAMESPACE
  # `S3method()` directive fills in. Each method is asserted to be in the table
  # belonging to its own generic, which for `coef()` and the rest is stats' and
  # for `print()` and `as.data.frame()` is base's.
  for (generic in c(
    "print",
    "coef",
    "vcov",
    "confint",
    "nobs",
    "as.data.frame"
  )) {
    expect_true(
      exists(
        paste0(generic, ".ipw_pooled"),
        envir = s3_methods_table(generic),
        inherits = FALSE
      ),
      info = generic
    )
  }
})

test_that("the pooled methods dispatch from outside the test frame", {
  # `UseMethod()` searches the frame the generic was called from before it
  # consults the method table, so a generic called directly inside `test_that()`
  # would find a method that is only a local function. These calls are made from
  # an environment that cannot see this frame, so they can only dispatch through
  # the table.
  res <- pooled_binary()

  expect_identical(dispatch_from_baseenv(coef, res), coef(res))
  expect_identical(dispatch_from_baseenv(vcov, res), vcov(res))
  expect_identical(dispatch_from_baseenv(nobs, res), 20L)
  expect_s3_class(dispatch_from_baseenv(as.data.frame, res), "data.frame")
})

test_that("the pooled result does not inherit the unpooled methods", {
  # `ipw_pooled` deliberately does not inherit from `ipw`. The two classes share
  # column names and differ in what the numbers under them mean, so a pooled
  # result reaching `confint.ipw()` would come back with normal bounds rather
  # than t ones and nothing would say so.
  res <- pooled_binary()

  expect_identical(class(res), "ipw_pooled")
  expect_false(inherits(res, "ipw"))
  # The fields `ipw` methods read are not there to be read.
  expect_null(res$outcome_mod)
  expect_null(res$wt_mod)
  expect_null(res$fit)

  # Each generic resolves to a method of its own for the pooled class rather
  # than to the unpooled one. The non-inheritance above is what makes that
  # necessary, and this is what says it was done: a pooled frame carries every
  # column `confint.ipw()` reads, so that method would answer for one without
  # erroring, with normal bounds where the pooled inference calls for t.
  for (generic in c("print", "coef", "vcov", "confint", "nobs")) {
    expect_false(
      identical(
        utils::getS3method(generic, "ipw_pooled"),
        utils::getS3method(generic, "ipw")
      ),
      info = generic
    )
  }
})
