# The accessors give an `ipw` result the interface a fitted model has:
# `coef()`, `vcov()`, `confint()`, `nobs()`, `df.residual()`, and `weights()`.
# They live here for the same reason `print.ipw()` does. Two packages each
# registering `coef.ipw()` would collide in the shared S3 method table, and a
# caller writing against a result would then get whichever package was installed
# last rather than the contract.
#
# Everything the accessors read is part of the `new_ipw()` contract: the
# `estimates` frame, the `ipw_vcov` attribute attached to it, `outcome_mod`, and
# `fit`. They never branch on `se_method`, and they reach into `fit` only
# through ordinary S3 dispatch. That is what lets a package whose variance
# object is a bare list use them, and it is asserted directly below rather than
# left to follow from the other tests.
#
# The estimates frames are written out literally, in the shape the `ipw()`
# return contract documents, rather than produced by fitting a model. Nothing
# here needs a real estimator, and this package holds no model-fitting
# dependency to build one with. The binary, categorical, and continuous frames
# hold the same numbers as their counterparts in `test-ipw-result.R`, so that
# the printed form and the accessors are asserted against one set of estimates.

# ---- fixtures ----------------------------------------------------------------

# A binary-exposure estimates frame: one row per effect measure, no `comparison`
# column, so the effect labels stand alone.
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

# A continuous-outcome estimates frame: a difference in means and nothing else,
# so every accessor has to work on a single row.
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

# A result whose stored interval is not the normal one: the bounds are
# asymmetric about the estimate, the way a bootstrap or profile interval is, and
# the stored level is 0.9 rather than the 0.95 that `confint()` defaults to.
# Both facts are what make "the stored bounds come back at the stored level"
# observable. No recomputation from `estimate` and `std.err` produces these
# numbers, and the level they are stored at is not the one a caller who names no
# level asks for.
stored_interval_estimates <- function() {
  data.frame(
    effect = c("rd", "log(rr)", "log(or)"),
    estimate = c(0.199882, 0.560414, 0.878313),
    std.err = c(0.092425, 0.273519, 0.418661),
    z = c(2.1626, 2.0489, 2.0979),
    ci.lower = c(0.061304, 0.128951, 0.207738),
    ci.upper = c(0.372015, 1.043772, 1.622904),
    conf.level = 0.9,
    p.value = c(0.030570, 0.040470, 0.035910)
  )
}

# The effect labels each frame's rows carry, written out rather than pasted
# together. The label rule is what these tests are for, so restating it with the
# same `paste()` the implementation uses would assert nothing.
binary_labels <- function() {
  c("rd", "log(rr)", "log(or)")
}

categorical_labels <- function() {
  c(
    "rd b vs a",
    "log(rr) b vs a",
    "log(or) b vs a",
    "rd c vs a",
    "log(rr) c vs a",
    "log(or) c vs a"
  )
}

# A covariance matrix of the reported effects, in the shape a fitting package
# attaches. The correlation falls off with the distance between rows, which is
# an ordinary shape for effects estimated from the same weighted means and, more
# to the point here, gives non-zero off-diagonal entries: a `vcov()` that built
# `diag(std.err^2)` from the frame instead of reading the attribute would agree
# on the diagonal and differ everywhere else.
#
# `outer()` makes the product symmetric to the bit, because multiplication is
# commutative in floating point even though it is not associative. That is what
# lets the symmetry assertion below use `expect_identical()`.
effects_vcov <- function(std_err, labels, rho = 0.9) {
  index <- seq_along(std_err)
  correlation <- rho^abs(outer(index, index, "-"))
  covariance <- correlation * outer(std_err, std_err)
  dimnames(covariance) <- list(labels, labels)
  covariance
}

binary_vcov <- function() {
  effects_vcov(binary_estimates()$std.err, binary_labels())
}

categorical_vcov <- function() {
  effects_vcov(categorical_estimates()$std.err, categorical_labels())
}

continuous_vcov <- function() {
  effects_vcov(continuous_estimates()$std.err, "diff")
}

# Fixed rather than simulated data, so the models never move between runs.
ipw_data <- function() {
  data.frame(
    x = rep(c(-1.5, -0.5, 0.5, 1.5), each = 5),
    z = rep(c(0, 1), 10),
    y = rep(c(0, 1, 1, 0, 1), 4)
  )
}

# Weights the outcome model is fitted with. They are named, which is the part
# that matters: `weights()` returns what the model frame stored, so the names
# coming back is what separates returning the stored vector from rebuilding one
# of the right length.
outcome_weights <- function() {
  setNames(rep(c(0.5, 1.5, 1, 2), each = 5), paste0("u", 1:20))
}

# A variance object that is not a fitted model: a named parameter vector and its
# covariance in a bare list, which is the shape an M-estimator from outside the
# ecosystem can take. `df.residual()` has nothing for it and `$vcov` is right
# there to be picked up by mistake, so it is the fixture that pins both of those
# contracts.
foreign_fit <- function() {
  theta <- c(mu0 = 0.352941, mu1 = 0.552823)
  list(
    theta = theta,
    vcov = matrix(
      c(0.004112, 0.000317, 0.000317, 0.005238),
      nrow = 2,
      dimnames = list(names(theta), names(theta))
    )
  )
}

# Assemble a result the way a method would, so that each test names only the
# part it is about. The covariance is attached to the `estimates` frame rather
# than passed to `new_ipw()`, which is where the contract puts it.
ipw_result <- function(
  estimates,
  vcov = NULL,
  fit = NULL,
  se_method = "mestimation",
  wts = outcome_weights()
) {
  if (!is.null(vcov)) {
    attr(estimates, "ipw_vcov") <- vcov
  }
  dat <- ipw_data()
  new_ipw(
    estimand = "ate",
    wt_mod = glm(z ~ x, family = binomial(), data = dat),
    outcome_mod = glm(
      y ~ z,
      family = quasibinomial(),
      data = dat,
      weights = wts
    ),
    estimates = estimates,
    se_method = se_method,
    fit = fit
  )
}

# Whether `print()` writes a row labelled exactly `label`.
#
# The estimate columns are numeric, so a row label is a line whose remainder
# after the label is spaces and then a number. Testing only that a line starts
# with the label would hold for a prefix of the real one: a `coef()` that
# labelled a categorical result `rd` would still match the line `print()` writes
# for `rd b vs a`, and the two surfaces would have drifted with nothing to say
# so. `printCoefmat()` also wraps significance stars onto lines of their own
# under the same labels when the frame is wide, and requiring a number is what
# keeps those out.
labels_a_printed_row <- function(out, label) {
  candidates <- out[startsWith(out, label)]
  remainder <- substring(candidates, nchar(label) + 1L)
  any(grepl("^ +-?[0-9]", remainder))
}

# ---- coef.ipw() --------------------------------------------------------------

test_that("coef() returns the estimates named by effect", {
  res <- ipw_result(binary_estimates())

  expect_identical(
    coef(res),
    c("rd" = 0.199882, "log(rr)" = 0.560414, "log(or)" = 0.878313)
  )
})

test_that("coef() keys a categorical result by effect and comparison", {
  # The effect labels repeat across contrasts, so `effect` alone would name
  # three of the six rows twice over and a caller could not tell which
  # comparison an estimate belonged to. The label is the two columns together.
  res <- ipw_result(categorical_estimates())

  expect_identical(
    coef(res),
    c(
      "rd b vs a" = 0.081945,
      "log(rr) b vs a" = 0.168870,
      "log(or) b vs a" = 0.328762,
      "rd c vs a" = 0.166939,
      "log(rr) c vs a" = 0.318293,
      "log(or) c vs a" = 0.676435
    )
  )
})

test_that("coef() returns a length-one vector for a single-row result", {
  res <- ipw_result(continuous_estimates())

  expect_identical(coef(res), c("diff" = 2.25255))
})

test_that("coef() names its result the way print() labels its rows", {
  # The two surfaces read the same rule, so a change to one that missed the
  # other would leave `coef(res)["rd b vs a"]` naming an estimate that the
  # printed table labels something else. This is the assertion that ties them
  # together; nothing else here would notice.
  for (estimates in list(
    binary_estimates(),
    categorical_estimates(),
    continuous_estimates()
  )) {
    res <- ipw_result(estimates)
    out <- capture.output(print(res))
    labels <- names(coef(res))

    expect_length(labels, nrow(estimates))

    unlabelled <- labels[
      !vapply(labels, function(l) labels_a_printed_row(out, l), logical(1))
    ]
    expect_identical(unlabelled, character())
  }
})

# ---- vcov.ipw() --------------------------------------------------------------

test_that("vcov() returns the covariance attached to the estimates", {
  covariance <- binary_vcov()
  res <- ipw_result(binary_estimates(), vcov = covariance)

  expect_identical(vcov(res), covariance)
  expect_identical(dimnames(vcov(res)), list(binary_labels(), binary_labels()))
  # A covariance matrix is symmetric, and the attribute is returned unaltered,
  # so the matrix that comes back is too. An implementation that reordered rows
  # without reordering columns would fail here and nowhere else.
  expect_identical(vcov(res), t(vcov(res)))
})

test_that("vcov() returns the covariance of a categorical result", {
  covariance <- categorical_vcov()
  res <- ipw_result(categorical_estimates(), vcov = covariance)

  expect_identical(vcov(res), covariance)
  expect_identical(dim(vcov(res)), c(6L, 6L))
  expect_identical(rownames(vcov(res)), categorical_labels())
})

test_that("vcov() returns a one by one matrix for a single-row result", {
  covariance <- continuous_vcov()
  res <- ipw_result(continuous_estimates(), vcov = covariance)

  expect_identical(vcov(res), covariance)
  expect_identical(dim(vcov(res)), c(1L, 1L))
})

test_that("coef(), vcov(), and confint() agree on the effect labels", {
  # The three surfaces are used together: a caller who reads a covariance out by
  # name expects the name to be the one `coef()` gave. Each is asserted against
  # a literal above; this is the assertion that they are the same literal.
  res <- ipw_result(categorical_estimates(), vcov = categorical_vcov())
  labels <- names(coef(res))

  expect_identical(rownames(vcov(res)), labels)
  expect_identical(colnames(vcov(res)), labels)
  expect_identical(rownames(confint(res)), labels)
})

test_that("vcov() errors when no covariance was attached", {
  # Older results and packages that have not adopted the contract carry no
  # `ipw_vcov` attribute. There is no honest fallback: the standard errors in
  # the frame give the diagonal and say nothing about the off-diagonal entries,
  # so a matrix built from them would report zero covariance between effects
  # that are estimated from the same data and are strongly correlated. The
  # caller is told instead.
  res <- ipw_result(binary_estimates())

  expect_error(vcov(res), class = "causalgenerics_no_vcov_ipw")
  expect_error(vcov(res), class = "causalgenerics_no_vcov")
  expect_snapshot(error = TRUE, vcov(res))
})

test_that("vcov() reads the attribute rather than the fit", {
  # The attribute is what the contract names, not the fit. Both halves are
  # needed to say so. A fit that has a covariance of its own must not satisfy a
  # missing attribute, and a result with no fit at all must still answer when
  # the attribute is there.
  #
  # The first fixture's fit is an `lm`, which `stats::vcov()` answers for, and
  # the second's is a bare list carrying a `vcov` element that `$` would hand
  # over. Neither is the covariance of the reported effects.
  dat <- ipw_data()

  from_lm <- ipw_result(binary_estimates(), fit = lm(y ~ x, data = dat))
  expect_error(vcov(from_lm), class = "causalgenerics_no_vcov")

  from_list <- ipw_result(binary_estimates(), fit = foreign_fit())
  expect_error(vcov(from_list), class = "causalgenerics_no_vcov")

  covariance <- binary_vcov()
  attached <- ipw_result(binary_estimates(), vcov = covariance, fit = NULL)
  expect_identical(vcov(attached), covariance)
})

# ---- confint.ipw() -----------------------------------------------------------

test_that("confint() returns the stored bounds at the stored level", {
  # The bounds in the frame are the ones the fitting package reported, and at
  # the level it reported them for they are returned rather than rebuilt. The
  # fixture's stored bounds are rounded to six decimals, so they differ from a
  # normal recomputation in the ninth; `expect_identical()` is what makes that
  # difference matter and separates the two implementations.
  estimates <- binary_estimates()
  res <- ipw_result(estimates)

  ci <- confint(res)

  expect_true(is.matrix(ci))
  expect_identical(dim(ci), c(3L, 2L))
  expect_identical(
    dimnames(ci),
    list(binary_labels(), c("2.5 %", "97.5 %"))
  )
  expect_identical(ci[, 1], setNames(estimates$ci.lower, binary_labels()))
  expect_identical(ci[, 2], setNames(estimates$ci.upper, binary_labels()))
})

test_that("confint() recomputes at a level other than the stored one", {
  # The half width is `qnorm(1 - (1 - level) / 2) * std.err`, and the spelling
  # of the quantile is part of that. `qnorm()` is not exactly antisymmetric:
  # `qnorm(0.05)` and `-qnorm(0.95)` differ in the last bit, so an
  # implementation that took the lower bound from the lower tail directly would
  # miss this assertion by one unit in the last place. The two bounds here come
  # from one half width, added and subtracted.
  estimates <- binary_estimates()
  res <- ipw_result(estimates)

  ci <- confint(res, level = 0.9)
  half_width <- qnorm(1 - (1 - 0.9) / 2) * estimates$std.err

  expect_identical(
    dimnames(ci),
    list(binary_labels(), c("5 %", "95 %"))
  )
  expect_identical(
    ci[, 1],
    setNames(estimates$estimate - half_width, binary_labels())
  )
  expect_identical(
    ci[, 2],
    setNames(estimates$estimate + half_width, binary_labels())
  )
  # The stored bounds are for 0.95 and have no business appearing here.
  expect_false(identical(
    ci[, 1],
    setNames(estimates$ci.lower, binary_labels())
  ))
})

test_that("confint() reads the stored level rather than assuming 0.95", {
  # A result stored at 0.9. Asking for 0.9 returns what it stored, and asking
  # for anything else recomputes, including the 0.95 that a caller who names no
  # level gets. The fixture's stored bounds are asymmetric about the estimate,
  # so no recomputation produces them and the first half cannot pass by
  # coincidence.
  estimates <- stored_interval_estimates()
  res <- ipw_result(estimates)

  stored <- confint(res, level = 0.9)

  expect_identical(
    dimnames(stored),
    list(binary_labels(), c("5 %", "95 %"))
  )
  expect_identical(stored[, 1], setNames(estimates$ci.lower, binary_labels()))
  expect_identical(stored[, 2], setNames(estimates$ci.upper, binary_labels()))

  default <- confint(res)
  half_width <- qnorm(1 - (1 - 0.95) / 2) * estimates$std.err

  expect_identical(
    dimnames(default),
    list(binary_labels(), c("2.5 %", "97.5 %"))
  )
  expect_identical(
    default[, 1],
    setNames(estimates$estimate - half_width, binary_labels())
  )
})

test_that("confint() labels its columns the way stats does", {
  # `confint()` methods across base R label the columns with the tail
  # probabilities as percentages, and a caller reading a bound out by column
  # name is entitled to the same labels here.
  res <- ipw_result(binary_estimates())

  expect_identical(colnames(confint(res, level = 0.99)), c("0.5 %", "99.5 %"))
  expect_identical(colnames(confint(res, level = 0.8)), c("10 %", "90 %"))
})

test_that("confint() selects rows by label", {
  estimates <- binary_estimates()
  res <- ipw_result(estimates)

  one <- confint(res, parm = "log(rr)")

  expect_identical(dim(one), c(1L, 2L))
  expect_identical(rownames(one), "log(rr)")
  expect_identical(one[1, 1], estimates$ci.lower[2])

  # The rows come back in the order they were asked for, not in the order the
  # frame stores them, which is what a caller who wrote the labels down expects.
  reordered <- confint(res, parm = c("log(or)", "rd"))

  expect_identical(rownames(reordered), c("log(or)", "rd"))
  expect_identical(
    reordered[, 1],
    setNames(estimates$ci.lower[c(3, 1)], c("log(or)", "rd"))
  )
})

test_that("confint() selects rows by position", {
  estimates <- binary_estimates()
  res <- ipw_result(estimates)

  expect_identical(rownames(confint(res, parm = 2)), "log(rr)")
  expect_identical(
    rownames(confint(res, parm = c(3, 1))),
    c("log(or)", "rd")
  )
  expect_identical(
    confint(res, parm = c(3, 1))[, 2],
    setNames(estimates$ci.upper[c(3, 1)], c("log(or)", "rd"))
  )
})

test_that("confint() selects a categorical row by its full label", {
  # The label is effect and comparison together, so `parm = "rd"` names nothing
  # in a categorical result even though `rd` is a value of the `effect` column.
  estimates <- categorical_estimates()
  res <- ipw_result(estimates)

  ci <- confint(res, parm = "rd c vs a")

  expect_identical(rownames(ci), "rd c vs a")
  expect_identical(ci[1, 1], estimates$ci.lower[4])

  expect_error(
    confint(res, parm = "rd"),
    class = "causalgenerics_invalid_argument_parm"
  )
})

test_that("confint() errors on labels the result does not have", {
  # A caller who mistypes a label gets an error naming the labels that are
  # available. Dropping the unmatched ones and returning the rest would answer a
  # question that was not asked, with a matrix of the wrong number of rows and
  # nothing to say why.
  res <- ipw_result(binary_estimates())

  expect_error(
    confint(res, parm = "rr"),
    class = "causalgenerics_invalid_argument_parm"
  )
  expect_error(
    confint(res, parm = "rr"),
    class = "causalgenerics_invalid_argument"
  )
  # A mix of matched and unmatched labels is refused as well: the matched ones
  # do not make the request answerable.
  expect_error(
    confint(res, parm = c("rd", "rr")),
    class = "causalgenerics_invalid_argument_parm"
  )

  expect_snapshot(error = TRUE, confint(res, parm = "rr"))
  expect_snapshot(error = TRUE, confint(res, parm = c("rd", "rr")))
})

# ---- nobs.ipw() and df.residual.ipw() ----------------------------------------

test_that("nobs() counts the outcome model's observations", {
  # The outcome model is fitted on sixteen of the twenty rows, the way it would
  # be after trimming, so a method that read `wt_mod` would report twenty and
  # fail here. The count is an integer, as it is for every model `nobs()`
  # answers for.
  dat <- ipw_data()
  res <- new_ipw(
    estimand = "ate",
    wt_mod = glm(z ~ x, family = binomial(), data = dat),
    outcome_mod = glm(y ~ z, family = quasibinomial(), data = dat[1:16, ]),
    estimates = binary_estimates(),
    se_method = "mestimation",
    fit = NULL
  )

  expect_identical(nobs(res), 16L)
  expect_type(nobs(res), "integer")
  expect_identical(nobs(ipw_result(binary_estimates())), 20L)
})

test_that("df.residual() reports the residual degrees of freedom of the fit", {
  # Twenty rows and two estimated parameters. The number differs from `nobs()`,
  # so a method that answered with the sample size would fail here.
  dat <- ipw_data()
  res <- ipw_result(binary_estimates(), fit = lm(y ~ x, data = dat))

  expect_identical(df.residual(res), 18L)
  expect_type(df.residual(res), "integer")
})

test_that("df.residual() is NA when there is no fit", {
  # The linearization path records `fit = NULL`, and a result with no fitted
  # variance object has no residual degrees of freedom to report. `NA_integer_`
  # rather than `NULL` keeps the return type the same whichever path produced
  # the result.
  res <- ipw_result(binary_estimates(), fit = NULL)

  expect_identical(df.residual(res), NA_integer_)
  expect_type(df.residual(res), "integer")
})

test_that("df.residual() is NA for a fit the generic has nothing for", {
  # A bare list reaches `df.residual.default()`, which returns `NULL` because
  # the list has no `df.residual` element. That is not an error and must not
  # become one: an M-estimator from outside the ecosystem is a legitimate value
  # for `fit`, and every other accessor works on the result carrying it.
  res <- ipw_result(binary_estimates(), fit = foreign_fit())

  expect_no_error(df.residual(res))
  expect_identical(df.residual(res), NA_integer_)
})

test_that("df.residual() converts a fit's answer to an integer safely", {
  # Two things a registered method can hand back that `as.integer()` alone
  # mishandles. A whole number stored as a double has to come back as an
  # integer, and `Inf` has to come back as `NA_integer_` without the "NAs
  # introduced by coercion" warning that `as.integer(Inf)` raises. A warning
  # would surface as noise in every run of a suite downstream, so the guard is
  # asserted here rather than left to a code reading.
  local_s3_method("df.residual", "cg_double_df", function(object, ...) 12)
  local_s3_method("df.residual", "cg_infinite_df", function(object, ...) Inf)
  local_s3_method("df.residual", "cg_empty_df", function(object, ...) numeric())

  as_double <- ipw_result(
    binary_estimates(),
    fit = structure(list(), class = "cg_double_df")
  )
  expect_identical(df.residual(as_double), 12L)

  infinite <- ipw_result(
    binary_estimates(),
    fit = structure(list(), class = "cg_infinite_df")
  )
  expect_no_warning(expect_identical(df.residual(infinite), NA_integer_))

  empty <- ipw_result(
    binary_estimates(),
    fit = structure(list(), class = "cg_empty_df")
  )
  expect_identical(df.residual(empty), NA_integer_)
})

# ---- weights.ipw() -----------------------------------------------------------

test_that("weights() returns the outcome model's weights as stored", {
  # The weights come back the way the model frame holds them, names and all. A
  # concrete weight class such as `psw` survives the same route, which is the
  # reason the contract is written as "as stored" rather than "as a double": a
  # caller who weighted with one gets it back and can read its estimand.
  wts <- outcome_weights()
  res <- ipw_result(binary_estimates(), wts = wts)

  expect_identical(weights(res), wts)
})

test_that("weights() reads the outcome model rather than the weighting model", {
  # Both models here are weighted, and with different vectors. Only the outcome
  # model's weights are the ones the estimate was computed under.
  dat <- ipw_data()
  wt_mod_weights <- rep(3, 20)
  outcome_mod_weights <- outcome_weights()

  res <- new_ipw(
    estimand = "ate",
    wt_mod = glm(
      z ~ x,
      family = binomial(),
      data = dat,
      weights = wt_mod_weights
    ),
    outcome_mod = glm(
      y ~ z,
      family = quasibinomial(),
      data = dat,
      weights = outcome_mod_weights
    ),
    estimates = binary_estimates(),
    se_method = "mestimation",
    fit = NULL
  )

  expect_identical(weights(res), outcome_mod_weights)
})

test_that("weights() is NULL for an unweighted outcome model", {
  # An unweighted model frame holds no weights column, so there is nothing to
  # return. This is the honest answer rather than a vector of ones, which would
  # claim the model was weighted.
  dat <- ipw_data()
  res <- new_ipw(
    estimand = "ate",
    wt_mod = glm(z ~ x, family = binomial(), data = dat),
    outcome_mod = glm(y ~ z, family = quasibinomial(), data = dat),
    estimates = binary_estimates(),
    se_method = "mestimation",
    fit = NULL
  )

  expect_null(weights(res))
})

# ---- the contract the accessors read -----------------------------------------

test_that("every accessor works on a result whose fit is not a model", {
  # The variance object a package records in `fit` is its own business. A bare
  # list is the shape that breaks an accessor which assumes a fitted model, and
  # apart from `df.residual()`, which has nothing to report for it, every
  # accessor answers exactly as it would with any other fit.
  covariance <- binary_vcov()
  res <- ipw_result(
    binary_estimates(),
    vcov = covariance,
    fit = foreign_fit()
  )

  expect_identical(
    coef(res),
    c("rd" = 0.199882, "log(rr)" = 0.560414, "log(or)" = 0.878313)
  )
  expect_identical(vcov(res), covariance)
  expect_identical(rownames(confint(res)), binary_labels())
  expect_identical(nobs(res), 20L)
  expect_identical(weights(res), outcome_weights())
  expect_identical(df.residual(res), NA_integer_)
})

test_that("the accessors do not branch on se_method", {
  # `se_method` records how the standard errors were computed, and the fields
  # the accessors read hold the result of that computation. Reading the label
  # would make an accessor answer differently for two results that carry the
  # same numbers, and would need updating every time a package added a method.
  covariance <- binary_vcov()
  estimates <- binary_estimates()

  mestimation <- ipw_result(
    estimates,
    vcov = covariance,
    se_method = "mestimation"
  )
  linearization <- ipw_result(
    estimates,
    vcov = covariance,
    se_method = "linearization"
  )

  expect_identical(coef(mestimation), coef(linearization))
  expect_identical(vcov(mestimation), vcov(linearization))
  expect_identical(confint(mestimation), confint(linearization))
  expect_identical(nobs(mestimation), nobs(linearization))
  expect_identical(weights(mestimation), weights(linearization))
})

# ---- registration ------------------------------------------------------------

test_that("the accessors are registered against ipw", {
  # Every assertion above reaches these methods from the test frame, and
  # `UseMethod()` searches the calling frame before it consults the method
  # table. All of them therefore pass whether or not the NAMESPACE carries an
  # `S3method()` directive. Downstream packages call `coef()` and the rest from
  # their own namespaces, where the table is the only route, so table membership
  # is the claim that matters and nothing else here makes it.
  #
  # All six generics live in stats, so their entries belong to stats' table
  # rather than to this package's.
  generics <- c("coef", "vcov", "confint", "nobs", "df.residual", "weights")

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
