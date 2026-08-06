# The accessors give an `ipw` result the interface a fitted model has:
# `coef()`, `vcov()`, `confint()`, `nobs()`, `df.residual()`, `weights()`, and
# `model.frame()`, alongside `estimand()`, which reports the estimand the
# weights targeted. They live here for the same reason `print.ipw()` does. Two
# packages each registering `coef.ipw()` would collide in the shared S3 method
# table, and a caller writing against a result would then get whichever package
# was installed last rather than the contract.
#
# Everything the accessors read is part of the `new_ipw()` contract: the
# `estimates` frame, the `ipw_vcov` attribute attached to it, `estimand`,
# `outcome_mod`, and `fit`. They never branch on `se_method`, and they reach
# into `fit` only through ordinary S3 dispatch. That is what lets a package
# whose variance object is a bare list use them, and it is asserted directly
# below rather than left to follow from the other tests.
#
# The estimates frames are written out literally, in the shape the `ipw()`
# return contract documents, rather than produced by fitting a model. Nothing
# here needs a real estimator, and this package holds no model-fitting
# dependency to build one with. The binary, categorical, and continuous frames
# hold the same numbers as their counterparts in `test-ipw-result.R`, so that
# the printed form and the accessors are asserted against one set of estimates.
#
# `coef()`, `vcov()`, and `confint()` also read the presentation mode the
# `effects` field records, and take an `effects` argument that names a reading
# for one call. The marginal reading is the one every assertion in the sections
# above is written against. The conditional reading reports the outcome model's
# coefficient surface, and its covariance is the corrected block a fitting
# package attached with `new_ipw_model()`, never the one the outcome model
# computed for itself. `nobs()`, `df.residual()`, `weights()`, `model.frame()`,
# and `estimand()` describe the fit rather than a surface of it, so they answer
# the same way in either mode.

# ---- fixtures ----------------------------------------------------------------

# A binary-exposure estimates frame: one row per effect measure, no `contrast`
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

# A categorical-exposure estimates frame as an earlier version of the contract
# stored one. The effect labels repeat across contrasts, so a column sits
# immediately after `effect` naming the non-reference and reference level of each
# one, and this vintage calls that column `comparison`. It is here to derive the
# canonical frame below from; the tests that are about the older name live in
# `test-ipw-result.R`, where the surfaces that read the column are.
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

# The same frame with the column naming its contrasts called `contrast`, which is
# what the constructor's contract calls it and the name `mice::pool()` and
# \pkg{marginaleffects} both read the tidier surface by. This is the categorical
# frame every assertion below is written against. It is derived from the frame
# above rather than written out again, so that the two differ in the column name
# and in nothing else.
contrast_estimates <- function() {
  estimates <- categorical_estimates()
  names(estimates)[names(estimates) == "comparison"] <- "contrast"
  estimates
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

# The names the outcome model's coefficient surface carries, which are the row
# labels of the conditional reading. Written out rather than read off the model
# for the reason the effect labels are: the claim is about which of the two
# surfaces comes back, and `names(coef(mod))` on both sides of an assertion
# would hold whichever one it was.
conditional_labels <- function() {
  c("(Intercept)", "z")
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
  effects_vcov(contrast_estimates()$std.err, categorical_labels())
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

# The outcome model every result below is built around. It is a function of its
# own rather than a call inside the result builder because the conditional
# reading is a statement about this model, and asserting it needs the model and
# the result at once. The fit is deterministic, so the model a test holds and the
# one the result carries are the same numbers.
outcome_model <- function(wts = outcome_weights()) {
  glm(y ~ z, family = quasibinomial(), data = ipw_data(), weights = wts)
}

# The corrected covariance of that model's coefficients, in the shape a fitting
# package hands to `new_ipw_model()`: the outcome block of the stacked sandwich,
# labelled with the coefficient names. It is built from the model's own
# covariance so that the dimnames agree without being restated, and scaled so
# that it differs from it. Accounting for the weights having been estimated
# ordinarily widens the block, and the scaling is what makes "the conditional
# covariance is the corrected one" separable from "the conditional covariance is
# whatever the outcome model computed for itself".
corrected_outcome_vcov <- function(mod = outcome_model()) {
  stats::vcov(mod) * 1.4
}

# The corrected block labelled in the reverse of the coefficient order, which is
# the shape a fitting package produces when its stacked system holds the outcome
# parameters the other way round. The labels are what say which coefficient each
# entry belongs to, so a block read by position rather than by name reports the
# intercept's variance for the slope and the slope's for the intercept. The two
# diagonal entries are far apart and their square roots are exact doubles, which
# is what lets the standard errors below be written down.
reversed_outcome_vcov <- function() {
  matrix(
    c(0.0625, 0.03125, 0.03125, 0.25),
    nrow = 2,
    dimnames = list(c("z", "(Intercept)"), c("z", "(Intercept)"))
  )
}

# That same block in coefficient order, written out rather than permuted from
# the fixture above. The reordering is the claim, so restating it with the
# indexing the implementation uses would assert nothing.
coefficient_order_vcov <- function() {
  matrix(
    c(0.25, 0.03125, 0.03125, 0.0625),
    nrow = 2,
    dimnames = list(c("(Intercept)", "z"), c("(Intercept)", "z"))
  )
}

# A block labelled with the parameter names of a stacked system rather than with
# the model's coefficient names. `new_ipw_model()` takes it, deliberately: only
# the fitting package knows which block of the sandwich belongs to which model,
# so the constructor checks that both margins are labelled rather than what the
# labels say. Where the block meets the coefficients is where the pairing has to
# be made, and there is none to make here.
theta_outcome_vcov <- function() {
  matrix(
    c(0.0625, 0.03125, 0.03125, 0.25),
    nrow = 2,
    dimnames = list(c("theta1", "theta2"), c("theta1", "theta2"))
  )
}

# A three by three block for a model with two coefficients, labelled with both
# coefficient names and one more. The extra term is what makes it discriminating:
# indexing it by the coefficient names alone gives a two by two matrix, so an
# implementation that only reordered would answer with the covariance of two
# coefficients of a different fit.
oversized_outcome_vcov <- function() {
  matrix(
    c(0.25, 0.03125, 0, 0.03125, 0.0625, 0, 0, 0, 0.5),
    nrow = 3,
    dimnames = list(
      c("(Intercept)", "z", "x"),
      c("(Intercept)", "z", "x")
    )
  )
}

# A conditional result whose outcome model reports its coefficients without
# names, carrying a block that is labelled. The two cannot be paired at all:
# there is nothing to match the labels against, and matching them by position is
# what the labels exist to prevent. The class carries no `coef()` method of its
# own, so the test that uses this registers one.
unnamed_coef_result <- function() {
  new_ipw(
    estimand = "ate",
    wt_mod = glm(z ~ x, family = binomial(), data = ipw_data()),
    outcome_mod = new_ipw_model(
      structure(list(), class = "cg_unnamed_coef"),
      coefficient_order_vcov()
    ),
    estimates = binary_estimates(),
    se_method = "mestimation",
    fit = NULL,
    effects = "conditional"
  )
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
# part it is about. The covariance of the effects is attached to the `estimates`
# frame rather than passed to `new_ipw()`, which is where the contract puts it,
# and the covariance of the outcome model's coefficients goes on the model
# itself, which is where `new_ipw_model()` puts it. A result built without
# `outcome_vcov` carries a plain `glm`, which is what a fitting package whose
# variance estimator has no stacked system to take the block from produces.
ipw_result <- function(
  estimates,
  vcov = NULL,
  fit = NULL,
  se_method = "mestimation",
  wts = outcome_weights(),
  effects = "marginal",
  outcome_vcov = NULL
) {
  if (!is.null(vcov)) {
    attr(estimates, "ipw_vcov") <- vcov
  }
  outcome_mod <- outcome_model(wts)
  if (!is.null(outcome_vcov)) {
    outcome_mod <- new_ipw_model(outcome_mod, outcome_vcov)
  }
  new_ipw(
    estimand = "ate",
    wt_mod = glm(z ~ x, family = binomial(), data = ipw_data()),
    outcome_mod = outcome_mod,
    estimates = estimates,
    se_method = se_method,
    fit = fit,
    effects = effects
  )
}

# A conditional result whose outcome model carries the wrapper class with no
# covariance behind it. `new_ipw_model()` cannot produce one: it validates the
# matrix before it prepends the class, so the two go on together. An object of
# this shape was assembled some other way, and it is a different case from a
# model that was never wrapped: that one is a fitting package which has not
# adopted the contract, and this one is a model that claims to carry a block and
# does not. The accessors have to tell the caller which of the two they met.
stripped_wrapper_result <- function(estimates = binary_estimates()) {
  res <- ipw_result(
    estimates,
    vcov = binary_vcov(),
    outcome_vcov = corrected_outcome_vcov(),
    effects = "conditional"
  )
  attr(res$outcome_mod, "ipw_vcov") <- NULL
  res
}

# A result of the shape the constructor built before the mode existed: six
# fields and no `effects`. Results stored from an earlier version of a fitting
# package have this shape, as does one built by hand against the earlier
# contract. It is derived from the seven-field result rather than written out
# again so that the two differ in the field and in nothing else.
legacy_result <- function(estimates = binary_estimates(), vcov = NULL) {
  fields <- unclass(ipw_result(estimates, vcov = vcov))
  fields$effects <- NULL
  structure(fields, class = "ipw")
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

# The columns the outcome model's frame holds once the weights column is gone:
# the response and the one term, as `ipw_data()` supplies them. Written out
# rather than read back off the model, for the reason the effect labels are.
# Which columns come back is the claim, and a comparison against the frame the
# model stored would hold whichever ones they were.
outcome_frame_columns <- function() {
  data.frame(y = ipw_data()$y, z = ipw_data()$z)
}

# A model frame with its `terms` attribute removed.
#
# Selecting columns from a data frame drops that attribute, so a frame that had
# no weights column to drop is not `identical()` to the one the model stored
# even though every column, name, and row name agrees. Taking it off both sides
# keeps the comparison an `expect_identical()` without asserting anything about
# the attribute in either direction.
without_terms <- function(frame) {
  attr(frame, "terms") <- NULL
  frame
}

# ---- coef.ipw() --------------------------------------------------------------

test_that("coef() returns the estimates named by effect", {
  res <- ipw_result(binary_estimates())

  expect_identical(
    coef(res),
    c("rd" = 0.199882, "log(rr)" = 0.560414, "log(or)" = 0.878313)
  )
})

test_that("coef() keys a categorical result by effect and contrast", {
  # The effect labels repeat across contrasts, so `effect` alone would name
  # three of the six rows twice over and a caller could not tell which contrast
  # an estimate belonged to. The label is the two columns together.
  res <- ipw_result(contrast_estimates())

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
    contrast_estimates(),
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
  res <- ipw_result(contrast_estimates(), vcov = covariance)

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
  res <- ipw_result(contrast_estimates(), vcov = categorical_vcov())
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
  # normal recomputation in the sixth; `expect_identical()` is what makes that
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

test_that("confint() gives a numeric parm its ordinary subscript meaning", {
  # Positions are indexed rather than matched, so a numeric `parm` reads the way
  # a subscript reads everywhere else. Zero selects nothing and a negative
  # position drops the row it names. An implementation that accepted only
  # positive positions in range would refuse both, and a caller who reached for
  # either would be told the surface does not have a row it does have.
  estimates <- binary_estimates()
  res <- ipw_result(estimates)

  none <- confint(res, parm = 0)

  expect_true(is.matrix(none))
  expect_identical(dim(none), c(0L, 2L))
  # An empty selection is still a labelled matrix, so a caller who built one
  # from a filter that matched nothing can read its columns or bind it to
  # another without treating the case apart.
  expect_identical(colnames(none), c("2.5 %", "97.5 %"))

  dropped <- confint(res, parm = -1)

  expect_identical(dim(dropped), c(2L, 2L))
  expect_identical(rownames(dropped), c("log(rr)", "log(or)"))
  expect_identical(
    dropped[, 1],
    setNames(estimates$ci.lower[c(2, 3)], c("log(rr)", "log(or)"))
  )
})

test_that("confint() selects a categorical row by its full label", {
  # The label is effect and contrast together, so `parm = "rd"` names nothing in
  # a categorical result even though `rd` is a value of the `effect` column.
  estimates <- contrast_estimates()
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

test_that("confint() errors on positions the result does not have", {
  # The other half of the same guard. A position past the last row indexes
  # nothing, and letting it through would give back a row of `NA` limits under
  # an `NA` label, which reads as an effect that was estimated and came out
  # unknown rather than as a subscript that named nothing. The message names the
  # range rather than the labels, since a caller who wrote a number is working
  # from the positions.
  res <- ipw_result(binary_estimates())

  expect_error(
    confint(res, parm = 99),
    class = "causalgenerics_invalid_argument_parm"
  )
  expect_error(
    confint(res, parm = 99),
    class = "causalgenerics_invalid_argument"
  )
  # A real position beside one past the end is refused as well, for the reason a
  # mix of labels is: the position that matched does not make the request
  # answerable.
  expect_error(
    confint(res, parm = c(1, 99)),
    class = "causalgenerics_invalid_argument_parm"
  )
  # `NA` reaches this branch too, and it is the subscript that would go furthest
  # unnoticed. A logical `NA` recycles to the length of the surface, so indexing
  # with it gives one `NA` row per effect: a matrix of exactly the shape a
  # caller expects, holding nothing.
  expect_error(
    confint(res, parm = NA),
    class = "causalgenerics_invalid_argument_parm"
  )
  expect_error(
    confint(res, parm = NA_integer_),
    class = "causalgenerics_invalid_argument_parm"
  )

  expect_snapshot(error = TRUE, confint(res, parm = 99))
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
  # The answers a registered method can hand back that `as.integer()` alone
  # mishandles. A whole number stored as a double has to come back as an
  # integer, and `Inf` has to come back as `NA_integer_` without the "NAs
  # introduced by coercion" warning that `as.integer(Inf)` raises. A whole
  # number past the largest integer raises that same warning, and it has to come
  # back as the double it is, since the count itself is a real one. A warning
  # would surface as noise in every run of a suite downstream, so the guard is
  # asserted here rather than left to a code reading.
  local_s3_method("df.residual", "cg_double_df", function(object, ...) 12)
  local_s3_method("df.residual", "cg_infinite_df", function(object, ...) Inf)
  local_s3_method("df.residual", "cg_empty_df", function(object, ...) numeric())
  local_s3_method("df.residual", "cg_huge_df", function(object, ...) 1e10)

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

  huge <- ipw_result(
    binary_estimates(),
    fit = structure(list(), class = "cg_huge_df")
  )
  expect_no_warning(expect_identical(df.residual(huge), 1e10))
  expect_type(df.residual(huge), "double")
})

test_that("df.residual() reports a fractional answer as it stands", {
  # Residual degrees of freedom are not always a whole number. A fit that spends
  # a fractional count on penalized or smooth terms reports one, and so does a
  # small-sample correction, so `12.4` is what the variance object has to say
  # rather than an imprecise way of saying `12`. Truncating it reports a fit that
  # spent 0.4 more parameters than it did, and nothing in the returned value
  # would say the number had been altered on the way out.
  #
  # The whole-number cases are unchanged and are pinned in the test above, where
  # `12` comes back as `12L` and `Inf` and `numeric()` come back as
  # `NA_integer_`. This is the one answer that has no integer to become.
  local_s3_method("df.residual", "cg_fractional_df", function(object, ...) 12.4)

  res <- ipw_result(
    binary_estimates(),
    fit = structure(list(), class = "cg_fractional_df")
  )

  expect_identical(df.residual(res), 12.4)
  expect_type(df.residual(res), "double")
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

# ---- model.frame.ipw() -------------------------------------------------------

test_that("model.frame() drops the weights column from the outcome frame", {
  # The frame a weighted model stores carries a `(weights)` column that the
  # formula never named. It is the fitting machinery's bookkeeping rather than
  # data the caller supplied, and `weights()` reports it already, so the frame
  # comes back holding the columns the formula named and nothing else.
  res <- ipw_result(binary_estimates())
  stored <- stats::model.frame(res$outcome_mod)
  mf <- model.frame(res)

  expect_true("(weights)" %in% colnames(stored))
  expect_false("(weights)" %in% colnames(mf))
  expect_identical(colnames(mf), c("y", "z"))
  expect_identical(nrow(mf), nrow(stored))
  expect_identical(mf, outcome_frame_columns())
})

test_that("model.frame() returns every column of an unweighted frame", {
  # An unweighted model frame has no weights column, so there is nothing to drop
  # and the whole frame comes back. The `terms` attribute does not: selecting
  # columns from a data frame drops it whether or not a column went with it, and
  # it is taken off both sides here rather than asserted either way, since the
  # claim is about the columns.
  res <- ipw_result(binary_estimates(), wts = NULL)
  stored <- stats::model.frame(res$outcome_mod)
  mf <- model.frame(res)

  expect_false("(weights)" %in% colnames(stored))
  expect_identical(without_terms(mf), without_terms(stored))
  expect_identical(mf, outcome_frame_columns())
  expect_identical(nrow(mf), 20L)
})

test_that("model.frame() reads the outcome model wrapped or not", {
  # The wrapper carries a corrected covariance and registers nothing else, so
  # `model.frame()` walks past it to the model's own method the way `coef()`
  # does, and the frame is the same either way.
  wrapped <- ipw_result(
    binary_estimates(),
    outcome_vcov = corrected_outcome_vcov()
  )
  bare <- ipw_result(binary_estimates())

  expect_s3_class(wrapped$outcome_mod, "ipw_model")
  expect_false(inherits(bare$outcome_mod, "ipw_model"))
  expect_identical(model.frame(wrapped), model.frame(bare))
  expect_identical(model.frame(wrapped), outcome_frame_columns())
})

test_that("model.frame() lets the outcome model's error through", {
  # Delegation and nothing else, with no guard for a model that has no frame to
  # give. A model whose own method refuses says why it refuses, and an error of
  # this package's own raised in front of it would replace the one explanation a
  # caller can act on with one about a result that is not the problem.
  local_s3_method("model.frame", "cg_no_frame", function(formula, ...) {
    stop(errorCondition(
      "this model records no model frame",
      class = "cg_no_frame_error"
    ))
  })

  res <- new_ipw(
    estimand = "ate",
    wt_mod = glm(z ~ x, family = binomial(), data = ipw_data()),
    outcome_mod = structure(list(), class = "cg_no_frame"),
    estimates = binary_estimates(),
    se_method = "mestimation",
    fit = NULL
  )

  expect_error(model.frame(res), class = "cg_no_frame_error")
  expect_error(
    model.frame(res),
    "this model records no model frame",
    fixed = TRUE
  )
})

# ---- estimand.ipw() ----------------------------------------------------------

test_that("estimand() returns the estimand the result records", {
  # The field rather than a constant. The two results here differ in it and in
  # nothing else, so a method that answered `"ate"` for every result would pass
  # the first assertion and fail the second.
  expect_identical(estimand(ipw_result(binary_estimates())), "ate")

  att <- new_ipw(
    estimand = "att",
    wt_mod = glm(z ~ x, family = binomial(), data = ipw_data()),
    outcome_mod = outcome_model(),
    estimates = binary_estimates(),
    se_method = "mestimation",
    fit = NULL
  )

  expect_identical(estimand(att), "att")
})

test_that("estimand<-() has no method for a result", {
  # The estimand a result records is a fact about the weights its method
  # targeted and the estimates computed under them. Assigning a new one would
  # relabel those numbers rather than recompute them, so there is no method for
  # the replacement generic here and a caller who writes the assignment is told
  # as much. The absence is deliberate, and this is what pins it.
  res <- ipw_result(binary_estimates())

  expect_error(
    estimand(res) <- "att",
    class = "causalgenerics_no_method_estimand<-"
  )
  expect_error(
    estimand(res) <- "att",
    class = "causalgenerics_no_method"
  )

  expect_snapshot(error = TRUE, estimand(res) <- "att")
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

# ---- the conditional reading -------------------------------------------------

test_that("coef() reports the outcome coefficients in the conditional mode", {
  # No argument at the call site. A caller handed a conditional result asks for
  # its coefficients the ordinary way, and that is what lets a package which
  # calls `coef()` and `vcov()` on whatever it is given work with one: the mode
  # travels on the object, so the caller does not have to know about it.
  res <- ipw_result(
    binary_estimates(),
    vcov = binary_vcov(),
    outcome_vcov = corrected_outcome_vcov(),
    effects = "conditional"
  )

  expect_identical(coef(res), coef(res$outcome_mod))
  expect_identical(names(coef(res)), conditional_labels())
  expect_length(coef(res), 2L)
  # The two readings are different vectors, so the assertions above are about
  # the mode rather than about a result whose surfaces happen to agree.
  expect_false(identical(coef(res), coef(as_marginal(res))))
})

test_that("coef() reads the outcome model wrapped or not", {
  # The wrapper is about the covariance. `coef()` walks past it to the model's
  # own method, so the coefficient surface is there either way and the guard
  # `vcov()` raises in the conditional mode has no counterpart here.
  wrapped <- ipw_result(
    binary_estimates(),
    outcome_vcov = corrected_outcome_vcov(),
    effects = "conditional"
  )
  bare <- ipw_result(binary_estimates(), effects = "conditional")

  expect_identical(names(coef(bare)), conditional_labels())
  expect_identical(coef(bare), coef(wrapped))
  expect_no_error(coef(bare))
})

test_that("vcov() reports the corrected block in the conditional mode", {
  # The conditional covariance is the one the fitting package attached to the
  # outcome model, which is the block of the stacked sandwich that accounts for
  # the weights having been estimated. It comes back as attached, dimnames
  # included, so a caller reading a variance out by coefficient name gets the
  # name `coef()` gave in the same mode.
  corrected <- corrected_outcome_vcov()
  res <- ipw_result(
    binary_estimates(),
    vcov = binary_vcov(),
    outcome_vcov = corrected,
    effects = "conditional"
  )

  expect_identical(vcov(res), corrected)
  expect_identical(vcov(res), vcov(res$outcome_mod))
  expect_identical(dim(vcov(res)), c(2L, 2L))
  expect_identical(
    dimnames(vcov(res)),
    list(conditional_labels(), conditional_labels())
  )
  # Not the covariance the outcome model computed for itself, and not the
  # covariance of the effects either. Both are reachable from this result, and
  # neither is the answer to this question.
  expect_false(identical(vcov(res), stats::vcov(outcome_model())))
  expect_false(identical(vcov(res), binary_vcov()))
})

test_that("vcov() reorders a conditional block into coefficient order", {
  # A fitting package labels the block from the stacked system it solved, and
  # the order the parameters sit in there is its own. The labels are what pair
  # the block with the coefficients, so the matrix comes back in the order
  # `coef()` reports and a caller reading the diagonal alongside the
  # coefficients gets each variance against the coefficient it belongs to.
  #
  # The test above is the other half of this: a block already in coefficient
  # order comes back as attached, so the reordering is skipped rather than
  # applied to a matrix that needs none.
  res <- ipw_result(
    binary_estimates(),
    vcov = binary_vcov(),
    outcome_vcov = reversed_outcome_vcov(),
    effects = "conditional"
  )

  expect_identical(vcov(res), coefficient_order_vcov())
  expect_identical(
    dimnames(vcov(res)),
    list(conditional_labels(), conditional_labels())
  )
  expect_identical(names(coef(res)), rownames(vcov(res)))
  # The variance each label names travels with it, which is what separates
  # reordering the matrix from relabelling it. Relabelling would produce a
  # matrix that satisfies every assertion about the dimnames while reporting
  # the intercept's variance for the slope.
  expect_identical(diag(vcov(res)), c("(Intercept)" = 0.25, "z" = 0.0625))
  # The block as attached is a different matrix, so the assertions above are
  # about the reordering rather than about a block that was already in order.
  expect_false(identical(vcov(res), reversed_outcome_vcov()))
})

test_that("vcov() refuses a conditional block labelled for other parameters", {
  # `new_ipw_model()` takes the block's labels on trust, since only the fitting
  # package knows which block of the sandwich belongs to which model. The
  # accessor is where they meet the coefficients, and labels that name a
  # stacked system's parameters cannot be paired with them at all. Pairing by
  # position instead would report a covariance of something else under this
  # model's coefficient names, which is the one answer a caller has no way to
  # check.
  res <- ipw_result(
    binary_estimates(),
    vcov = binary_vcov(),
    outcome_vcov = theta_outcome_vcov(),
    effects = "conditional"
  )

  expect_error(vcov(res), class = "causalgenerics_conditional_vcov_mismatch")
  expect_error(vcov(res), class = "causalgenerics_no_vcov")
  expect_error(
    confint(res),
    class = "causalgenerics_conditional_vcov_mismatch"
  )

  # Not the refusal the reading raises for a model that was never wrapped,
  # which says the fitting package attached no block and is answered by
  # wrapping the model. This model is wrapped and carries a block. The two
  # conditions share the general class, so telling them apart is a claim about
  # the specific one.
  cnd <- tryCatch(vcov(res), error = identity)
  expect_false(inherits(cnd, "causalgenerics_no_conditional_vcov"))

  # The coefficients are the model's own and need no block to report, so the
  # pairing is not something `coef()` has to make and this result is still one
  # a caller can read the estimates of.
  expect_identical(coef(res), coef(res$outcome_mod))

  expect_snapshot(error = TRUE, vcov(res))
})

test_that("vcov() refuses a conditional block of the wrong size", {
  # Three labels for two coefficients: the block of a model with another term
  # in it. Both coefficient names are among its labels, so an implementation
  # that indexed by name and stopped there would take the two by two corner and
  # answer with it, which is the covariance of two coefficients of a different
  # fit rather than of these.
  res <- ipw_result(
    binary_estimates(),
    vcov = binary_vcov(),
    outcome_vcov = oversized_outcome_vcov(),
    effects = "conditional"
  )

  expect_error(vcov(res), class = "causalgenerics_conditional_vcov_mismatch")
  expect_error(vcov(res), class = "causalgenerics_no_vcov")
  expect_error(
    confint(res),
    class = "causalgenerics_conditional_vcov_mismatch"
  )
})

test_that("vcov() refuses a conditional block it cannot pair by name", {
  # A model whose coefficients have no names leaves the labels nothing to match
  # against. The pairing is unverifiable rather than wrong, and assuming the
  # positions are it would be the assumption the labels exist to replace.
  local_s3_method("coef", "cg_unnamed_coef", function(object, ...) {
    c(-0.847298, 1.694596)
  })
  res <- unnamed_coef_result()

  expect_null(names(coef(res)))
  expect_error(vcov(res), class = "causalgenerics_conditional_vcov_mismatch")
  expect_error(vcov(res), class = "causalgenerics_no_vcov")
  expect_error(
    confint(res),
    class = "causalgenerics_conditional_vcov_mismatch"
  )
})

test_that("vcov() refuses the conditional mode without a corrected block", {
  # An outcome model that was not wrapped carries no covariance the two-step
  # estimation implies, and there is nothing to fall back on. The model has a
  # `vcov()` of its own, and it is the wrong answer rather than an approximate
  # one: it treats the estimated weights as fixed and reports an uncertainty the
  # estimate does not have. The caller is told instead.
  res <- ipw_result(
    binary_estimates(),
    vcov = binary_vcov(),
    effects = "conditional"
  )

  expect_error(vcov(res), class = "causalgenerics_no_conditional_vcov")
  expect_error(vcov(res), class = "causalgenerics_no_vcov")
  expect_error(
    vcov(as_marginal(res), effects = "conditional"),
    class = "causalgenerics_no_conditional_vcov"
  )
  # The fixture is discriminating: the naive matrix is right there to be
  # returned by an implementation that let the guard fall through.
  expect_no_error(stats::vcov(res$outcome_mod))

  expect_snapshot(error = TRUE, vcov(res))
})

test_that("vcov() refuses a conditional block the outcome model lost", {
  # A wrapper carrying nothing is refused for its own reason. The guard on the
  # class is passed, since the class is there, and the model itself is the one
  # that has nothing to report.
  res <- stripped_wrapper_result()

  expect_identical(class(res$outcome_mod)[[1]], "ipw_model")

  expect_error(vcov(res), class = "causalgenerics_no_vcov_ipw_model")
  expect_error(vcov(res), class = "causalgenerics_no_vcov")

  # Not the conditional error, which says the fitting package attached no block
  # and is answered by wrapping the model. This model was wrapped, so that
  # advice has nothing left to tell the caller to do, and the object is what is
  # wrong. The two conditions share the general class, so telling them apart is
  # a claim about the specific one.
  cnd <- tryCatch(vcov(res), error = identity)
  expect_false(inherits(cnd, "causalgenerics_no_conditional_vcov"))
})

test_that("confint() bounds the outcome coefficients in the conditional mode", {
  # The interval is the normal one built from the corrected block: the estimate
  # plus and minus one half width, taken from the square root of the diagonal.
  # `qnorm()` is not exactly antisymmetric, so a lower limit taken from the
  # lower tail instead would miss this by one unit in the last place.
  corrected <- corrected_outcome_vcov()
  res <- ipw_result(
    binary_estimates(),
    vcov = binary_vcov(),
    outcome_vcov = corrected,
    effects = "conditional"
  )

  ci <- confint(res)
  half_width <- qnorm(1 - (1 - 0.95) / 2) * sqrt(diag(corrected))
  estimate <- coef(res$outcome_mod)

  expect_true(is.matrix(ci))
  expect_identical(dim(ci), c(2L, 2L))
  expect_identical(
    dimnames(ci),
    list(conditional_labels(), c("2.5 %", "97.5 %"))
  )
  expect_identical(
    ci[, 1],
    setNames(estimate - half_width, conditional_labels())
  )
  expect_identical(
    ci[, 2],
    setNames(estimate + half_width, conditional_labels())
  )
})

test_that("confint() pairs a reordered block with its own coefficients", {
  # The half width of a coefficient's interval is taken from that coefficient's
  # variance, which is the entry the labels name rather than the entry in its
  # position. The block's diagonal here is 0.25 for the intercept and 0.0625
  # for the slope, so the standard errors are 0.5 and 0.25 exactly and the
  # limits can be written down.
  res <- ipw_result(
    binary_estimates(),
    vcov = binary_vcov(),
    outcome_vcov = reversed_outcome_vcov(),
    effects = "conditional"
  )

  ci <- confint(res)
  estimate <- coef(res$outcome_mod)
  half_width <- qnorm(1 - (1 - 0.95) / 2) * c(0.5, 0.25)

  expect_identical(
    dimnames(ci),
    list(conditional_labels(), c("2.5 %", "97.5 %"))
  )
  expect_identical(
    ci[, 1],
    setNames(estimate - half_width, conditional_labels())
  )
  expect_identical(
    ci[, 2],
    setNames(estimate + half_width, conditional_labels())
  )

  # The pairing the positions give is a matrix of the same shape carrying the
  # same labels, so nothing but the numbers tells the two apart. Here the
  # intervals it gives are twice and half the width of the right ones.
  crossed <- qnorm(1 - (1 - 0.95) / 2) * c(0.25, 0.5)
  expect_false(identical(
    ci[, 1],
    setNames(estimate - crossed, conditional_labels())
  ))
})

test_that("confint() recomputes in the conditional mode at every level", {
  # The stored limits belong to the marginal reading, and the level they are
  # stored at is the 0.95 asked for above, so the fast path has a row to offer
  # at every level a caller is likely to name. It has nothing to do here: the
  # conditional limits are built from the outcome coefficients whatever the
  # level, and the two-row matrix above is already a different shape from the
  # three stored effects.
  corrected <- corrected_outcome_vcov()
  res <- ipw_result(
    stored_interval_estimates(),
    vcov = binary_vcov(),
    outcome_vcov = corrected,
    effects = "conditional"
  )

  ci <- confint(res, level = 0.9)
  half_width <- qnorm(1 - (1 - 0.9) / 2) * sqrt(diag(corrected))
  estimate <- coef(res$outcome_mod)

  expect_identical(
    dimnames(ci),
    list(conditional_labels(), c("5 %", "95 %"))
  )
  expect_identical(
    ci[, 1],
    setNames(estimate - half_width, conditional_labels())
  )
  expect_identical(
    ci[, 2],
    setNames(estimate + half_width, conditional_labels())
  )
  # The columns are labelled the way the marginal reading labels them, since a
  # caller reading a limit out by column name is entitled to the same names in
  # either mode.
  expect_identical(
    colnames(confint(res, level = 0.99)),
    c("0.5 %", "99.5 %")
  )
})

test_that("confint() selects conditional rows by position and by name", {
  # `parm` matches the surface being reported, so in this mode it indexes and
  # names the outcome coefficients. The rows come back in the order they were
  # asked for rather than the order the model stores them.
  corrected <- corrected_outcome_vcov()
  res <- ipw_result(
    binary_estimates(),
    vcov = binary_vcov(),
    outcome_vcov = corrected,
    effects = "conditional"
  )
  half_width <- qnorm(1 - (1 - 0.95) / 2) * sqrt(diag(corrected))
  estimate <- coef(res$outcome_mod)

  expect_identical(rownames(confint(res, parm = 2)), "z")
  expect_identical(dim(confint(res, parm = 2)), c(1L, 2L))
  expect_identical(
    confint(res, parm = 2)[1, 1],
    unname(estimate - half_width)[2]
  )

  by_name <- confint(res, parm = c("z", "(Intercept)"))

  expect_identical(rownames(by_name), c("z", "(Intercept)"))
  expect_identical(
    by_name[, 2],
    setNames(
      unname(estimate + half_width)[c(2, 1)],
      c("z", "(Intercept)")
    )
  )
})

test_that("confint() empties a conditional selection the same way", {
  # A subscript means the same thing in either reading. The two build their
  # limits in different places in the method, so an empty selection that came
  # back as a labelled matrix in one and as something else in the other would be
  # a difference a caller who switched modes would meet with nothing to warn
  # them.
  res <- ipw_result(
    binary_estimates(),
    vcov = binary_vcov(),
    outcome_vcov = corrected_outcome_vcov(),
    effects = "conditional"
  )

  none <- confint(res, parm = 0)

  expect_true(is.matrix(none))
  expect_identical(dim(none), c(0L, 2L))
  expect_identical(colnames(none), c("2.5 %", "97.5 %"))
  # The surface it emptied is the coefficient one: two rows rather than the
  # three effects the other reading reports, so the assertions above are about
  # this reading rather than about a result that answered in the other.
  expect_identical(dim(confint(res)), c(2L, 2L))
})

test_that("confint() refuses a conditional parm the outcome model lacks", {
  # `rd` is an effect the result reports in its other reading, which is exactly
  # the mistake worth catching: the labels a caller has in hand belong to the
  # surface they were reading, and asking the other one for them names nothing.
  # The condition is the one the marginal reading raises for an unmatched label,
  # so a caller handles one class in either mode.
  res <- ipw_result(
    binary_estimates(),
    vcov = binary_vcov(),
    outcome_vcov = corrected_outcome_vcov(),
    effects = "conditional"
  )

  expect_error(
    confint(res, parm = "rd"),
    class = "causalgenerics_invalid_argument_parm"
  )
  expect_error(
    confint(res, parm = "rd"),
    class = "causalgenerics_invalid_argument"
  )
  expect_error(
    confint(res, parm = c("z", "rd")),
    class = "causalgenerics_invalid_argument_parm"
  )
  expect_error(
    confint(res, parm = 3),
    class = "causalgenerics_invalid_argument_parm"
  )
  # A position is matched against the surface the call names, so the reading
  # asked for at the call site is refused the same way as the one the result
  # records. Two coefficients are all this model has, and `99` indexes none of
  # them either way.
  expect_error(
    confint(res, parm = 99, effects = "conditional"),
    class = "causalgenerics_invalid_argument_parm"
  )

  expect_snapshot(error = TRUE, confint(res, parm = "rd"))
})

test_that("confint() refuses the conditional mode without a corrected block", {
  # The limits are built from the covariance, so an interval in this mode needs
  # the block `vcov()` needs. Widths taken from the model's own covariance would
  # be the naive ones, and an interval is where a reader would be least likely
  # to notice.
  res <- ipw_result(
    binary_estimates(),
    vcov = binary_vcov(),
    effects = "conditional"
  )

  expect_error(confint(res), class = "causalgenerics_no_conditional_vcov")
  expect_error(confint(res), class = "causalgenerics_no_vcov")
  expect_error(
    confint(res, parm = 1, level = 0.9),
    class = "causalgenerics_no_conditional_vcov"
  )
  expect_no_error(stats::vcov(res$outcome_mod))

  expect_snapshot(error = TRUE, confint(res))
})

test_that("confint() refuses a conditional block the outcome model lost", {
  # The limits are built from the covariance, so this reading needs the block
  # `vcov()` needs, whatever the reason it is missing. A wrapper carrying
  # nothing is where an interval goes wrong most quietly: the widths taken from
  # it are an empty vector, and the limits come back as `NA` at every level with
  # nothing said.
  res <- stripped_wrapper_result()

  expect_error(confint(res), class = "causalgenerics_no_vcov_ipw_model")
  expect_error(confint(res), class = "causalgenerics_no_vcov")
  expect_error(
    confint(res, parm = 1, level = 0.9),
    class = "causalgenerics_no_vcov_ipw_model"
  )

  # The covariance is asked for before `parm` is matched, the way the reading's
  # own refusal documents. A label the surface does not have is beside the point
  # when no interval can be built for any of them.
  expect_error(
    confint(res, parm = "nope", effects = "conditional"),
    class = "causalgenerics_no_vcov_ipw_model"
  )

  # Not the conditional error, for the reason `vcov()` does not raise it here:
  # the model was wrapped, and wrapping it again is no answer.
  cnd <- tryCatch(confint(res), error = identity)
  expect_false(inherits(cnd, "causalgenerics_no_conditional_vcov"))
})

# ---- naming a reading for one call -------------------------------------------

test_that("the effects argument overrides the mode a result records", {
  # The stored mode is the default and the argument is the override, in both
  # directions. A caller who holds a marginal result and wants the coefficient
  # surface for one line does not have to flip the object to get it.
  corrected <- corrected_outcome_vcov()
  marginal <- ipw_result(
    binary_estimates(),
    vcov = binary_vcov(),
    outcome_vcov = corrected
  )
  conditional <- as_conditional(marginal)
  effects <- setNames(binary_estimates()$estimate, binary_labels())

  expect_identical(coef(marginal, effects = "conditional"), coef(conditional))
  expect_identical(vcov(marginal, effects = "conditional"), corrected)
  expect_identical(
    confint(marginal, effects = "conditional"),
    confint(conditional)
  )

  expect_identical(coef(conditional, effects = "marginal"), effects)
  expect_identical(vcov(conditional, effects = "marginal"), binary_vcov())
  expect_identical(
    confint(conditional, effects = "marginal"),
    confint(marginal)
  )
})

test_that("the effects argument reaches parm and level in the same call", {
  # The three arguments are independent, and a caller who names all of them
  # names the surface, the rows of it, and the width at once. The mode reaches
  # `confint()` through an argument here rather than through the field, so an
  # implementation that resolved it after selecting the rows would select from
  # the wrong surface and the level would be applied to the wrong estimates.
  corrected <- corrected_outcome_vcov()
  marginal <- ipw_result(
    binary_estimates(),
    vcov = binary_vcov(),
    outcome_vcov = corrected
  )
  half_width <- qnorm(1 - (1 - 0.9) / 2) * sqrt(diag(corrected))
  estimate <- coef(marginal$outcome_mod)

  ci <- confint(marginal, parm = "z", level = 0.9, effects = "conditional")

  expect_identical(dim(ci), c(1L, 2L))
  expect_identical(dimnames(ci), list("z", c("5 %", "95 %")))
  expect_identical(ci[1, 1], unname(estimate - half_width)[2])
  expect_identical(ci[1, 2], unname(estimate + half_width)[2])
})

test_that("a NULL effects argument declines to override", {
  # `NULL` is the default, and it means the caller named no reading rather than
  # naming a third one. Passing it through from a wrapper that computes an
  # argument therefore reads the result the way omitting it does.
  res <- ipw_result(
    binary_estimates(),
    vcov = binary_vcov(),
    outcome_vcov = corrected_outcome_vcov(),
    effects = "conditional"
  )

  expect_identical(coef(res, effects = NULL), coef(res))
  expect_identical(coef(res, effects = NULL), coef(res$outcome_mod))
  expect_identical(vcov(res, effects = NULL), vcov(res))
  expect_identical(confint(res, effects = NULL), confint(res))
})

test_that("the marginal reading is the one that was there before the mode", {
  # Everything the sections above assert of a result holds of a marginal result
  # and of a conditional one asked for the marginal reading. The mode adds a
  # surface; it does not move the one that was already reported.
  estimates <- binary_estimates()
  covariance <- binary_vcov()
  marginal <- ipw_result(
    estimates,
    vcov = covariance,
    outcome_vcov = corrected_outcome_vcov()
  )
  conditional <- as_conditional(marginal)

  for (res in list(marginal, conditional)) {
    expect_identical(
      coef(res, effects = "marginal"),
      setNames(estimates$estimate, binary_labels())
    )
    expect_identical(vcov(res, effects = "marginal"), covariance)

    ci <- confint(res, effects = "marginal")

    expect_identical(
      dimnames(ci),
      list(binary_labels(), c("2.5 %", "97.5 %"))
    )
    # The stored limits at the stored level, which is the fast path the marginal
    # reading takes and the conditional one has no use for.
    expect_identical(ci[, 1], setNames(estimates$ci.lower, binary_labels()))
    expect_identical(ci[, 2], setNames(estimates$ci.upper, binary_labels()))
    expect_identical(
      rownames(confint(res, parm = "log(rr)", effects = "marginal")),
      "log(rr)"
    )
  }

  # Without this the loop would hold of accessors that ignored the argument and
  # the mode alike, which is what they did before this contract existed.
  expect_false(identical(coef(conditional), coef(marginal)))
  expect_false(identical(vcov(conditional), vcov(marginal)))
})

test_that("the accessors refuse an effects argument that is not a reading", {
  # There are two readings and no third, so anything else is a misspelling or a
  # wrong argument, and answering it with either surface would give the caller
  # one they did not ask for. The condition is the one the constructor raises
  # for the same value, so a caller handles a single class wherever it came
  # from.
  res <- ipw_result(binary_estimates(), vcov = binary_vcov())

  expect_error(
    coef(res, effects = "banana"),
    class = "causalgenerics_invalid_argument_effects"
  )
  expect_error(
    coef(res, effects = "banana"),
    class = "causalgenerics_invalid_argument"
  )
  # A character vector of length one that names nothing, which `==` alone would
  # answer with `NA` rather than with a refusal.
  expect_error(
    coef(res, effects = NA_character_),
    class = "causalgenerics_invalid_argument_effects"
  )
  expect_error(
    vcov(res, effects = 1),
    class = "causalgenerics_invalid_argument_effects"
  )
  expect_error(
    confint(res, effects = c("marginal", "conditional")),
    class = "causalgenerics_invalid_argument_effects"
  )

  expect_snapshot(error = TRUE, coef(res, effects = "banana"))
  expect_snapshot(error = TRUE, vcov(res, effects = 1))
  expect_snapshot(
    error = TRUE,
    confint(res, effects = c("marginal", "conditional"))
  )
})

test_that("the accessors refuse a stored mode that is not a reading", {
  # The field is set through `as_marginal()` and `as_conditional()`, which can
  # only put one of the two readings there, and the constructor checks the value
  # it is handed. A field that says something else was assigned to directly,
  # which is outside the contract, and there is no reading to fall back on:
  # reporting the marginal surface would answer for a result that named neither
  # and leave the corruption in place for the next caller to inherit.
  res <- ipw_result(binary_estimates(), vcov = binary_vcov())
  res$effects <- "banana"

  expect_error(coef(res), class = "causalgenerics_invalid_argument_effects")
  expect_error(coef(res), class = "causalgenerics_invalid_argument")
  expect_error(vcov(res), class = "causalgenerics_invalid_argument_effects")
  expect_error(confint(res), class = "causalgenerics_invalid_argument_effects")
  # Naming a reading is still answerable. The argument says which surface to
  # report, so the stored field is not consulted at all.
  expect_identical(vcov(res, effects = "marginal"), binary_vcov())

  expect_snapshot(error = TRUE, coef(res))
})

# ---- what the mode does not reach --------------------------------------------

test_that("nobs(), df.residual(), and weights() do not read the mode", {
  # These three describe the fit rather than a surface of it. The number of
  # observations the outcome model saw, the residual degrees of freedom of the
  # variance object, and the weights it was fitted with are the same facts
  # whichever reading the result presents, so a flip cannot move them.
  marginal <- ipw_result(
    binary_estimates(),
    vcov = binary_vcov(),
    fit = lm(y ~ x, data = ipw_data()),
    outcome_vcov = corrected_outcome_vcov()
  )
  conditional <- as_conditional(marginal)

  expect_identical(nobs(conditional), nobs(marginal))
  expect_identical(nobs(conditional), 20L)
  expect_identical(df.residual(conditional), df.residual(marginal))
  expect_identical(df.residual(conditional), 18L)
  expect_identical(weights(conditional), weights(marginal))
  expect_identical(weights(conditional), outcome_weights())

  # The flip is one the other three accessors do act on, so the assertions above
  # are about these three rather than about a result nothing reads the mode of.
  expect_false(identical(coef(conditional), coef(marginal)))
  expect_false(identical(vcov(conditional), vcov(marginal)))
})

test_that("a result with no mode reads as marginal", {
  # A result stored before the field existed carries six fields. Every accessor
  # reports the reading it reported then, and a caller who names the other one
  # gets it, so an older result needs no migration to be read either way.
  estimates <- binary_estimates()
  covariance <- binary_vcov()
  legacy <- legacy_result(estimates, vcov = covariance)

  expect_length(legacy, 6L)
  expect_null(legacy$effects)

  expect_identical(
    coef(legacy),
    setNames(estimates$estimate, binary_labels())
  )
  expect_identical(vcov(legacy), covariance)
  expect_identical(rownames(confint(legacy)), binary_labels())
  expect_identical(
    confint(legacy)[, 1],
    setNames(estimates$ci.lower, binary_labels())
  )

  expect_identical(coef(legacy, effects = "marginal"), coef(legacy))
  expect_identical(
    coef(legacy, effects = "conditional"),
    coef(legacy$outcome_mod)
  )
})

test_that("model.frame() does not read the mode", {
  # The frame is the outcome model's, and which surface a result presents says
  # nothing about the data that model was fitted on. A result stored before the
  # field existed answers the same way, so an older one needs no migration to be
  # asked for its frame.
  res <- ipw_result(binary_estimates())
  legacy <- legacy_result()

  expect_length(legacy, 6L)
  expect_identical(model.frame(as_conditional(res)), model.frame(res))
  expect_identical(model.frame(legacy), model.frame(res))
  expect_identical(model.frame(res), outcome_frame_columns())
})

test_that("estimand() does not read the mode", {
  # The estimand is what the weights targeted, and both readings of a result are
  # readings of estimates computed under those weights. A flip cannot move it,
  # and a result recording no mode reports it unchanged.
  res <- ipw_result(binary_estimates())
  legacy <- legacy_result()

  expect_length(legacy, 6L)
  expect_identical(estimand(as_conditional(res)), estimand(res))
  expect_identical(estimand(legacy), estimand(res))
  expect_identical(estimand(res), "ate")
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

test_that("model.frame() is registered against ipw", {
  # `model.frame()` lives in stats as the six accessors above do, so its entry
  # belongs to stats' table rather than to this package's. The claim is the one
  # the test above makes: every assertion written for this method reaches it
  # from the test frame, and a downstream package calling `model.frame()` on a
  # result from its own namespace has only the table.
  table <- s3_methods_table("model.frame")

  expect_false(is.null(table))
  expect_true(exists("model.frame.ipw", envir = table, inherits = FALSE))
})

test_that("estimand() is registered against ipw", {
  # `estimand()` is this package's own generic, so its entry belongs to this
  # package's table rather than to stats'. `s3_methods_table()` reads the
  # generic's own environment and finds either one.
  table <- s3_methods_table("estimand")

  expect_false(is.null(table))
  expect_true(exists("estimand.ipw", envir = table, inherits = FALSE))
})
