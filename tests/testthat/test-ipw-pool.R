# `pool_ipw()` combines the `ipw` results fitted to a set of multiply imputed
# datasets into one result, by Rubin's rules, and returns it under the class
# `ipw_pooled`. It lives here for the reason the rest of the result layer does:
# every package that supplies an `ipw()` method returns the same class, so the
# pooling of those results is written once against the contract rather than once
# per package.
#
# The function reads a plain list of results or a `mira`, which is what
# `mice::with()` returns. A `mira` is recognized by its class and read through
# its `analyses` element, and nothing else about mice is assumed: this package
# holds no dependency on it, and the fixtures below build the class by hand.
#
# The estimates frames are written out literally, in the shape the `ipw()`
# return contract documents, rather than produced by fitting a model. The model
# slots hold real fits, since the complete-data degrees of freedom, the
# observation count, and the outcome link are all read off them.
#
# The pooled arithmetic is asserted twice over. The `rd` label of the binary
# fixture carries numbers chosen so that every pooled quantity is an exact
# fraction that can be written down, and the comments on those tests work each
# one out. Every other label is checked against `rubin_rules()` below, which is
# transcribed from Rubin's rules and the Barnard-Rubin small-sample adjustment
# as `mice::pool()` implements them.

# ---- fixtures ----------------------------------------------------------------

# Fixed rather than simulated data, so the models never move between runs. The
# row count is an argument because three steps of the complete-data degrees of
# freedom chain are about fits that differ in how much data they saw.
pool_data <- function(rows = 20) {
  dat <- data.frame(
    x = rep(c(-1.5, -0.5, 0.5, 1.5), each = 5),
    z = rep(c(0, 1), 10),
    y = rep(c(0, 1, 1, 0, 1), 4)
  )
  dat[seq_len(rows), , drop = FALSE]
}

# One imputation of that dataset: a single cell filled in differently, which is
# what separates two imputations of the same incomplete data. The coefficients
# move with it, so the between-imputation variance the conditional reading
# reports is not zero.
pool_imputed_data <- function(i) {
  dat <- pool_data()
  dat$y[i] <- 1 - dat$y[i]
  dat
}

# The outcome model each result is built around. Two coefficients on twenty
# rows, so it reports eighteen residual degrees of freedom and twenty
# observations, and its link is `logit`.
pool_outcome_model <- function(rows = 20) {
  glm(y ~ z, family = quasibinomial(), data = pool_data(rows))
}

pool_wt_model <- function(rows = 20) {
  glm(z ~ x, family = binomial(), data = pool_data(rows))
}

# One imputation's binary estimates frame. The per-fit numbers are arguments
# rather than part of the body, so the arithmetic the pooling tests assert
# against is written at the call site where it can be read.
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

# The three imputations' estimates, one row per fit and one column per effect.
# The `rd` column is what makes the literal test below possible: the three
# estimates are 0.20, 0.30, and 0.40 against standard errors of 0.10, 0.20, and
# 0.30, and every pooled quantity those give is an exact fraction.
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

# The effect labels the binary and categorical frames carry, written out rather
# than pasted together. The label rule is part of what these tests assert, so
# restating it with the same `paste()` the implementation uses would assert
# nothing.
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

# One imputation's categorical estimates frame: the effect measures repeat
# across contrasts, so a column sits immediately after `effect` naming each one.
pool_categorical_estimates <- function(
  estimate,
  std.err,
  contrast = c("b vs a", "c vs a")
) {
  half_width <- stats::qnorm(0.975) * std.err
  data.frame(
    effect = rep(c("rd", "log(rr)", "log(or)"), times = 2),
    contrast = rep(contrast, each = 3),
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

# A covariance of the reported effects, in the shape a fitting package attaches.
# The diagonal is the squared standard errors the frame reports, so the pooled
# block's diagonal has to agree with the pooled standard errors, and the
# off-diagonal entries fall off with the distance between rows. Non-zero
# off-diagonals are the point: a pooling that combined the diagonals alone would
# agree with this one there and differ everywhere else.
pool_effects_vcov <- function(std_err, labels, rho = 0.9) {
  index <- seq_along(std_err)
  covariance <- rho^abs(outer(index, index, "-")) * outer(std_err, std_err)
  dimnames(covariance) <- list(labels, labels)
  covariance
}

# Assemble one imputation's result the way a fitting package would, so that each
# test names only the part it is about.
pool_fit <- function(
  estimates,
  vcov = NULL,
  fit = NULL,
  estimand = "ate",
  se_method = "linearization",
  effects = "marginal",
  readings = c("marginal", "conditional"),
  outcome_mod = pool_outcome_model(),
  wt_mod = pool_wt_model()
) {
  if (!is.null(vcov)) {
    attr(estimates, "ipw_vcov") <- vcov
  }
  new_ipw(
    estimand = estimand,
    wt_mod = wt_mod,
    outcome_mod = outcome_mod,
    estimates = estimates,
    se_method = se_method,
    fit = fit,
    effects = effects,
    readings = readings
  )
}

# The three binary results a pooling reads. `vcov = TRUE` attaches the
# covariance of the effects to each one, which is what the pooled block is built
# from; the results a fitting package with no covariance to report produces
# carry none.
pool_binary_fits <- function(
  estimate = pool_binary_estimate(),
  std_err = pool_binary_std_err(),
  vcov = FALSE,
  ...
) {
  lapply(seq_len(nrow(estimate)), function(i) {
    covariance <- if (vcov) {
      pool_effects_vcov(std_err[i, ], pool_binary_labels())
    } else {
      NULL
    }
    pool_fit(
      pool_binary_estimates(estimate[i, ], std_err[i, ]),
      vcov = covariance,
      ...
    )
  })
}

pool_categorical_fits <- function(
  estimate = pool_categorical_estimate(),
  std_err = pool_categorical_std_err(),
  vcov = FALSE,
  ...
) {
  lapply(seq_len(nrow(estimate)), function(i) {
    covariance <- if (vcov) {
      pool_effects_vcov(std_err[i, ], pool_categorical_labels())
    } else {
      NULL
    }
    pool_fit(
      pool_categorical_estimates(estimate[i, ], std_err[i, ]),
      vcov = covariance,
      ...
    )
  })
}

# One imputation's estimates frame for a result reported once per level of a
# grouping variable. The effect measures repeat across the subgroups, so a
# column sits immediately after `effect` saying which subgroup each row
# describes, written as the `"var = value"` string the contract spells a group
# with. The numbers are the categorical fixture's, since what differs between
# the two shapes is which column completes a row's identity rather than what the
# other columns hold.
pool_group_estimates <- function(
  estimate,
  std.err,
  group = c("sex = 0", "sex = 1")
) {
  half_width <- stats::qnorm(0.975) * std.err
  data.frame(
    effect = rep(c("rd", "log(rr)", "log(or)"), times = 2),
    group = rep(group, each = 3),
    estimate = estimate,
    std.err = std.err,
    z = estimate / std.err,
    ci.lower = estimate - half_width,
    ci.upper = estimate + half_width,
    conf.level = 0.95,
    p.value = 2 * stats::pnorm(-abs(estimate / std.err))
  )
}

pool_group_labels <- function() {
  c(
    "rd sex = 0",
    "log(rr) sex = 0",
    "log(or) sex = 0",
    "rd sex = 1",
    "log(rr) sex = 1",
    "log(or) sex = 1"
  )
}

pool_group_fits <- function(
  estimate = pool_categorical_estimate(),
  std_err = pool_categorical_std_err(),
  vcov = FALSE,
  group = c("sex = 0", "sex = 1"),
  ...
) {
  lapply(seq_len(nrow(estimate)), function(i) {
    covariance <- if (vcov) {
      pool_effects_vcov(std_err[i, ], pool_group_labels())
    } else {
      NULL
    }
    pool_fit(
      pool_group_estimates(estimate[i, ], std_err[i, ], group = group),
      vcov = covariance,
      ...
    )
  })
}

# One imputation's estimates frame crossing the contrast with the group: two
# effect measures, two contrasts, and two subgroups. Eight rows, and no one of
# the three key columns names one of them, nor any two of them together, which
# is what makes this the fixture that says a row's identity is all three.
pool_contrast_group_estimates <- function(
  estimate,
  std.err,
  group = c("sex = 0", "sex = 1")
) {
  half_width <- stats::qnorm(0.975) * std.err
  data.frame(
    effect = rep(c("rd", "log(rr)"), times = 4),
    contrast = rep(rep(c("b vs a", "c vs a"), each = 2), times = 2),
    group = rep(group, each = 4),
    estimate = estimate,
    std.err = std.err,
    z = estimate / std.err,
    ci.lower = estimate - half_width,
    ci.upper = estimate + half_width,
    conf.level = 0.95,
    p.value = 2 * stats::pnorm(-abs(estimate / std.err))
  )
}

pool_contrast_group_estimate <- function() {
  rbind(
    c(0.08, 0.17, 0.16, 0.31, 0.06, 0.13, 0.12, 0.25),
    c(0.10, 0.20, 0.20, 0.35, 0.08, 0.16, 0.15, 0.29),
    c(0.12, 0.23, 0.24, 0.39, 0.10, 0.19, 0.18, 0.33)
  )
}

pool_contrast_group_std_err <- function() {
  rbind(
    c(0.05, 0.10, 0.04, 0.09, 0.06, 0.11, 0.05, 0.10),
    c(0.06, 0.11, 0.05, 0.10, 0.07, 0.12, 0.06, 0.11),
    c(0.07, 0.12, 0.06, 0.11, 0.08, 0.13, 0.07, 0.12)
  )
}

pool_contrast_group_labels <- function() {
  c(
    "rd b vs a sex = 0",
    "log(rr) b vs a sex = 0",
    "rd c vs a sex = 0",
    "log(rr) c vs a sex = 0",
    "rd b vs a sex = 1",
    "log(rr) b vs a sex = 1",
    "rd c vs a sex = 1",
    "log(rr) c vs a sex = 1"
  )
}

pool_contrast_group_fits <- function(
  estimate = pool_contrast_group_estimate(),
  std_err = pool_contrast_group_std_err(),
  vcov = FALSE,
  ...
) {
  lapply(seq_len(nrow(estimate)), function(i) {
    covariance <- if (vcov) {
      pool_effects_vcov(std_err[i, ], pool_contrast_group_labels())
    } else {
      NULL
    }
    pool_fit(
      pool_contrast_group_estimates(estimate[i, ], std_err[i, ]),
      vcov = covariance,
      ...
    )
  })
}

# The same three results with one field replaced in each. Each step of the
# complete-data degrees of freedom chain reads a different field, so each of
# those tests needs fits that differ in one field and agree in the rest.
pool_fits_varying <- function(field, values, ...) {
  fits <- pool_binary_fits(...)
  Map(
    function(fit, value) {
      fit[[field]] <- value
      fit
    },
    fits,
    values
  )
}

# Three results that disagree about the mode they present. With no mode named
# the stored field decides which surface is pooled and this set names none;
# naming one at the call site answers the question. The two tests that read this
# fixture are the two halves of that one contract, so they share it rather than
# building the same disagreement twice.
pool_mixed_mode_fits <- function() {
  pool_fits_varying("effects", list("marginal", "conditional", "marginal"))
}

# The outcome models of three imputations, each wrapped with the corrected
# covariance the joint estimation of the weights and the outcome implies. The
# conditional reading pools that surface: the coefficients are the estimates and
# the square roots of the block's diagonal are their standard errors. The
# scaling differs by imputation so that the within-imputation variance is not
# constant either.
pool_conditional_models <- function() {
  lapply(1:3, function(i) {
    mod <- glm(y ~ z, family = quasibinomial(), data = pool_imputed_data(i))
    new_ipw_model(mod, stats::vcov(mod) * (1 + i / 10))
  })
}

pool_conditional_fits <- function(
  effects = "conditional",
  readings = c("marginal", "conditional")
) {
  estimate <- pool_binary_estimate()
  std_err <- pool_binary_std_err()
  Map(
    function(i, mod) {
      pool_fit(
        pool_binary_estimates(estimate[i, ], std_err[i, ]),
        outcome_mod = mod,
        effects = effects,
        readings = readings
      )
    },
    1:3,
    pool_conditional_models()
  )
}

# Three results whose two surfaces are both poolable: the outcome models carry
# the corrected block the conditional reading is pooled from, and the estimates
# frames carry the covariance of the effects the marginal reading pools. A
# fixture holding only one of the two would leave the other surface's
# covariance out of whatever a two-sided test asserted.
#
# The stored level is named one per imputation. The marginal surface is the
# only one that records a level, so it is the surface whose level the two
# readings resolve differently, and a set that disagrees about it is one of the
# ways the reading that was not asked for cannot be pooled.
pool_both_surface_fits <- function(
  effects = "marginal",
  conf.level = c(0.95, 0.95, 0.95)
) {
  estimate <- pool_binary_estimate()
  std_err <- pool_binary_std_err()
  Map(
    function(i, mod) {
      pool_fit(
        pool_binary_estimates(estimate[i, ], std_err[i, ], conf.level[[i]]),
        vcov = pool_effects_vcov(std_err[i, ], pool_binary_labels()),
        outcome_mod = mod,
        effects = effects
      )
    },
    1:3,
    pool_conditional_models()
  )
}

# The outcome models of three imputations of a three-level exposure, wrapped
# the way the binary ones above are. Their coefficients are named for the
# levels rather than for the contrasts the marginal surface reports, and there
# are three of them against six reported effects, so the two surfaces of a
# categorical result are separable by their rows alone.
pool_categorical_models <- function() {
  lapply(1:3, function(i) {
    dat <- pool_imputed_data(i)
    dat$g <- factor(rep(c("a", "b", "c", "a"), each = 5))
    mod <- glm(y ~ g, family = quasibinomial(), data = dat)
    new_ipw_model(mod, stats::vcov(mod) * (1 + i / 10))
  })
}

pool_categorical_both_surface_fits <- function(effects = "marginal") {
  estimate <- pool_categorical_estimate()
  std_err <- pool_categorical_std_err()
  Map(
    function(i, mod) {
      pool_fit(
        pool_categorical_estimates(estimate[i, ], std_err[i, ]),
        outcome_mod = mod,
        effects = effects
      )
    },
    1:3,
    pool_categorical_models()
  )
}

# A carrier for an outcome model that reports no residual degrees of freedom: a
# bare list whose class has no `df.residual()` method, so the default one looks
# for an element that is not there and finds nothing. What it does answer for is
# registered by the tests that use it, since a registration has to be undone
# when the test that made it exits.
#
# The coefficients sit in the list rather than in the method, so three carriers
# built here report three different coefficient surfaces through one
# registration.
pool_bare_model <- function(coefficients = NULL) {
  structure(list(coefficients = coefficients), class = "cg_pool_bare")
}

# Three outcome models whose conditional surface is poolable and which report no
# residual degrees of freedom. `new_ipw_model()` registers nothing but `vcov()`,
# so wrapping the carrier above hands the conditional reading the corrected
# block it pools while `df.residual()` still finds nothing on the model. A
# fitted model holds those two properties apart, since one that reports
# coefficients reports a residual count as well, which is why the carrier is
# built rather than fitted.
#
# The coefficients and the block move by imputation, for the reason the fitted
# models above do: the between-imputation variance of the surface is then not
# zero.
pool_dfless_models <- function() {
  lapply(1:3, function(i) {
    coefficients <- c("(Intercept)" = -0.4 + i / 10, z = 0.5 + i / 20)
    new_ipw_model(
      pool_bare_model(coefficients),
      pool_effects_vcov(c(0.2, 0.3) * (1 + i / 10), names(coefficients))
    )
  })
}

# Three results built around those models, recording the marginal reading, so
# the conditional one is the alternate. Both readings of the set can be pooled
# and nothing in it reports a complete-data count, which is the set the
# large-sample fallback has to be exercised against: where the other reading
# cannot be pooled, the count is never wanted on that side, whether or not it is
# resolved once for the call.
pool_dfless_fits <- function() {
  estimate <- pool_binary_estimate()
  std_err <- pool_binary_std_err()
  Map(
    function(i, mod) {
      pool_fit(
        pool_binary_estimates(estimate[i, ], std_err[i, ]),
        outcome_mod = mod
      )
    },
    1:3,
    pool_dfless_models()
  )
}

# The three binary results as a fitting package stored them before the mode
# existed: six fields, with neither `effects` nor `readings` after them. The
# absent mode reads as marginal, which is the reading every method produced
# then, and the outcome models carry no corrected block either, which is the
# pair of properties such a set has.
pool_legacy_fits <- function() {
  lapply(pool_binary_fits(), function(fit) {
    fields <- unclass(fit)
    structure(
      fields[!names(fields) %in% c("effects", "readings")],
      class = "ipw"
    )
  })
}

# The same three results with the middle one supporting the conditional reading
# alone, which is the set where a single result decides what the pooling can
# report. Every result records the conditional mode, so the set agrees about
# the reading it presents, and replacing the field on one of them leaves a
# state the constructor would build: the mode that result records is one of the
# readings it names.
pool_mixed_reading_fits <- function() {
  fits <- pool_both_surface_fits(effects = "conditional")
  fits[[2]]$readings <- "conditional"
  fits
}

# Three results whose marginal surfaces disagree about the effects they report.
# A set that disagrees that way cannot be pooled on the marginal reading, and a
# set that also supports the conditional reading alone cannot be pooled on it
# for a second reason: which of the two a caller is told about says which
# question the pooling asked first.
pool_mislabeled_fits <- function(readings = c("marginal", "conditional")) {
  fits <- pool_conditional_fits(readings = readings)
  fits[[2]]$estimates$effect <- c("rd", "log(rr)", "log(hr)")
  fits
}

# Rubin's rules and the Barnard-Rubin small-sample adjustment for one effect,
# transcribed from the formulas `mice::pool()` implements. The tests recompute
# with this from the fixture inputs rather than reading the pooled frame back
# into itself, which would assert nothing.
#
# `t` here is the total variance, whose square root is the pooled standard
# error. The `t` column of the pooled estimates frame is the test statistic,
# which is a different quantity under the same letter.
rubin_rules <- function(estimate, std.err, dfcom) {
  m <- length(estimate)
  qbar <- mean(estimate)
  ubar <- mean(std.err^2)
  b <- stats::var(estimate)
  t <- ubar + (1 + 1 / m) * b

  lambda <- (1 + 1 / m) * b / t
  dfold <- (m - 1) / lambda^2
  df <- if (is.infinite(dfcom)) {
    dfold
  } else {
    tmp <- (1 - lambda) * (1 + dfcom) * dfcom
    (m - 1) * tmp / ((dfcom + 3) * (m - 1) + lambda^2 * tmp)
  }

  riv <- (1 + 1 / m) * b / ubar
  list(
    estimate = qbar,
    std.err = sqrt(t),
    df = df,
    ubar = ubar,
    b = b,
    riv = riv,
    lambda = lambda,
    fmi = (riv + 2 / (df + 3)) / (riv + 1)
  )
}

# The pooled quantities of every column of a fixture, as one list per effect in
# the fixture's own column order.
rubin_rules_by_column <- function(estimate, std_err, dfcom) {
  lapply(
    seq_len(ncol(estimate)),
    function(j) rubin_rules(estimate[, j], std_err[, j], dfcom)
  )
}

# One field of those, gathered into the vector the pooled frame holds.
rubin_column <- function(pooled, field) {
  vapply(pooled, function(x) x[[field]], numeric(1))
}

# Every warning one call raises, with the call's value beside them. Counting
# them takes a calling handler rather than `expect_warning()`, which muffles the
# first warning it matches and leaves a test unable to say whether a second one
# followed.
pool_warnings <- function(expr) {
  raised <- character()
  value <- withCallingHandlers(
    expr,
    warning = function(w) {
      raised <<- c(raised, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  list(value = value, warnings = raised)
}

# ---- what `fits` may be ------------------------------------------------------

test_that("pool_ipw() refuses a `fits` argument that is not a list", {
  # The argument is the set of results one analysis produced across the
  # imputations, and there is no reading of a number or a string as one.
  expect_error(pool_ipw(1:3), class = "causalgenerics_invalid_argument_fits")
  expect_error(pool_ipw(1:3), class = "causalgenerics_invalid_argument")
  expect_error(pool_ipw("a"), class = "causalgenerics_invalid_argument_fits")
  expect_error(pool_ipw(NULL), class = "causalgenerics_invalid_argument_fits")

  expect_snapshot(error = TRUE, pool_ipw(1:3))
})

test_that("pool_ipw() refuses elements that are not results", {
  # A result is itself a list, so a caller who passes one rather than the set of
  # them reaches the same guard from the other side: the fields of an `ipw`
  # object are not `ipw` objects. The positions are named, since a set of twenty
  # imputations is not something a caller reads through to find the odd one.
  fits <- pool_binary_fits()

  expect_error(
    pool_ipw(fits[[1]]),
    class = "causalgenerics_invalid_argument_fits"
  )

  one_bad <- list(fits[[1]], "nope", fits[[3]])
  expect_error(
    pool_ipw(one_bad),
    class = "causalgenerics_invalid_argument_fits"
  )
  expect_error(pool_ipw(one_bad), regexp = "\\b2\\b")

  two_bad <- list(fits[[1]], "nope", 3)
  cnd <- tryCatch(pool_ipw(two_bad), error = identity)
  expect_match(conditionMessage(cnd), "\\b2\\b")
  expect_match(conditionMessage(cnd), "\\b3\\b")

  expect_snapshot(error = TRUE, pool_ipw(fits[[1]]))
  expect_snapshot(error = TRUE, pool_ipw(one_bad))
  expect_snapshot(error = TRUE, pool_ipw(two_bad))
})

test_that("pool_ipw() refuses fewer than two results", {
  # Rubin's rules have no between-imputation variance to estimate from one
  # analysis, and none at all from zero. `mice::pool()` warns and hands the
  # single fit back; that returns an object of a different class than the call
  # asked for, and a caller who did not read the warning carries on with it.
  fits <- pool_binary_fits()

  expect_error(pool_ipw(list()), class = "causalgenerics_invalid_argument_fits")
  expect_error(
    pool_ipw(fits[1]),
    class = "causalgenerics_invalid_argument_fits"
  )
  expect_error(
    pool_ipw(structure(list(analyses = fits[1]), class = "mira")),
    class = "causalgenerics_invalid_argument_fits"
  )

  expect_snapshot(error = TRUE, pool_ipw(list()))
  expect_snapshot(error = TRUE, pool_ipw(fits[1]))
})

test_that("pool_ipw() reads a mira the way it reads a list", {
  # `mice::with()` returns a `mira`, which is a list of fits under an
  # `analyses` element. The class is all that is read, so this package needs no
  # dependency on mice and a hand-built object of that class is answered the
  # same way. The two results are the same object rather than merely the same
  # numbers, which is what says the wrapper is unwrapped and nothing else.
  fits <- pool_binary_fits(vcov = TRUE)
  mira <- structure(list(analyses = fits), class = "mira")

  expect_identical(pool_ipw(mira), pool_ipw(fits))
  expect_identical(
    pool_ipw(mira, conf_level = 0.8),
    pool_ipw(fits, conf_level = 0.8)
  )
})

# ---- the results have to agree -----------------------------------------------

test_that("pool_ipw() refuses results targeting different estimands", {
  # The pooled estimate is an average of the per-imputation ones, and averaging
  # an ATE with an ATT gives a number that estimates neither. The differing
  # values travel on the condition, so a caller reports them without reading
  # them back out of the sentence.
  fits <- pool_fits_varying("estimand", list("ate", "att", "ate"))

  expect_error(pool_ipw(fits), class = "causalgenerics_pool_mismatch_estimand")
  expect_error(pool_ipw(fits), class = "causalgenerics_pool_mismatch")

  cnd <- tryCatch(pool_ipw(fits), error = identity)
  expect_identical(cnd$field, "estimand")
  expect_setequal(unlist(cnd$values), c("ate", "att"))

  expect_snapshot(error = TRUE, pool_ipw(fits))
})

test_that("pool_ipw() refuses results whose standard errors differ in kind", {
  # The within-imputation variance is the average of the per-imputation ones, so
  # the standard errors have to be the same quantity across the set. A
  # linearization standard error and an M-estimation one answer different
  # questions about the same estimate.
  fits <- pool_fits_varying(
    "se_method",
    list("linearization", "mestimation", "linearization")
  )

  expect_error(pool_ipw(fits), class = "causalgenerics_pool_mismatch_se_method")
  expect_error(pool_ipw(fits), class = "causalgenerics_pool_mismatch")

  cnd <- tryCatch(pool_ipw(fits), error = identity)
  expect_identical(cnd$field, "se_method")
  expect_setequal(unlist(cnd$values), c("linearization", "mestimation"))

  expect_snapshot(error = TRUE, pool_ipw(fits))
})

test_that("pool_ipw() refuses results recording different presentation modes", {
  # With no mode named at the call site the stored one decides which surface is
  # pooled, and a set that disagrees about it names no surface at all.
  fits <- pool_mixed_mode_fits()

  expect_error(pool_ipw(fits), class = "causalgenerics_pool_mismatch_effects")
  expect_error(pool_ipw(fits), class = "causalgenerics_pool_mismatch")

  cnd <- tryCatch(pool_ipw(fits), error = identity)
  expect_identical(cnd$field, "effects")
  expect_setequal(unlist(cnd$values), c("marginal", "conditional"))

  expect_snapshot(error = TRUE, pool_ipw(fits))
})

test_that("pool_ipw() names the mode a mismatched set does not", {
  # The other half of the refusal above, read from the same fixture. A mode
  # named at the call site says which surface to pool, so the stored field is
  # not consulted and the set that named none is poolable again. This is how
  # `resolve_ipw_effects()` reads an accessor's `effects` argument, and the same
  # request has the same answer here.
  fits <- pool_mixed_mode_fits()

  expect_error(pool_ipw(fits), class = "causalgenerics_pool_mismatch_effects")

  res <- pool_ipw(fits, effects = "marginal", dfcom = 17)

  expect_identical(res$effects, "marginal")
  expect_identical(res$estimates$effect, pool_binary_labels())
  # These results differ from a set that agrees on the mode in the stored field
  # and in nothing else, so naming the mode has to give back exactly what that
  # set gives. An implementation that refused first and read the argument second
  # would never reach this, and one that let the stored field color the result
  # would answer with something else.
  expect_identical(res, pool_ipw(pool_binary_fits(), dfcom = 17))
})

test_that("pool_ipw() refuses results reporting different effects", {
  # Two imputations of a categorical exposure whose contrast sets differ, which
  # is what a factor with a level observed in one imputation and not in another
  # produces. There is no row in the second set to average the `c vs a` row of
  # the first with, and pooling by position would combine `c vs a` with
  # `d vs a`.
  fits <- pool_categorical_fits()
  estimate <- pool_categorical_estimate()
  std_err <- pool_categorical_std_err()
  fits[[2]]$estimates <- pool_categorical_estimates(
    estimate[2, ],
    std_err[2, ],
    contrast = c("b vs a", "d vs a")
  )

  expect_error(pool_ipw(fits), class = "causalgenerics_pool_mismatch_labels")
  expect_error(pool_ipw(fits), class = "causalgenerics_pool_mismatch")

  cnd <- tryCatch(pool_ipw(fits), error = identity)
  expect_identical(cnd$field, "labels")
  expect_true("log(or) d vs a" %in% unlist(cnd$values))

  expect_snapshot(error = TRUE, pool_ipw(fits))
})

test_that("pool_ipw() refuses results reporting effects in different orders", {
  # The same six contrasts in the other order. The set matches, so a check that
  # compared sets would take these as poolable and then average the `b vs a`
  # rows of one imputation with the `c vs a` rows of another. The labels are
  # what say which row is which, and the pooled frame reports them in an order,
  # so the order is part of what has to agree.
  fits <- pool_categorical_fits()
  estimate <- pool_categorical_estimate()
  std_err <- pool_categorical_std_err()
  fits[[3]]$estimates <- pool_categorical_estimates(
    estimate[3, ],
    std_err[3, ],
    contrast = c("c vs a", "b vs a")
  )

  expect_error(pool_ipw(fits), class = "causalgenerics_pool_mismatch_labels")
  expect_error(pool_ipw(fits), class = "causalgenerics_pool_mismatch")
})

test_that("pool_ipw() refuses stored levels that disagree with no level named", {
  # With `conf_level = NULL` the level the pooled bounds report is the one the
  # results were reported at, and a set that disagrees names none. Saying which
  # level to use answers the question, so the same set pools once a level is
  # named.
  estimate <- pool_binary_estimate()
  std_err <- pool_binary_std_err()
  levels <- c(0.95, 0.9, 0.95)
  fits <- lapply(1:3, function(i) {
    pool_fit(
      pool_binary_estimates(estimate[i, ], std_err[i, ], conf.level = levels[i])
    )
  })

  expect_error(
    pool_ipw(fits),
    class = "causalgenerics_pool_mismatch_conf_level"
  )
  expect_error(pool_ipw(fits), class = "causalgenerics_pool_mismatch")

  cnd <- tryCatch(pool_ipw(fits), error = identity)
  expect_identical(cnd$field, "conf_level")
  expect_setequal(unlist(cnd$values), c(0.95, 0.9))

  named <- pool_ipw(fits, conf_level = 0.95)
  expect_true(all(named$estimates$conf.level == 0.95))

  expect_snapshot(error = TRUE, pool_ipw(fits))
})

test_that("pool_ipw() refuses outcome models fitted on different links", {
  # The effects are on the scale the outcome model's link puts them on, so a
  # `log(or)` row from a logit fit and a `log(rr)` row from a log one are not
  # the same quantity even when they are labeled the same way. The link is
  # capturable from both of these models, which is what makes them separable.
  logit <- pool_outcome_model()
  identity_link <- glm(y ~ z, family = gaussian(), data = pool_data())
  fits <- pool_fits_varying(
    "outcome_mod",
    list(logit, identity_link, logit)
  )

  expect_error(
    pool_ipw(fits),
    class = "causalgenerics_pool_mismatch_outcome_link"
  )
  expect_error(pool_ipw(fits), class = "causalgenerics_pool_mismatch")

  cnd <- tryCatch(pool_ipw(fits), error = identity)
  expect_identical(cnd$field, "outcome_link")
  expect_setequal(unlist(cnd$values), c("logit", "identity"))

  expect_snapshot(error = TRUE, pool_ipw(fits))
})

# ---- the arguments -----------------------------------------------------------

test_that("pool_ipw() refuses an `effects` value that names no reading", {
  # The same two readings the rest of the result layer has, checked the same
  # way. `"fixed"` is the word a caller who has used a mixed-model pooling
  # reaches for, and it names neither surface here.
  fits <- pool_binary_fits()

  expect_error(
    pool_ipw(fits, effects = "fixed"),
    class = "causalgenerics_invalid_argument_effects"
  )
  expect_error(
    pool_ipw(fits, effects = "fixed"),
    class = "causalgenerics_invalid_argument"
  )
  expect_error(
    pool_ipw(fits, effects = "both"),
    class = "causalgenerics_invalid_argument_effects"
  )
  expect_error(
    pool_ipw(fits, effects = c("marginal", "conditional")),
    class = "causalgenerics_invalid_argument_effects"
  )
  expect_error(
    pool_ipw(fits, effects = 1),
    class = "causalgenerics_invalid_argument_effects"
  )
  expect_error(
    pool_ipw(fits, effects = NA_character_),
    class = "causalgenerics_invalid_argument_effects"
  )

  expect_snapshot(error = TRUE, pool_ipw(fits, effects = "fixed"))
})

test_that("pool_ipw() refuses a `conf_level` that is not a probability", {
  # The same refusals `as.data.frame()` makes for the level it reports bounds
  # at. The bounds of the open interval go with everything outside it, since no
  # finite pair of limits exists at level 0 or 1, and a vector of levels names
  # several intervals where the frame reports one.
  fits <- pool_binary_fits()

  for (level in list(0, 1, -0.5, 1.5, NA_real_, c(0.9, 0.95), "0.95")) {
    expect_error(
      pool_ipw(fits, conf_level = level),
      class = "causalgenerics_invalid_argument"
    )
  }

  # Keyed to the argument as well as carrying the general class, the way every
  # other refused argument in this package is. The spelling of the key is the
  # implementation's to choose between the argument's own name and the
  # `conf.level` the frame stores it under, so the assertion is on the prefix.
  cnd <- tryCatch(pool_ipw(fits, conf_level = 2), error = identity)
  expect_true(any(startsWith(
    class(cnd),
    "causalgenerics_invalid_argument_conf"
  )))

  expect_no_error(pool_ipw(fits, conf_level = 0.8))

  expect_snapshot(error = TRUE, pool_ipw(fits, conf_level = 2))
})

test_that("pool_ipw() refuses a dfcom that is not a count", {
  # A missing count is the one worth refusing loudest. It is not an error
  # anywhere it is used: it travels through the small-sample adjustment into an
  # `NaN` degrees of freedom, and from there into an `NaN` bound and an `NaN`
  # p-value for every effect, with nothing in the returned result to say where
  # it came from.
  fits <- pool_binary_fits()

  for (value in list("18", c(18, 20), NA_real_, NA, list(18), character())) {
    expect_error(
      pool_ipw(fits, dfcom = value),
      class = "causalgenerics_invalid_argument_dfcom"
    )
  }
  expect_error(
    pool_ipw(fits, dfcom = NA_real_),
    class = "causalgenerics_invalid_argument"
  )

  # `Inf` is the large-sample assumption written out and the value the fallback
  # itself reaches, so it is a count this argument takes rather than one of the
  # non-finite values refused above.
  expect_no_error(pool_ipw(fits, dfcom = Inf))
  expect_no_error(pool_ipw(fits, dfcom = 18L))

  expect_snapshot(error = TRUE, pool_ipw(fits, dfcom = "18"))
})

test_that("pool_ipw() refuses arguments it has no name for", {
  # Everything after `fits` is matched by name, so a positional argument and a
  # misspelled one both land in the dots. Ignoring them would run the pooling at
  # the default the caller was trying to change and return a result that looks
  # like the one asked for.
  fits <- pool_binary_fits()

  expect_error(
    pool_ipw(fits, "conditional"),
    class = "causalgenerics_invalid_argument"
  )
  expect_error(
    pool_ipw(fits, effect = "conditional"),
    class = "causalgenerics_invalid_argument"
  )
  expect_error(
    pool_ipw(fits, dfcomm = 12),
    class = "causalgenerics_invalid_argument"
  )

  expect_no_error(
    pool_ipw(fits, effects = "marginal", dfcom = 12, conf_level = 0.9)
  )

  expect_snapshot(error = TRUE, pool_ipw(fits, "conditional"))
  expect_snapshot(error = TRUE, pool_ipw(fits, effect = "conditional"))
})

# ---- Rubin's rules -----------------------------------------------------------

test_that("pool_ipw() pools one effect by numbers that can be written down", {
  # The `rd` row of the fixture: estimates of 0.20, 0.30, and 0.40 against
  # standard errors of 0.10, 0.20, and 0.30, over m = 3 imputations.
  #
  #   qbar  = (0.20 + 0.30 + 0.40) / 3                     = 0.30
  #   ubar  = (0.01 + 0.04 + 0.09) / 3                     = 0.14 / 3
  #   b     = (0.01 + 0.00 + 0.01) / 2                     = 0.01
  #   t     = ubar + (1 + 1/3) * b = 0.14/3 + 0.04/3       = 0.06
  #
  # so the pooled standard error is sqrt(0.06) exactly. None of the three
  # standard errors is that number and neither is their mean, which is what
  # separates the total variance from the within-imputation part of it.
  fits <- pool_binary_fits()
  res <- pool_ipw(fits, dfcom = 17)
  rd <- res$estimates[res$estimates$effect == "rd", ]

  expect_identical(rd$estimate, 0.30)
  expect_equal(rd$std.err, sqrt(0.06))
  expect_false(isTRUE(all.equal(rd$std.err, sqrt(0.14 / 3))))
})

test_that("pool_ipw() pools every effect by Rubin's rules", {
  # The other two labels, recomputed from the fixture inputs. The three effects
  # carry different amounts of between-imputation variance, so a pooling that
  # dropped the `(1 + 1/m) * b` term would miss all three by different amounts.
  estimate <- pool_binary_estimate()
  std_err <- pool_binary_std_err()
  fits <- pool_binary_fits()
  res <- pool_ipw(fits, dfcom = 17)
  expected <- rubin_rules_by_column(estimate, std_err, dfcom = 17)

  expect_identical(res$estimates$effect, pool_binary_labels())
  expect_equal(res$estimates$estimate, rubin_column(expected, "estimate"))
  expect_equal(res$estimates$std.err, rubin_column(expected, "std.err"))
  # The pooled estimate is the mean of the stored ones and nothing else, which
  # is the half of Rubin's rules a weighted average would get wrong.
  expect_equal(res$estimates$estimate, colMeans(estimate))
})

test_that("pool_ipw() adjusts the degrees of freedom for the sample size", {
  # The Barnard-Rubin small-sample adjustment for the `rd` row, worked out from
  # the quantities above with dfcom = 17.
  #
  #   lambda = (1 + 1/3) * 0.01 / 0.06                         = 2 / 9
  #   dfold  = (m - 1) / lambda^2 = 2 / (4/81)                 = 40.5
  #   tmp    = (1 - 2/9) * (1 + 17) * 17 = (7/9) * 306         = 238
  #   df     = 2 * 238 / (20 * 2 + (4/81) * 238)
  #          = 476 / (4192/81) = 38556 / 4192                  = 9639 / 1048
  #
  # which is 9.1975..., far below both the eighteen residual degrees of freedom
  # of a single fit and the large-sample limit. An implementation that reported
  # `dfold` here would be out by a factor of four.
  fits <- pool_binary_fits()
  res <- pool_ipw(fits, dfcom = 17)
  rd <- res$estimates[res$estimates$effect == "rd", ]

  expect_equal(rd$df, 9639 / 1048)
  expect_false(isTRUE(all.equal(rd$df, 40.5)))
  expect_false(isTRUE(all.equal(rd$df, 17)))
})

test_that("pool_ipw() reports the large-sample degrees of freedom at infinity", {
  # With no complete-data degrees of freedom to spend, the adjustment has
  # nothing to adjust and the answer is the unadjusted `(m - 1) / lambda^2`.
  # For the `rd` row that is 2 / (2/9)^2 = 40.5 exactly.
  fits <- pool_binary_fits()
  res <- pool_ipw(fits, dfcom = Inf)
  rd <- res$estimates[res$estimates$effect == "rd", ]

  expect_equal(rd$df, 40.5)

  estimate <- pool_binary_estimate()
  std_err <- pool_binary_std_err()
  expected <- rubin_rules_by_column(estimate, std_err, dfcom = Inf)
  expect_equal(res$estimates$df, rubin_column(expected, "df"))
  # The adjusted numbers are smaller everywhere, so the assertions above are
  # about the limit rather than about a fixture where the two agree.
  adjusted <- pool_ipw(fits, dfcom = 17)$estimates$df
  expect_true(all(adjusted < res$estimates$df))
})

# ---- the pooled estimates frame ----------------------------------------------

test_that("pool_ipw() reports the pooled estimates in the documented shape", {
  # One row per effect and the columns in this order. The frame is the surface a
  # caller reads and a tidier is written against, so the names and their order
  # are the contract rather than a presentation detail.
  fits <- pool_binary_fits()
  res <- pool_ipw(fits, dfcom = 17)

  expect_s3_class(res$estimates, "data.frame")
  expect_identical(nrow(res$estimates), 3L)
  expect_identical(
    names(res$estimates),
    c(
      "effect",
      "estimate",
      "std.err",
      "t",
      "df",
      "ci.lower",
      "ci.upper",
      "conf.level",
      "p.value"
    )
  )
  expect_identical(res$estimates$effect, pool_binary_labels())
})

test_that("pool_ipw() reports a t statistic and a p-value on its own df", {
  # The statistic is the pooled estimate over the pooled standard error, and it
  # is referred to the pooled degrees of freedom rather than to the normal.
  # Those degrees of freedom are in single figures here, so a normal p-value
  # would be visibly smaller than the right one.
  fits <- pool_binary_fits()
  res <- pool_ipw(fits, dfcom = 17)
  estimates <- res$estimates

  expect_equal(estimates$t, estimates$estimate / estimates$std.err)
  expect_equal(
    estimates$p.value,
    2 * stats::pt(-abs(estimates$t), estimates$df)
  )
  expect_false(isTRUE(all.equal(
    estimates$p.value,
    2 * stats::pnorm(-abs(estimates$t))
  )))
})

test_that("pool_ipw() bounds the pooled estimates on its own df", {
  # The limits are one half width added and subtracted, taken from the t
  # quantile at the pooled degrees of freedom. `qt()` is not exactly
  # antisymmetric, so a lower limit taken from the lower tail would differ in
  # the last bit, and a normal quantile would be too narrow throughout.
  fits <- pool_binary_fits()
  res <- pool_ipw(fits, dfcom = 17)
  estimates <- res$estimates
  half_width <- stats::qt(0.975, estimates$df) * estimates$std.err

  expect_equal(estimates$ci.lower, estimates$estimate - half_width)
  expect_equal(estimates$ci.upper, estimates$estimate + half_width)
  expect_identical(estimates$conf.level, rep(0.95, 3))
})

test_that("pool_ipw() reports the bounds at the level asked for", {
  # A level named at the call site is used for the bounds and recorded in the
  # column, whatever the results were reported at.
  fits <- pool_binary_fits()
  res <- pool_ipw(fits, dfcom = 17, conf_level = 0.8)
  estimates <- res$estimates
  half_width <- stats::qt(0.9, estimates$df) * estimates$std.err

  expect_identical(estimates$conf.level, rep(0.8, 3))
  expect_equal(estimates$ci.lower, estimates$estimate - half_width)
  expect_equal(estimates$ci.upper, estimates$estimate + half_width)
})

test_that("pool_ipw() takes the level the results agree on when none is named", {
  # With no level named the stored one is used rather than a default of 0.95,
  # which is what makes the mismatch above worth refusing. The fixture stores
  # 0.9, so a hardcoded default would fail here.
  estimate <- pool_binary_estimate()
  std_err <- pool_binary_std_err()
  fits <- lapply(1:3, function(i) {
    pool_fit(
      pool_binary_estimates(estimate[i, ], std_err[i, ], conf.level = 0.9)
    )
  })

  res <- pool_ipw(fits, dfcom = 17)
  half_width <- stats::qt(0.95, res$estimates$df) * res$estimates$std.err

  expect_identical(res$estimates$conf.level, rep(0.9, 3))
  expect_equal(res$estimates$ci.lower, res$estimates$estimate - half_width)
})

test_that("pool_ipw() keys a categorical result by effect and contrast", {
  # Six rows, the contrast column in the position the `ipw()` contract puts it,
  # and the labels in the order the results reported them. Pooling groups the
  # rows by label, so a categorical result is where a grouping that read only
  # the `effect` column would combine the `b vs a` and `c vs a` rows of each
  # measure into one.
  fits <- pool_categorical_fits()
  estimate <- pool_categorical_estimate()
  std_err <- pool_categorical_std_err()
  res <- pool_ipw(fits, dfcom = 17)
  expected <- rubin_rules_by_column(estimate, std_err, dfcom = 17)

  expect_identical(nrow(res$estimates), 6L)
  expect_identical(
    names(res$estimates),
    c(
      "effect",
      "contrast",
      "estimate",
      "std.err",
      "t",
      "df",
      "ci.lower",
      "ci.upper",
      "conf.level",
      "p.value"
    )
  )
  expect_identical(
    res$estimates$effect,
    rep(c("rd", "log(rr)", "log(or)"), times = 2)
  )
  expect_identical(
    res$estimates$contrast,
    rep(c("b vs a", "c vs a"), each = 3)
  )
  expect_equal(res$estimates$estimate, rubin_column(expected, "estimate"))
  expect_equal(res$estimates$std.err, rubin_column(expected, "std.err"))
})

test_that("pool_ipw() pools a legacy contrast column as the canonical one", {
  # A result stored against an earlier version of the contract names the column
  # holding its contrasts `comparison`. The labels the pooling groups by are
  # read through the same helper the accessors read them with, so such a result
  # is keyed exactly as one built now is: it pools with the same rows in the
  # same order and reports them under the canonical heading. The two results are
  # the same object rather than merely the same numbers, which is what says the
  # older spelling is read as an alias rather than as a shape of its own.
  canonical <- pool_categorical_fits(vcov = TRUE)
  legacy <- lapply(canonical, function(fit) {
    names(fit$estimates)[names(fit$estimates) == "contrast"] <- "comparison"
    fit
  })

  expect_identical(names(legacy[[1]]$estimates)[[2]], "comparison")
  # The rename is the only difference, covariance attribute included, so the
  # comparison below is about the column name and not about two fixtures that
  # drifted apart.
  expect_identical(
    attr(legacy[[1]]$estimates, "ipw_vcov", exact = TRUE),
    attr(canonical[[1]]$estimates, "ipw_vcov", exact = TRUE)
  )

  res <- pool_ipw(legacy, dfcom = 17)

  expect_identical(res, pool_ipw(canonical, dfcom = 17))
  expect_true("contrast" %in% names(res$estimates))
  expect_false("comparison" %in% names(res$estimates))
  expect_identical(res$estimates$contrast, rep(c("b vs a", "c vs a"), each = 3))
  expect_identical(res$pooling$contrast, res$estimates$contrast)
})

test_that("pool_ipw() pools the stored values on the scale they were stored", {
  # The ratio rows are on the log scale, which is where the inference is done
  # and where an average of three estimates means something. Pooling the
  # exponentials instead would give a different number under the same label, and
  # nothing in the frame would say the scale had changed.
  estimate <- pool_binary_estimate()
  fits <- pool_binary_fits()
  res <- pool_ipw(fits, dfcom = 17)

  expect_identical(res$estimates$effect, c("rd", "log(rr)", "log(or)"))
  expect_false(any(res$estimates$effect %in% c("rr", "or")))
  expect_equal(res$estimates$estimate[2], mean(estimate[, 2]))
  expect_false(isTRUE(all.equal(
    res$estimates$estimate[2],
    log(mean(exp(estimate[, 2])))
  )))
})

# ---- the pooling diagnostics -------------------------------------------------

test_that("pool_ipw() reports the pooling diagnostics keyed the same way", {
  # The quantities the pooling itself produces, which say how much of the
  # uncertainty came from the imputation rather than from the data. They sit in
  # a frame of their own rather than beside the estimates, keyed by the same
  # labels so that a row of one is found from a row of the other.
  estimate <- pool_binary_estimate()
  std_err <- pool_binary_std_err()
  fits <- pool_binary_fits()
  res <- pool_ipw(fits, dfcom = 17)
  expected <- rubin_rules_by_column(estimate, std_err, dfcom = 17)

  expect_s3_class(res$pooling, "data.frame")
  expect_identical(
    names(res$pooling),
    c("effect", "ubar", "b", "riv", "lambda", "fmi")
  )
  expect_identical(res$pooling$effect, pool_binary_labels())
  expect_equal(res$pooling$ubar, rubin_column(expected, "ubar"))
  expect_equal(res$pooling$b, rubin_column(expected, "b"))
  expect_equal(res$pooling$riv, rubin_column(expected, "riv"))
  expect_equal(res$pooling$lambda, rubin_column(expected, "lambda"))
  expect_equal(res$pooling$fmi, rubin_column(expected, "fmi"))
})

test_that("pool_ipw() reports diagnostics that can be written down", {
  # The `rd` row again, from the quantities the literal test above works out.
  #
  #   ubar   = 0.14 / 3
  #   b      = 0.01
  #   riv    = (1 + 1/3) * 0.01 / (0.14/3) = 0.04 / 0.14   = 2 / 7
  #   lambda = (1 + 1/3) * 0.01 / 0.06                     = 2 / 9
  #
  # `riv` is measured against the within-imputation variance and `lambda`
  # against the total, so the two differ. An implementation that computed one
  # and reported it twice would agree with only one of these.
  fits <- pool_binary_fits()
  res <- pool_ipw(fits, dfcom = 17)
  rd <- res$pooling[res$pooling$effect == "rd", ]

  expect_equal(rd$ubar, 0.14 / 3)
  expect_equal(rd$b, 0.01)
  expect_equal(rd$riv, 2 / 7)
  expect_equal(rd$lambda, 2 / 9)
  # The fraction of missing information is built from `riv` and the pooled
  # degrees of freedom, so it moves with the small-sample adjustment.
  expect_equal(rd$fmi, (2 / 7 + 2 / (9639 / 1048 + 3)) / (2 / 7 + 1))
})

test_that("pool_ipw() keys the diagnostics of a categorical result too", {
  fits <- pool_categorical_fits()
  estimate <- pool_categorical_estimate()
  std_err <- pool_categorical_std_err()
  res <- pool_ipw(fits, dfcom = 17)
  expected <- rubin_rules_by_column(estimate, std_err, dfcom = 17)

  expect_identical(
    names(res$pooling),
    c("effect", "contrast", "ubar", "b", "riv", "lambda", "fmi")
  )
  expect_identical(res$pooling$effect, res$estimates$effect)
  expect_identical(res$pooling$contrast, res$estimates$contrast)
  expect_equal(res$pooling$fmi, rubin_column(expected, "fmi"))
})

# ---- the pooled covariance ---------------------------------------------------

test_that("pool_ipw() pools the covariance of the effects", {
  # The same rule the standard errors follow, applied to the whole matrix: the
  # average of the within-imputation covariances plus the between-imputation one
  # inflated by `1 + 1/m`. The off-diagonal entries are what the matrix is for,
  # and they are not recoverable from the pooled standard errors.
  estimate <- pool_binary_estimate()
  std_err <- pool_binary_std_err()
  fits <- pool_binary_fits(vcov = TRUE)
  res <- pool_ipw(fits, dfcom = 17)

  within <- lapply(
    seq_len(nrow(std_err)),
    function(i) pool_effects_vcov(std_err[i, ], pool_binary_labels())
  )
  ubar_mat <- Reduce(`+`, within) / length(within)
  between <- stats::cov(estimate)
  expected <- ubar_mat + (1 + 1 / 3) * between
  dimnames(expected) <- list(pool_binary_labels(), pool_binary_labels())

  covariance <- attr(res$estimates, "ipw_vcov", exact = TRUE)

  expect_equal(covariance, expected)
  expect_identical(
    dimnames(covariance),
    list(pool_binary_labels(), pool_binary_labels())
  )
  # The diagonal is the square of the pooled standard errors, since the two are
  # the same combination of the same numbers. A block built any other way would
  # report a variance beside a standard error that is not its square root.
  expect_equal(diag(covariance), res$estimates$std.err^2, ignore_attr = TRUE)
  # Non-zero off-diagonals, so the assertions above are about the covariance
  # rather than about a fixture where the effects happen to be uncorrelated.
  expect_true(all(covariance[upper.tri(covariance)] > 0))
})

test_that("pool_ipw() names the pooled covariance for a categorical result", {
  fits <- pool_categorical_fits(vcov = TRUE)
  res <- pool_ipw(fits, dfcom = 17)
  covariance <- attr(res$estimates, "ipw_vcov", exact = TRUE)

  expect_identical(dim(covariance), c(6L, 6L))
  expect_identical(
    dimnames(covariance),
    list(pool_categorical_labels(), pool_categorical_labels())
  )
})

test_that("pool_ipw() attaches no covariance when a result records none", {
  # The average of the within-imputation covariances needs one from every
  # imputation. A pooling over the results that have one would report a matrix
  # built from a subset of the imputations beside estimates built from all of
  # them, and nothing on the object would say so. The estimates are pooled
  # either way, since the standard errors are in the frame.
  fits <- pool_binary_fits(vcov = TRUE)
  fits[[2]]$estimates <- pool_binary_estimates(
    pool_binary_estimate()[2, ],
    pool_binary_std_err()[2, ]
  )

  res <- pool_ipw(fits, dfcom = 17)

  expect_null(attr(res$estimates, "ipw_vcov", exact = TRUE))
  expect_identical(nrow(res$estimates), 3L)

  none <- pool_ipw(pool_binary_fits(), dfcom = 17)
  expect_null(attr(none$estimates, "ipw_vcov", exact = TRUE))
})

# ---- the complete-data degrees of freedom ------------------------------------

test_that("pool_ipw() uses the complete-data df it is given", {
  # A number at the call site settles it, whatever the results record. The
  # fixture's own chain would give eighteen, so a supplied twelve is separable
  # from it.
  fits <- pool_binary_fits()
  res <- pool_ipw(fits, dfcom = 12)

  expect_equal(res$dfcom, 12)
  expect_equal(
    res$estimates$df,
    rubin_column(
      rubin_rules_by_column(
        pool_binary_estimate(),
        pool_binary_std_err(),
        dfcom = 12
      ),
      "df"
    )
  )
})

test_that("pool_ipw() floors the complete-data df at one", {
  # Zero and below name no fit at all, and the adjustment is undefined there.
  # One is the smallest count a fit can have left, so that is what is used.
  fits <- pool_binary_fits()

  expect_equal(pool_ipw(fits, dfcom = 0)$dfcom, 1)
  expect_equal(pool_ipw(fits, dfcom = -3)$dfcom, 1)
  expect_equal(
    pool_ipw(fits, dfcom = 0)$estimates$df,
    pool_ipw(fits, dfcom = 1)$estimates$df
  )
})

test_that("pool_ipw() reads the smallest df the variance objects report", {
  # The fitted variance object is the first place the complete-data degrees of
  # freedom are looked for, through `df.residual()` on the result. The smallest
  # is taken rather than the first, since the pooled inference is no stronger
  # than the weakest imputation supports and `mice::pool()`'s habit of reading
  # only the first fit hides a result fitted on less data.
  fits <- pool_fits_varying(
    "fit",
    lapply(
      c(25, 17, 21),
      function(df) structure(list(df.residual = df), class = "cg_pool_fit")
    )
  )

  expect_identical(df.residual(fits[[2]]), 17L)

  res <- pool_ipw(fits)

  expect_equal(res$dfcom, 17)
  expect_equal(res$estimates$df[[1]], 9639 / 1048)
})

test_that("pool_ipw() reads past the results that report no df", {
  # One result records a fitted variance object and the other two record none,
  # which is the shape a set takes when the imputations did not all reach the
  # same estimator. The count comes from the results that report one.
  #
  # A minimum taken across the missing ones as well is `NA`, and `NA` travels:
  # the small-sample adjustment gives `NaN` degrees of freedom, every bound and
  # every p-value in the frame comes back `NaN`, and nothing on the object says
  # where that came from. The outcome models here report eighteen, so the
  # fallback below is reachable and is not what answers.
  fits <- pool_binary_fits()
  fits[[1]]$fit <- structure(list(df.residual = 25), class = "cg_pool_fit")

  expect_identical(df.residual(fits[[1]]), 25L)
  expect_identical(df.residual(fits[[2]]), NA_integer_)
  expect_identical(stats::df.residual(fits[[2]]$outcome_mod), 18L)

  res <- pool_ipw(fits)

  expect_equal(res$dfcom, 25)
  expect_false(is.na(res$dfcom))
  expect_equal(
    res$estimates$df,
    rubin_column(
      rubin_rules_by_column(
        pool_binary_estimate(),
        pool_binary_std_err(),
        dfcom = 25
      ),
      "df"
    )
  )
  expect_false(anyNA(res$estimates$df))
  expect_false(anyNA(res$estimates$p.value))
  # Not the outcome models' eighteen, which is the answer the step below this
  # one would give and is close enough to pass an assertion that only checked
  # for a number.
  expect_false(isTRUE(all.equal(res$dfcom, 18)))
})

test_that("pool_ipw() falls back on the outcome models' residual df", {
  # The linearization path records no fitted variance object, so `df.residual()`
  # on the result is `NA` and the outcome model is where the count comes from.
  # Three models on different amounts of data, so the smallest is separable from
  # the first and from the last.
  fits <- pool_fits_varying(
    "outcome_mod",
    lapply(c(20, 16, 18), pool_outcome_model)
  )

  expect_identical(df.residual(fits[[1]]), NA_integer_)
  expect_identical(stats::df.residual(fits[[2]]$outcome_mod), 14L)

  res <- pool_ipw(fits)

  expect_equal(res$dfcom, 14)
  expect_equal(
    res$estimates$df,
    rubin_column(
      rubin_rules_by_column(
        pool_binary_estimate(),
        pool_binary_std_err(),
        dfcom = 14
      ),
      "df"
    )
  )
})

test_that("pool_ipw() assumes a large sample when nothing reports a df", {
  # A result whose variance object and outcome model both report nothing leaves
  # no count to adjust with. The large-sample limit is the honest answer, and it
  # is the widest one, so it is said out loud rather than taken silently: the
  # intervals it gives are narrower than the ones a real complete-data count
  # would give, which is the direction a reader would not question.
  # The suite is run under `options(warn = 2)` so that a stray warning fails it,
  # and this is the one test whose subject is a warning. `expect_warning()`
  # muffles before that option converts anything, but `expect_snapshot()` does
  # not, so the snapshot below would record an error rather than the warning it
  # is there to pin. The option is restored to its default for this block only;
  # the warning is still asserted, by class and by message.
  withr::local_options(warn = 0)

  local_s3_method("nobs", "cg_pool_bare", function(object, ...) 20L)
  bare <- pool_bare_model()
  fits <- pool_fits_varying("outcome_mod", list(bare, bare, bare))

  expect_null(stats::df.residual(bare))

  expect_warning(
    pool_ipw(fits),
    class = "causalgenerics_pool_large_sample"
  )

  res <- suppressWarnings(pool_ipw(fits))

  expect_identical(res$dfcom, Inf)
  expect_equal(res$estimates$df[[1]], 40.5)
  expect_equal(
    res$estimates$df,
    rubin_column(
      rubin_rules_by_column(
        pool_binary_estimate(),
        pool_binary_std_err(),
        dfcom = Inf
      ),
      "df"
    )
  )
  # Naming a count answers the question, so the warning belongs to the fallback
  # rather than to the fits.
  expect_no_warning(pool_ipw(fits, dfcom = 17))

  expect_snapshot(pooled <- pool_ipw(fits))
})

# ---- the pooled result -------------------------------------------------------

test_that("pool_ipw() returns an object of its own class", {
  # A pooled result is not an `ipw` result. It reports several imputations
  # rather than one fit, holds no component models, and reports its inference
  # against t rather than z, so the methods registered for `ipw` would answer
  # for fields it does not carry.
  fits <- pool_binary_fits()
  res <- pool_ipw(fits, dfcom = 17)

  expect_s3_class(res, "ipw_pooled")
  expect_identical(class(res), "ipw_pooled")
  expect_false(inherits(res, "ipw"))
})

test_that("pool_ipw() records the documented fields in order", {
  # The field names and their order are the contract, since callers read fields
  # by name and a printed form writes them positionally.
  fits <- pool_binary_fits()
  res <- pool_ipw(fits, dfcom = 17)

  expect_identical(
    names(res),
    c(
      "estimand",
      "estimates",
      "pooling",
      "se_method",
      "effects",
      "m",
      "dfcom",
      "nobs",
      "outcome_link",
      "alternate"
    )
  )
  expect_identical(res$estimand, "ate")
  expect_identical(res$se_method, "linearization")
  expect_identical(res$effects, "marginal")
  expect_identical(res$m, 3L)
  expect_identical(res$outcome_link, "logit")
})

test_that("pool_ipw() records the fields the results agree on", {
  # The three scalar fields are carried across rather than defaulted, so a set
  # of results recording something else reports that instead.
  fits <- pool_binary_fits(
    estimand = "att",
    se_method = "mestimation"
  )
  res <- pool_ipw(fits, dfcom = 17)

  expect_identical(res$estimand, "att")
  expect_identical(res$se_method, "mestimation")
})

test_that("pool_ipw() counts the imputations", {
  fits <- pool_binary_fits()

  expect_identical(pool_ipw(fits, dfcom = 17)$m, 3L)
  expect_identical(pool_ipw(fits[1:2], dfcom = 17)$m, 2L)
})

test_that("pool_ipw() reports the smallest observation count", {
  # Imputations of the same incomplete data are the same size, but a result
  # whose method dropped rows, such as one fitted after trimming, is smaller.
  # The pooled result reports the count the weakest imputation contributed, for
  # the reason the complete-data degrees of freedom are the smallest rather than
  # the first.
  fits <- pool_fits_varying(
    "outcome_mod",
    lapply(c(20, 16, 18), pool_outcome_model)
  )

  expect_identical(nobs(fits[[2]]), 16L)

  res <- pool_ipw(fits, dfcom = 17)

  expect_identical(res$nobs, 16L)
  expect_identical(pool_ipw(pool_binary_fits(), dfcom = 17)$nobs, 20L)
})

test_that("pool_ipw() records the link the outcome models share", {
  # The link says what scale the pooled effects are on, which is what makes the
  # refusal of a set that disagrees about it meaningful. A model that is not a
  # generalized linear one reports no link of its own, and the identity is what
  # its coefficients are on.
  logit <- pool_ipw(pool_binary_fits(), dfcom = 17)
  expect_identical(logit$outcome_link, "logit")

  linear <- pool_ipw(
    pool_fits_varying(
      "outcome_mod",
      rep(list(lm(y ~ z, data = pool_data())), 3)
    ),
    dfcom = 17
  )
  expect_identical(linear$outcome_link, "identity")

  gaussian_glm <- pool_ipw(
    pool_fits_varying(
      "outcome_mod",
      rep(list(glm(y ~ z, family = gaussian(), data = pool_data())), 3)
    ),
    dfcom = 17
  )
  expect_identical(gaussian_glm$outcome_link, "identity")
})

# ---- the conditional reading -------------------------------------------------

test_that("pool_ipw() pools the outcome coefficients in the conditional mode", {
  # The conditional reading presents the outcome model's coefficient surface, so
  # that is what the pooling combines: the coefficients are the estimates and
  # the square roots of the corrected block's diagonal are their standard
  # errors. The rules are the same ones the marginal reading uses.
  models <- pool_conditional_models()
  fits <- pool_conditional_fits()
  labels <- names(coef(models[[1]]))
  estimate <- t(vapply(models, coef, numeric(2)))
  std_err <- t(vapply(models, function(m) sqrt(diag(vcov(m))), numeric(2)))
  expected <- rubin_rules_by_column(estimate, std_err, dfcom = 18)

  res <- pool_ipw(fits, dfcom = 18)

  expect_identical(res$effects, "conditional")
  expect_identical(res$estimates$effect, labels)
  expect_false("contrast" %in% names(res$estimates))
  expect_identical(nrow(res$estimates), 2L)
  expect_equal(res$estimates$estimate, rubin_column(expected, "estimate"))
  expect_equal(res$estimates$std.err, rubin_column(expected, "std.err"))
  expect_equal(res$pooling$b, rubin_column(expected, "b"))
  # The between-imputation variance is real here, so the pooled standard errors
  # are wider than the average within-imputation one.
  expect_true(all(res$pooling$b > 0))
  # Not the marginal surface, which these same results also carry. The two
  # tables differ in their rows as well as their numbers, so a pooling that
  # reported the wrong one could not pass this.
  expect_false(any(res$estimates$effect %in% pool_binary_labels()))
})

test_that("pool_ipw() pools the corrected blocks in the conditional mode", {
  # Every result in this reading carries a covariance by construction, since the
  # standard errors were taken from its diagonal, so the pooled block is
  # attached rather than withheld. It combines the same way the marginal one
  # does, and its margins are named by coefficient rather than by effect: a
  # caller reading a variance out by the name `coef()` gave in this mode has to
  # get that coefficient's.
  models <- pool_conditional_models()
  fits <- pool_conditional_fits()
  labels <- names(coef(models[[1]]))
  estimate <- t(vapply(models, coef, numeric(2)))
  within <- lapply(models, vcov)
  expected <- Reduce(`+`, within) / 3 + (1 + 1 / 3) * stats::cov(estimate)
  dimnames(expected) <- list(labels, labels)

  res <- pool_ipw(fits, dfcom = 18)
  covariance <- attr(res$estimates, "ipw_vcov", exact = TRUE)

  expect_equal(covariance, expected)
  expect_identical(dimnames(covariance), list(labels, labels))
  expect_identical(dim(covariance), c(2L, 2L))
  # The diagonal is the square of the pooled standard errors, since the two are
  # the same combination of the same numbers.
  expect_equal(diag(covariance), res$estimates$std.err^2, ignore_attr = TRUE)
  # Not the effect labels, which is what the other reading would have named it
  # with, and not the average of the corrected blocks on its own either.
  expect_false(any(labels %in% pool_binary_labels()))
  expect_false(isTRUE(all.equal(covariance, Reduce(`+`, within) / 3)))
})

test_that("pool_ipw() reads the mode the results record", {
  # With no mode named at the call site the stored one decides, the same way it
  # does for the accessors. Naming the mode a set already records is the same
  # request, so the two results are the same object.
  stored <- pool_ipw(pool_conditional_fits(), dfcom = 18)
  named <- pool_ipw(
    pool_conditional_fits(effects = "marginal"),
    dfcom = 18,
    effects = "conditional"
  )

  expect_identical(stored$effects, "conditional")
  expect_identical(named, stored)
})

test_that("pool_ipw() refuses the conditional mode without corrected blocks", {
  # The standard errors of the conditional reading are the corrected ones, and
  # the covariance an unwrapped outcome model computed for itself is not a
  # substitute: it treats the estimated weights as fixed. The refusal is the one
  # the accessors raise for the same model, so a caller handles one class
  # whichever surface asked.
  fits <- pool_binary_fits()

  expect_error(
    pool_ipw(fits, effects = "conditional"),
    class = "causalgenerics_no_conditional_vcov"
  )
  expect_error(
    pool_ipw(fits, effects = "conditional"),
    class = "causalgenerics_no_vcov"
  )
  # The fixture is discriminating: the naive covariance is right there to be
  # taken by an implementation that let the guard fall through.
  expect_no_error(stats::vcov(fits[[1]]$outcome_mod))
  # The marginal reading of the same results is unaffected.
  expect_no_error(pool_ipw(fits))

  expect_snapshot(error = TRUE, pool_ipw(fits, effects = "conditional"))
})

# ---- the alternate surface ---------------------------------------------------

test_that("pool_ipw() pools the reading it was not asked for as well", {
  # Both readings are pooled from one set of results, so which one the call
  # named decides which pair of frames is the active one rather than what gets
  # computed at all. The other pair is recorded under `alternate`, in the shape
  # the active fields have for that reading, which is what lets a caller move a
  # pooled result between the two readings afterwards.
  fits <- pool_both_surface_fits()
  labels <- names(coef(pool_conditional_models()[[1]]))

  res <- pool_ipw(fits, dfcom = 18)

  expect_identical(res$effects, "marginal")
  expect_identical(names(res$alternate), c("effects", "estimates", "pooling"))
  expect_identical(res$alternate$effects, "conditional")

  # The same two frames the active reading holds, keyed by coefficient rather
  # than by effect label.
  expect_s3_class(res$alternate$estimates, "data.frame")
  expect_s3_class(res$alternate$pooling, "data.frame")
  expect_identical(names(res$alternate$estimates), names(res$estimates))
  expect_identical(names(res$alternate$pooling), names(res$pooling))
  expect_identical(res$alternate$estimates$effect, labels)
  expect_identical(res$alternate$pooling$effect, labels)
  expect_identical(nrow(res$alternate$estimates), 2L)

  # The pooled block travels with the frame it describes and is named the way
  # that frame is named. A block carrying the active reading's labels would
  # report one surface's covariance under the other surface's names.
  covariance <- attr(res$alternate$estimates, "ipw_vcov", exact = TRUE)
  expect_identical(dimnames(covariance), list(labels, labels))
  expect_equal(
    diag(covariance),
    res$alternate$estimates$std.err^2,
    ignore_attr = TRUE
  )
  # The two surfaces are separable by their labels alone, so nothing asserted
  # above could be met by the active frames.
  expect_false(any(labels %in% pool_binary_labels()))
})

test_that("pool_ipw() gives the same four frames whichever reading it pools", {
  # The property the two-surface pooling is for: pooling the conditional
  # reading and pooling the marginal one over the same results compute the same
  # four frames, and differ in which pair of them is the active one and in
  # nothing else. Every shared field is settled by the results rather than by
  # the reading, so those come back the same as well.
  fits <- pool_both_surface_fits()
  shared <- c("estimand", "se_method", "m", "dfcom", "nobs", "outcome_link")

  marginal <- pool_ipw(fits, effects = "marginal", dfcom = 18)
  conditional <- pool_ipw(fits, effects = "conditional", dfcom = 18)

  expect_identical(conditional$alternate$estimates, marginal$estimates)
  expect_identical(conditional$alternate$pooling, marginal$pooling)
  expect_identical(marginal$alternate$estimates, conditional$estimates)
  expect_identical(marginal$alternate$pooling, conditional$pooling)
  expect_identical(conditional$alternate$effects, marginal$effects)
  expect_identical(marginal$alternate$effects, conditional$effects)
  expect_identical(marginal[shared], conditional[shared])

  # The two active frames report different rows, so the identities above are
  # about a swap rather than about a fixture whose surfaces agree.
  expect_false(identical(marginal$estimates, conditional$estimates))
})

test_that("pool_ipw() swaps the same frames at a level named for both", {
  # A level named at the call site settles the bounds of both readings, so the
  # swap is between frames built at that level on either side.
  fits <- pool_both_surface_fits()

  marginal <- pool_ipw(
    fits,
    effects = "marginal",
    dfcom = 18,
    conf_level = 0.8
  )
  conditional <- pool_ipw(
    fits,
    effects = "conditional",
    dfcom = 18,
    conf_level = 0.8
  )

  expect_identical(conditional$alternate$estimates, marginal$estimates)
  expect_identical(conditional$alternate$pooling, marginal$pooling)
  expect_identical(marginal$alternate$estimates, conditional$estimates)
  expect_identical(marginal$alternate$pooling, conditional$pooling)
  expect_identical(marginal$estimates$conf.level, rep(0.8, 3))
  expect_identical(marginal$alternate$estimates$conf.level, rep(0.8, 2))
})

test_that("pool_ipw() resolves the level of each surface on its own", {
  # The level is stored on the marginal surface alone: an estimates frame
  # records the level its bounds were reported at, and a coefficient surface
  # records none. With no level named the two readings therefore resolve
  # different ones, each from what its own surface says, and which reading is
  # the active one does not decide the other's bounds.
  fits <- pool_both_surface_fits(conf.level = rep(0.9, 3))

  res <- pool_ipw(fits, effects = "conditional", dfcom = 18)
  alternate <- res$alternate$estimates
  # The same two levels the other way round when the marginal reading is the
  # active one.
  flipped <- pool_ipw(fits, effects = "marginal", dfcom = 18)

  expect_identical(res$estimates$conf.level, rep(0.95, 2))
  expect_identical(alternate$conf.level, rep(0.9, 3))
  expect_identical(flipped$estimates$conf.level, rep(0.9, 3))
  expect_identical(flipped$alternate$estimates$conf.level, rep(0.95, 2))

  # The bounds are built at the level the frame records rather than merely
  # labeled with it.
  half_width <- stats::qt(0.95, alternate$df) * alternate$std.err

  expect_equal(alternate$ci.lower, alternate$estimate - half_width)
  expect_equal(alternate$ci.upper, alternate$estimate + half_width)
})

test_that("pool_ipw() records why the other reading could not be pooled", {
  # The conditional reading is pooled from the corrected block a fitting
  # package attaches with `new_ipw_model()`, and results from the linearization
  # path carry none. The pooling the caller asked for is done regardless, and
  # what stopped the other reading is recorded rather than raised, since
  # nothing about the result that was asked for is wrong.
  fits <- pool_binary_fits()
  cnd <- tryCatch(
    pool_ipw(fits, effects = "conditional", dfcom = 17),
    error = identity
  )

  res <- pool_ipw(fits, dfcom = 17)

  expect_s3_class(cnd, "causalgenerics_no_conditional_vcov")
  expect_identical(names(res$alternate), c("effects", "reason"))
  expect_identical(res$alternate$effects, "conditional")
  expect_type(res$alternate$reason, "character")
  expect_length(res$alternate$reason, 1L)
  # The reason is what the same refusal said when the caller asked for that
  # reading, rather than a second wording of it.
  expect_identical(res$alternate$reason, conditionMessage(cnd))
  # The reading that was asked for is pooled exactly as it is when the other
  # one is available: the `rd` row is the 0.30 the literal test works out.
  expect_identical(res$estimates$effect, pool_binary_labels())
  expect_identical(res$estimates$estimate[[1]], 0.30)
})

test_that("pool_ipw() records a reason when only some results carry a block", {
  # An average of the per-imputation coefficients needs a corrected block from
  # every imputation, so a set where one outcome model was never wrapped cannot
  # be pooled on that reading at all. What is recorded is the refusal that one
  # result raises, since it is the whole of why the reading is unavailable.
  fits <- pool_both_surface_fits()
  fits[[2]]$outcome_mod <- pool_outcome_model()
  cnd <- tryCatch(
    pool_ipw(fits, effects = "conditional", dfcom = 18),
    error = identity
  )

  res <- pool_ipw(fits, dfcom = 18)

  expect_s3_class(cnd, "causalgenerics_no_conditional_vcov")
  expect_identical(names(res$alternate), c("effects", "reason"))
  expect_identical(res$alternate$effects, "conditional")
  expect_identical(res$alternate$reason, conditionMessage(cnd))
  # The results that do carry a block are not enough on their own, and the
  # reading the caller asked for is unaffected either way.
  expect_no_error(conditional_vcov(fits[[1]]$outcome_mod))
  expect_identical(res$estimates$effect, pool_binary_labels())
})

test_that("pool_ipw() pools results stored before the mode existed", {
  # A result from an earlier version of a fitting package carries six fields,
  # with no `effects` among them, which reads as marginal, and its outcome
  # model carries no corrected block either. Both halves of that are answered here: the marginal
  # reading is pooled and the conditional one is recorded as unavailable.
  fits <- pool_legacy_fits()

  expect_length(fits[[1]], 6L)
  expect_null(fits[[1]]$effects)

  res <- pool_ipw(fits, dfcom = 17)

  expect_identical(res$effects, "marginal")
  expect_identical(names(res$alternate), c("effects", "reason"))
  expect_identical(res$alternate$effects, "conditional")
  # The absent field is the only difference from the same results carrying it,
  # so the pooled result has to be the one those give, alternate included.
  expect_identical(res, pool_ipw(pool_binary_fits(), dfcom = 17))
})

test_that("pool_ipw() records a level disagreement on the other reading", {
  # The levels the results stored have to agree for the marginal reading to
  # report bounds at one of them, and that disagreement is a refusal when the
  # marginal reading is the one asked for. On the other side it is one more
  # reason the reading could not be pooled, since what the estimates frames
  # stored says nothing about the conditional pooling that was asked for.
  fits <- pool_both_surface_fits(conf.level = c(0.95, 0.9, 0.95))
  cnd <- tryCatch(
    pool_ipw(fits, effects = "marginal", dfcom = 18),
    error = identity
  )

  res <- pool_ipw(fits, effects = "conditional", dfcom = 18)

  expect_s3_class(cnd, "causalgenerics_pool_mismatch_conf_level")
  expect_identical(names(res$alternate), c("effects", "reason"))
  expect_identical(res$alternate$effects, "marginal")
  expect_identical(res$alternate$reason, conditionMessage(cnd))

  # A level named at the call site settles the question for both readings, so
  # the same set has both surfaces again.
  named <- pool_ipw(
    fits,
    effects = "conditional",
    dfcom = 18,
    conf_level = 0.9
  )

  expect_identical(names(named$alternate), c("effects", "estimates", "pooling"))
  expect_identical(named$alternate$estimates$conf.level, rep(0.9, 3))
})

test_that("pool_ipw() appends the alternate after the documented nine", {
  # The nine components the return contract documents keep their names and
  # their order, since callers read fields by name and a printed form writes
  # them positionally. The alternate is a tenth after them rather than a field
  # among them.
  fits <- pool_binary_fits()
  res <- pool_ipw(fits, dfcom = 17)

  expect_length(res, 10L)
  expect_identical(
    names(res)[seq_len(9)],
    c(
      "estimand",
      "estimates",
      "pooling",
      "se_method",
      "effects",
      "m",
      "dfcom",
      "nobs",
      "outcome_link"
    )
  )
  # Everything after those nine, so the assertion is that the alternate is the
  # last field as well as the tenth.
  expect_identical(names(res)[-seq_len(9)], "alternate")
})

test_that("pool_ipw() keys the two surfaces of a categorical result apart", {
  # A categorical exposure reports one row per measure per contrast on the
  # marginal reading and one row per coefficient on the conditional one, so the
  # two surfaces of the same result differ in their rows as well as in their
  # numbers. Each frame is keyed the way its own surface is: the active one
  # carries the contrast column the `ipw()` contract puts after `effect`, and
  # the alternate names coefficients and has no contrast to report.
  fits <- pool_categorical_both_surface_fits()
  labels <- names(coef(pool_categorical_models()[[1]]))

  res <- pool_ipw(fits, dfcom = 17)

  expect_identical(names(res$alternate), c("effects", "estimates", "pooling"))
  expect_identical(nrow(res$estimates), 6L)
  expect_identical(res$estimates$contrast, rep(c("b vs a", "c vs a"), each = 3))

  expect_identical(nrow(res$alternate$estimates), 3L)
  expect_identical(res$alternate$estimates$effect, labels)
  expect_identical(res$alternate$pooling$effect, labels)
  expect_false("contrast" %in% names(res$alternate$estimates))
  expect_false("contrast" %in% names(res$alternate$pooling))
  expect_identical(
    dimnames(attr(res$alternate$estimates, "ipw_vcov", exact = TRUE)),
    list(labels, labels)
  )
})

test_that("pool_ipw() pools the other reading without saying anything", {
  # A caller who asked for one reading is told nothing about the other: a set
  # whose conditional surface cannot be pooled records why rather than warning
  # about it, and the case where both readings are poolable says nothing
  # either.
  available <- pool_both_surface_fits()

  expect_silent(pool_ipw(pool_binary_fits(), dfcom = 17))
  expect_silent(pool_ipw(pool_legacy_fits(), dfcom = 17))
  expect_silent(pool_ipw(available, dfcom = 18))
  expect_silent(pool_ipw(available, effects = "conditional", dfcom = 18))
})

test_that("pool_ipw() warns once when nothing reports a complete-data df", {
  # The fallback that assumes a large sample belongs to the pooling rather than
  # to a reading of it, since the complete-data count is one of the fields both
  # readings share. It is settled once, before either surface is pooled, so
  # pooling the other reading beside the one that was asked for does not raise
  # that warning a second time, which a caller reading two of them would take as
  # two separate assumptions.
  #
  # The option is restored to its default for this block for the reason the
  # fallback's own test above restores it: the suite treats a stray warning as
  # a failure, and the subject here is the warning itself.
  withr::local_options(warn = 0)

  local_s3_method("nobs", "cg_pool_bare", function(object, ...) 20L)
  local_s3_method("coef", "cg_pool_bare", function(object, ...) {
    object$coefficients
  })
  fits <- pool_dfless_fits()

  # What makes the set discriminating: both readings can be pooled, and neither
  # the results nor their outcome models report a count. A count resolved once
  # per surface would warn twice for this set, and a set whose other reading is
  # refused before the count is wanted cannot tell the two apart.
  expect_no_error(conditional_vcov(fits[[1]]$outcome_mod))
  expect_identical(df.residual(fits[[1]]), NA_integer_)
  expect_null(stats::df.residual(fits[[1]]$outcome_mod))

  pooled <- pool_warnings(pool_ipw(fits))

  expect_length(pooled$warnings, 1L)
  expect_match(pooled$warnings, "a large sample is assumed")
  expect_identical(pooled$value$dfcom, Inf)
  # The other reading was pooled rather than passed over, which is what the
  # single warning is asserted against.
  expect_identical(
    names(pooled$value$alternate),
    c("effects", "estimates", "pooling")
  )
  expect_identical(pooled$value$alternate$effects, "conditional")
  expect_identical(
    pooled$value$alternate$estimates$effect,
    c("(Intercept)", "z")
  )

  # A set whose other reading cannot be pooled at all is answered the same way.
  # The count was settled before either surface was read, so a reading that is
  # never pooled neither adds a warning nor takes one away.
  unavailable <- pool_fits_varying(
    "outcome_mod",
    rep(list(pool_bare_model()), 3)
  )
  refused <- pool_warnings(pool_ipw(unavailable))

  expect_length(refused$warnings, 1L)
  expect_identical(names(refused$value$alternate), c("effects", "reason"))
})

# ---- the group column --------------------------------------------------------

# A result may report each effect measure once per level of a grouping variable,
# and the column naming that level is `group`. It is the third component of a
# row's identity, after the effect and the column naming the contrast, and its
# values are written as `"var = value"` strings such as `"sex = 0"`.
#
# Pooling reads that identity twice over. The rows of the per-imputation frames
# are aligned by it, so the `sex = 0` row of one imputation is averaged with the
# `sex = 0` row of the next rather than with whatever sits in the same position;
# and the columns that name a row are carried onto both pooled frames, so a
# pooled result is keyed the way the results it pooled were. A pooling that read
# only the effect would combine the two subgroups of each measure into one row
# and report an average of estimates that answer different questions, which is
# the failure this section guards against and which does not error.
#
# What the results have to agree on is unchanged: the effect labels, as an
# ordered vector. The group is part of a label, so a set whose subgroups differ
# disagrees about its labels and is refused through the condition that
# disagreement already raises, rather than through one of its own.

test_that("pool_ipw() keys a grouped result by effect and group", {
  # Six rows, the group column in the position the identity puts it, and the
  # labels in the order the results reported them. The numbers are the pooled
  # ones for each column of the fixture, so a pooling that aligned the rows any
  # other way would report a different average under the same label.
  fits <- pool_group_fits()
  estimate <- pool_categorical_estimate()
  std_err <- pool_categorical_std_err()
  res <- pool_ipw(fits, dfcom = 17)
  expected <- rubin_rules_by_column(estimate, std_err, dfcom = 17)

  expect_identical(nrow(res$estimates), 6L)
  expect_identical(
    names(res$estimates),
    c(
      "effect",
      "group",
      "estimate",
      "std.err",
      "t",
      "df",
      "ci.lower",
      "ci.upper",
      "conf.level",
      "p.value"
    )
  )
  expect_false("contrast" %in% names(res$estimates))
  expect_identical(
    res$estimates$effect,
    rep(c("rd", "log(rr)", "log(or)"), times = 2)
  )
  expect_identical(res$estimates$group, rep(c("sex = 0", "sex = 1"), each = 3))
  expect_equal(res$estimates$estimate, rubin_column(expected, "estimate"))
  expect_equal(res$estimates$std.err, rubin_column(expected, "std.err"))

  # The diagnostics frame is keyed the same way, so a row of one is found from a
  # row of the other.
  expect_identical(
    names(res$pooling),
    c("effect", "group", "ubar", "b", "riv", "lambda", "fmi")
  )
  expect_identical(res$pooling$effect, res$estimates$effect)
  expect_identical(res$pooling$group, res$estimates$group)
  expect_equal(res$pooling$fmi, rubin_column(expected, "fmi"))
})

test_that("pool_ipw() keys a crossed result by effect, contrast, and group", {
  # All three columns at once, in the order the identity puts them. No two of
  # them name a row here, so a pooled frame missing any one would carry rows a
  # caller cannot tell apart.
  fits <- pool_contrast_group_fits()
  estimate <- pool_contrast_group_estimate()
  std_err <- pool_contrast_group_std_err()
  res <- pool_ipw(fits, dfcom = 17)
  expected <- rubin_rules_by_column(estimate, std_err, dfcom = 17)

  expect_identical(nrow(res$estimates), 8L)
  expect_identical(names(res$estimates)[1:3], c("effect", "contrast", "group"))
  expect_identical(
    names(res$pooling),
    c("effect", "contrast", "group", "ubar", "b", "riv", "lambda", "fmi")
  )

  expect_identical(res$estimates$effect, rep(c("rd", "log(rr)"), times = 4))
  expect_identical(
    res$estimates$contrast,
    rep(rep(c("b vs a", "c vs a"), each = 2), times = 2)
  )
  expect_identical(res$estimates$group, rep(c("sex = 0", "sex = 1"), each = 4))
  expect_equal(res$estimates$estimate, rubin_column(expected, "estimate"))
  expect_equal(res$estimates$std.err, rubin_column(expected, "std.err"))
})

test_that("pool_ipw() names the pooled covariance by the grouped labels", {
  # The dimnames are the effect labels, which carry the group, so a caller who
  # reads an entry out by the name `coef()` gives gets the variance of the row
  # that name belongs to.
  res <- pool_ipw(pool_group_fits(vcov = TRUE), dfcom = 17)
  covariance <- attr(res$estimates, "ipw_vcov", exact = TRUE)

  expect_identical(dim(covariance), c(6L, 6L))
  expect_identical(
    dimnames(covariance),
    list(pool_group_labels(), pool_group_labels())
  )
  expect_equal(diag(covariance), res$estimates$std.err^2, ignore_attr = TRUE)

  crossed <- pool_ipw(pool_contrast_group_fits(vcov = TRUE), dfcom = 17)
  expect_identical(
    dimnames(attr(crossed$estimates, "ipw_vcov", exact = TRUE)),
    list(pool_contrast_group_labels(), pool_contrast_group_labels())
  )
})

test_that("pool_ipw() refuses results whose subgroups disagree", {
  # Two imputations reporting different levels of the grouping variable, which
  # is what a subgroup observed in one imputation and not in another produces.
  # There is no row in the second set to average the `sex = 1` row of the first
  # with, and pooling by position would combine `sex = 1` with `sex = 2`. The
  # group is part of a label, so the refusal is the one a disagreement about the
  # effects reported already raises rather than a second condition beside it.
  fits <- pool_group_fits()
  estimate <- pool_categorical_estimate()
  std_err <- pool_categorical_std_err()
  fits[[2]]$estimates <- pool_group_estimates(
    estimate[2, ],
    std_err[2, ],
    group = c("sex = 0", "sex = 2")
  )

  expect_error(pool_ipw(fits), class = "causalgenerics_pool_mismatch_labels")
  expect_error(pool_ipw(fits), class = "causalgenerics_pool_mismatch")

  cnd <- tryCatch(pool_ipw(fits), error = identity)
  expect_identical(cnd$field, "labels")
  expect_true("log(or) sex = 2" %in% unlist(cnd$values))

  expect_snapshot(error = TRUE, pool_ipw(fits))
})

test_that("pool_ipw() refuses subgroups reported in different orders", {
  # The same two subgroups the other way round. The set of labels matches, so a
  # check that compared sets would take these as poolable and then average the
  # `sex = 0` rows of one imputation with the `sex = 1` rows of another. The
  # labels are what say which row is which, and the pooled frame reports them in
  # an order, so the order is part of what has to agree. Nothing about this set
  # errors on its own: the rows line up by position and the numbers combine, so
  # the wrong answer is one a caller has no way to see.
  fits <- pool_group_fits()
  estimate <- pool_categorical_estimate()
  std_err <- pool_categorical_std_err()
  fits[[3]]$estimates <- pool_group_estimates(
    estimate[3, ],
    std_err[3, ],
    group = c("sex = 1", "sex = 0")
  )

  expect_error(pool_ipw(fits), class = "causalgenerics_pool_mismatch_labels")
  expect_error(pool_ipw(fits), class = "causalgenerics_pool_mismatch")
})

test_that("pool_ipw() refuses a set where only some results name subgroups", {
  # One result reporting its effects for the whole sample and the others
  # reporting them by subgroup. The rows do not correspond at all: three of
  # these results name six and one names three. The row counts already tell
  # them apart, and what the labels have to say is which rows the counts belong
  # to, so the values the condition carries are the grouped labels rather than
  # the bare effects on both sides.
  fits <- pool_group_fits()
  fits[[3]]$estimates <- pool_binary_estimates(
    pool_binary_estimate()[3, ],
    pool_binary_std_err()[3, ]
  )

  expect_error(pool_ipw(fits), class = "causalgenerics_pool_mismatch_labels")
  expect_error(pool_ipw(fits), class = "causalgenerics_pool_mismatch")

  cnd <- tryCatch(pool_ipw(fits), error = identity)
  expect_identical(cnd$field, "labels")
  expect_true("rd" %in% unlist(cnd$values))
  expect_true("rd sex = 0" %in% unlist(cnd$values))
})

test_that("pool_ipw() keys an ungrouped result the way it always did", {
  # The backward-compatible half. Every set a fitting package pools today names
  # no subgroups, and reading a column that is not there must leave both pooled
  # frames keyed by the effect, or by the effect and the contrast, and nothing
  # else.
  binary <- pool_ipw(pool_binary_fits(vcov = TRUE), dfcom = 17)
  categorical <- pool_ipw(pool_categorical_fits(vcov = TRUE), dfcom = 17)

  expect_false("group" %in% names(binary$estimates))
  expect_false("group" %in% names(binary$pooling))
  expect_false("group" %in% names(categorical$estimates))
  expect_false("group" %in% names(categorical$pooling))

  expect_identical(names(binary$estimates)[[1]], "effect")
  expect_identical(
    names(binary$pooling),
    c("effect", "ubar", "b", "riv", "lambda", "fmi")
  )
  expect_identical(names(categorical$estimates)[1:2], c("effect", "contrast"))
  expect_identical(
    dimnames(attr(binary$estimates, "ipw_vcov", exact = TRUE)),
    list(pool_binary_labels(), pool_binary_labels())
  )
  expect_identical(
    dimnames(attr(categorical$estimates, "ipw_vcov", exact = TRUE)),
    list(pool_categorical_labels(), pool_categorical_labels())
  )
})

# ---- the readings the results support ----------------------------------------

# A result records which of the two readings its analysis has, and a fitting
# package names one alone when the other has no meaning for what it estimated:
# an exposure entering the outcome model through several columns has no single
# coefficient to read as the conditional effect, and an analysis reported on the
# coefficient scale alone has no contrast to read as the marginal one. Pooling a
# set of such results has to honor what they declare, since a pooled surface
# built from a reading no result supports would report numbers none of them
# claims.

test_that("pool_ipw() pools the one reading a set of results supports", {
  # With no mode named the stored one decides, and a set supporting one reading
  # records that reading as its mode. What is pooled is therefore the surface
  # the results have, and it is pooled the way it is for results that support
  # both: the declared set says which requests are answerable rather than
  # changing the answer to the ones that are.
  fits <- pool_conditional_fits(readings = "conditional")

  res <- pool_ipw(fits, dfcom = 18)
  both <- pool_ipw(pool_conditional_fits(), dfcom = 18)
  shared <- setdiff(names(res), "alternate")

  expect_identical(res$effects, "conditional")
  expect_identical(res[shared], both[shared])
  # The two results differ in the alternate alone, which is where the readings
  # the results declared show up at all.
  expect_false(identical(res$alternate, both$alternate))
})

test_that("pool_ipw() records no readings field for a one-reading set", {
  # The readings are a property of the results rather than of the pooled object.
  # A pooled result presents one reading and records the other or the reason it
  # has none, which is what a caller asks it about, so the ten components the
  # return contract documents are what it carries either way.
  res <- pool_ipw(pool_conditional_fits(readings = "conditional"), dfcom = 18)

  expect_length(res, 10L)
  expect_identical(
    names(res),
    c(
      "estimand",
      "estimates",
      "pooling",
      "se_method",
      "effects",
      "m",
      "dfcom",
      "nobs",
      "outcome_link",
      "alternate"
    )
  )
  expect_false("readings" %in% names(res))
})

test_that("pool_ipw() refuses a reading the results do not support", {
  # Naming a reading no result has is refused rather than answered from the
  # surface the results happen to carry beside it. The refusal is the one a
  # result raises when the same reading is asked of it directly, so a caller
  # handles one class whether the request went to a result or to the pooling of
  # a set of them.
  fits <- pool_conditional_fits(readings = "conditional")
  cnd <- tryCatch(
    pool_ipw(fits, effects = "marginal", dfcom = 18),
    error = identity
  )
  direct <- tryCatch(
    as.data.frame(fits[[1]], effects = "marginal"),
    error = identity
  )

  expect_error(
    pool_ipw(fits, effects = "marginal", dfcom = 18),
    class = "causalgenerics_unsupported_reading_marginal"
  )
  expect_error(
    pool_ipw(fits, effects = "marginal", dfcom = 18),
    class = "causalgenerics_unsupported_reading"
  )
  # Both facts travel as fields, so a handler reports them without parsing the
  # sentence for them.
  expect_identical(cnd$effects, "marginal")
  expect_identical(cnd$readings, "conditional")
  expect_identical(conditionMessage(cnd), conditionMessage(direct))
  # The reading they do support is unaffected, named or not.
  expect_no_error(pool_ipw(fits, effects = "conditional", dfcom = 18))
  expect_no_error(pool_ipw(fits, dfcom = 18))

  expect_snapshot(
    error = TRUE,
    pool_ipw(fits, effects = "marginal", dfcom = 18)
  )
})

test_that("pool_ipw() refuses a reading one result of the set lacks", {
  # The pooled estimate of an effect is an average over the imputations, so the
  # reading has to be there in every one of them. A set where a single result
  # has no reading of the one asked for is refused for that result, since it is
  # the whole of why the reading cannot be pooled, and the ones that do support
  # it are not enough on their own.
  fits <- pool_mixed_reading_fits()
  cnd <- tryCatch(
    pool_ipw(fits, effects = "marginal", dfcom = 18),
    error = identity
  )

  expect_error(
    pool_ipw(fits, effects = "marginal", dfcom = 18),
    class = "causalgenerics_unsupported_reading_marginal"
  )
  expect_error(
    pool_ipw(fits, effects = "marginal", dfcom = 18),
    class = "causalgenerics_unsupported_reading"
  )
  expect_identical(cnd$readings, "conditional")
  # The fixture is discriminating: two of the three results answer the request
  # the pooling refuses.
  expect_no_error(as.data.frame(fits[[1]], effects = "marginal"))
  expect_no_error(as.data.frame(fits[[3]], effects = "marginal"))
  # The reading every result supports is pooled as it always was.
  expect_no_error(pool_ipw(fits, effects = "conditional", dfcom = 18))

  expect_snapshot(
    error = TRUE,
    pool_ipw(fits, effects = "marginal", dfcom = 18)
  )
})

test_that("pool_ipw() refuses the conditional reading of a marginal-only set", {
  # The mirror of the case above. Neither reading is the one the guard is
  # written for, so a set supporting the marginal reading alone is refused the
  # conditional one in the same words with the reading names the other way
  # round.
  fits <- pool_binary_fits(readings = "marginal")
  cnd <- tryCatch(
    pool_ipw(fits, effects = "conditional", dfcom = 17),
    error = identity
  )

  expect_error(
    pool_ipw(fits, effects = "conditional", dfcom = 17),
    class = "causalgenerics_unsupported_reading_conditional"
  )
  expect_error(
    pool_ipw(fits, effects = "conditional", dfcom = 17),
    class = "causalgenerics_unsupported_reading"
  )
  expect_identical(cnd$effects, "conditional")
  expect_identical(cnd$readings, "marginal")
  expect_no_error(pool_ipw(fits, dfcom = 17))

  expect_snapshot(
    error = TRUE,
    pool_ipw(fits, effects = "conditional", dfcom = 17)
  )
})

test_that("pool_ipw() reads the declared readings before the surfaces", {
  # Two things stand between each of these sets and the reading asked for: the
  # results declare no such reading, and the surface it would be pooled from
  # cannot be read either, because the estimates frames disagree about their
  # labels on one side and the outcome models carry no corrected block on the
  # other. What a caller is told is that the analysis has no such reading, since
  # that question is settled before a frame or an outcome model is read at all.
  mislabeled <- pool_mislabeled_fits(readings = "conditional")
  unwrapped <- pool_binary_fits(readings = "marginal")

  marginal <- tryCatch(
    pool_ipw(mislabeled, effects = "marginal", dfcom = 18),
    error = identity
  )
  conditional <- tryCatch(
    pool_ipw(unwrapped, effects = "conditional", dfcom = 17),
    error = identity
  )

  expect_s3_class(marginal, "causalgenerics_unsupported_reading_marginal")
  expect_false(inherits(marginal, "causalgenerics_pool_mismatch"))
  expect_s3_class(conditional, "causalgenerics_unsupported_reading_conditional")
  expect_false(inherits(conditional, "causalgenerics_no_vcov"))

  # The fixtures are discriminating: the same two sets declaring both readings
  # are refused for the other reason instead, so the assertions above are about
  # which question was asked first rather than about sets with one fault each.
  expect_error(
    pool_ipw(pool_mislabeled_fits(), effects = "marginal", dfcom = 18),
    class = "causalgenerics_pool_mismatch"
  )
  expect_error(
    pool_ipw(pool_binary_fits(), effects = "conditional", dfcom = 17),
    class = "causalgenerics_no_conditional_vcov"
  )
})

test_that("pool_ipw() records an unsupported reading on the alternate", {
  # The reading the caller did not ask for is pooled under a guard, and a
  # reading the results do not support is one more thing that guard catches. The
  # pooled result therefore records why it has no such reading, in the words the
  # refusal used, rather than refusing the call the results could answer.
  fits <- pool_conditional_fits(readings = "conditional")
  cnd <- tryCatch(
    pool_ipw(fits, effects = "marginal", dfcom = 18),
    error = identity
  )
  labels <- names(coef(pool_conditional_models()[[1]]))

  res <- pool_ipw(fits, dfcom = 18)

  expect_identical(names(res$alternate), c("effects", "reason"))
  expect_identical(res$alternate$effects, "marginal")
  expect_type(res$alternate$reason, "character")
  expect_length(res$alternate$reason, 1L)
  expect_identical(res$alternate$reason, conditionMessage(cnd))
  # The reading the results do support is the active one and is unaffected: the
  # coefficient surface is what the frames report, and the marginal labels of
  # the same results are nowhere in them.
  expect_identical(res$estimates$effect, labels)
  expect_false(any(res$estimates$effect %in% pool_binary_labels()))
  # Nothing is said out loud about the reading that could not be pooled.
  expect_silent(pool_ipw(fits, dfcom = 18))
})
