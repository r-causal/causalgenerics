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
    effects = effects
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

pool_conditional_fits <- function(effects = "conditional") {
  estimate <- pool_binary_estimate()
  std_err <- pool_binary_std_err()
  Map(
    function(i, mod) {
      pool_fit(
        pool_binary_estimates(estimate[i, ], std_err[i, ]),
        outcome_mod = mod,
        effects = effects
      )
    },
    1:3,
    pool_conditional_models()
  )
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

# ---- what `fits` may be ------------------------------------------------------

test_that("pool_ipw() refuses a `fits` argument that is not a list", {
  # The argument is the set of results one analysis produced across the
  # imputations, and there is no reading of a number or a string as one.
  expect_error(pool_ipw(1:3), class = "causalgenerics_invalid_argument_fits")
  expect_error(pool_ipw(1:3), class = "causalgenerics_invalid_argument")
  expect_error(pool_ipw("a"), class = "causalgenerics_invalid_argument_fits")
  expect_error(pool_ipw(NULL), class = "causalgenerics_invalid_argument_fits")
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
  # would never reach this, and one that let the stored field colour the result
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
})

test_that("pool_ipw() refuses outcome models fitted on different links", {
  # The effects are on the scale the outcome model's link puts them on, so a
  # `log(or)` row from a logit fit and a `log(rr)` row from a log one are not
  # the same quantity even when they are labelled the same way. The link is
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
  local_s3_method("nobs", "cg_pool_bare", function(object, ...) 20L)
  bare <- structure(list(), class = "cg_pool_bare")
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
      "outcome_link"
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
})
