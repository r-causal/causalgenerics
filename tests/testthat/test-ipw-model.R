# `new_ipw_model()` wraps a component model of an IPW fit so that the model
# carries the covariance the two-step estimation implies rather than the one it
# computed for itself. A weighted outcome model fitted as if its weights were
# fixed understates its own uncertainty, because the weights were estimated from
# the same data; the corrected block comes out of the stacked estimating
# equations the fitting package solved.
#
# The wrapper is a class prepended to the model's own class vector and an
# attribute alongside it, so everything else about the model keeps working
# through inheritance. That is the whole design: a caller writes
# `vcov(fit$outcome_mod)` and gets the corrected matrix, while `predict()`,
# `coef()`, and `model.frame()` behave as they did before.

# ---- fixtures ----------------------------------------------------------------

# Fixed rather than simulated data, so the fits never move between runs.
ipw_model_data <- function() {
  data.frame(
    x = rep(c(-1.5, -0.5, 0.5, 1.5), each = 5),
    z = rep(c(0, 1), 10),
    y = rep(c(0, 1, 1, 0, 1), 4)
  )
}

ipw_model_weights <- function() {
  setNames(rep(c(0.5, 1.5, 1, 2), each = 5), paste0("u", 1:20))
}

outcome_model <- function() {
  glm(
    y ~ z,
    family = quasibinomial(),
    data = ipw_model_data(),
    weights = ipw_model_weights()
  )
}

# A covariance in the shape a fitting package hands over: the outcome block of
# the stacked sandwich, whose dimnames are the model's coefficient names. It is
# built from the model's own covariance so that the dimnames agree without being
# restated, and scaled so that it differs from it. Accounting for the weights
# having been estimated ordinarily widens the interval, and the scaling is what
# makes "the wrapped model reports the corrected covariance" separable from "the
# wrapped model reports the covariance it always did".
corrected_vcov <- function(mod) {
  stats::vcov(mod) * 1.4
}

# ---- new_ipw_model() ---------------------------------------------------------

test_that("new_ipw_model() prepends the class and attaches the covariance", {
  mod <- outcome_model()
  corrected <- corrected_vcov(mod)

  wrapped <- new_ipw_model(mod, corrected)

  # Prepended rather than replaced. The model's own classes have to stay, and
  # stay in order, or every method it inherits is lost.
  expect_identical(class(wrapped), c("ipw_model", "glm", "lm"))
  expect_identical(attr(wrapped, "ipw_vcov"), corrected)
})

test_that("new_ipw_model() wraps a model whose class vector has one element", {
  mod <- lm(y ~ z, data = ipw_model_data())

  wrapped <- new_ipw_model(mod, corrected_vcov(mod))

  expect_identical(class(wrapped), c("ipw_model", "lm"))
})

test_that("new_ipw_model() leaves the rest of the model alone", {
  # The wrapper is exactly two changes. Taking both of them back has to give the
  # model that went in, which is a stronger claim than checking a few components
  # by name: a constructor that dropped a component, reordered the list, or
  # rebuilt the fit would fail here.
  mod <- outcome_model()

  wrapped <- new_ipw_model(mod, corrected_vcov(mod))

  stripped <- wrapped
  attr(stripped, "ipw_vcov") <- NULL
  class(stripped) <- class(mod)

  expect_identical(stripped, mod)
})

test_that("new_ipw_model() does not check the covariance against the model", {
  # A `new_*()` constructor is the low-level developer interface and validates
  # only what its own contract names: a square numeric matrix with dimnames. The
  # fitting package is the one that knows which block of the sandwich belongs to
  # which model, and it labels the block itself. Checking the dimnames against
  # `coef()` here would reject a legitimate wrapper whose labels are the
  # estimating-equation parameter names rather than the coefficient names.
  mod <- outcome_model()
  unrelated <- matrix(
    c(1, 0, 0, 1),
    nrow = 2,
    dimnames = list(c("theta1", "theta2"), c("theta1", "theta2"))
  )

  expect_no_error(new_ipw_model(mod, unrelated))
  expect_identical(vcov(new_ipw_model(mod, unrelated)), unrelated)
})

# ---- vcov.ipw_model() --------------------------------------------------------

test_that("vcov() on a wrapped model returns the corrected covariance", {
  mod <- outcome_model()
  corrected <- corrected_vcov(mod)

  wrapped <- new_ipw_model(mod, corrected)

  expect_identical(vcov(wrapped), corrected)
  # The fixture is discriminating: the model's own covariance is not the
  # corrected one, so the assertion above cannot pass through inheritance.
  expect_false(identical(corrected, stats::vcov(mod)))
  expect_false(identical(vcov(wrapped), stats::vcov(mod)))
})

test_that("vcov() errors when a wrapped model carries no covariance", {
  # The class and the attribute go on together, and `new_ipw_model()` cannot
  # produce one without the other: it validates the matrix before it prepends
  # the class. An object carrying the class with nothing behind it came from
  # somewhere else, such as a class vector assembled by hand, or a model that
  # passed through code which dropped its attributes and kept its class, and the
  # method is where that is found out. Handing the attribute back as it was
  # found answers `NULL`, and `NULL` travels: the standard errors taken from it
  # are an empty vector, and the limits built from those come back as `NA` with
  # nothing said about why.
  mod <- outcome_model()
  stripped <- new_ipw_model(mod, corrected_vcov(mod))
  attr(stripped, "ipw_vcov") <- NULL

  # The class is still there, so dispatch reaches this method, which is what
  # makes the case worth guarding: the method runs and has nothing to report.
  expect_identical(class(stripped), c("ipw_model", "glm", "lm"))

  expect_error(vcov(stripped), class = "causalgenerics_no_vcov_ipw_model")
  expect_error(vcov(stripped), class = "causalgenerics_no_vcov")

  # Keyed to the component model rather than to the result, so a handler can
  # tell a model that carries no block from a result that records none.
  cnd <- tryCatch(vcov(stripped), error = identity)
  expect_identical(cnd$result, "ipw_model")

  # The fixture is discriminating: the model's own covariance is reachable
  # through inheritance, and it is the answer this class exists to replace, so
  # a method that fell back rather than refused would have one to give.
  expect_no_error(stats::vcov(mod))

  expect_snapshot(error = TRUE, vcov(stripped))
})

test_that("a wrapped model still works through inheritance", {
  # Nothing is registered for `ipw_model` beyond `vcov()`, so every other
  # generic walks past it to the model's own methods. This is what the class
  # being prepended rather than replacing buys, and it is the reason the
  # corrected covariance can travel with the model instead of beside it.
  mod <- outcome_model()

  wrapped <- new_ipw_model(mod, corrected_vcov(mod))

  expect_identical(coef(wrapped), coef(mod))
  expect_identical(predict(wrapped), predict(mod))
  expect_identical(fitted(wrapped), fitted(mod))
  expect_identical(model.frame(wrapped), model.frame(mod))
  expect_identical(formula(wrapped), formula(mod))
  expect_identical(nobs(wrapped), nobs(mod))
})

test_that("an ipw result carrying a wrapped outcome model keeps its accessors", {
  # The wrapper goes on at the fitting package's construction site, so the
  # accessors on the result meet a wrapped model rather than a bare one. They
  # read the model through the same generics as before, and the corrected
  # covariance is reachable through the component.
  mod <- outcome_model()
  corrected <- corrected_vcov(mod)
  dat <- ipw_model_data()

  res <- new_ipw(
    estimand = "ate",
    wt_mod = glm(z ~ x, family = binomial(), data = dat),
    outcome_mod = new_ipw_model(mod, corrected),
    estimates = data.frame(
      effect = "diff",
      estimate = 2.25255,
      std.err = 0.17524,
      z = 12.854,
      ci.lower = 1.909083,
      ci.upper = 2.596017,
      conf.level = 0.95,
      p.value = 4.5e-38
    ),
    se_method = "mestimation",
    fit = NULL
  )

  expect_identical(nobs(res), 20L)
  expect_identical(weights(res), ipw_model_weights())
  expect_identical(vcov(res$outcome_mod), corrected)
})

# ---- validation --------------------------------------------------------------

test_that("new_ipw_model() rejects a covariance that is not a matrix", {
  # A variance estimator that returned a vector, or a package that handed over a
  # data frame of the same numbers, would otherwise be stored and reported as a
  # covariance matrix by `vcov()`, and the caller would find out at whatever
  # matrix algebra ran next.
  mod <- outcome_model()

  expect_error(
    new_ipw_model(mod, c(theta1 = 1, theta2 = 1)),
    class = "causalgenerics_invalid_argument_vcov"
  )
  expect_error(
    new_ipw_model(mod, as.data.frame(corrected_vcov(mod))),
    class = "causalgenerics_invalid_argument_vcov"
  )
  expect_error(
    new_ipw_model(mod, c(theta1 = 1, theta2 = 1)),
    class = "causalgenerics_invalid_argument"
  )
})

test_that("new_ipw_model() rejects a covariance that is not numeric", {
  mod <- outcome_model()
  labels <- list(c("a", "b"), c("a", "b"))

  expect_error(
    new_ipw_model(mod, matrix(letters[1:4], nrow = 2, dimnames = labels)),
    class = "causalgenerics_invalid_argument_vcov"
  )
})

test_that("new_ipw_model() rejects a covariance that is not square", {
  # A covariance matrix has one row and one column per parameter. A rectangular
  # block is some other quantity, most likely a cross-covariance between two
  # blocks of the stacked system, and reporting it as `vcov()` would be wrong in
  # a way that no later call announces.
  mod <- outcome_model()

  expect_error(
    new_ipw_model(
      mod,
      matrix(
        c(1, 0, 0, 1, 0, 0),
        nrow = 2,
        dimnames = list(c("a", "b"), c("a", "b", "c"))
      )
    ),
    class = "causalgenerics_invalid_argument_vcov"
  )
})

test_that("new_ipw_model() rejects a covariance with no dimnames", {
  # The labels are what make the matrix readable: a caller pulls a variance out
  # by coefficient name, and a summary lines the diagonal up with `coef()` by
  # name rather than by position. Both margins are needed, because a matrix
  # labelled down one side only cannot be indexed by name in both directions
  # and, on a square matrix, gives no way to notice a transposed block.
  mod <- outcome_model()

  expect_error(
    new_ipw_model(mod, matrix(c(1, 0, 0, 1), nrow = 2)),
    class = "causalgenerics_invalid_argument_vcov"
  )
  expect_error(
    new_ipw_model(
      mod,
      matrix(c(1, 0, 0, 1), nrow = 2, dimnames = list(c("a", "b"), NULL))
    ),
    class = "causalgenerics_invalid_argument_vcov"
  )
})

test_that("the covariance error states the contract", {
  mod <- outcome_model()

  expect_snapshot(error = TRUE, new_ipw_model(mod, c(theta1 = 1, theta2 = 1)))
  expect_snapshot(
    error = TRUE,
    new_ipw_model(mod, matrix(c(1, 0, 0, 1), nrow = 2))
  )
})

# ---- registration and export -------------------------------------------------

test_that("vcov.ipw_model is registered and new_ipw_model() is exported", {
  # `vcov()` is a stats generic, so the method entry belongs to stats' table
  # rather than to this package's. Every assertion above reaches the method from
  # the test frame, which `UseMethod()` searches before the table, so table
  # membership is the only claim here about the NAMESPACE directive downstream
  # packages depend on.
  #
  # The constructor is exported because the fitting packages are the ones that
  # call it. It is the second half of the contract: this package owns the class
  # and the method, and they own the covariance that goes into it.
  table <- s3_methods_table("vcov")

  expect_true(!is.null(table))
  expect_true(exists("vcov.ipw_model", envir = table, inherits = FALSE))
  expect_true("new_ipw_model" %in% getNamespaceExports("causalgenerics"))
})
