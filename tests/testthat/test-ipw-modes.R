# `as_marginal()` and `as_conditional()` set the presentation mode an `ipw`
# result reports its effects in. The mode is the `effects` field of the
# `new_ipw()` contract, and these generics are the supported way to move a
# result between the two readings: a caller writes `as_conditional(res)` rather
# than reaching into the field and setting it.
#
# The generics live here for the reason `print.ipw()` does. Two packages each
# registering `as_conditional.ipw()` would collide in the shared S3 method
# table, and a caller writing against a result would then get whichever package
# was installed last rather than the contract.
#
# The `ipw` methods are total. Every result has one of the two modes, so asking
# for either is always answerable: the methods never error, they are idempotent,
# and they round-trip. They set one field and read nothing else, which is what
# the assertions below say by comparing the whole object rather than a few
# components of it.
#
# A result produced before the field existed carries six fields rather than
# seven. `ipw_effects()` reads that shape as marginal, so an older result and a
# result built today behave the same way.
#
# The estimates frames are written out literally, in the shape the `ipw()`
# return contract documents, rather than produced by fitting a model. Nothing
# here needs a real estimator, and this package holds no model-fitting
# dependency to build one with.

# ---- fixtures ----------------------------------------------------------------

# A binary-exposure estimates frame: one row per effect measure, no `contrast`
# column. It holds the same numbers as its counterpart in `test-ipw-result.R`.
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

# The covariance of those effects, in the shape a fitting package attaches it.
# It rides along on the `estimates` frame as an attribute, which is the part of
# the result a flip is most likely to lose: rebuilding the frame rather than
# carrying it over drops the attribute and leaves `vcov()` with nothing.
effects_vcov <- function() {
  std_err <- binary_estimates()$std.err
  labels <- c("rd", "log(rr)", "log(or)")
  index <- seq_along(std_err)
  covariance <- 0.9^abs(outer(index, index, "-")) * outer(std_err, std_err)
  dimnames(covariance) <- list(labels, labels)
  covariance
}

# Fixed rather than simulated data, so the fits never move between runs.
ipw_data <- function() {
  data.frame(
    x = rep(c(-1.5, -0.5, 0.5, 1.5), each = 5),
    z = rep(c(0, 1), 10),
    y = rep(c(0, 1, 1, 0, 1), 4)
  )
}

ipw_models <- function() {
  dat <- ipw_data()
  list(
    wt_mod = glm(z ~ x, family = binomial(), data = dat),
    outcome_mod = glm(y ~ z, family = quasibinomial(), data = dat)
  )
}

# A result in the mode a method that names none produces. The mode is left at
# its default here rather than passed, because the marginal result these tests
# flip is the one every method downstream builds today.
marginal_result <- function(estimates = binary_estimates()) {
  mods <- ipw_models()
  new_ipw(
    estimand = "ate",
    wt_mod = mods$wt_mod,
    outcome_mod = mods$outcome_mod,
    estimates = estimates,
    se_method = "mestimation",
    fit = NULL
  )
}

# A result a method built in the conditional mode, rather than one flipped into
# it. A test of `as_marginal()` that started from `as_conditional()` would hold
# whatever the pair did to each other and say nothing about a stored mode.
conditional_result <- function(estimates = binary_estimates()) {
  mods <- ipw_models()
  new_ipw(
    estimand = "ate",
    wt_mod = mods$wt_mod,
    outcome_mod = mods$outcome_mod,
    estimates = estimates,
    se_method = "mestimation",
    fit = NULL,
    effects = "conditional"
  )
}

# A result with the covariance attached, which is where the `new_ipw()` contract
# puts it.
result_with_vcov <- function() {
  estimates <- binary_estimates()
  attr(estimates, "ipw_vcov") <- effects_vcov()
  marginal_result(estimates)
}

# A result of the shape the constructor built before the mode existed: six
# fields and no `effects`. Results stored from an earlier version of a fitting
# package have this shape, as does one built by hand against the earlier
# contract.
legacy_result <- function() {
  mods <- ipw_models()
  structure(
    list(
      estimand = "ate",
      wt_mod = mods$wt_mod,
      outcome_mod = mods$outcome_mod,
      estimates = binary_estimates(),
      se_method = "mestimation",
      fit = NULL
    ),
    class = "ipw"
  )
}

# The other shape an absent mode takes: the field is there and holds `NULL`,
# which is what building the list with `effects = NULL` gives. A reader that
# tested for the name rather than the value would take this one for a mode.
null_effects_result <- function() {
  legacy <- legacy_result()
  fields <- unclass(legacy)
  fields["effects"] <- list(NULL)
  structure(fields, class = "ipw")
}

# ---- as_conditional() --------------------------------------------------------

test_that("as_conditional() sets the mode on a marginal result", {
  res <- marginal_result()

  conditional <- as_conditional(res)

  expect_s3_class(conditional, "ipw", exact = TRUE)
  expect_identical(conditional$effects, "conditional")
  expect_identical(names(conditional), names(res))
})

test_that("as_conditional() leaves every other field alone", {
  # The method sets one field. Setting that field back has to give the object
  # that went in, which is a stronger claim than comparing a few components by
  # name: a method that rebuilt the result, reordered the list, or dropped the
  # covariance attached to the estimates would fail here and nowhere else.
  res <- result_with_vcov()

  conditional <- as_conditional(res)

  restored <- conditional
  restored$effects <- res$effects

  expect_identical(restored, res)
  expect_identical(
    attr(conditional$estimates, "ipw_vcov", exact = TRUE),
    effects_vcov()
  )
})

test_that("as_conditional() leaves a conditional result as it is", {
  # Asking a result for the mode it already reports is a legitimate call, not a
  # mistake, so it answers with the result rather than refusing.
  res <- conditional_result()

  expect_no_error(as_conditional(res))
  expect_identical(as_conditional(res), res)
})

# ---- as_marginal() -----------------------------------------------------------

test_that("as_marginal() sets the mode on a conditional result", {
  res <- conditional_result()

  marginal <- as_marginal(res)

  expect_s3_class(marginal, "ipw", exact = TRUE)
  expect_identical(marginal$effects, "marginal")
  expect_identical(names(marginal), names(res))
})

test_that("as_marginal() leaves every other field alone", {
  estimates <- binary_estimates()
  attr(estimates, "ipw_vcov") <- effects_vcov()
  res <- conditional_result(estimates)

  marginal <- as_marginal(res)

  restored <- marginal
  restored$effects <- res$effects

  expect_identical(restored, res)
  expect_identical(
    attr(marginal$estimates, "ipw_vcov", exact = TRUE),
    effects_vcov()
  )
})

test_that("as_marginal() leaves a marginal result as it is", {
  res <- marginal_result()

  expect_no_error(as_marginal(res))
  expect_identical(as_marginal(res), res)
})

# ---- idempotence and the round trip ------------------------------------------

test_that("both modes are idempotent", {
  # The mode is a state rather than a toggle, so applying either generic twice
  # says the same thing as applying it once.
  conditional <- as_conditional(marginal_result())
  marginal <- as_marginal(conditional_result())

  expect_identical(as_conditional(conditional), conditional)
  expect_identical(as_conditional(as_conditional(conditional)), conditional)
  expect_identical(as_marginal(marginal), marginal)
  expect_identical(as_marginal(as_marginal(marginal)), marginal)
})

test_that("the two modes round-trip", {
  # A caller who looks at a result the other way round and then goes back is
  # holding what they started with, whichever mode they started in. Nothing is
  # dropped on the way out and nothing is invented on the way back.
  marginal <- marginal_result()
  conditional <- conditional_result()

  expect_identical(as_marginal(as_conditional(marginal)), marginal)
  expect_identical(as_conditional(as_marginal(conditional)), conditional)
})

test_that("the mode methods read no field but the mode", {
  # A result whose components are the awkward ones: a weighting object that
  # cannot be subset, an outcome model that records no call, and an estimates
  # frame with no rows. None of them is anything a flip has business reading,
  # and a method that validated the result rather than setting a field would
  # fail on one of them.
  res <- new_ipw(
    estimand = "ate",
    wt_mod = structure(1:3, class = "cg_unsubsettable"),
    outcome_mod = structure(list(), class = "cg_no_call"),
    estimates = binary_estimates()[0, ],
    se_method = "linearization",
    fit = NULL
  )

  expect_no_error(as_conditional(res))
  expect_no_error(as_marginal(as_conditional(res)))
  expect_identical(as_conditional(res)$effects, "conditional")
  expect_identical(as_marginal(as_conditional(res)), res)
})

# ---- results that predate the mode -------------------------------------------

test_that("ipw_effects() reads the stored mode", {
  expect_identical(ipw_effects(marginal_result()), "marginal")
  expect_identical(ipw_effects(conditional_result()), "conditional")
})

test_that("ipw_effects() reads a result with no mode as marginal", {
  # Both shapes an absent mode takes. A result stored before the field existed
  # has six fields; a list built with `effects = NULL` has seven, one of them
  # empty. Neither records a mode, and the marginal reading is the one every
  # method produced when neither did.
  legacy <- legacy_result()
  null_effects <- null_effects_result()

  expect_length(legacy, 6L)
  expect_null(legacy$effects)
  expect_identical(ipw_effects(legacy), "marginal")

  expect_length(null_effects, 7L)
  expect_true("effects" %in% names(null_effects))
  expect_null(null_effects$effects)
  expect_identical(ipw_effects(null_effects), "marginal")
})

test_that("as_conditional() gives a result with no mode the field", {
  # The flip is from the reading the older result already had, so a result that
  # records no mode is flipped as the marginal one it is. Every field it does
  # carry survives the trip.
  legacy <- legacy_result()

  conditional <- as_conditional(legacy)

  expect_s3_class(conditional, "ipw", exact = TRUE)
  expect_identical(conditional$effects, "conditional")
  expect_identical(ipw_effects(conditional), "conditional")
  expect_identical(unclass(conditional)[names(legacy)], unclass(legacy))
})

test_that("as_marginal() leaves a result with no mode reading as marginal", {
  # The older result already reads as marginal, so asking for that reading
  # cannot change what it reports, whether or not the field is written out.
  legacy <- legacy_result()

  marginal <- as_marginal(legacy)

  expect_s3_class(marginal, "ipw", exact = TRUE)
  expect_identical(ipw_effects(marginal), "marginal")
  expect_identical(unclass(marginal)[names(legacy)], unclass(legacy))
})

# ---- the default methods -----------------------------------------------------

test_that("the mode generics have no default", {
  # There is no marginal or conditional reading of an object that is not an IPW
  # result, so the default refuses rather than inventing one. The condition
  # carries the classes `stop_no_method()` builds, specific first, which is what
  # a downstream handler matches on.
  x <- structure(list(), class = "cg_unregistered")

  expect_error(
    dispatch_from_baseenv(as_marginal, x),
    class = "causalgenerics_no_method_as_marginal"
  )
  expect_error(
    dispatch_from_baseenv(as_marginal, x),
    class = "causalgenerics_no_method"
  )
  expect_error(
    dispatch_from_baseenv(as_conditional, x),
    class = "causalgenerics_no_method_as_conditional"
  )
  expect_error(
    dispatch_from_baseenv(as_conditional, x),
    class = "causalgenerics_no_method"
  )
})

test_that("the mode generics refuse an object with no class of its own", {
  # A bare vector reaches the default the same way, which is the call a caller
  # who passed the wrong argument makes.
  expect_error(
    dispatch_from_baseenv(as_marginal, 1:10),
    class = "causalgenerics_no_method"
  )
  expect_error(
    dispatch_from_baseenv(as_conditional, 1:10),
    class = "causalgenerics_no_method"
  )
})

test_that("the no-method error names the generic and the class", {
  x <- structure(list(), class = "cg_unregistered")

  expect_snapshot(error = TRUE, as_marginal(x))
  expect_snapshot(error = TRUE, as_conditional(x))
})

# ---- registration and export -------------------------------------------------

test_that("the mode generics dispatch through the method table", {
  # Every assertion above reaches these methods from the test frame, and
  # `UseMethod()` searches the calling frame before it consults the method
  # table. All of them therefore pass whether or not the NAMESPACE carries an
  # `S3method()` directive. `baseenv()` cannot see the test frame, so a call
  # made from there reaches the registered method or nothing at all, which is
  # the route a downstream package's own namespace takes.
  res <- marginal_result()

  expect_identical(
    dispatch_from_baseenv(as_conditional, res)$effects,
    "conditional"
  )
  expect_identical(
    dispatch_from_baseenv(as_marginal, conditional_result())$effects,
    "marginal"
  )
})

test_that("the mode generics and their methods are registered", {
  # Both generics belong to this package, so their entries are in this package's
  # own method table rather than base's or stats'.
  generics <- c("as_marginal", "as_conditional")

  registered <- vapply(
    generics,
    function(generic) {
      table <- s3_methods_table(generic)
      !is.null(table) &&
        exists(paste0(generic, ".ipw"), envir = table, inherits = FALSE) &&
        exists(paste0(generic, ".default"), envir = table, inherits = FALSE)
    },
    logical(1)
  )

  expect_identical(generics[!registered], character())
})

test_that("the mode generics are exported and ipw_effects() is not", {
  # The generics are the supported way to move a result between the two
  # readings, so a caller and a fitting package both need them. `ipw_effects()`
  # is how the methods here read a field that an older result may not carry, and
  # it stays internal.
  #
  # A namespace's parent chain runs out to the search path, where
  # `pkgload::load_all()` has attached every object in the package, so
  # `inherits = FALSE` is what makes the second assertion about the namespace
  # itself rather than about the attached copy.
  exports <- getNamespaceExports("causalgenerics")

  expect_true("as_marginal" %in% exports)
  expect_true("as_conditional" %in% exports)

  expect_true(
    exists(
      "ipw_effects",
      envir = asNamespace("causalgenerics"),
      inherits = FALSE
    )
  )
  expect_false("ipw_effects" %in% exports)
})
