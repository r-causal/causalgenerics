# The exposure-type resolver the r-causal packages share. Whether an exposure
# reads as binary, categorical, or continuous decides which weight formula
# propensity applies, which balance statistic halfmoon computes, and which
# diagnostic positively runs, so the reading has to be the same in all of them.
# Four packages carried a copy of it, and copies drift: two of them counted a
# missing value as a level of its own, which takes a binary exposure with any
# missingness off the binary path entirely. The canonical reading is
# propensity's, and these tests pin it here so the other packages can delete
# their copies and call this one.
#
# The machinery is three layers, and the tests follow them. `observed_values()`,
# `has_two_levels()`, and `is_categorical()` read the data; the middle layer,
# `detect_exposure_type()`, classifies from them and announces what it read; the
# resolution layer, `match_exposure_type()` and `check_forced_type()`, decides
# what a caller's declared type means and refuses the declarations that cannot
# be honored.
#
# Two properties are asserted everywhere rather than once. A missing value is
# not a level, so every count is taken over the observed values. And the
# announcement is output, which makes it part of the contract: the tests pin its
# wording, the argument it names, and both routes that silence it.

# Fixtures. Fixed vectors rather than draws, so a classification never depends
# on the state of the random number generator. `rep_len()` builds the counted
# cases exactly: the categorical numeric takes 15 values over 300 observations,
# a ratio of 0.05, and the continuous one takes 300 values over 300.
cg_binary_numeric <- function() {
  c(0, 1, 1, 0, 1)
}

cg_binary_missing <- function() {
  c(0, 1, NA, 1)
}

cg_one_level_factor <- function() {
  factor(c("treated", "treated", "treated"))
}

cg_two_level_character <- function() {
  c("treated", "untreated", "treated")
}

cg_three_level_factor <- function() {
  factor(c("low", "medium", "high", "low"))
}

cg_categorical_numeric <- function() {
  rep_len(1:15, 300)
}

cg_continuous_numeric <- function() {
  seq_len(300) / 300
}

# The unique-value heuristic on either side of its threshold: 19 values over 100
# observations is 0.19, and 20 over 100 is exactly 0.20.
cg_ratio_below_threshold <- function() {
  rep_len(1:19, 100)
}

cg_ratio_at_threshold <- function() {
  rep_len(1:20, 100)
}

# 19 observed values over 96 non-missing observations, a ratio of 0.198, plus
# two missing values. Counting the missing values as a value of their own gives
# 20 over 96, or 0.208, so the two readings fall on opposite sides of the
# threshold and this vector separates them.
cg_ratio_missing_discriminating <- function() {
  c(rep_len(1:19, 96), NA, NA)
}

# Wrappers whose only job is to be named in a refusal. Each takes the default
# `call`, so the frame the condition reports is the wrapper's rather than the
# resolver's, which is what a consumer package relies on to have its own
# function named.
cg_wrap_match <- function(exposure_type, .exposure, valid_types) {
  match_exposure_type(exposure_type, .exposure, valid_types)
}

cg_wrap_forced <- function(forced, .exposure) {
  check_forced_type(forced, .exposure)
}

# ---- observed_values() -------------------------------------------------

test_that("observed_values() drops missing values and keeps the rest in order", {
  # First-appearance order rather than sorted order, because the value is only
  # ever counted, and sorting it would impose an ordering on a categorical
  # exposure that nothing else in the resolver honors.
  expect_identical(observed_values(cg_binary_missing()), c(0, 1))
  expect_identical(observed_values(c(3, 1, 2, 1)), c(3, 1, 2))
})

test_that("observed_values() returns a bare vector with no attributes", {
  # `unique(x[!is.na(x)])` rather than `unique(stats::na.omit(x))`, which
  # records the dropped positions in an `na.action` attribute that survives
  # `unique()`. The counts are the same either way, but the attribute would
  # travel out to every caller of an exported function.
  expect_null(attributes(observed_values(cg_binary_missing())))
})

test_that("observed_values() keeps the declared levels of a factor", {
  # A factor's declared levels are part of its type, so the observed values of
  # one are a factor over the same levels. Only the values are subset, which is
  # what lets the count separate declared levels from observed ones.
  x <- factor(c("a", "b", NA), levels = c("a", "b", "c"))

  expect_identical(observed_values(x), factor(c("a", "b"), levels = levels(x)))
})

test_that("observed_values() returns nothing when there is nothing observed", {
  expect_identical(observed_values(numeric(0)), numeric(0))
  expect_identical(observed_values(c(NA_real_, NA_real_)), numeric(0))
})

# ---- has_two_levels() --------------------------------------------------

test_that("has_two_levels() counts the observed values", {
  expect_true(has_two_levels(cg_binary_numeric()))
  expect_true(has_two_levels(cg_two_level_character()))
  expect_false(has_two_levels(cg_one_level_factor()))
  expect_false(has_two_levels(cg_three_level_factor()))
  expect_false(has_two_levels(numeric(0)))
})

test_that("has_two_levels() does not count a missing value as a level", {
  # The divergence that made hoisting this worth doing. Counting `NA` as a level
  # reads a binary exposure with any missingness as a three-level one, which
  # takes it off the binary path in every package that consumes the reading.
  expect_true(has_two_levels(cg_binary_missing()))
  expect_false(has_two_levels(c(0, NA, NA)))
})

test_that("has_two_levels() counts observed levels rather than declared ones", {
  # A factor carrying a level nothing takes is still a two-level exposure in
  # the data, and the weight formulas are fitted from the data.
  x <- factor(c("a", "b", "a"), levels = c("a", "b", "c"))

  expect_true(has_two_levels(x))
})

# ---- is_categorical() --------------------------------------------------

test_that("is_categorical() applies the 20 percent unique-value heuristic", {
  expect_true(is_categorical(cg_categorical_numeric()))
  expect_true(is_categorical(cg_ratio_below_threshold()))
  expect_false(is_categorical(cg_continuous_numeric()))
})

test_that("is_categorical() takes the threshold strictly", {
  # A ratio of exactly 0.2 is not categorical. The boundary is pinned so that a
  # later rewrite of the comparison cannot move it without failing here.
  expect_false(is_categorical(cg_ratio_at_threshold()))
})

test_that("is_categorical() counts observed values over non-missing ones", {
  # Both halves of the ratio ignore the missing values, and this vector is the
  # one where that matters: 19 observed values over 96 non-missing observations
  # is below the threshold, and counting the missing values as a value of their
  # own would put it above.
  expect_true(is_categorical(cg_ratio_missing_discriminating()))
})

test_that("is_categorical() is FALSE when nothing is observed", {
  # There is no ratio to take, so the answer is the one that leaves the
  # classification on the continuous branch rather than a `NaN` that travels.
  expect_false(is_categorical(numeric(0)))
  expect_false(is_categorical(c(NA_real_, NA_real_)))
})

# ---- detect_exposure_type(): classification ----------------------------

test_that("detect_exposure_type() reads two observed values as binary", {
  expect_identical(
    detect_exposure_type(cg_binary_numeric(), announce = FALSE),
    "binary"
  )
  expect_identical(
    detect_exposure_type(cg_two_level_character(), announce = FALSE),
    "binary"
  )
  expect_identical(
    detect_exposure_type(c(TRUE, FALSE, TRUE), announce = FALSE),
    "binary"
  )
})

test_that("detect_exposure_type() reads a binary exposure with NAs as binary", {
  expect_identical(
    detect_exposure_type(cg_binary_missing(), announce = FALSE),
    "binary"
  )
})

test_that("detect_exposure_type() reads a one-level factor as binary", {
  # The factor branch asks whether there are more than two observed values, so
  # one value falls to binary rather than to categorical. A single-arm exposure
  # is a degenerate binary exposure, and refusing it belongs to the function
  # that goes on to fit something, which knows what it needs.
  expect_identical(
    detect_exposure_type(cg_one_level_factor(), announce = FALSE),
    "binary"
  )
  expect_identical(
    detect_exposure_type(c("treated", "treated"), announce = FALSE),
    "binary"
  )
})

test_that("detect_exposure_type() counts observed factor levels, not declared", {
  x <- factor(c("a", "b", "a"), levels = c("a", "b", "c"))

  expect_identical(detect_exposure_type(x, announce = FALSE), "binary")
})

test_that("detect_exposure_type() reads more than two levels as categorical", {
  expect_identical(
    detect_exposure_type(cg_three_level_factor(), announce = FALSE),
    "categorical"
  )
  expect_identical(
    detect_exposure_type(c("low", "medium", "high"), announce = FALSE),
    "categorical"
  )
})

test_that("detect_exposure_type() reads a sparse numeric as categorical", {
  expect_identical(
    detect_exposure_type(cg_categorical_numeric(), announce = FALSE),
    "categorical"
  )
  expect_identical(
    detect_exposure_type(cg_ratio_below_threshold(), announce = FALSE),
    "categorical"
  )
})

test_that("detect_exposure_type() reads a varied numeric as continuous", {
  expect_identical(
    detect_exposure_type(cg_continuous_numeric(), announce = FALSE),
    "continuous"
  )
  expect_identical(
    detect_exposure_type(cg_ratio_at_threshold(), announce = FALSE),
    "continuous"
  )
})

test_that("detect_exposure_type() reads a constant numeric as categorical", {
  # Not a case anyone declares, but the branch order decides it and the answer
  # should not change silently. One value is not two, so the binary branch is
  # not taken; the vector is not a factor; and one value over ten observations
  # is below the threshold, so the heuristic claims it.
  expect_identical(
    detect_exposure_type(rep(1, 10), announce = FALSE),
    "categorical"
  )
})

test_that("detect_exposure_type() reads an unobserved exposure as continuous", {
  # Nothing observed means no ratio to take, so `is_categorical()` is `FALSE`
  # and the classification falls through to continuous. Detection classifies
  # whatever it is handed; a function that cannot fit an empty exposure refuses
  # it on its own terms.
  expect_identical(
    detect_exposure_type(numeric(0), announce = FALSE),
    "continuous"
  )
  expect_identical(
    detect_exposure_type(c(NA_real_, NA_real_), announce = FALSE),
    "continuous"
  )
})

# ---- detect_exposure_type(): announcement ------------------------------

test_that("detect_exposure_type() announces the type it read", {
  expect_snapshot(detect_exposure_type(cg_binary_numeric()))
  expect_snapshot(detect_exposure_type(cg_three_level_factor()))
  expect_snapshot(detect_exposure_type(cg_continuous_numeric()))
})

test_that("detect_exposure_type() announces the argument it was told to name", {
  # A package whose exposure argument is not called `.exposure` passes its own
  # name, so the announcement never points at an argument the calling function
  # does not have.
  expect_snapshot(
    detect_exposure_type(cg_binary_numeric(), arg = ".exposures")
  )
})

test_that("the announcement ends without a period", {
  # The wording is shared across the ecosystem, so it is asserted directly as
  # well as through the snapshot: the type is the last thing on the line.
  cnd <- rlang::catch_cnd(
    detect_exposure_type(cg_binary_numeric()),
    classes = "message"
  )

  expect_match(conditionMessage(cnd), "Treating `\\.exposure` as binary\\n$")
})

test_that("announce = FALSE silences the announcement", {
  expect_no_message(
    expect_identical(
      detect_exposure_type(cg_binary_numeric(), announce = FALSE),
      "binary"
    )
  )
})

test_that("the causalgenerics.quiet option silences the announcement", {
  # The option is the user's route and `announce` is the caller's, and either
  # one is enough on its own.
  withr::local_options(causalgenerics.quiet = TRUE)

  expect_no_message(
    expect_identical(detect_exposure_type(cg_binary_numeric()), "binary")
  )
})

test_that("the announcement is made when the quiet option is FALSE", {
  withr::local_options(causalgenerics.quiet = FALSE)

  expect_message(detect_exposure_type(cg_binary_numeric()), "as binary")
})

# ---- match_exposure_type(): resolution ---------------------------------

test_that("match_exposure_type() returns an explicit type without the data", {
  # The strongest statement of the rule that a declaration wins: `.exposure` is
  # not supplied at all, and the call still resolves. An implementation that
  # consulted the data on this path would error on the missing argument.
  expect_no_message(
    expect_identical(match_exposure_type("binary"), "binary")
  )
  expect_no_message(
    expect_identical(match_exposure_type("continuous"), "continuous")
  )
})

test_that("match_exposure_type() does not check an explicit type against data", {
  # Detection would read each of these the other way. The caller's declaration
  # is how the exposure is modeled, and a heuristic over unique values does not
  # get to overrule it.
  expect_no_message(
    expect_identical(
      match_exposure_type("continuous", cg_three_level_factor()),
      "continuous"
    )
  )
  expect_no_message(
    expect_identical(
      match_exposure_type("categorical", cg_continuous_numeric()),
      "categorical"
    )
  )
})

test_that("match_exposure_type() detects from the data under auto", {
  expect_identical(
    match_exposure_type("auto", cg_binary_numeric(), announce = FALSE),
    "binary"
  )
  expect_identical(
    match_exposure_type("auto", cg_three_level_factor(), announce = FALSE),
    "categorical"
  )
  expect_identical(
    match_exposure_type("auto", cg_continuous_numeric(), announce = FALSE),
    "continuous"
  )
})

test_that("match_exposure_type() detects when exposure_type is not supplied", {
  # The default is the full set of choices, which resolves to `"auto"` the way
  # a matched argument does.
  expect_identical(
    match_exposure_type(.exposure = cg_binary_numeric(), announce = FALSE),
    "binary"
  )
})

test_that("match_exposure_type() announces the detected type", {
  expect_snapshot(match_exposure_type("auto", cg_categorical_numeric()))
})

test_that("match_exposure_type() threads the argument name into detection", {
  expect_snapshot(
    match_exposure_type("auto", cg_binary_numeric(), arg = ".exposures")
  )
})

test_that("match_exposure_type() threads announce into detection", {
  expect_no_message(
    expect_identical(
      match_exposure_type("auto", cg_binary_numeric(), announce = FALSE),
      "binary"
    )
  )
})

# ---- match_exposure_type(): refusals -----------------------------------

test_that("match_exposure_type() refuses a string that is not a type", {
  expect_error(
    match_exposure_type("bimary", cg_binary_numeric()),
    class = "rlang_error"
  )
  expect_snapshot(
    error = TRUE,
    match_exposure_type("bimary", cg_binary_numeric())
  )
})

test_that("match_exposure_type() refuses a known type the caller cannot fit", {
  # The two-tier resolution. A type this package knows but this function does
  # not answer is not an unrecognized value, and the refusal says which types
  # the function does answer rather than listing the ones it recognizes.
  expect_error(
    match_exposure_type(
      "continuous",
      cg_binary_numeric(),
      valid_types = c("auto", "binary", "categorical")
    ),
    class = "causalgenerics_unsupported_exposure_type"
  )
  expect_error(
    match_exposure_type(
      "continuous",
      cg_binary_numeric(),
      valid_types = c("auto", "binary", "categorical")
    ),
    class = "causalgenerics_error"
  )
  expect_snapshot(
    error = TRUE,
    match_exposure_type(
      "continuous",
      cg_binary_numeric(),
      valid_types = c("auto", "binary", "categorical")
    )
  )
})

test_that("the unsupported refusal covers a detected type as well", {
  # Naming a type and letting it be detected reach the same refusal, so a
  # caller cannot get further by leaving the type off. Detection still
  # announces what it read before the refusal, which is what makes the refusal
  # legible: the reading is the reason for it. The snapshot below is where that
  # announcement is asserted; the two class assertions silence it, since a
  # message escaping a test is output nobody asked for.
  expect_error(
    suppressMessages(match_exposure_type(
      "auto",
      cg_continuous_numeric(),
      valid_types = c("auto", "binary", "categorical")
    )),
    class = "causalgenerics_unsupported_exposure_type"
  )
  expect_error(
    suppressMessages(match_exposure_type(
      "auto",
      cg_continuous_numeric(),
      valid_types = c("auto", "binary", "categorical")
    )),
    class = "causalgenerics_error"
  )
  expect_snapshot(
    error = TRUE,
    match_exposure_type(
      "auto",
      cg_continuous_numeric(),
      valid_types = c("auto", "binary", "categorical")
    )
  )
})

test_that("the unsupported refusal does not offer auto as a type", {
  # `"auto"` is how a caller asks for detection, not an answer detection can
  # give, so the supported set names the three readable types only.
  expect_snapshot(
    error = TRUE,
    match_exposure_type(
      "categorical",
      cg_binary_numeric(),
      valid_types = c("auto", "binary")
    )
  )
})

test_that("both refusals name the function the caller called", {
  # The condition reports the wrapper's call rather than
  # `match_exposure_type()`'s, which is what a consumer package needs: the user
  # called its function, and never called this one.
  expect_snapshot(
    error = TRUE,
    cg_wrap_match("continuous", cg_binary_numeric(), c("auto", "binary"))
  )
  expect_snapshot(
    error = TRUE,
    cg_wrap_match("bimary", cg_binary_numeric(), c("auto", "binary"))
  )
})

test_that("a caller narrowing valid_types narrows its own default too", {
  # Leaving `exposure_type` at the full set of choices while supporting fewer
  # of them is a mismatch in the calling function's signature, and it is
  # reported as one rather than silently read as `"auto"`. A consumer keeps the
  # two in step, as propensity's weight functions do.
  expect_snapshot(
    error = TRUE,
    match_exposure_type(
      .exposure = cg_binary_numeric(),
      valid_types = c("auto", "binary")
    )
  )
})

# ---- check_forced_type() -----------------------------------------------

test_that("check_forced_type() passes a declaration the data can carry", {
  expect_invisible(check_forced_type("binary", cg_binary_numeric()))
  expect_null(check_forced_type("binary", cg_binary_numeric()))
  expect_null(check_forced_type("continuous", cg_continuous_numeric()))
})

test_that("check_forced_type() counts observed values for a binary claim", {
  # The same missing-value rule the rest of the machinery uses, applied to the
  # structural check: a binary exposure with missing values is still binary.
  expect_null(check_forced_type("binary", cg_binary_missing()))
})

test_that("check_forced_type() always passes a categorical claim", {
  # Every vector can be read as the set of levels it takes, so no data
  # contradicts a categorical declaration.
  expect_null(check_forced_type("categorical", cg_continuous_numeric()))
  expect_null(check_forced_type("categorical", cg_one_level_factor()))
  expect_null(check_forced_type("categorical", cg_binary_numeric()))
})

test_that("check_forced_type() passes a continuous claim over a sparse numeric", {
  # Detection reads this one as categorical, and the declaration still stands:
  # a dose taking fifteen values is a dose. Only data that cannot carry the
  # declaration at all is refused.
  expect_null(check_forced_type("continuous", cg_categorical_numeric()))
})

test_that("check_forced_type() refuses a binary claim over more than two levels", {
  expect_error(
    check_forced_type("binary", cg_three_level_factor()),
    class = "causalgenerics_forced_exposure_type"
  )
  expect_error(
    check_forced_type("binary", cg_three_level_factor()),
    class = "causalgenerics_error"
  )
  expect_snapshot(
    error = TRUE,
    check_forced_type("binary", cg_three_level_factor())
  )
})

test_that("check_forced_type() refuses a continuous claim over a factor", {
  expect_error(
    check_forced_type("continuous", cg_three_level_factor()),
    class = "causalgenerics_forced_exposure_type"
  )
  expect_error(
    check_forced_type("continuous", cg_three_level_factor()),
    class = "causalgenerics_error"
  )
  expect_snapshot(
    error = TRUE,
    check_forced_type("continuous", cg_two_level_character())
  )
})

test_that("check_forced_type() reports what detection would have read", {
  # The refusal names the route out of it, and the route is dropping the
  # declaration, so it says what dropping it would give. The reading is
  # computed for the message and never announced: the caller declared a type,
  # so an alert reporting a different one would read as the answer rather than
  # as the alternative.
  expect_snapshot(error = TRUE, check_forced_type("binary", rep(1, 10)))
})

test_that("check_forced_type() announces nothing on either path", {
  expect_no_message(check_forced_type("categorical", cg_binary_numeric()))
  expect_no_message(
    expect_error(check_forced_type("binary", cg_three_level_factor()))
  )
})

test_that("check_forced_type() names the argument it was told to name", {
  expect_snapshot(
    error = TRUE,
    check_forced_type("binary", cg_three_level_factor(), arg = ".exposures")
  )
})

test_that("check_forced_type() names the function the caller called", {
  expect_snapshot(
    error = TRUE,
    cg_wrap_forced("binary", cg_three_level_factor())
  )
})

test_that("check_forced_type() refuses a type it does not know", {
  # An exported check answers for the three readable types. `"auto"` is a
  # request for detection rather than a declaration, and anything else is a
  # typo; either would otherwise pass silently and leave the caller believing a
  # structural check had run.
  expect_error(check_forced_type("auto", cg_binary_numeric()))
  expect_snapshot(error = TRUE, check_forced_type("auto", cg_binary_numeric()))
  expect_snapshot(
    error = TRUE,
    check_forced_type("bimary", cg_binary_numeric())
  )
})
