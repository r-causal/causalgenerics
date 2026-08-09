# `joint_exposure()` is a data-level declaration: it says that one exposure
# variable is the crossing of two treatments, and it carries which two, in what
# order, and which cell is the reference. That declaration is what turns a
# discrete two-treatment interaction into an ordinary categorical-exposure
# problem downstream, so every package that touches the analysis has to read it
# the same way. The package that builds the weights, the package that checks
# balance over the cells, and the package that reports the effects all need the
# same answer to which cells exist and which one is the reference, and none of
# them can own the declaration without the other two depending on it. The
# shared-vocabulary package is where it belongs.
#
# The class is a vctrs vector over integer codes built the way
# `new_causal_wts()` builds a weight vector, with one addition: `"factor"` is in
# the class vector, ahead of `"vctrs_vctr"`, so that `model.matrix()` and the
# rest of the formula machinery treat a joint exposure as the factor it is
# rather than as an opaque S3 object. That position is also the hazard the
# subsetting and degradation sections are written around. `[.factor` and
# `levels.vctrs_vctr` both run before anything registered on `vctrs_vctr` or on
# nothing at all, and each of them answers for a plain factor in a way that is
# wrong for this one: the first re-attaches the levels and the class and says
# nothing about the crossing metadata, and the second reports no levels at all.
#
# Fixtures are file-local and built from literal vectors. Where a test compares
# the class against a plain factor, the plain factor is written out from the
# same literals rather than derived from the joint vector, so that the
# comparison cannot be satisfied by two wrong answers agreeing.

# ---- fixtures ----------------------------------------------------------------

# The class vector downstream packages will pin. `"factor"` precedes
# `"vctrs_vctr"` deliberately, and `exact = TRUE` everywhere below is what holds
# the order in place.
joint_class <- function() {
  c("joint_exposure", "factor", "vctrs_vctr", "integer")
}

# Ten observations over two binary treatments, with a different number of
# observations in each of the four cells. Equal counts would let a cell ordering
# that is wrong in any way still satisfy `table()`.
smoking_qsmk <- function() {
  c(0, 1, 0, 1, 0, 0, 1, 1, 1, 0)
}

smoking_exercise <- function() {
  c("no", "no", "yes", "yes", "no", "yes", "yes", "yes", "yes", "no")
}

joint_smoking <- function() {
  joint_exposure(qsmk = smoking_qsmk(), exercise = smoking_exercise())
}

# The crossing, labelled and ordered as the contract states: the first component
# varies fastest within the second, and the cell where both components sit at
# their reference level comes first.
joint_smoking_levels <- function() {
  c(
    "qsmk = 0, exercise = no",
    "qsmk = 1, exercise = no",
    "qsmk = 0, exercise = yes",
    "qsmk = 1, exercise = yes"
  )
}

# Which cell each of the ten observations falls in, as an index into the levels.
joint_smoking_codes <- function() {
  c(1L, 2L, 3L, 4L, 1L, 3L, 4L, 4L, 4L, 1L)
}

joint_smoking_cells <- function() {
  joint_smoking_levels()[joint_smoking_codes()]
}

joint_smoking_components <- function() {
  list(qsmk = c("0", "1"), exercise = c("no", "yes"))
}

# The plain factor a joint exposure has to behave like wherever the formula
# machinery is concerned. Built from the same literals rather than from the
# joint vector, so that a shared mistake cannot make the two agree.
joint_smoking_baseline <- function() {
  factor(joint_smoking_cells(), levels = joint_smoking_levels())
}

# A second joint over the same two components with the same levels in the same
# order, and different data. Combining this with `joint_smoking()` is the
# compatible path.
joint_smoking_other <- function() {
  joint_exposure(
    qsmk = c(1, 0, 1, 0),
    exercise = c("yes", "yes", "no", "no")
  )
}

# Three ways the metadata of two joints can disagree, one per clause of the
# coercion rule: the components are named differently, their levels are
# different values, and their levels are the same values in a different order.
joint_renamed <- function() {
  joint_exposure(
    smoke = c(0, 1, 0, 1),
    exercise = c("no", "no", "yes", "yes")
  )
}

joint_relevelled <- function() {
  joint_exposure(
    qsmk = c(0, 1, 0, 1),
    exercise = c("none", "none", "some", "some")
  )
}

joint_reordered <- function() {
  joint_exposure(
    qsmk = c(0, 1, 0, 1),
    exercise = factor(
      c("no", "no", "yes", "yes"),
      levels = c("yes", "no")
    )
  )
}

# Every attribute `input` carried has to come back on `result` with the same
# value. The names are compared as a set first so that a failure names the
# fields that went missing rather than only reporting that two objects differ,
# and the whole set is compared rather than the fields this file can name, so
# that carrying the levels forward and dropping the crossing cannot satisfy it.
expect_joint_metadata_preserved <- function(result, input) {
  metadata <- attributes(input)

  expect_identical(
    setdiff(names(metadata), names(attributes(result))),
    character()
  )
  expect_identical(attributes(result)[names(metadata)], metadata)
}

# ---- construction ------------------------------------------------------------

test_that("joint_exposure() builds the documented class vector", {
  x <- joint_smoking()

  expect_s3_class(x, joint_class(), exact = TRUE)
  # `is.factor()` is what `model.frame()` and `model.matrix()` branch on, and it
  # reads the class vector rather than the class of the object, so it is the
  # assertion that says the formula machinery will treat this natively.
  expect_true(is.factor(x))
  expect_length(x, 10L)
})

test_that("joint_exposure() crosses the two components into labelled levels", {
  x <- joint_smoking()

  expect_identical(levels(x), joint_smoking_levels())
  expect_identical(nlevels(x), 4L)
  expect_identical(as.character(x), joint_smoking_cells())
})

test_that("the first component varies fastest within the second", {
  # This is the ordering `interaction()` produces with its default
  # `lex.order = FALSE`, and it is what makes the reference cell first: the
  # cell where both components sit at their reference level is the only one that
  # can come before every other. An ordering that varied the second component
  # fastest would put the same four labels in a different sequence, and every
  # contrast the formula machinery builds downstream would name a different
  # comparison.
  baseline <- interaction(
    factor(smoking_qsmk()),
    factor(smoking_exercise()),
    sep = ", "
  )
  expect_identical(
    levels(baseline),
    c("0, no", "1, no", "0, yes", "1, yes")
  )

  expect_identical(levels(joint_smoking()), joint_smoking_levels())
})

test_that("the reference cell crosses the two components' reference levels", {
  x <- joint_smoking()

  expect_identical(joint_reference(x), "qsmk = 0, exercise = no")
  expect_identical(levels(x)[[1]], joint_reference(x))
})

test_that("a factor component's own level order settles its reference level", {
  # A factor already states which of its levels is the reference, and a crossing
  # that re-sorted them would silently change which comparison the analysis
  # reports. Here `"yes"` is the reference of the second component even though
  # `"no"` sorts before it.
  x <- joint_reordered()

  expect_identical(
    levels(x),
    c(
      "qsmk = 0, exercise = yes",
      "qsmk = 1, exercise = yes",
      "qsmk = 0, exercise = no",
      "qsmk = 1, exercise = no"
    )
  )
  expect_identical(joint_reference(x), "qsmk = 0, exercise = yes")
  expect_identical(
    joint_components(x),
    list(qsmk = c("0", "1"), exercise = c("yes", "no"))
  )
})

test_that("zero is the reference level of a binary numeric component", {
  # A 0/1 numeric exposure has no level order of its own to read, and the
  # untreated cell is the one an effect is reported against. Writing the levels
  # as the characters `"0"` and `"1"` is what puts the label in the crossing.
  x <- joint_smoking()

  expect_identical(joint_components(x)$qsmk, c("0", "1"))
  expect_identical(joint_reference(x), "qsmk = 0, exercise = no")
})

test_that("the component types accepted agree with one another", {
  # Factor, character, and binary numeric are three spellings of the same
  # declaration, and a caller should not have to know which one the package
  # prefers. The whole object is compared, so the codes, the levels, and the
  # crossing metadata all have to agree, not only the labels.
  from_literals <- joint_smoking()
  from_factors <- joint_exposure(
    qsmk = factor(smoking_qsmk()),
    exercise = factor(smoking_exercise())
  )

  expect_identical(from_factors, from_literals)
})

test_that("constructing a joint exposure signals nothing", {
  expect_no_condition(joint_smoking())
})

# ---- validation --------------------------------------------------------------

test_that("joint_exposure() refuses anything but two components", {
  # The class declares a crossing of exactly two treatments. One component is
  # not a crossing, and three is a design the categorical-exposure route
  # downstream has no reference cell rule for.
  expect_error(
    joint_exposure(qsmk = smoking_qsmk()),
    class = "causalgenerics_joint_exposure_two_components"
  )
  expect_error(
    joint_exposure(qsmk = smoking_qsmk()),
    class = "causalgenerics_invalid_joint_exposure"
  )
  expect_error(
    joint_exposure(),
    class = "causalgenerics_joint_exposure_two_components"
  )
  expect_error(
    joint_exposure(
      qsmk = c(0, 1, 0, 1),
      exercise = c("no", "no", "yes", "yes"),
      sex = c("f", "m", "f", "m")
    ),
    class = "causalgenerics_joint_exposure_two_components"
  )
})

test_that("the two-component refusal states the contract", {
  expect_snapshot(error = TRUE, joint_exposure(qsmk = smoking_qsmk()))
})

test_that("joint_exposure() refuses an unnamed component", {
  # The names are the exposure names the cell labels are written from, so a
  # component with none leaves a label that names no variable. Both spellings
  # of the fault are refused: no names at all, and one name missing.
  expect_error(
    joint_exposure(smoking_qsmk(), smoking_exercise()),
    class = "causalgenerics_joint_exposure_unnamed_component"
  )
  expect_error(
    joint_exposure(smoking_qsmk(), smoking_exercise()),
    class = "causalgenerics_invalid_joint_exposure"
  )
  expect_error(
    joint_exposure(smoking_qsmk(), exercise = smoking_exercise()),
    class = "causalgenerics_joint_exposure_unnamed_component"
  )

  # The second position unnamed rather than the first. R has no surface syntax
  # for an empty argument name, so the call is assembled instead.
  components <- stats::setNames(
    list(smoking_qsmk(), smoking_exercise()),
    c("qsmk", "")
  )
  expect_error(
    do.call(joint_exposure, components),
    class = "causalgenerics_joint_exposure_unnamed_component"
  )
})

test_that("the unnamed-component refusal states the contract", {
  expect_snapshot(
    error = TRUE,
    joint_exposure(smoking_qsmk(), exercise = smoking_exercise())
  )
})

test_that("joint_exposure() refuses a component with one observed level", {
  # A component that never varies is not a treatment in these data, and the
  # crossing it produces has half its cells empty. The refusal is keyed to the
  # component rather than to the cells so that the message points at the
  # variable the caller can do something about.
  constant <- c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1)

  expect_error(
    joint_exposure(qsmk = constant, exercise = smoking_exercise()),
    class = "causalgenerics_joint_exposure_constant_component"
  )
  expect_error(
    joint_exposure(qsmk = constant, exercise = smoking_exercise()),
    class = "causalgenerics_invalid_joint_exposure"
  )

  # A factor that declares two levels and shows one is the same fault. The
  # declared levels say a crossing is possible and the data say it was never
  # observed, which is the case a check written against `levels()` alone would
  # let through.
  expect_error(
    joint_exposure(
      qsmk = smoking_qsmk(),
      exercise = factor(rep("yes", 10), levels = c("no", "yes"))
    ),
    class = "causalgenerics_joint_exposure_constant_component"
  )
})

test_that("the constant-component refusal names the component", {
  constant <- c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
  cnd <- tryCatch(
    joint_exposure(qsmk = constant, exercise = smoking_exercise()),
    error = identity
  )

  expect_identical(cnd$component, "qsmk")
  expect_snapshot(
    error = TRUE,
    joint_exposure(qsmk = constant, exercise = smoking_exercise())
  )
})

test_that("joint_exposure() refuses a missing value in a component", {
  # A missing exposure has no cell, and the crossing has no way to record one.
  # Refusing at construction keeps the decision with the caller, who is the only
  # one who knows whether the row should be dropped or the value recovered.
  with_na <- c("no", NA, "yes", "yes", "no", "yes", "yes", "yes", "yes", "no")

  expect_error(
    joint_exposure(qsmk = smoking_qsmk(), exercise = with_na),
    class = "causalgenerics_joint_exposure_missing_value"
  )
  expect_error(
    joint_exposure(qsmk = smoking_qsmk(), exercise = with_na),
    class = "causalgenerics_invalid_joint_exposure"
  )

  cnd <- tryCatch(
    joint_exposure(qsmk = smoking_qsmk(), exercise = with_na),
    error = identity
  )
  expect_identical(cnd$component, "exercise")
})

test_that("the missing-value refusal states the contract", {
  with_na <- c("no", NA, "yes", "yes", "no", "yes", "yes", "yes", "yes", "no")

  expect_snapshot(
    error = TRUE,
    joint_exposure(qsmk = smoking_qsmk(), exercise = with_na)
  )
})

test_that("joint_exposure() refuses a factor that declares NA as a level", {
  # `anyNA()` reads the values, and a factor that declares `NA` as one of its
  # levels gives its missing observations an ordinary code, so they are not
  # missing by that test. The exposure is still unknown for those rows, and the
  # cell the crossing would write for them names an absence rather than a
  # treatment, which is the same fault an `NA` value is and takes the same
  # refusal.
  declared_na <- addNA(factor(c("no", NA, "yes", "yes")))

  expect_error(
    joint_exposure(qsmk = c(0, 1, 0, 1), exercise = declared_na),
    class = "causalgenerics_joint_exposure_missing_value"
  )
  expect_error(
    joint_exposure(qsmk = c(0, 1, 0, 1), exercise = declared_na),
    class = "causalgenerics_invalid_joint_exposure"
  )

  cnd <- tryCatch(
    joint_exposure(qsmk = c(0, 1, 0, 1), exercise = declared_na),
    error = identity
  )
  expect_identical(cnd$component, "exercise")

  # `exclude = NULL` is the other spelling of the same declaration, and a caller
  # should not have to know which one produced the level to know it is refused.
  excluded_none <- factor(c("no", NA, "yes", "yes"), exclude = NULL)

  expect_error(
    joint_exposure(qsmk = c(0, 1, 0, 1), exercise = excluded_none),
    class = "causalgenerics_joint_exposure_missing_value"
  )
  expect_error(
    joint_exposure(qsmk = c(0, 1, 0, 1), exercise = excluded_none),
    class = "causalgenerics_invalid_joint_exposure"
  )
})

test_that("the declared-NA-level refusal states the contract", {
  declared_na <- addNA(factor(c("no", NA, "yes", "yes")))

  expect_snapshot(
    error = TRUE,
    joint_exposure(qsmk = c(0, 1, 0, 1), exercise = declared_na)
  )
})

test_that("an unpopulated NA level refuses as a missing value", {
  # Here the `NA` level is declared and nothing falls in it, so the crossing does
  # have empty cells and the positivity refusal would fire on its own terms. It
  # would be the wrong answer twice over: the cells are empty because the level
  # should not be there, and a caller told to coarsen a component or restrict the
  # analysis has been sent to fix something that is not broken.
  declared_only <- factor(
    c("no", "no", "yes", "yes"),
    levels = c("no", "yes", NA),
    exclude = NULL
  )

  expect_error(
    joint_exposure(qsmk = c(0, 1, 0, 1), exercise = declared_only),
    class = "causalgenerics_joint_exposure_missing_value"
  )
  expect_error(
    joint_exposure(qsmk = c(0, 1, 0, 1), exercise = declared_only),
    class = "causalgenerics_invalid_joint_exposure"
  )
})

test_that("joint_exposure() refuses a crossing with an empty cell", {
  # Both components vary, so nothing about either one on its own is wrong. What
  # is wrong is the crossing: nobody in these data both quit and exercises, so
  # the joint effect of the two exposures is not identified here and no weight,
  # no contrast, and no balance check over the crossing can be computed. This is
  # a positivity violation, and the refusal has to say so rather than report a
  # level count.
  never_both <- c(0, 1, 0, 0, 1, 0)
  exercise <- c("no", "no", "yes", "yes", "no", "yes")

  expect_error(
    joint_exposure(qsmk = never_both, exercise = exercise),
    class = "causalgenerics_joint_exposure_empty_cell"
  )
  expect_error(
    joint_exposure(qsmk = never_both, exercise = exercise),
    class = "causalgenerics_invalid_joint_exposure"
  )
})

test_that("the empty-cell refusal names the cell and the positivity failure", {
  never_both <- c(0, 1, 0, 0, 1, 0)
  exercise <- c("no", "no", "yes", "yes", "no", "yes")
  cnd <- tryCatch(
    joint_exposure(qsmk = never_both, exercise = exercise),
    error = identity
  )

  # The cells are a field as well as part of the sentence, so a handler reports
  # them without parsing the message for them.
  expect_identical(cnd$cells, "qsmk = 1, exercise = yes")
  expect_true(grepl("qsmk = 1, exercise = yes", conditionMessage(cnd)))
  expect_true(grepl("identified", conditionMessage(cnd)))

  expect_snapshot(
    error = TRUE,
    joint_exposure(qsmk = never_both, exercise = exercise)
  )
})

# ---- accessors ---------------------------------------------------------------

test_that("is_joint_exposure() separates the class from a plain factor", {
  expect_true(is_joint_exposure(joint_smoking()))
  expect_false(is_joint_exposure(joint_smoking_baseline()))
  expect_false(is_joint_exposure(factor(c("a", "b"))))
  expect_false(is_joint_exposure(joint_smoking_cells()))
  expect_false(is_joint_exposure(NULL))
})

test_that("joint_components() reports the component names and their levels", {
  # The names are the order the crossing was declared in, and each component's
  # levels are in the order the crossing varies them. Comparing the whole list
  # pins both at once; comparing `names()` alone would let a class that recorded
  # the levels in sorted order pass.
  expect_identical(
    joint_components(joint_smoking()),
    joint_smoking_components()
  )
  expect_identical(
    names(joint_components(joint_smoking())),
    c("qsmk", "exercise")
  )
})

test_that("joint_reference() reports the reference cell label", {
  expect_identical(joint_reference(joint_smoking()), "qsmk = 0, exercise = no")
  expect_identical(joint_reference(joint_renamed()), "smoke = 0, exercise = no")
  expect_identical(
    joint_reference(joint_relevelled()),
    "qsmk = 0, exercise = none"
  )
})

# ---- coercion ----------------------------------------------------------------

test_that("identical metadata makes the class the common type of two joints", {
  x <- joint_smoking()
  ptype <- vctrs::vec_ptype2(x, joint_smoking_other())

  expect_s3_class(ptype, joint_class(), exact = TRUE)
  expect_length(ptype, 0L)
  expect_identical(levels(ptype), joint_smoking_levels())
  expect_identical(joint_components(ptype), joint_smoking_components())
})

test_that("combining two joints with identical metadata keeps the class", {
  x <- joint_smoking()
  combined <- vctrs::vec_c(x, joint_smoking_other())

  expect_s3_class(combined, joint_class(), exact = TRUE)
  expect_length(combined, 14L)
  expect_identical(levels(combined), joint_smoking_levels())
  expect_identical(joint_components(combined), joint_smoking_components())
  expect_identical(joint_reference(combined), joint_reference(x))
  expect_identical(
    as.character(combined),
    c(
      joint_smoking_cells(),
      "qsmk = 1, exercise = yes",
      "qsmk = 0, exercise = yes",
      "qsmk = 1, exercise = no",
      "qsmk = 0, exercise = no"
    )
  )
})

test_that("joints whose metadata disagrees warn and degrade to a bare factor", {
  # Three ways to disagree, one per clause of the rule. Each one means the two
  # vectors declare different crossings, and a combine that kept the class would
  # hand downstream code one declaration for data drawn from two.
  for (other in list(joint_renamed(), joint_relevelled(), joint_reordered())) {
    expect_warning(
      vctrs::vec_c(joint_smoking(), other),
      class = "causalgenerics_joint_exposure_incompatible_metadata"
    )
    expect_warning(
      vctrs::vec_c(joint_smoking(), other),
      class = "causalgenerics_joint_exposure_downgrade"
    )

    combined <- suppressWarnings(vctrs::vec_c(joint_smoking(), other))
    expect_s3_class(combined, "factor", exact = TRUE)
    expect_false(is_joint_exposure(combined))
    expect_identical(
      levels(combined),
      union(joint_smoking_levels(), levels(other))
    )
  }
})

test_that("the downgrade between two joints warns once per combine", {
  # A caller sizes its own warning output on one combine producing one warning,
  # and `expect_warning()` proves only that at least one was signalled.
  signalled <- character()
  combined <- withCallingHandlers(
    vctrs::vec_c(joint_smoking(), joint_renamed()),
    warning = function(w) {
      signalled <<- c(signalled, class(w)[[1]])
      invokeRestart("muffleWarning")
    }
  )

  expect_identical(
    signalled,
    "causalgenerics_joint_exposure_incompatible_metadata"
  )
  expect_s3_class(combined, "factor", exact = TRUE)
})

test_that("combining a joint with a bare factor warns and degrades", {
  # The bare factor carries no declaration at all, so the crossing cannot be
  # claimed for the combined vector even when the levels happen to line up.
  # Both argument orders take the same route.
  x <- joint_smoking()
  plain <- joint_smoking_baseline()

  expect_warning(
    vctrs::vec_c(x, plain),
    class = "causalgenerics_joint_exposure_foreign_type"
  )
  expect_warning(
    vctrs::vec_c(x, plain),
    class = "causalgenerics_joint_exposure_downgrade"
  )
  expect_warning(
    vctrs::vec_c(plain, x),
    class = "causalgenerics_joint_exposure_foreign_type"
  )

  combined <- suppressWarnings(vctrs::vec_c(x, plain))
  expect_s3_class(combined, "factor", exact = TRUE)
  expect_false(is_joint_exposure(combined))
  expect_identical(levels(combined), joint_smoking_levels())
  expect_length(combined, 20L)
})

test_that("combining a joint with character warns and degrades to character", {
  x <- joint_smoking()

  expect_warning(
    vctrs::vec_c(x, "qsmk = 1, exercise = no"),
    class = "causalgenerics_joint_exposure_foreign_type"
  )
  expect_warning(
    vctrs::vec_c("qsmk = 1, exercise = no", x),
    class = "causalgenerics_joint_exposure_foreign_type"
  )

  combined <- suppressWarnings(vctrs::vec_c(x, "qsmk = 1, exercise = no"))
  expect_identical(
    combined,
    c(joint_smoking_cells(), "qsmk = 1, exercise = no")
  )
})

test_that("casting a factor or a character vector into the class is refused", {
  # A crossing cannot be recovered from labels. Anything that reproduced one by
  # parsing the level strings would accept a factor whose labels merely look
  # like a crossing and record a declaration nobody made.
  ptype <- vctrs::vec_ptype2(joint_smoking(), joint_smoking_other())

  expect_error(
    vctrs::vec_cast(joint_smoking_baseline(), ptype),
    class = "vctrs_error_incompatible_type"
  )
  expect_error(
    vctrs::vec_cast(joint_smoking_cells(), ptype),
    class = "vctrs_error_incompatible_type"
  )
})

test_that("casting out of the class gives the labels and the codes", {
  # `table()` reaches the codes through `as.integer()`, and every report of the
  # cells reaches the labels through `as.character()`, so both directions out of
  # the class have to answer.
  x <- joint_smoking()

  expect_identical(vctrs::vec_cast(x, character()), joint_smoking_cells())
  expect_identical(vctrs::vec_cast(x, integer()), joint_smoking_codes())
  expect_identical(as.integer(x), joint_smoking_codes())
})

test_that("base c() keeps the declaration where vec_c() does", {
  # `"factor"` precedes `"vctrs_vctr"` in the class vector, so `c.factor()` runs
  # before anything vctrs would reach and combines the codes into a bare factor
  # without a word. That is the one degradation that would happen silently, and
  # it happens on the spelling a caller is most likely to write. `c()` has to
  # answer what `vec_c()` answers, declaration and warnings alike.
  x <- joint_smoking()

  doubled <- c(x, x)
  expect_s3_class(doubled, joint_class(), exact = TRUE)
  expect_length(doubled, 20L)
  expect_joint_metadata_preserved(doubled, x)
  expect_identical(as.character(doubled), rep(joint_smoking_cells(), 2L))

  expect_identical(
    c(x, joint_smoking_other()),
    vctrs::vec_c(x, joint_smoking_other())
  )
})

test_that("base c() warns and degrades where vec_c() does", {
  x <- joint_smoking()

  expect_warning(
    c(x, joint_renamed()),
    class = "causalgenerics_joint_exposure_incompatible_metadata"
  )
  expect_warning(
    c(x, joint_renamed()),
    class = "causalgenerics_joint_exposure_downgrade"
  )
  expect_s3_class(
    suppressWarnings(c(x, joint_renamed())),
    "factor",
    exact = TRUE
  )

  expect_warning(
    c(x, joint_smoking_baseline()),
    class = "causalgenerics_joint_exposure_foreign_type"
  )
  expect_warning(
    c(x, joint_smoking_baseline()),
    class = "causalgenerics_joint_exposure_downgrade"
  )

  degraded <- suppressWarnings(c(x, joint_smoking_baseline()))
  expect_s3_class(degraded, "factor", exact = TRUE)
  expect_false(is_joint_exposure(degraded))
  expect_identical(levels(degraded), joint_smoking_levels())
})

# ---- subsetting --------------------------------------------------------------

test_that("subsetting keeps the class and the whole declaration", {
  # `[.factor` runs before anything registered against `vctrs_vctr`, because
  # `"factor"` precedes it in the class vector, and it re-attaches the levels
  # and the class and nothing else. A slice that came back with the right class,
  # the right levels, and no crossing metadata would satisfy every assertion
  # written against the levels alone, so the whole attribute set is compared.
  x <- joint_smoking()

  expect_s3_class(x[2:4], joint_class(), exact = TRUE)
  expect_joint_metadata_preserved(x[2:4], x)
  expect_identical(as.character(x[2:4]), joint_smoking_cells()[2:4])

  # Every cell stays a level of the slice even when no observation in it falls
  # there. A contrast over the crossing is defined by the declaration, not by
  # what a subset happens to contain.
  expect_identical(levels(x[c(1L, 5L)]), joint_smoking_levels())
  expect_identical(joint_components(x[c(1L, 5L)]), joint_smoking_components())
  expect_identical(joint_reference(x[c(1L, 5L)]), joint_reference(x))

  expect_s3_class(
    x[c(TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE)],
    joint_class(),
    exact = TRUE
  )
  expect_joint_metadata_preserved(
    x[c(TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE)],
    x
  )
})

test_that("vec_slice() keeps the class and the whole declaration", {
  # tibble and dplyr never reach a column through `[`. They slice it through
  # `vec_proxy()`, `vec_slice()`, and `vec_restore()`, which is a separate path
  # with its own way of going wrong: a proxy that drops to the bare integer
  # codes and a restore that rebuilds around them gives back a vector with the
  # right class and the right values and no crossing on it, and a
  # `dplyr::filter()` or an `arrange()` would lose the declaration without a
  # word. Nothing asserted against `[` reaches this path, so it is pinned on its
  # own terms rather than as a corollary.
  x <- joint_smoking()

  sliced <- vctrs::vec_slice(x, 2:4)
  expect_s3_class(sliced, joint_class(), exact = TRUE)
  expect_joint_metadata_preserved(sliced, x)
  expect_identical(levels(sliced), joint_smoking_levels())
  expect_identical(as.character(sliced), joint_smoking_cells()[2:4])

  # A slice holding two of the four cells keeps all four as levels, for the same
  # reason a `[` slice does: the crossing is what was declared, not what the
  # rows that survived a filter happen to contain.
  partial <- vctrs::vec_slice(x, c(1L, 5L, 2L))
  expect_s3_class(partial, joint_class(), exact = TRUE)
  expect_joint_metadata_preserved(partial, x)
  expect_identical(levels(partial), joint_smoking_levels())
  expect_identical(joint_reference(partial), joint_reference(x))

  # A logical subscript is how a filter arrives, and a reordering subscript is
  # how a sort does. Both take the same route and neither may lose anything.
  filtered <- vctrs::vec_slice(
    x,
    c(TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE)
  )
  expect_s3_class(filtered, joint_class(), exact = TRUE)
  expect_joint_metadata_preserved(filtered, x)
  expect_identical(
    as.character(filtered),
    joint_smoking_cells()[c(1L, 3L, 5L, 7L, 9L)]
  )

  arranged <- vctrs::vec_slice(x, c(10L, 1L, 5L, 2L))
  expect_s3_class(arranged, joint_class(), exact = TRUE)
  expect_joint_metadata_preserved(arranged, x)
  expect_identical(
    as.character(arranged),
    joint_smoking_cells()[c(10L, 1L, 5L, 2L)]
  )
})

test_that("a bare subscript returns the joint exposure unchanged", {
  x <- joint_smoking()

  expect_identical(x[], x)
})

test_that("a single-element subscript keeps the class and the declaration", {
  x <- joint_smoking()

  expect_s3_class(x[[3]], joint_class(), exact = TRUE)
  expect_length(x[[3]], 1L)
  expect_joint_metadata_preserved(x[[3]], x)
  expect_identical(as.character(x[[3]]), "qsmk = 0, exercise = yes")
})

test_that("reordering keeps the class and the declaration", {
  # `rev()` and `sort()` both reach the vector through `[`, so they are the two
  # ways a caller meets the subsetting hazard without writing a subscript. The
  # sorted vector is also what says the codes order by the declared cell order
  # rather than by label.
  x <- joint_smoking()

  reversed <- rev(x)
  expect_s3_class(reversed, joint_class(), exact = TRUE)
  expect_joint_metadata_preserved(reversed, x)
  expect_identical(as.character(reversed), rev(joint_smoking_cells()))

  sorted <- sort(x)
  expect_s3_class(sorted, joint_class(), exact = TRUE)
  expect_joint_metadata_preserved(sorted, x)
  expect_identical(
    as.character(sorted),
    joint_smoking_levels()[sort(joint_smoking_codes())]
  )

  manual <- x[c(10L, 1L, 5L, 2L)]
  expect_s3_class(manual, joint_class(), exact = TRUE)
  expect_joint_metadata_preserved(manual, x)
})

# ---- formula machinery -------------------------------------------------------

test_that("model.matrix() gives what a plain factor of the same levels gives", {
  # The whole reason `"factor"` is in the class vector. The baseline half runs
  # first and on its own terms: a plain factor built from the same literals,
  # asserted to produce one column per non-reference cell and none for the
  # reference. The joint vector then has to produce that matrix exactly,
  # including the `assign` and `contrasts` attributes a downstream `predict()`
  # reads.
  baseline <- model.matrix(~x, data.frame(x = joint_smoking_baseline()))

  expect_identical(
    colnames(baseline),
    c(
      "(Intercept)",
      "xqsmk = 1, exercise = no",
      "xqsmk = 0, exercise = yes",
      "xqsmk = 1, exercise = yes"
    )
  )
  expect_false("xqsmk = 0, exercise = no" %in% colnames(baseline))

  joint <- model.matrix(~x, data.frame(x = joint_smoking()))

  expect_identical(joint, baseline)
})

test_that("table() counts the four cells in the declared order", {
  baseline <- table(joint_smoking_baseline())

  expect_identical(names(baseline), joint_smoking_levels())
  expect_identical(as.integer(baseline), c(3L, 1L, 2L, 4L))

  counted <- table(joint_smoking())

  expect_identical(names(counted), joint_smoking_levels())
  expect_identical(as.integer(counted), as.integer(baseline))
})

test_that("a joint exposure survives a round trip through a data frame", {
  # A fitting package holds the declaration in the model frame it was handed,
  # and reads it back off the column to decide which contrasts to report. A
  # column that came back a bare factor would still fit and would report the
  # crossing as though nobody had declared one.
  x <- joint_smoking()
  data <- data.frame(y = seq_len(10), x = x)

  expect_s3_class(data$x, joint_class(), exact = TRUE)
  expect_identical(data$x, x)
  expect_joint_metadata_preserved(data[3:6, "x"], x)
})

# ---- degradation -------------------------------------------------------------

test_that("droplevels() warns and hands back a bare factor", {
  # A joint exposure is defined by the full crossing, so a vector with cells
  # missing from its levels is not one. `droplevels()` is a request for exactly
  # that, and the answer is to give the caller what they asked for as a plain
  # factor and say that the declaration was dropped with it. Base R already
  # returns a bare factor here; the warning is what stops the loss being silent.
  x <- joint_smoking()
  partial <- x[c(1L, 5L, 2L)]

  expect_warning(
    droplevels(partial),
    class = "causalgenerics_joint_exposure_dropped_levels"
  )
  expect_warning(
    droplevels(partial),
    class = "causalgenerics_joint_exposure_downgrade"
  )

  dropped <- suppressWarnings(droplevels(partial))
  expect_s3_class(dropped, "factor", exact = TRUE)
  expect_false(is_joint_exposure(dropped))
  expect_identical(
    levels(dropped),
    c("qsmk = 0, exercise = no", "qsmk = 1, exercise = no")
  )

  # The rule does not turn on whether anything would actually be dropped. A
  # vector with every cell populated goes the same way, so a caller never has to
  # know which case they are in to know what class comes back.
  expect_warning(
    droplevels(x),
    class = "causalgenerics_joint_exposure_dropped_levels"
  )
  expect_s3_class(suppressWarnings(droplevels(x)), "factor", exact = TRUE)
})

test_that("a dropping subscript warns and hands back a bare factor", {
  # `x[i, drop = TRUE]` is `droplevels()` written as a subscript, and it reaches
  # the same place through `[.factor`. Left alone it either errors inside vctrs
  # or produces a factor with cells missing, so the same rule is stated here.
  x <- joint_smoking()

  expect_warning(
    x[c(1L, 5L, 2L), drop = TRUE],
    class = "causalgenerics_joint_exposure_dropped_levels"
  )

  dropped <- suppressWarnings(x[c(1L, 5L, 2L), drop = TRUE])
  expect_s3_class(dropped, "factor", exact = TRUE)
  expect_identical(
    levels(dropped),
    c("qsmk = 0, exercise = no", "qsmk = 1, exercise = no")
  )

  # `drop = FALSE` is the ordinary slice and stays a joint exposure.
  expect_no_condition(x[c(1L, 5L, 2L), drop = FALSE])
  expect_s3_class(x[c(1L, 5L, 2L), drop = FALSE], joint_class(), exact = TRUE)
})

test_that("replacing the levels warns and hands back a bare factor", {
  # `levels<-.factor` writes the levels attribute and leaves the class and the
  # crossing metadata where they were, so left alone this produces a vector that
  # claims a declaration its labels no longer match. Nothing downstream could
  # detect that, which is why the operation degrades rather than being refused
  # in place.
  x <- joint_smoking()

  expect_warning(
    levels(x) <- c("a", "b", "c", "d"),
    class = "causalgenerics_joint_exposure_replaced_levels"
  )
  expect_warning(
    {
      y <- joint_smoking()
      levels(y) <- c("a", "b", "c", "d")
    },
    class = "causalgenerics_joint_exposure_downgrade"
  )

  expect_s3_class(x, "factor", exact = TRUE)
  expect_false(is_joint_exposure(x))
  expect_identical(levels(x), c("a", "b", "c", "d"))
  expect_identical(
    as.character(x),
    c("a", "b", "c", "d", "a", "c", "d", "d", "d", "a")
  )
  expect_null(attr(x, "components"))
})

# ---- print and format --------------------------------------------------------

test_that("printing names the class and the two components", {
  # The snapshot records the exact form. These assertions are what the form has
  # to carry whatever it looks like: a reader has to be able to tell a joint
  # exposure from a factor of pasted labels, and to see which two variables were
  # crossed without reading the levels.
  printed <- paste(capture.output(print(joint_smoking())), collapse = "\n")

  expect_true(grepl("joint_exposure", printed, fixed = TRUE))
  expect_true(grepl("qsmk", printed, fixed = TRUE))
  expect_true(grepl("exercise", printed, fixed = TRUE))
})

test_that("the printed and formatted forms are recorded", {
  x <- joint_exposure(
    qsmk = c(0, 1, 0, 1),
    exercise = c("no", "no", "yes", "yes")
  )

  expect_snapshot(print(x))
  expect_snapshot(format(x))
  expect_snapshot(vctrs::vec_ptype_full(x))
  expect_snapshot(vctrs::vec_ptype_abbr(x))
})

test_that("format() writes one cell label per observation", {
  # Whether the labels are padded to a common width is the printed form's
  # business and the snapshot's to record. What `format()` owes every caller is
  # one entry per observation, in order, carrying that observation's cell.
  x <- joint_smoking()
  formatted <- format(x)

  expect_type(formatted, "character")
  expect_length(formatted, 10L)
  expect_null(names(formatted))
  expect_identical(trimws(formatted), joint_smoking_cells())
})

# ---- registration ------------------------------------------------------------

test_that("the base methods the class overrides are registered", {
  # Downstream packages reach these through the S3 method table a NAMESPACE
  # `S3method()` directive fills in, and the table an entry goes to belongs to
  # the environment of the generic rather than to this package. Every generic
  # here is base's, and every one of them already resolves to something that is
  # wrong for a joint exposure: `levels.vctrs_vctr` reports no levels at all,
  # `[.factor` and `[[.factor` say nothing about the crossing metadata,
  # `format.factor` and `print.factor` report a plain factor, and
  # `droplevels.factor` and `levels<-.factor` change the level set while leaving
  # the declaration in place, and `c.factor` combines the codes into a bare
  # factor without a word.
  generics <- c(
    "levels",
    "levels<-",
    "[",
    "[[",
    "format",
    "print",
    "droplevels",
    "c"
  )

  registered <- vapply(
    generics,
    function(generic) {
      table <- s3_methods_table(generic)
      !is.null(table) &&
        exists(
          paste0(generic, ".joint_exposure"),
          envir = table,
          inherits = FALSE
        )
    },
    logical(1)
  )

  expect_identical(generics[!registered], character())
})

test_that("the vctrs coercion methods are registered", {
  # vctrs resolves `vec_ptype2()` and `vec_cast()` by an exact lookup on the
  # first element of each class vector, so every pairing needs its own entry and
  # each one is stored under its full two-class name. The table belongs to
  # vctrs, and `s3_methods_table()` resolves a generic from this package's
  # namespace, so it reaches that table only when the generic is imported. A
  # failure to find `vec_ptype2` or `vec_cast` at all reports a missing
  # `importFrom` rather than a missing method.
  methods <- c(
    "vec_ptype2.joint_exposure.joint_exposure",
    "vec_ptype2.joint_exposure.factor",
    "vec_ptype2.factor.joint_exposure",
    "vec_ptype2.joint_exposure.character",
    "vec_ptype2.character.joint_exposure",
    "vec_cast.joint_exposure.joint_exposure",
    "vec_cast.factor.joint_exposure",
    "vec_cast.character.joint_exposure",
    "vec_cast.integer.joint_exposure"
  )

  registered <- vapply(
    methods,
    function(method) {
      generic <- if (startsWith(method, "vec_ptype2")) {
        "vec_ptype2"
      } else {
        "vec_cast"
      }
      table <- s3_methods_table(generic)
      !is.null(table) && exists(method, envir = table, inherits = FALSE)
    },
    logical(1)
  )

  expect_identical(methods[!registered], character())
})

test_that("the vctrs printing methods are registered", {
  # Without these the class falls back to the vctrs defaults, which write an
  # abbreviated class name where the crossing should appear, in the place a data
  # frame or a tibble names the column's type.
  generics <- c("vec_ptype_abbr", "vec_ptype_full")

  registered <- vapply(
    generics,
    function(generic) {
      table <- s3_methods_table(generic)
      !is.null(table) &&
        exists(
          paste0(generic, ".joint_exposure"),
          envir = table,
          inherits = FALSE
        )
    },
    logical(1)
  )

  expect_identical(generics[!registered], character())
})

test_that("the joint_exposure methods dispatch from outside the test frame", {
  # `UseMethod()` searches the frame the generic was called from before it
  # consults the method table, so every assertion above would pass on a method
  # that is only a local function here. `baseenv()` cannot see this frame, so a
  # call made from there reaches the registered method or nothing at all, which
  # is the route a downstream package's own namespace takes.
  x <- joint_smoking()

  expect_identical(dispatch_from_baseenv(levels, x), joint_smoking_levels())
  expect_identical(
    trimws(dispatch_from_baseenv(format, x)),
    joint_smoking_cells()
  )
  expect_s3_class(
    dispatch_from_baseenv(`[`, x, 2:4),
    joint_class(),
    exact = TRUE
  )
})
