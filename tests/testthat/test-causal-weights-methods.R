# The read-only unwrappers are registered against the abstract class, so every
# assertion below builds a `cg_probe_wts` vector and registers nothing at all
# for it. Anything that works here does so because dispatch walked past
# `cg_probe_wts` and found a method on `causal_wts`.
#
# The fixture also supplies no `vec_ptype2()` or `vec_cast()` method, which is
# deliberate on two counts. It is what a concrete class is allowed to look like
# partway through being built, and it is the instrument the comparison
# assertions rely on: combining this fixture with a double is a hard error, so a
# comparison that returns a logical vector proves it never reached
# `vec_cast_common()`.

cg_probe_wts_class <- c("cg_probe_wts", "causal_wts", "vctrs_vctr", "double")

cg_probe_wts <- function(x = c(2, 1, 3)) {
  new_causal_wts(x, subclass = "cg_probe_wts")
}

# The fixture for every operation that has to hand its result back intact. An
# operation of that kind must delegate, and delegation is observable only
# through metadata this package never names: `new_causal_wts()` documents that
# it "neither names nor types" the fields in the dots, and `estimand` is the one
# field the abstract layer knows exists. So an implementation that rebuilds the
# result around fresh data and copies `estimand` onto it produces the right
# class and the right values, and an assertion written against named fields
# alone cannot see the difference. `stabilized` and `groups` are shaped like
# what the concrete weight classes really carry, a scalar flag and a list of
# per-level indices, and no method here could know to preserve either.
cg_probe_wts_meta <- function(x = c(2, 1, 3)) {
  new_causal_wts(
    x,
    subclass = "cg_probe_wts",
    estimand = "ate",
    stabilized = TRUE,
    groups = list(exposed = c(1L, 3L), unexposed = 2L)
  )
}

# Every attribute `input` carried has to come back on `result` with the same
# value. The names are compared as a set first so that a failure names the
# fields that went missing rather than only reporting that two lists differ, and
# the comparison is written against the whole set rather than field by field so
# that carrying any particular subset forward cannot satisfy it.
expect_metadata_preserved <- function(result, input) {
  metadata <- attributes(input)

  expect_identical(
    setdiff(names(metadata), names(attributes(result))),
    character()
  )
  expect_identical(attributes(result)[names(metadata)], metadata)
}

test_that("vec_math() keeps the class for the cumulative functions", {
  # The class vector alone does not separate restoring the result from `.x` from
  # rebuilding a fresh weight vector around the cumulated data, and neither do
  # the values. Only the metadata does, so the whole attribute set is what gets
  # asserted: a cumulative sum of weights that quietly forgot what it was
  # carrying is the failure worth pinning, whichever field it forgot.
  wts <- cg_probe_wts_meta()

  expect_s3_class(cumsum(wts), cg_probe_wts_class, exact = TRUE)
  expect_identical(as.vector(cumsum(wts)), c(2, 3, 6))
  expect_metadata_preserved(cumsum(wts), wts)

  expect_s3_class(cumprod(wts), cg_probe_wts_class, exact = TRUE)
  expect_identical(as.vector(cumprod(wts)), c(2, 2, 6))
  expect_metadata_preserved(cumprod(wts), wts)

  expect_s3_class(cummin(wts), cg_probe_wts_class, exact = TRUE)
  expect_identical(as.vector(cummin(wts)), c(2, 1, 1))
  expect_metadata_preserved(cummin(wts), wts)

  expect_s3_class(cummax(wts), cg_probe_wts_class, exact = TRUE)
  expect_identical(as.vector(cummax(wts)), c(2, 2, 3))
  expect_metadata_preserved(cummax(wts), wts)
})

test_that("vec_math() unwraps every other function to a bare double", {
  wts <- cg_probe_wts(c(4, 1, 9))

  expect_identical(sqrt(wts), c(2, 1, 3))
  expect_identical(log(wts), log(c(4, 1, 9)))
  expect_identical(exp(wts), exp(c(4, 1, 9)))

  # `expect_identical()` already compares attributes, but a weight vector whose
  # class survived a transformation that is no longer a weight is the exact
  # failure this pins, so it is worth naming.
  expect_null(attributes(sqrt(wts)))
  expect_null(attributes(log(wts)))
})

test_that("the Summary group members return bare numeric", {
  wts <- cg_probe_wts()

  expect_identical(sum(wts), 6)
  expect_identical(prod(wts), 6)
})

test_that("the Summary group members cover all() and any()", {
  # `all()` and `any()` are group members too, so they reach the same method,
  # and both coerce a double to logical with a warning of base R's own. Every
  # nonzero weight becomes `TRUE` under that coercion, so the first fixture
  # below mixes a zero in with nonzero weights to keep the two answers from
  # agreeing by construction. The other two pin the ends that mixed fixture
  # cannot reach: weights that are all nonzero, where `all()` is `TRUE`, and
  # weights that are all zero, where `any()` is `FALSE`.
  wts <- cg_probe_wts(c(2, 0, 3))
  coercion <- "coercing argument of type 'double' to logical"

  expect_warning(all_wts <- all(wts), coercion)
  expect_identical(all_wts, FALSE)

  expect_warning(any_wts <- any(wts), coercion)
  expect_identical(any_wts, TRUE)

  expect_warning(all_nonzero <- all(cg_probe_wts()), coercion)
  expect_identical(all_nonzero, TRUE)

  expect_warning(any_zero <- any(cg_probe_wts(c(0, 0, 0))), coercion)
  expect_identical(any_zero, FALSE)
})

test_that("the Summary group members unwrap every argument in the dots", {
  wts <- cg_probe_wts()

  expect_identical(sum(wts, cg_probe_wts(c(1, 1, 1))), 9)
  expect_identical(sum(wts, 4), 10)
  expect_identical(prod(wts, 2), 12)
})

test_that("min(), max(), and range() return bare numeric", {
  wts <- cg_probe_wts()

  expect_identical(min(wts), 1)
  expect_identical(max(wts), 3)
  expect_identical(range(wts), c(1, 3))
})

test_that("min(), max(), and range() unwrap every argument in the dots", {
  wts <- cg_probe_wts()

  expect_identical(min(wts, cg_probe_wts(c(0.5, 4))), 0.5)
  expect_identical(max(wts, cg_probe_wts(c(0.5, 4))), 4)
  expect_identical(range(wts, cg_probe_wts(c(0.5, 4))), c(0.5, 4))
  expect_identical(min(wts, 0.25), 0.25)
})

test_that("the summary methods drop missing values when na.rm is TRUE", {
  wts <- cg_probe_wts(c(2, NA, 3))

  expect_identical(sum(wts, na.rm = TRUE), 5)
  expect_identical(min(wts, na.rm = TRUE), 2)
  expect_identical(max(wts, na.rm = TRUE), 3)
  expect_identical(range(wts, na.rm = TRUE), c(2, 3))
  expect_identical(median(wts, na.rm = TRUE), 2.5)
})

test_that("the summary methods keep missing values by default", {
  wts <- cg_probe_wts(c(2, NA, 3))

  expect_identical(sum(wts), NA_real_)
  expect_identical(min(wts), NA_real_)
  expect_identical(max(wts), NA_real_)
  expect_identical(range(wts), c(NA_real_, NA_real_))
  expect_identical(median(wts), NA_real_)
})

test_that("median() returns a bare numeric summary", {
  expect_identical(median(cg_probe_wts()), 2)
  expect_identical(median(cg_probe_wts(c(1, 2, 3, 4))), 2.5)
})

test_that("quantile() returns a bare named numeric rather than a weight vector", {
  wts <- cg_probe_wts()

  expect_identical(
    quantile(wts),
    c(`0%` = 1, `25%` = 1.5, `50%` = 2, `75%` = 2.5, `100%` = 3)
  )
  expect_identical(quantile(wts, probs = 0.5), c(`50%` = 2))
  expect_identical(
    quantile(cg_probe_wts(c(2, NA, 3)), probs = 0.5, na.rm = TRUE),
    c(`50%` = 2.5)
  )
})

test_that("summary() returns base R's summary of the underlying double", {
  # A concrete weight class is not obliged to supply a `vec_cast()` method, and
  # the fixture here has none, so the payload has to be reached without one.
  wts <- cg_probe_wts()

  expect_identical(summary(wts), summary(c(2, 1, 3)))
  expect_s3_class(summary(wts), "summaryDefault")
})

test_that("anyDuplicated() returns the index of the first duplicate", {
  # `anyDuplicated.vctrs_vctr()` answers with a logical, so the integer index
  # that callers of `anyDuplicated()` expect, and that the base default already
  # returns for a bare double, is the point of the method.
  expect_identical(anyDuplicated(cg_probe_wts(c(1, 1, 2))), 2L)
  expect_identical(anyDuplicated(cg_probe_wts(c(1, 2, 3))), 0L)
})

test_that("diff() returns a bare numeric difference", {
  expect_identical(diff(cg_probe_wts()), c(-1, 2))
  expect_identical(diff(cg_probe_wts(c(1, 2, 4, 8)), lag = 2L), c(3, 6))
  expect_identical(diff(cg_probe_wts(c(1, 2, 4, 8)), differences = 2L), c(1, 2))
})

test_that("a bare subscript returns the weight vector unchanged", {
  # Reading `i` at all in the method would force the missing argument and error,
  # so this branch exists to be a no-op. The fixture carries metadata for the
  # same reason the subsetting one below does: unchanged means every attribute,
  # not only the ones this package can name. One `expect_identical()` covers
  # that here, since it compares the class vector and the whole attribute set
  # along with the values.
  wts <- cg_probe_wts_meta()

  expect_identical(wts[], wts)
})

test_that("ordinary subsetting keeps the class", {
  # Handing the subscript on to the `vctrs_vctr` method and building a fresh
  # weight vector around the sliced data both produce the right class vector and
  # the right values. The metadata is the only thing that tells the two apart,
  # and a slice must carry all of it forward, not just the estimand: propensity's
  # `psw` carries a stabilization flag that its own tests assert across a subset.
  #
  # Carried is all this asserts. A field such as balancing's per-level indices
  # comes back holding indices into the pre-slice vector, which is wrong for the
  # sliced one, and nothing here says otherwise. Reindexing a field the abstract
  # layer cannot name belongs to the concrete class's own `vec_restore()`.
  wts <- cg_probe_wts_meta()

  expect_s3_class(wts[2:3], cg_probe_wts_class, exact = TRUE)
  expect_identical(as.vector(wts[2:3]), c(1, 3))
  expect_metadata_preserved(wts[2:3], wts)

  expect_s3_class(wts[c(TRUE, FALSE, TRUE)], cg_probe_wts_class, exact = TRUE)
  expect_identical(as.vector(wts[c(TRUE, FALSE, TRUE)]), c(2, 3))
  expect_metadata_preserved(wts[c(TRUE, FALSE, TRUE)], wts)
})

test_that("a matrix or array subscript returns the bare data", {
  # vctrs rejects a matrix subscript outright, so a caller that indexes weights
  # with one gets an error unless the method hands the subscript to the
  # underlying double.
  wts <- cg_probe_wts()

  expect_identical(wts[matrix(c(1L, 3L), ncol = 1)], c(2, 3))
  expect_identical(wts[matrix(c(TRUE, FALSE, TRUE), ncol = 1)], c(2, 3))
  expect_identical(wts[array(2L)], 1)
})

test_that("the comparison operators return logical", {
  wts <- cg_probe_wts()

  expect_identical(wts == 1, c(FALSE, TRUE, FALSE))
  expect_identical(wts != 1, c(TRUE, FALSE, TRUE))
  expect_identical(wts < 2, c(FALSE, TRUE, FALSE))
  expect_identical(wts > 2, c(FALSE, FALSE, TRUE))
  expect_identical(wts <= 2, c(TRUE, TRUE, FALSE))
  expect_identical(wts >= 2, c(TRUE, FALSE, TRUE))
  expect_identical(wts == cg_probe_wts(c(2, 2, 2)), c(TRUE, FALSE, FALSE))
})

test_that("the comparison operators work with the weight vector on the right", {
  wts <- cg_probe_wts()

  expect_identical(1 == wts, c(FALSE, TRUE, FALSE))
  expect_identical(1 != wts, c(TRUE, FALSE, TRUE))
  expect_identical(2 > wts, c(FALSE, TRUE, FALSE))
  expect_identical(2 < wts, c(FALSE, FALSE, TRUE))
  expect_identical(2 >= wts, c(TRUE, TRUE, FALSE))
  expect_identical(2 <= wts, c(TRUE, FALSE, TRUE))
})

test_that("comparison never combines the weight vector with the other operand", {
  # `glm.fit()` evaluates `weights == 0` and `weights > 0` repeatedly inside
  # `profile.glm()`. A comparison that routed through `vec_equal()` would ask
  # `vec_cast_common()` to combine the weight vector with a double on every one
  # of those calls, which for a concrete class that warns on the downgrade means
  # one warning per call and over a hundred from a single `tidy()`. This fixture
  # has no `vec_ptype2()` method at all, so the combine would raise
  # `vctrs_error_incompatible_type` here rather than warn.
  wts <- cg_probe_wts()

  expect_no_condition(wts == 0)
  expect_no_condition(wts > 0)
  expect_no_condition(wts != 0)
})

test_that("comparison enforces vctrs N-or-1 recycling rather than base's", {
  # The size rule is the one thing the short circuit must not relax. Base R
  # would recycle a length-two operand against a length-three weight vector and
  # only warn about it, so this is what separates the two behaviors.
  wts <- cg_probe_wts()

  expect_identical(wts == 2, c(TRUE, FALSE, FALSE))
  expect_identical(wts == c(2, 1, 1), c(TRUE, TRUE, FALSE))

  expect_error(wts == c(1, 2), class = "vctrs_error_incompatible_size")
  expect_error(wts > c(1, 2), class = "vctrs_error_incompatible_size")
  expect_error(c(1, 2) < wts, class = "vctrs_error_incompatible_size")
  expect_error(
    wts != cg_probe_wts(c(1, 2)),
    class = "vctrs_error_incompatible_size"
  )
})

test_that("the unwrapper methods are registered against causal_wts", {
  # Downstream packages reach these through the S3 method table a NAMESPACE
  # `S3method()` directive fills in, and the table an entry goes to belongs to
  # the environment of the generic rather than to this package. Every one of
  # these generics is foreign, so this asserts both that the methods exist and
  # that they went somewhere dispatch will look.
  #
  # It is also the only assertion that separates this package's `vec_math()` and
  # `[` methods from the `vctrs_vctr` defaults for the cases where the two agree.
  #
  # `vec_math` goes through the same lookup as everything else, which makes this
  # assert one more thing on its own: `s3_methods_table()` resolves a generic
  # from this package's namespace, so it can only reach vctrs' table for
  # `vec_math` if the generic is imported. The package cannot be loaded without
  # that import, so a failure here would report the cause rather than the
  # symptom.
  generics <- c(
    "vec_math",
    "Summary",
    "min",
    "max",
    "range",
    "median",
    "quantile",
    "summary",
    "anyDuplicated",
    "diff",
    "[",
    "==",
    "!=",
    "<",
    ">",
    "<=",
    ">="
  )

  registered <- vapply(
    generics,
    function(generic) {
      table <- s3_methods_table(generic)
      !is.null(table) &&
        exists(paste0(generic, ".causal_wts"), envir = table, inherits = FALSE)
    },
    logical(1)
  )

  expect_identical(generics[!registered], character())
})
