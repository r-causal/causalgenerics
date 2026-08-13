# The exposure-type resolver the r-causal ecosystem shares. Whether an exposure
# reads as binary, categorical, or continuous decides which weight formula is
# applied, which balance statistic is computed, and which diagnostic is run, so
# the reading has to be the same wherever it is taken. Four packages carried a
# copy of it, and copies drift: two of them counted a missing value as a level
# of its own, which reads a binary exposure recorded with any missingness as a
# three-level one and takes it off the binary path entirely. The reading here is
# the one propensity established, and it is the reading the ecosystem inherits.
#
# The machinery is three layers. `observed_values()`, `has_two_levels()`, and
# `is_categorical()` read the data. `detect_exposure_type()` classifies from
# those readings and announces what it read. `match_exposure_type()` resolves a
# caller's declared type against the set its function supports, and
# `check_forced_type()` refuses a declaration the data cannot carry.
#
# The division between the last two is the design. A declared type is how the
# caller models the exposure and wins over the heuristics outright, so
# `match_exposure_type()` returns it without looking at the data at all: a
# ten-level dose declared continuous is fitted as a dose, and no unique-value
# ratio gets to overrule that. The only declarations that are impossible are
# structural, and a function that wants them refused asks for that separately
# through `check_forced_type()`.

#' Detect the type of an exposure
#'
#' @description
#' `detect_exposure_type()` reads an exposure vector and reports whether the
#' analysis it feeds should treat it as binary, categorical, or continuous. It
#' is the reading the r-causal packages take when a caller declares no type, so
#' the classification lives in one place and every package that needs it gives
#' the same answer for the same column.
#'
#' @details
#' The reading is taken in this order. A vector with exactly two observed values
#' is binary, whatever its type. A factor or character vector taking more than
#' two observed values is categorical, and one taking fewer is binary, so a
#' single-level factor is read as a degenerate binary exposure rather than as a
#' categorical one. Any other vector is categorical when `is_categorical()`
#' finds few enough distinct values in it, and continuous otherwise.
#'
#' A missing value is not a level. Every count is taken over the observed
#' values, so an exposure recorded with missing values is read as the exposure
#' it is: two observed values are binary however many observations are missing,
#' and a categorical exposure reports the categories it has rather than one
#' more.
#'
#' Detection classifies what it is handed and refuses nothing. An exposure with
#' one observed value, and one with none, each leave the branches above with a
#' type rather than an error, because what a degenerate exposure means depends
#' on what the calling function is about to fit with it. That refusal belongs to
#' the calling function. [check_forced_type()] is the structural check for a
#' type the caller declared.
#'
#' @param .exposure The exposure vector to classify.
#' @param arg The name the announcement gives the exposure argument. A function
#'   whose exposure argument is named something else passes its own name, so the
#'   announcement never points at an argument the caller does not have.
#' @param announce Should the detected type be announced? The announcement
#'   reports which branch the analysis took, so a caller passes `FALSE` where
#'   that is not worth reporting, such as one detecting several exposures at
#'   once. It is suppressed for every call, whatever this argument says, when
#'   `options(causalgenerics.quiet = TRUE)`.
#' @param call The environment a refusal reports as its calling context.
#'   Detection raises none of its own; the argument is here so that a caller
#'   threading its own call through the resolver passes it the same way to every
#'   part of it.
#'
#' @return A single string: `"binary"`, `"categorical"`, or `"continuous"`.
#'
#' @seealso [match_exposure_type()], which resolves a declared type against the
#'   types a function supports and calls this when the caller asked for
#'   detection.
#'
#' @export
#'
#' @examples
#' # Two observed values are binary, and the missing one is not a third value.
#' detect_exposure_type(c(0, 1, NA, 1))
#'
#' # More than two levels are categorical.
#' detect_exposure_type(factor(c("low", "medium", "high")))
#'
#' # A numeric exposure taking many distinct values is continuous.
#' detect_exposure_type(seq_len(300) / 300)
#'
#' # A function whose exposure argument is named otherwise says so.
#' detect_exposure_type(c(0, 1), arg = ".exposures")
#'
#' # `announce = FALSE` returns the reading without reporting it.
#' detect_exposure_type(c(0, 1), announce = FALSE)
detect_exposure_type <- function(
  .exposure,
  arg = ".exposure",
  announce = TRUE,
  call = rlang::caller_env()
) {
  exposure_type <- if (has_two_levels(.exposure)) {
    "binary"
  } else if (is.factor(.exposure) || is.character(.exposure)) {
    if (length(observed_values(.exposure)) > 2) {
      "categorical"
    } else {
      "binary"
    }
  } else if (is_categorical(.exposure)) {
    "categorical"
  } else {
    "continuous"
  }

  if (announce) {
    announce_exposure_type(exposure_type, arg)
  }

  exposure_type
}

#' Resolve the exposure type a function will work on
#'
#' @description
#' `match_exposure_type()` turns the `exposure_type` argument a function
#' exposes into the type it will actually work on. A caller either declares one
#' of the three types or asks for `"auto"`, and a function declares through
#' `valid_types` which of them it can answer.
#'
#' @details
#' Resolution has two tiers, and which one a request falls into is what the
#' refusal says. A type this package knows but the calling function does not
#' support is refused as unsupported, with an error of class
#' `causalgenerics_unsupported_exposure_type`, and the message names the types
#' the function does support. A string that is no exposure type at all is
#' refused as an unrecognized value instead, by [rlang::arg_match()], which
#' names the choices and suggests the nearest one. Both refusals report the
#' environment given in `call`, so a consumer package names its own function
#' rather than this one.
#'
#' A detected type outside `valid_types` reaches the same unsupported refusal as
#' a declared one. Naming a type and letting it be detected therefore give the
#' same answer, and a caller cannot get further into a function by leaving the
#' type off.
#'
#' A declared type that the function supports is returned as it stands. The data
#' is not consulted on that path and nothing is announced: the declaration is
#' how the caller models the exposure, and a heuristic over unique values does
#' not overrule it. A function that also wants a declaration the data cannot
#' carry refused calls [check_forced_type()], which is a separate decision.
#'
#' A function that narrows `valid_types` narrows the default of its own
#' `exposure_type` argument to match. The two are matched against each other, so
#' a function offering all four choices while supporting three of them reports
#' the mismatch on the first call that leaves the argument at its default.
#'
#' @inheritParams detect_exposure_type
#' @param exposure_type The type the caller declared, or `"auto"` to detect it
#'   from `.exposure`. A calling function passes its own argument through.
#' @param valid_types The types the calling function can answer, including
#'   `"auto"` where it offers detection. The unsupported refusal names this set
#'   without `"auto"`, which is a request for detection rather than a type an
#'   exposure can have.
#' @param call The environment both refusals report as their calling context.
#'
#' @return A single string: `"binary"`, `"categorical"`, or `"continuous"`.
#'
#' @seealso [detect_exposure_type()], which supplies the reading under
#'   `"auto"`, and [check_forced_type()], which refuses a declaration the data
#'   cannot carry.
#'
#' @export
#'
#' @examples
#' # A declared type is returned as it stands, without consulting the data.
#' match_exposure_type("continuous", c(0, 1, 2, 3))
#'
#' # `"auto"` reads the exposure and announces what it read.
#' match_exposure_type("auto", c(0, 1, 1, 0))
#'
#' # A function that works on two of the types says so.
#' try(match_exposure_type(
#'   "continuous",
#'   c(0, 1),
#'   valid_types = c("auto", "binary", "categorical")
#' ))
match_exposure_type <- function(
  exposure_type = c("auto", "binary", "categorical", "continuous"),
  .exposure,
  valid_types = c("auto", "binary", "categorical", "continuous"),
  arg = ".exposure",
  announce = TRUE,
  call = rlang::caller_env()
) {
  # A known type the function cannot answer is refused before the match, so that
  # it is reported as unsupported rather than as a value nobody recognizes. The
  # two say different things: one is a function that does not go this way, and
  # the other is a typo.
  if (
    rlang::is_string(exposure_type) &&
      exposure_type %in% exposure_types &&
      !exposure_type %in% valid_types
  ) {
    stop_unsupported_exposure_type(exposure_type, valid_types, call = call)
  }

  exposure_type <- rlang::arg_match(
    exposure_type,
    valid_types,
    error_call = call
  )

  if (exposure_type != "auto") {
    return(exposure_type)
  }

  exposure_type <- detect_exposure_type(
    .exposure,
    arg = arg,
    announce = announce,
    call = call
  )

  if (!exposure_type %in% valid_types) {
    stop_unsupported_exposure_type(exposure_type, valid_types, call = call)
  }

  exposure_type
}

#' Refuse an exposure type the data cannot carry
#'
#' @description
#' `check_forced_type()` is the structural check on a declared exposure type. A
#' declaration wins over the detection heuristics wherever the data can carry
#' it, so this refuses only the declarations that no data of this shape could
#' support.
#'
#' @details
#' Two declarations can be contradicted. A binary exposure takes exactly two
#' observed values, and a continuous one is numeric. Every vector can be read as
#' the set of levels it takes, so a categorical declaration is always possible
#' and always passes.
#'
#' Everything else stands. A numeric dose taking fifteen values declared
#' continuous is fitted as a dose even though detection would read it as
#' categorical, and a numeric exposure taking two values declared continuous is
#' the caller's decision. A check that refused those would put the unique-value
#' heuristic in charge of a reading the caller had already made.
#'
#' The refusal quotes what detection would have read, since dropping the
#' declaration is the route out of it. That reading is computed and not
#' announced: the caller declared a type, so an alert naming a different one
#' would read as the answer rather than as the alternative.
#'
#' The refusal names `exposure_type` as the argument the declaration came from,
#' which is what the ecosystem calls it.
#'
#' @inheritParams detect_exposure_type
#' @param forced The type the caller declared: `"binary"`, `"categorical"`, or
#'   `"continuous"`. `"auto"` is a request for detection rather than a
#'   declaration, and is refused here along with anything else that is not one
#'   of the three.
#' @param arg The name the refusal gives the exposure argument.
#' @param call The environment the refusal reports as its calling context.
#'
#' @return `NULL`, invisibly. The function is called for the refusal.
#'
#' @seealso [match_exposure_type()], which resolves the declaration this checks,
#'   and [detect_exposure_type()], which supplies the reading it quotes.
#'
#' @export
#'
#' @examples
#' # A declaration the data can carry passes silently.
#' check_forced_type("binary", c(0, 1, NA, 1))
#'
#' # So does one detection would have read the other way.
#' check_forced_type("continuous", rep_len(1:15, 300))
#'
#' # A binary exposure takes two values, and this one takes three.
#' try(check_forced_type("binary", factor(c("low", "medium", "high"))))
#'
#' # A continuous exposure is numeric.
#' try(check_forced_type("continuous", c("a", "b", "c")))
check_forced_type <- function(
  forced,
  .exposure,
  arg = ".exposure",
  call = rlang::caller_env()
) {
  forced <- rlang::arg_match0(forced, exposure_types, error_call = call)

  # The problem is written as a template and interpolated at the refusal, so
  # that the branch that found it is the branch that words it.
  problem <- switch(
    forced,
    binary = if (!has_two_levels(.exposure)) {
      paste0(
        "A {.val binary} exposure takes exactly two observed values, and ",
        "{.arg {arg}} takes {n_observed}."
      )
    },
    continuous = if (!is.numeric(.exposure)) {
      paste0(
        "A {.val continuous} exposure is numeric, and {.arg {arg}} is ",
        "{.obj_type_friendly {exposure_vec}}."
      )
    },
    categorical = NULL
  )

  if (is.null(problem)) {
    return(invisible(NULL))
  }

  # The facts the refusal words itself from, read once the refusal is certain.
  # `exposure_vec` is the exposure under a name cli can interpolate: a leading
  # dot reads as the start of inline markup rather than as a value.
  exposure_vec <- .exposure
  n_observed <- length(observed_values(.exposure))
  detected <- detect_exposure_type(
    .exposure,
    arg = arg,
    announce = FALSE,
    call = call
  )

  stop_exposure_type(
    c(
      paste0(
        "{.arg exposure_type} was set to {.val {forced}}, but {.arg {arg}} ",
        "cannot be treated that way."
      ),
      x = problem,
      i = paste0(
        "Drop {.arg exposure_type} to detect the type from the data, which ",
        "reads {.arg {arg}} as {.val {detected}}."
      )
    ),
    error_class = "causalgenerics_forced_exposure_type",
    call = call
  )
}

#' Read the values an exposure takes
#'
#' @description
#' The three readings [detect_exposure_type()] classifies from, exported so that
#' a package taking its own decision from the same facts reads them the same
#' way.
#'
#' `observed_values()` returns the distinct values a vector takes.
#' `has_two_levels()` reports whether there are exactly two of them.
#' `is_categorical()` applies the unique-value heuristic that separates a
#' categorical numeric exposure from a continuous one.
#'
#' @details
#' A missing value is not a level. `observed_values()` drops the missing values
#' before taking the distinct ones, and that is what the other two count, so an
#' exposure recorded with missing values is read as the exposure it is. Counting
#' the missing values as a value of their own would give a binary exposure three
#' levels, which takes it off the binary path entirely, and would make a
#' categorical one look as though it had one category more than it has.
#'
#' `is_categorical()` treats a vector as categorical when the number of values
#' it takes is less than 20 percent of the number of observations that are not
#' missing. The threshold is strict, so a ratio of exactly 0.2 is not
#' categorical. An exposure with nothing observed has no ratio to take and is
#' not categorical, which leaves it on the continuous branch of detection rather
#' than returning a missing value into the classification.
#'
#' The distinct values are returned in the order they first appear rather than
#' sorted, since they are only ever counted and sorting them would impose an
#' ordering on a categorical exposure that nothing else here honors. The values
#' keep the type they came in as, so the observed values of a factor are a
#' factor over the same declared levels.
#'
#' @param .x A vector.
#' @param .exposure An exposure vector.
#'
#' @return `observed_values()` returns the distinct observed values, as a vector
#'   of the same type as its input. `has_two_levels()` and `is_categorical()`
#'   each return a single logical value.
#'
#' @seealso [detect_exposure_type()], which combines these into a type.
#'
#' @name exposure-type-helpers
#'
#' @examples
#' observed_values(c(0, 1, NA, 1))
#' has_two_levels(c(0, 1, NA, 1))
#'
#' # 15 values over 300 observations is categorical; 300 over 300 is not.
#' is_categorical(rep_len(1:15, 300))
#' is_categorical(seq_len(300) / 300)
NULL

#' @rdname exposure-type-helpers
#' @export
observed_values <- function(.x) {
  unique(.x[!is.na(.x)])
}

#' @rdname exposure-type-helpers
#' @export
has_two_levels <- function(.x) {
  length(observed_values(.x)) == 2
}

#' @rdname exposure-type-helpers
#' @export
is_categorical <- function(.exposure) {
  n_non_na <- sum(!is.na(.exposure))
  if (n_non_na == 0) {
    return(FALSE)
  }

  ratio <- length(observed_values(.exposure)) / n_non_na
  if (is.nan(ratio)) {
    return(FALSE)
  }

  ratio < 0.2
}

# The types an exposure can be read as. `"auto"` is not among them: it is how a
# caller asks for detection, not an answer detection can give, which is why the
# unsupported refusal names the supported set without it.
exposure_types <- c("binary", "categorical", "continuous")

# Report the reading. The option is the user's route to silence and `announce`
# is the calling function's, and either one is enough on its own, so a user who
# has asked for quiet gets it whatever the function passes. `isTRUE()` rather
# than the option's own truthiness, so a value that is not a flag leaves the
# announcement in place rather than silencing it by accident.
announce_exposure_type <- function(exposure_type, arg) {
  if (isTRUE(getOption("causalgenerics.quiet", FALSE))) {
    return(invisible(NULL))
  }

  cli::cli_alert_info("Treating {.arg {arg}} as {exposure_type}")

  invisible(NULL)
}

# The condition shape of the resolver. The classes follow `stop_no_method()` in
# `R/utils.R`: one keyed to the fault, for tests and for handlers that care
# about this one refusal, and one general class for callers that want any
# refusal this package raised.
#
# These are the only conditions in the package built with cli rather than with
# base `errorCondition()`. The messages are the ones the ecosystem already
# shows, down to the `{.val}` and `{.arg}` markup the consuming packages format
# their own refusals with, and rewriting them in plain text would have made this
# package's copy the odd one out in every transcript. The rest of the package's
# conditions are unaffected and stay as they are.
#
# `.envir` defaults to the caller's frame, so the facts a refusal words itself
# from are read where they were computed.
stop_exposure_type <- function(
  message,
  error_class,
  call,
  .envir = parent.frame()
) {
  cli::cli_abort(
    message,
    class = c(error_class, "causalgenerics_error"),
    call = call,
    .envir = .envir
  )
}

# Signal that the type is one this package reads but the calling function cannot
# work on. The supported set is named without `"auto"`, which is a request
# rather than a type, and it is what the caller can do something about: the
# reading itself was not wrong.
stop_unsupported_exposure_type <- function(exposure_type, valid_types, call) {
  supported <- setdiff(valid_types, "auto")
  stop_exposure_type(
    c(
      "Exposure type {.val {exposure_type}} is not supported.",
      i = "Supported exposure types: {.val {supported}}."
    ),
    error_class = "causalgenerics_unsupported_exposure_type",
    call = call
  )
}
