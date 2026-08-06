# Signal that a generic was called on an object none of its registered methods
# handle. The condition carries a distinct class so downstream packages and
# tests can match it precisely, and it stays dependency-free by building on
# base R's `errorCondition()` rather than a messaging package.
stop_no_method <- function(generic, x, call = sys.call(-1)) {
  message <- paste0(
    "No `",
    generic,
    "()` method for an object of class <",
    class(x)[[1]],
    ">."
  )
  stop(errorCondition(
    message,
    generic = generic,
    class = c(
      paste0("causalgenerics_no_method_", generic),
      "causalgenerics_no_method"
    ),
    call = call
  ))
}

# Signal that an argument does not meet the contract the function documents.
# The classes follow `stop_no_method()`: one keyed to the argument, for tests and
# for handlers that care about this one argument, and one general class for
# callers that want any argument this package rejected. `must` completes the
# sentence "`<arg>` must ...".
stop_invalid_argument <- function(arg, must, call = sys.call(-1)) {
  message <- paste0("`", arg, "` must ", must, ".")
  stop(errorCondition(
    message,
    arg = arg,
    class = c(
      paste0("causalgenerics_invalid_argument_", arg),
      "causalgenerics_invalid_argument"
    ),
    call = call
  ))
}

# Signal that the results being pooled disagree about something they have to
# agree on. The classes follow `stop_invalid_argument()`: one keyed to the field
# they disagree about, for tests and for handlers that care about one kind of
# disagreement, and one general class for callers that want any refusal to pool.
# The distinct values are a field as well as part of the message, so a handler
# reports them without parsing the sentence for them, and they are a list rather
# than a vector because the effect labels a result reports are themselves a
# vector.
stop_pool_mismatch <- function(field, values, call = sys.call(-1)) {
  message <- paste0(
    pool_field_label(field),
    " must be the same in every result, but they report ",
    format_series(vapply(values, format_pool_value, character(1))),
    "."
  )
  stop(errorCondition(
    message,
    field = field,
    values = values,
    class = c(
      paste0("causalgenerics_pool_mismatch_", field),
      "causalgenerics_pool_mismatch"
    ),
    call = call
  ))
}

# What the results disagreed about, as the message names it. The field names are
# the pooled result's own, and four of the six are a component of it, but a
# message reading "`effects` must be the same" would name a field where the
# reader is thinking about a property of their analysis.
pool_field_label <- function(field) {
  switch(
    field,
    estimand = "The estimand",
    se_method = "The standard error method",
    effects = "The presentation mode",
    labels = "The effects reported",
    conf_level = "The confidence level",
    outcome_link = "The outcome model's link"
  )
}

# One of those values, written the way the message reports it. Strings are
# quoted, since an estimand and a link are read back as the labels they are, and
# a set of several values is parenthesized so that the commas separating one
# result's effects are not read as separating one result from the next.
#
# Recording nothing is one of the values a result can disagree on: an estimates
# frame with no `conf.level` column, or one whose rows disagree about the level,
# names no level for the pooled bounds to be reported at. It is written as a
# word, since an empty pair of parentheses in the middle of the sentence reads
# as a formatting fault rather than as the absence it reports.
format_pool_value <- function(value) {
  if (length(value) == 0L) {
    return("none")
  }

  written <- if (is.character(value)) {
    encodeString(value, quote = '"')
  } else {
    format(value, trim = TRUE)
  }
  if (length(written) == 1L) {
    return(written)
  }
  paste0("(", toString(written), ")")
}

# Join written items into a phrase rather than a list. `toString()` alone ends a
# sentence on "they report "ate", "att"", which reads as though it had been cut
# off; the conjunction is what makes the last item the last one.
format_series <- function(x) {
  if (length(x) < 2L) {
    return(paste0(x, collapse = ""))
  }
  if (length(x) == 2L) {
    return(paste(x, collapse = " and "))
  }
  paste0(toString(x[-length(x)]), ", and ", x[[length(x)]])
}

# Signal that nothing in the results says how much data the complete-data
# analysis had, so the pooled degrees of freedom are the large-sample ones. This
# is a warning rather than an error because the pooling is still done and the
# answer is still the honest one for a large sample; it is said out loud because
# the intervals it gives are the narrowest the adjustment can produce, which is
# the direction a reader is least likely to question. The class is keyed to the
# assumption, and there is no general warning class in this package for it to
# join.
warn_pool_large_sample <- function(call = sys.call(-1)) {
  warning(warningCondition(
    paste0(
      "No result reports the residual degrees of freedom of its complete-data ",
      "analysis, so a large sample is assumed; pass `dfcom` to name the count ",
      "the analyses were fitted with."
    ),
    class = "causalgenerics_pool_large_sample",
    call = call
  ))
}

# Signal that an object records no covariance to report. The classes follow
# `stop_no_method()`: one keyed to the class of the object, for tests and for
# handlers that care about one kind of object, and one general class for callers
# that want any missing covariance. The helper answers for a result and for a
# component model of one alike, so the message names the object rather than the
# result and points at the `ipw_vcov` contract both of them carry a covariance
# under. There is no fallback to offer, since the standard errors a result
# stores give the diagonal of that matrix and say nothing about the rest of it.
stop_no_vcov <- function(result, call = sys.call(-1)) {
  message <- paste0(
    "This `",
    result,
    "` object records no covariance to report; the package that produced it ",
    "attaches one when it supports the `ipw_vcov` contract."
  )
  stop(errorCondition(
    message,
    result = result,
    class = c(
      paste0("causalgenerics_no_vcov_", result),
      "causalgenerics_no_vcov"
    ),
    call = call
  ))
}

# Signal that the conditional reading has no covariance to report, because the
# outcome model carries none from the joint estimation. The classes follow
# `stop_no_vcov()`, whose general class this one also carries: a handler written
# for any covariance this package cannot report catches both, and a handler
# written for the conditional reading tells them apart. The specific class is
# keyed to the reading rather than to the result class, since it is the reading
# that cannot be answered. There is no fallback to offer: the covariance the
# outcome model computed for itself treats the estimated weights as fixed and
# reports an uncertainty the coefficients do not have.
stop_no_conditional_vcov <- function(call = sys.call(-1)) {
  message <- paste0(
    "The conditional reading reports the covariance the joint estimation of ",
    "the weights and the outcome implies, and this result's outcome model ",
    "records none; the package that produced the result attaches one by ",
    "wrapping the model with `new_ipw_model()`."
  )
  stop(errorCondition(
    message,
    class = c(
      "causalgenerics_no_conditional_vcov",
      "causalgenerics_no_vcov"
    ),
    call = call
  ))
}

# Signal that the corrected covariance an outcome model carries cannot be paired
# with the coefficients it is meant to report. The classes follow
# `stop_no_conditional_vcov()`, whose general class this one also carries, and
# deliberately not its specific one: that condition says the fitting package
# attached no block and is answered by wrapping the model, and this one is
# raised for a model that is wrapped already. The label sets are fields as well
# as parts of the message, so a handler reports them without parsing the
# sentence for them. A model whose coefficients carry no names is one of the
# ways the pairing fails, and the clause that would list them says so rather
# than listing an empty set, which would read as a model reporting no
# coefficients at all.
stop_conditional_vcov_mismatch <- function(
  block_labels,
  coef_labels,
  call = sys.call(-1)
) {
  reported <- if (is.null(coef_labels)) {
    "reports unnamed coefficients"
  } else {
    paste0(
      "reports coefficients named ",
      toString(encodeString(coef_labels, quote = '"'))
    )
  }
  message <- paste0(
    "The conditional covariance is labelled ",
    toString(encodeString(block_labels, quote = '"')),
    " and the outcome model ",
    reported,
    "; the package that produced the result attaches the block labelled by ",
    "coefficient name with `new_ipw_model()`."
  )
  stop(errorCondition(
    message,
    block_labels = block_labels,
    coef_labels = coef_labels,
    class = c(
      "causalgenerics_conditional_vcov_mismatch",
      "causalgenerics_no_vcov"
    ),
    call = call
  ))
}
