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

# Signal that a conditional table was asked to exponentiate coefficients that
# are not on a scale an exponential undoes. The marginal reading picks the rows
# to move out by label, so there is always a well-defined subset of them; a
# conditional reading has no such labels, and the outcome model's link settles
# the question for the whole table at once.
#
# Three classes rather than the usual two. The keyed class is the one
# `check_flag()` already raises for an `exponentiate` that is not a flag, so a
# caller who wants only this refusal has nothing to match on among those two;
# the specific class in front of them is that. The keyed class stays so that a
# handler written for anything wrong with this argument catches both, which is
# what the specific-then-general shape means everywhere else here.
stop_exponentiate_link <- function(link, exponentiable, call = sys.call(-1)) {
  message <- paste0(
    "`exponentiate` needs coefficients on a scale an exponential undoes, and ",
    "the outcome models were fitted with the ",
    encodeString(link, quote = '"'),
    " link, whose coefficients are not on one; only the ",
    format_series(encodeString(exponentiable, quote = '"')),
    " links exponentiate in the conditional reading."
  )
  stop(errorCondition(
    message,
    link = link,
    exponentiable = exponentiable,
    class = c(
      "causalgenerics_exponentiate_link",
      "causalgenerics_invalid_argument_exponentiate",
      "causalgenerics_invalid_argument"
    ),
    call = call
  ))
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

# Signal that a result was asked for a reading it does not support. A result
# records the readings it can present, and a fitting package records one of them
# when the other has no meaning for the analysis it ran: an exposure entering the
# outcome model through several columns has no single coefficient to read as the
# conditional effect. The classes follow `stop_no_method()`: one keyed to the
# reading that was asked for, since a caller who wants only that one has nothing
# to match on otherwise, and one general class for callers that care that a
# reading is absent rather than which one it was.
#
# The message says the result supports one reading only, which is what the set
# always holds here: a result supporting both refuses neither, so the reading
# asked for is absent only when the other one is the whole set. Both facts are
# fields as well as parts of the sentence, so a handler reports them without
# parsing the message for them. There is nothing for the caller to do to the
# result, since the set is written where the result is built, so the sentence
# says where it comes from rather than offering a route.
#
# `position` is how the pooling says which result of the set the refusal is
# about, and it is `NULL` everywhere a single result was asked directly, where
# "this result" is the only result there is. The position is a sentence of its
# own after the one the direct refusal ends with, so the two refusals read the
# same as far as the direct one goes, and it is a field as well, so a caller
# handling the condition goes back to the result without counting through the
# sentence.
stop_unsupported_reading <- function(
  effects,
  readings,
  position = NULL,
  call = sys.call(-1)
) {
  message <- paste0(
    "This result supports the ",
    readings,
    " reading only, so there is no ",
    effects,
    " reading of it to report; the package that produced it records the ",
    "readings it supports when it builds the result."
  )
  if (!is.null(position)) {
    message <- paste0(
      message,
      " It is the result at position ",
      position,
      " of `fits`."
    )
  }
  stop(errorCondition(
    message,
    effects = effects,
    readings = readings,
    position = position,
    class = c(
      paste0("causalgenerics_unsupported_reading_", effects),
      "causalgenerics_unsupported_reading"
    ),
    call = call
  ))
}

# Signal that a pooled result was asked to present a reading it does not carry.
# The classes follow `stop_no_method()`: one keyed to the reading that was asked
# for, since a caller who wants only the conditional reading has nothing to match
# on otherwise, and one general class for callers that care that the reading is
# absent rather than which one it was.
#
# Two facts get the same pair of classes and different sentences. A reading that
# was tried and could not be pooled has a reason recorded, and it is repeated
# here verbatim rather than paraphrased: it is the refusal the same request would
# have raised at pooling time, and a second wording of it would be a second thing
# to keep true. A result stored before both readings were kept has no reason to
# report, and what its caller has to do is pool the results again rather than
# change how they were fitted.
stop_pool_missing_surface <- function(effects, reason, call = sys.call(-1)) {
  message <- if (is.null(reason)) {
    paste0(
      "This pooled result carries no ",
      effects,
      " reading: it was pooled before both readings were kept, so there is ",
      "nothing recorded to present; pooling the results again gives a result ",
      "that carries both."
    )
  } else {
    paste0(
      "This pooled result carries no ",
      effects,
      " reading, since pooling that reading over the same results was ",
      "refused. ",
      reason
    )
  }
  stop(errorCondition(
    message,
    effects = effects,
    reason = reason,
    class = c(
      paste0("causalgenerics_pool_missing_surface_", effects),
      "causalgenerics_pool_missing_surface"
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
    "The conditional covariance is labeled ",
    toString(encodeString(block_labels, quote = '"')),
    " and the outcome model ",
    reported,
    "; the package that produced the result attaches the block labeled by ",
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

# Signal that the components handed to `joint_exposure()` are not a crossing of
# two treatments. The classes follow `stop_invalid_argument()`: one keyed to the
# fault, and one general class shared by every refusal to declare a joint
# exposure, so a caller who cares only that the declaration was rejected has one
# class to catch. The count is written into the sentence rather than carried as a
# field, since a caller who wants it has the call that supplied it.
stop_joint_exposure_two_components <- function(n, call = sys.call(-1)) {
  supplied <- if (n == 1L) "1 was supplied" else paste0(n, " were supplied")
  message <- paste0(
    "A joint exposure is the crossing of exactly two treatments, and ",
    supplied,
    "; pass one named vector per treatment."
  )
  stop(errorCondition(
    message,
    class = c(
      "causalgenerics_joint_exposure_two_components",
      "causalgenerics_invalid_joint_exposure"
    ),
    call = call
  ))
}

# Signal that a component of a joint exposure arrived without a name. The names
# are what the cell labels are written from, so an unnamed component leaves a
# label that names no variable and a declaration nothing downstream can read.
stop_joint_exposure_unnamed_component <- function(call = sys.call(-1)) {
  message <- paste0(
    "Every component of a joint exposure must be named, because the names are ",
    "the treatment names its cell labels are written from; pass each treatment ",
    "as `name = value`."
  )
  stop(errorCondition(
    message,
    class = c(
      "causalgenerics_joint_exposure_unnamed_component",
      "causalgenerics_invalid_joint_exposure"
    ),
    call = call
  ))
}

# Signal that both components of a joint exposure arrived under one name. The
# names are what separate the two treatments in everything the declaration is
# read for: each cell label names a variable on either side of it, and
# `joint_components()` reports one component per name. One name used twice
# leaves a crossing nothing downstream can read, and which component a cell
# varies cannot be recovered from anything the vector carries. The name is a
# field as well as part of the sentence, so a handler reports it without parsing
# the message for it.
stop_joint_exposure_shared_name <- function(component, call = sys.call(-1)) {
  message <- paste0(
    "Both components of a joint exposure are named `",
    component,
    "`, so its cells would name that one treatment twice and nothing reading ",
    "the declaration could tell which component a cell varies; give the two ",
    "treatments distinct names."
  )
  stop(errorCondition(
    message,
    component = component,
    class = c(
      "causalgenerics_joint_exposure_shared_name",
      "causalgenerics_invalid_joint_exposure"
    ),
    call = call
  ))
}

# Signal that a component takes one value in the data it was given. The refusal
# is keyed to the component rather than to the cells the crossing would leave
# empty, so that the message points at the variable the caller can do something
# about. The component is a field as well as part of the sentence, so a handler
# reports it without parsing the message for it.
stop_joint_exposure_constant_component <- function(
  component,
  call = sys.call(-1)
) {
  message <- paste0(
    "`",
    component,
    "` takes one value in these data, so it is not a treatment the crossing ",
    "can vary; a joint exposure needs both of its components observed at two ",
    "or more levels."
  )
  stop(errorCondition(
    message,
    component = component,
    class = c(
      "causalgenerics_joint_exposure_constant_component",
      "causalgenerics_invalid_joint_exposure"
    ),
    call = call
  ))
}

# Signal that a component admits a missing exposure. An observation with no
# exposure falls in no cell, and the crossing has nowhere to record that, so the
# decision stays with the caller, who is the only one who knows whether the row
# should be dropped or the value recovered. A factor that declares `NA` as one of
# its levels is the same fault even when no observation currently takes it, so
# the sentence names the level as well as the values.
stop_joint_exposure_missing_value <- function(component, call = sys.call(-1)) {
  message <- paste0(
    "`",
    component,
    "` admits a missing exposure, and an observation whose exposure is unknown ",
    "falls in no cell of the crossing; drop or recover those observations, and ",
    "any `NA` the component declares as a level, before declaring the joint ",
    "exposure."
  )
  stop(errorCondition(
    message,
    component = component,
    class = c(
      "causalgenerics_joint_exposure_missing_value",
      "causalgenerics_invalid_joint_exposure"
    ),
    call = call
  ))
}

# Signal that the crossing declares a cell nothing in the data falls in. Both
# components may vary and the fault still be here: what is wrong is the joint
# distribution, not either margin. That is a positivity violation, so the message
# says so rather than reporting a count of levels. The cells are a field as well
# as part of the sentence, so a handler reports them without parsing it.
stop_joint_exposure_empty_cell <- function(cells, call = sys.call(-1)) {
  message <- paste0(
    "Nothing in these data falls in the ",
    if (length(cells) == 1L) "cell " else "cells ",
    format_series(encodeString(cells, quote = '"')),
    ", so the crossing violates positivity and the joint effect of the two ",
    "treatments is not identified here; coarsen a component or restrict the ",
    "analysis to the cells that are populated."
  )
  stop(errorCondition(
    message,
    cells = cells,
    class = c(
      "causalgenerics_joint_exposure_empty_cell",
      "causalgenerics_invalid_joint_exposure"
    ),
    call = call
  ))
}

# Signal that two joint exposures declare different crossings, so the vector
# they combine into carries no declaration. The classes follow the refusals
# above: one keyed to the reason the declaration was given up, and one general
# class shared by every operation that gives it up, so a caller who cares only
# that a joint exposure degraded has one class to catch.
warn_joint_exposure_incompatible_metadata <- function(call = sys.call(-1)) {
  warning(warningCondition(
    paste0(
      "These joint exposures declare different crossings, so the result is a ",
      "plain factor over the cells of both; a combined vector carries a ",
      "crossing only when both sides declare the same one."
    ),
    class = c(
      "causalgenerics_joint_exposure_incompatible_metadata",
      "causalgenerics_joint_exposure_downgrade"
    ),
    call = call
  ))
}

# Signal that a joint exposure was combined with a vector that declares no
# crossing. The type is a field as well as part of the sentence, since a handler
# that reports the combine reports what it was combined with.
warn_joint_exposure_foreign_type <- function(type, call = sys.call(-1)) {
  warning(warningCondition(
    paste0(
      "A joint exposure combines with a ",
      type,
      " vector only by giving up its crossing, because the other side declares ",
      "none; the result carries the cell labels and nothing about the two ",
      "treatments they cross."
    ),
    type = type,
    class = c(
      "causalgenerics_joint_exposure_foreign_type",
      "causalgenerics_joint_exposure_downgrade"
    ),
    call = call
  ))
}

# Signal that an operation asked for a level set narrower than the crossing. A
# joint exposure is defined by the full crossing, so a vector with cells missing
# from its levels is not one, and the caller gets what they asked for as a plain
# factor. The warning fires whether or not a level would actually go, so that a
# caller never has to know which case they are in to know what class comes back.
warn_joint_exposure_dropped_levels <- function(call = sys.call(-1)) {
  warning(warningCondition(
    paste0(
      "A joint exposure is defined by the full crossing of its two treatments, ",
      "so dropping unused cells gives up the declaration; the result is a ",
      "plain factor over the cells that remain."
    ),
    class = c(
      "causalgenerics_joint_exposure_dropped_levels",
      "causalgenerics_joint_exposure_downgrade"
    ),
    call = call
  ))
}

# Signal that the cell labels were rewritten. Replacing the levels relabels the
# cells without touching the crossing they were written from, which would leave a
# vector claiming a declaration its labels no longer match, and nothing
# downstream could detect that. The operation degrades rather than being refused
# so that relabeling a factor stays available to the caller who wants it.
warn_joint_exposure_replaced_levels <- function(call = sys.call(-1)) {
  warning(warningCondition(
    paste0(
      "Replacing the levels of a joint exposure relabels its cells without ",
      "changing the crossing they were written from, so the declaration is ",
      "given up; the result is a plain factor with the new levels."
    ),
    class = c(
      "causalgenerics_joint_exposure_replaced_levels",
      "causalgenerics_joint_exposure_downgrade"
    ),
    call = call
  ))
}
