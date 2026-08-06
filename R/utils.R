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
