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

# Signal that a result records no covariance of the effects it reports. The
# classes follow `stop_no_method()`: one keyed to the result class, for tests and
# for handlers that care about one kind of result, and one general class for
# callers that want any missing covariance. There is no fallback to offer, since
# the standard errors a result stores give the diagonal of that matrix and say
# nothing about the rest of it.
stop_no_vcov <- function(result, call = sys.call(-1)) {
  message <- paste0(
    "This `",
    result,
    "` result records no covariance of the effects it reports; refit it with a ",
    "current version of the package that produced it."
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
