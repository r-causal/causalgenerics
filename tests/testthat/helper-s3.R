# Helpers for exercising real S3 method registration.
#
# Downstream packages reach these generics through the S3 method table that a
# NAMESPACE `S3method()` directive fills in. A test that defines its method as a
# local function inside `test_that()` never touches that table, so it passes
# whether or not registration works. The helpers below register into the table
# for real and assert dispatch from an environment that cannot see the test
# frame.

# The S3 method table that `registerS3method()` writes to for `generic`.
#
# The `envir` argument of `registerS3method()` says where to look the generic
# up, not where to store the method. The method goes into the method table of
# the environment the generic itself lives in, so a table hardcoded to this
# package is right only for the generics this package owns. A method for
# `median()` goes to stats' table and a method for `[` goes to base's, and a
# lookup that misses that both warns on cleanup and leaves the entry behind for
# the rest of the session. The rules below mirror `registerS3method()`: group
# generics belong to base, as does anything that is not a closure, which covers
# every primitive. The one rule not mirrored is base's S4 branch, which unwraps
# an S4 generic through `methods::finalDefaultMethod()` before taking its
# environment; that branch is reachable only if an attached package masks one of
# these generics with a `setGeneric()` promotion, which none of them do.
s3_methods_table <- function(generic) {
  group_generics <- c("Math", "Ops", "matrixOps", "Summary", "Complex")
  home <- if (generic %in% group_generics) {
    asNamespace("base")
  } else {
    genfun <- get(generic, envir = asNamespace("causalgenerics"))
    if (typeof(genfun) == "closure") {
      environment(genfun)
    } else {
      asNamespace("base")
    }
  }
  home$.__S3MethodsTable__.
}

# Register `fn` as the method for `generic` and `class`, then put the method
# table back the way it was when `env` exits.
#
# The cleanup restores rather than removes, because a registration may replace a
# method that is already there. Deleting the entry for, say, `ess.default` would
# take this package's own method out of the table for the rest of the session,
# and every later test of the default behavior would then be asserting base R's
# dispatch failure rather than ours. Restoring also keeps registrations from
# leaking between tests, so a later test that registers the same class sees its
# own method and never a stale one. No path emits a warning as long as `fn` is a
# function, which is the only form used here: the entry is removed only after
# checking that it is still there. `registerS3method()` itself warns when `fn` is
# the character name of a method that `envir` does not hold, so a character path
# would need a check of its own.
local_s3_method <- function(generic, class, fn, env = parent.frame()) {
  name <- paste0(generic, ".", class)
  table <- s3_methods_table(generic)
  replaced <- !is.null(table) && exists(name, envir = table, inherits = FALSE)
  previous <- if (replaced) get(name, envir = table, inherits = FALSE) else NULL

  registerS3method(generic, class, fn, envir = asNamespace("causalgenerics"))

  withr::defer(
    {
      table <- s3_methods_table(generic)
      if (replaced) {
        assign(name, previous, envir = table)
      } else if (exists(name, envir = table, inherits = FALSE)) {
        rm(list = name, envir = table)
      }
    },
    envir = env
  )
  invisible(fn)
}

# Call `generic` from `baseenv()` instead of from the test frame.
#
# `UseMethod()` looks for a method in the environment the generic was called
# from before it consults the method table, so a generic called directly inside
# `test_that()` finds a method that is only a local function in the test. An
# assertion made that way says nothing about registration. `baseenv()` cannot
# see the test frame and its parent is the empty environment, so a call made
# here can only dispatch through the method table. That indirection is the
# point of these assertions; do not fold it back into a direct call.
#
# `generic` is the function itself here, not its name, which is the opposite of
# `local_s3_method()`, where `generic` is a character name. `do.call()` looks a
# character `what` up in `envir`, and `baseenv()` holds none of this package's
# generics. `do.call()` also runs with `quote = FALSE`, so it evaluates each
# element of `...` again inside `baseenv()`; a language object therefore cannot
# be passed through to the method unevaluated.
dispatch_from_baseenv <- function(generic, ...) {
  do.call(generic, list(...), envir = baseenv())
}
