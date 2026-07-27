# Helpers for exercising real S3 method registration.
#
# Downstream packages reach these generics through the S3 method table that a
# NAMESPACE `S3method()` directive fills in. A test that defines its method as a
# local function inside `test_that()` never touches that table, so it passes
# whether or not registration works. The helpers below register into the table
# for real and assert dispatch from an environment that cannot see the test
# frame.

# causalgenerics' S3 method table. `registerS3method()` writes here for any
# generic defined in this package, because it stores a method in the environment
# of the generic it is registering against.
s3_methods_table <- function() {
  asNamespace("causalgenerics")$.__S3MethodsTable__.
}

# Register `fn` as the method for `generic` and `class` in causalgenerics' S3
# method table, then remove the entry again when `env` exits. The removal keeps
# registrations from leaking between tests, so a later test that registers the
# same class sees its own method and never a stale one.
local_s3_method <- function(generic, class, fn, env = parent.frame()) {
  registerS3method(generic, class, fn, envir = asNamespace("causalgenerics"))
  withr::defer(
    rm(list = paste0(generic, ".", class), envir = s3_methods_table()),
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
dispatch_from_baseenv <- function(generic, ...) {
  do.call(generic, list(...), envir = baseenv())
}
