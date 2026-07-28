# Helper for exercising a deprecation that warns once per session.

# Clear the flag that holds the `ps_mod` deprecation warning to one warning per
# session, then put it back the way it was when `env` exits.
#
# The flag lives in the package's own environment and so outlives any one test.
# Without a reset, whether a call warns depends on which tests ran before it,
# which makes an assertion about the warning pass or fail on test order alone.
# The restore matters for the same reason in the other direction: a test that
# leaves the flag cleared hands the next one a fresh warning it did not expect.
#
# A flag that was never set restores as `NULL`, which reads as not yet warned
# exactly as `FALSE` does. Nothing observes the binding itself.
local_ps_mod_warning_reset <- function(env = parent.frame()) {
  previous <- deprecation_state$ipw_ps_mod

  withr::defer(
    assign("ipw_ps_mod", previous, envir = deprecation_state),
    envir = env
  )

  deprecation_state$ipw_ps_mod <- FALSE
  invisible(NULL)
}
