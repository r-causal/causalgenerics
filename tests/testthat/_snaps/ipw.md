# check_ipw_reentry() refuses an argument list that would recurse

    Code
      check_ipw_reentry(list(outcome_mod = NULL))
    Condition
      Error:
      ! The re-entered `ipw()` call must supply `wt_mod` and drop `ps_mod`.

# the `ps_mod` deprecation warning reads as intended

    The `ps_mod` argument of `ipw()` is deprecated; name it `wt_mod` instead.

