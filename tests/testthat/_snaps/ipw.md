# ipw() errors when the weighting object is named both ways

    Code
      ipw(wt_mod = wt_mod, ps_mod = other)
    Condition
      Error in `ipw()`:
      ! `ps_mod` must not be supplied together with `wt_mod`, which names the same argument.

# check_ipw_reentry() refuses an argument list that would recurse

    Code
      check_ipw_reentry(list(outcome_mod = NULL))
    Condition
      Error:
      ! The re-entered `ipw()` call must supply `wt_mod` and drop `ps_mod`.

# the `ps_mod` deprecation warning reads as intended

    The `ps_mod` argument of `ipw()` is deprecated; name it `wt_mod` instead.

