# Changelog

## causalgenerics 0.0.0.9000

- Added
  [`new_causal_wts()`](https://r-causal.github.io/causalgenerics/reference/new_causal_wts.md),
  the low-level constructor for the abstract `causal_wts` weight class
  that concrete weight classes in the ecosystem are built on.

- Added `causal_wts` methods for
  [`is_causal_wt()`](https://r-causal.github.io/causalgenerics/reference/causal-weights.md),
  [`estimand()`](https://r-causal.github.io/causalgenerics/reference/causal-weights.md),
  and `estimand<-()`, so a concrete class inherits all three.

- Added `causal_wts` methods for the read-only operations that name no
  metadata field and so give the answer the underlying double would:
  `vec_math()`, the `Summary` group generic,
  [`min()`](https://rdrr.io/r/base/Extremes.html),
  [`max()`](https://rdrr.io/r/base/Extremes.html),
  [`range()`](https://rdrr.io/r/base/range.html),
  [`median()`](https://rdrr.io/r/stats/median.html),
  [`quantile()`](https://rdrr.io/r/stats/quantile.html),
  [`summary()`](https://rdrr.io/r/base/summary.html),
  [`anyDuplicated()`](https://rdrr.io/r/base/duplicated.html),
  [`diff()`](https://rdrr.io/r/base/diff.html), `[`, and the six
  comparison operators. Most of them unwrap the weights and hand the
  work to the base or vctrs implementation; `[` instead delegates to the
  next method, except when the index is a matrix or an array. Either
  way, the concrete class’s own `vec_restore()` stays in charge of the
  metadata wherever the result is still a weight vector.

- Added
  [`ess.default()`](https://r-causal.github.io/causalgenerics/reference/ess.md),
  which computes the Kish effective sample size for any numeric vector.
  Weight vectors are numeric underneath, so a concrete weight class
  inherits a working
  [`ess()`](https://r-causal.github.io/causalgenerics/reference/ess.md)
  without writing a method. The default rejects input that is not
  numeric, including `NULL`, and returns `NaN` for numeric input with no
  information to summarize, such as a zero-length vector or weights that
  are all zero.

- Added
  [`causal_wts_ptype2()`](https://r-causal.github.io/causalgenerics/reference/causal_wts_ptype2.md),
  the coercion rule that a concrete weight class calls from its own
  `vec_ptype2()` method for two vectors of that class: they combine only
  when their estimands are identical, and otherwise the common type is a
  plain double. The condition signalled on the downgrade path stays with
  the caller.

- Added
  [`new_ipw()`](https://r-causal.github.io/causalgenerics/reference/new_ipw.md),
  the low-level constructor for the object every
  [`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.md)
  method returns, together with
  [`print()`](https://rdrr.io/r/base/print.html) and
  [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) methods
  for the `ipw` class. The field names and their order are now a
  cross-package contract, so an IPW estimate reads the same way
  whichever package produced it and a package supplying an
  [`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.md)
  method inherits both methods rather than writing its own.
  [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) gains
  an `exponentiate` argument that moves the `log(rr)` and `log(or)` rows
  to their natural scale, exponentiating the point estimate and the
  confidence limits only. Standard errors, z statistics, and p-values
  stay on the log scale, where the inference is done.

## causalgenerics 0.0.0.9000

- Initial release with the shared generics
  [`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.md),
  [`ess()`](https://r-causal.github.io/causalgenerics/reference/ess.md),
  [`is_causal_wt()`](https://r-causal.github.io/causalgenerics/reference/causal-weights.md),
  [`estimand()`](https://r-causal.github.io/causalgenerics/reference/causal-weights.md),
  and `estimand<-()`.
