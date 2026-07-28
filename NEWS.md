# causalgenerics 0.0.0.9001

* The weighting argument of `ipw()` is now named `wt_mod`, which reads for a
  weighting object of any kind rather than only a propensity score model. A call
  that names it `ps_mod`, the spelling propensity 0.1.0 released, still works: it
  warns once per session and then behaves in every respect as though the object
  had been passed as `wt_mod`, dispatch included. `ps_mod` is not an argument of
  the generic, and new code should use `wt_mod`.

* Added `new_causal_wts()`, the low-level constructor for the abstract
  `causal_wts` weight class that concrete weight classes in the ecosystem are
  built on.

* Added `causal_wts` methods for `is_causal_wt()`, `estimand()`, and
  `estimand<-()`, so a concrete class inherits all three.

* Added `causal_wts` methods for the read-only operations that name no metadata
  field and so give the answer the underlying double would: `vec_math()`, the
  `Summary` group generic, `min()`, `max()`, `range()`, `median()`, `quantile()`,
  `summary()`, `anyDuplicated()`, `diff()`, `[`, and the six comparison
  operators. Most of them unwrap the weights and hand the work to the base or
  vctrs implementation; `[` instead delegates to the next method, except when
  the index is a matrix or an array. Either way, the concrete class's own
  `vec_restore()` stays in charge of the metadata wherever the result is still a
  weight vector.

* Added `ess.default()`, which computes the Kish effective sample size for any
  numeric vector. Weight vectors are numeric underneath, so a concrete weight
  class inherits a working `ess()` without writing a method. The default rejects
  input that is not numeric, including `NULL`, and returns `NaN` for numeric
  input with no information to summarize, such as a zero-length vector or
  weights that are all zero.

* Added `causal_wts_ptype2()`, the coercion rule that a concrete weight class
  calls from its own `vec_ptype2()` method for two vectors of that class: they
  combine only when their estimands are identical, and otherwise the common type
  is a plain double. The condition signalled on the downgrade path stays with
  the caller.

* Added `new_ipw()`, the low-level constructor for the object every `ipw()`
  method returns, together with `print()` and `as.data.frame()` methods for the
  `ipw` class. The field names and their order are now a cross-package contract,
  so an IPW estimate reads the same way whichever package produced it and a
  package supplying an `ipw()` method inherits both methods rather than writing
  its own. `as.data.frame()` gains an `exponentiate` argument that moves the
  `log(rr)` and `log(or)` rows to their natural scale, exponentiating the point
  estimate and the confidence limits only. Standard errors, z statistics, and
  p-values stay on the log scale, where the inference is done.

# causalgenerics 0.0.0.9000

* Initial release with the shared generics `ipw()`, `ess()`, `is_causal_wt()`, `estimand()`, and `estimand<-()`.
