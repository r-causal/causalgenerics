# causalgenerics 0.0.0.9001

* Added `new_causal_wts()`, the low-level constructor for the abstract
  `causal_wts` weight class that concrete weight classes in the ecosystem are
  built on.

* Added `causal_wts` methods for `is_causal_wt()`, `estimand()`, and
  `estimand<-()`, so a concrete class inherits all three.

* Added `causal_wts` methods for the read-only operations that name no metadata
  field and so delegate to the underlying double: `vec_math()`, the `Summary`
  group generic, `min()`, `max()`, `range()`, `median()`, `quantile()`,
  `summary()`, `anyDuplicated()`, `diff()`, `[`, and the six comparison
  operators. Each unwraps the weights and hands the work to the base or vctrs
  implementation, which leaves the concrete class's own `vec_restore()` in
  charge of the metadata wherever the result is still a weight vector.

* Added `causal_wts_ptype2()`, the coercion rule that a concrete weight class
  calls from its own `vec_ptype2()` method for two vectors of that class: they
  combine only when their estimands are identical, and otherwise the common type
  is a plain double. The condition signalled on the downgrade path stays with
  the caller.

# causalgenerics 0.0.0.9000

* Initial release with the shared generics `ipw()`, `ess()`, `is_causal_wt()`, `estimand()`, and `estimand<-()`.
