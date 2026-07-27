# causalgenerics 0.0.0.9001

* Added `new_causal_wts()`, the low-level constructor for the abstract
  `causal_wts` weight class that concrete weight classes in the ecosystem are
  built on.

* Added `causal_wts` methods for `is_causal_wt()`, `estimand()`, and
  `estimand<-()`, so a concrete class inherits all three.

* Added `causal_wts` methods for the read-only operations that unwrap to the
  underlying double: `vec_math()`, the `Summary` group generic, `min()`,
  `max()`, `range()`, `median()`, `quantile()`, `summary()`, `anyDuplicated()`,
  `diff()`, `[`, and the six comparison operators.

# causalgenerics 0.0.0.9000

* Initial release with the shared generics `ipw()`, `ess()`, `is_causal_wt()`, `estimand()`, and `estimand<-()`.
