# causalgenerics (development version)

* `new_ipw()` gains an `effects` argument, stored as a seventh field, recording
  which reading a result presents: `"marginal"` for the causal contrast
  estimates or `"conditional"` for the outcome model's coefficient surface. A
  method that names no mode reports marginal effects, as every method did
  before.

* New generics `as_marginal()` and `as_conditional()` move a result between the
  two readings. Both surfaces exist on every result, so the methods on `ipw`
  set the mode and leave every other field alone.

* `ipw` results gain the accessors a fitted model answers to: `coef()`,
  `vcov()`, `confint()`, `nobs()`, `df.residual()`, and `weights()`. A method
  that can compute the covariance of the effects it reports now attaches it to
  its `estimates` data frame as the `ipw_vcov` attribute, which is what `vcov()`
  reads.

* New `new_ipw_model()` wraps a component model of an IPW fit with the
  covariance the two-step estimation implies, so that `vcov()` on the model
  accounts for the weights having been estimated rather than fixed.

# causalgenerics 0.1.0

* Initial CRAN submission.
