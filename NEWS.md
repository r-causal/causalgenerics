# causalgenerics (development version)

* `print()` on an `ipw` result names the reading it is showing, on an `Effects:`
  line beside the estimand and in the heading of the table. The marginal
  reading tabulates the effect estimates as before, under `Marginal estimates:`.
  The conditional reading tabulates the outcome model's coefficients with the
  standard errors the corrected covariance implies, and prints the coefficients
  on their own, with a note, when the outcome model records no such covariance.

* `coef()`, `vcov()`, and `confint()` on an `ipw` result report the reading the
  `effects` field records, and each gains an `effects` argument that names a
  reading for one call without changing the result. The conditional reading
  reports the outcome model's coefficient surface: its coefficients, the
  covariance the joint estimation of the weights and the outcome implies, and
  the normal limits built from that covariance. `nobs()`, `df.residual()`, and
  `weights()` describe the fit rather than a surface of it and answer the same
  way in either reading.

* `vcov()` and `confint()` refuse the conditional reading when the outcome model
  carries no corrected covariance, with an error of class
  `causalgenerics_no_conditional_vcov`. The covariance a weighted model computes
  for itself treats the estimated weights as fixed, so it is not a substitute; a
  fitting package supplies the corrected block by wrapping the outcome model
  with `new_ipw_model()`.

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
