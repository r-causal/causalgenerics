# causalgenerics (development version)

* `df.residual()` on an `ipw` result reports a fractional degrees of freedom as
  the double the fit reports, rather than truncating it to an integer. A
  penalized or smooth fit spends a fractional count of parameters, and
  truncating it would report a fit that spent more of them than it did. A whole
  number still comes back as an integer.

* The conditional reading pairs the corrected covariance with the outcome
  model's coefficients by name. A block a fitting package attached in the order
  of its own stacked system is reported in coefficient order, so the variance
  read beside a coefficient is that coefficient's, and `print()` writes each
  standard error against the coefficient it belongs to. A block whose labels
  cannot be paired with the coefficients is refused with an error of class
  `causalgenerics_conditional_vcov_mismatch`. A block of another size, one
  labelled with the parameter names of a stacked system, and a model whose
  coefficients carry no names are the three ways the pairing fails, and reading
  such a block by position would report the covariance of other parameters under
  this model's coefficient names.

* `print()` on a conditional result whose outcome model reports no coefficients
  says so in place of the table, rather than writing an empty one. The
  covariance is not asked for on that path, since a table with no rows has no
  use for one.

* `vcov()` on a model wrapped by `new_ipw_model()` raises an error of class
  `causalgenerics_no_vcov_ipw_model` when the wrapper no longer carries a
  covariance, rather than returning `NULL`. The class and the attribute go on
  together, so a model carrying one without the other passed through code that
  dropped its attributes and kept its class. `NULL` travelled from there:
  standard errors taken from it are an empty vector, and the limits built from
  those come back as `NA` with nothing said about why. `print()` on a
  conditional result whose wrapper lost its covariance refuses for the same
  reason, rather than printing the coefficients under a note telling the reader
  to wrap a model that is already wrapped.

* The error a missing covariance raises now says that the package which produced
  the object attaches one when it supports the `ipw_vcov` contract, in place of
  telling the reader to refit with a current version. The same helper answers
  for a result and for a component model of one, and a package that never
  attached a covariance is not one an upgrade would fix.

* `model.frame()` on an `ipw` result returns the outcome model's model frame, so
  prediction and averaging tooling can recover the data an estimate was computed
  from out of the result itself. The `(weights)` column a weighted frame carries
  is dropped, since a package that reads the columns of a frame as the variables
  a model was fitted on would otherwise treat the estimation weights as one of
  them. `weights()` reports those.

* `estimand()` on an `ipw` result returns the estimand the weights it was
  computed under targeted. There is no replacement method for a result:
  assigning a new estimand would relabel the estimates rather than recompute
  them.

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
