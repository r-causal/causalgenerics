# Changelog

## causalgenerics (development version)

- [`print()`](https://rdrr.io/r/base/print.html) on an `ipw` result
  names the reading it is showing, on an `Effects:` line beside the
  estimand and in the heading of the table. The marginal reading
  tabulates the effect estimates as before, under `Marginal estimates:`.
  The conditional reading tabulates the outcome model’s coefficients
  with the standard errors the corrected covariance implies, and prints
  the coefficients on their own, with a note, when the outcome model
  records no such covariance.

- [`coef()`](https://rdrr.io/r/stats/coef.html),
  [`vcov()`](https://rdrr.io/r/stats/vcov.html), and
  [`confint()`](https://rdrr.io/r/stats/confint.html) on an `ipw` result
  report the reading the `effects` field records, and each gains an
  `effects` argument that names a reading for one call without changing
  the result. The conditional reading reports the outcome model’s
  coefficient surface: its coefficients, the covariance the joint
  estimation of the weights and the outcome implies, and the normal
  limits built from that covariance.
  [`nobs()`](https://rdrr.io/r/stats/nobs.html),
  [`df.residual()`](https://rdrr.io/r/stats/df.residual.html), and
  [`weights()`](https://rdrr.io/r/stats/weights.html) describe the fit
  rather than a surface of it and answer the same way in either reading.

- [`vcov()`](https://rdrr.io/r/stats/vcov.html) and
  [`confint()`](https://rdrr.io/r/stats/confint.html) refuse the
  conditional reading when the outcome model carries no corrected
  covariance, with an error of class
  `causalgenerics_no_conditional_vcov`. The covariance a weighted model
  computes for itself treats the estimated weights as fixed, so it is
  not a substitute; a fitting package supplies the corrected block by
  wrapping the outcome model with
  [`new_ipw_model()`](https://r-causal.github.io/causalgenerics/reference/new_ipw_model.md).

- [`new_ipw()`](https://r-causal.github.io/causalgenerics/reference/new_ipw.md)
  gains an `effects` argument, stored as a seventh field, recording
  which reading a result presents: `"marginal"` for the causal contrast
  estimates or `"conditional"` for the outcome model’s coefficient
  surface. A method that names no mode reports marginal effects, as
  every method did before.

- New generics
  [`as_marginal()`](https://r-causal.github.io/causalgenerics/reference/ipw-modes.md)
  and
  [`as_conditional()`](https://r-causal.github.io/causalgenerics/reference/ipw-modes.md)
  move a result between the two readings. Both surfaces exist on every
  result, so the methods on `ipw` set the mode and leave every other
  field alone.

- `ipw` results gain the accessors a fitted model answers to:
  [`coef()`](https://rdrr.io/r/stats/coef.html),
  [`vcov()`](https://rdrr.io/r/stats/vcov.html),
  [`confint()`](https://rdrr.io/r/stats/confint.html),
  [`nobs()`](https://rdrr.io/r/stats/nobs.html),
  [`df.residual()`](https://rdrr.io/r/stats/df.residual.html), and
  [`weights()`](https://rdrr.io/r/stats/weights.html). A method that can
  compute the covariance of the effects it reports now attaches it to
  its `estimates` data frame as the `ipw_vcov` attribute, which is what
  [`vcov()`](https://rdrr.io/r/stats/vcov.html) reads.

- New
  [`new_ipw_model()`](https://r-causal.github.io/causalgenerics/reference/new_ipw_model.md)
  wraps a component model of an IPW fit with the covariance the two-step
  estimation implies, so that
  [`vcov()`](https://rdrr.io/r/stats/vcov.html) on the model accounts
  for the weights having been estimated rather than fixed.

## causalgenerics 0.1.0

- Initial CRAN submission.
