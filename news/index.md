# Changelog

## causalgenerics (development version)

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
