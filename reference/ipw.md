# Inverse probability weighted estimation

`ipw()` is the generic for bring-your-own-model inverse probability
weighted estimation of causal effects. A method takes a fitted weighting
or propensity score model together with a fitted weighted outcome model
and returns causal effect estimates with standard errors that account
for the two-step estimation process.

## Usage

``` r
ipw(wt_mod, outcome_mod, ...)
```

## Arguments

- wt_mod:

  The weighting object that produced the weights, for example a fitted
  propensity score model. `ipw()` dispatches on this argument.

- outcome_mod:

  A fitted weighted outcome model.

- ...:

  Arguments passed to methods.

## Value

An object of class `ipw`, holding the causal effect estimates and their
standard errors alongside the models they came from. See
[`new_ipw()`](https://r-causal.github.io/causalgenerics/reference/new_ipw.md)
for the components and their order.

## Details

This package defines the generic and the shared result class its methods
return. The methods themselves live in the packages that own the
relevant model classes, such as `propensity` for propensity score models
and `balancing` for balancing weight fits. Dispatch is on `wt_mod`, the
weighting object; each method documents the outcome model classes and
further arguments it accepts.

A method builds its return value with
[`new_ipw()`](https://r-causal.github.io/causalgenerics/reference/new_ipw.md)
rather than a result object of its own. The field names and their order
are a cross-package contract, so that an IPW estimate reads the same way
whichever package produced it, and constructing through
[`new_ipw()`](https://r-causal.github.io/causalgenerics/reference/new_ipw.md)
is also what gives a method the shared
[`print()`](https://rdrr.io/r/base/print.html) and
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) methods.

## See also

[`new_ipw()`](https://r-causal.github.io/causalgenerics/reference/new_ipw.md),
the constructor every method returns through, and the `propensity` and
`balancing` packages for methods.

## Examples

``` r
# A confounded toy data set: the confounder `x` raises both the chance of
# exposure `z` and the risk of outcome `y`. Within either level of `x` the
# risk difference is 0.25.
dat <- data.frame(
  x = rep(c(0, 1), each = 10),
  z = c(0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0),
  y = c(1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1, 0)
)

# The crude comparison ignores `x` and overstates that risk difference.
coef(lm(y ~ z, data = dat))[["z"]]
#> [1] 0.4

# A weighting object of a toy class, standing in for the fitted weighting
# models that the methods in other packages dispatch on. Weighting by the
# inverse probability of the exposure received breaks the dependence of `z`
# on `x`, so the weighted outcome model recovers the 0.25 risk difference.
ps <- fitted(glm(z ~ x, family = binomial(), data = dat))
wt_mod <- structure(
  list(wts = ifelse(dat$z == 1, 1 / ps, 1 / (1 - ps))),
  class = "cg_toy_model"
)

# A method for that class. It reads the exposure coefficient of the weighted
# outcome model as a risk difference and returns through `new_ipw()`, which is
# what every `ipw()` method does. Inference is normal-based throughout, as the
# `z` column of the estimates contract calls for. The standard error is the
# model-based one the weighted fit reports; a real method uses a variance
# estimator that also accounts for the weights having been estimated.
ipw.cg_toy_model <- function(wt_mod, outcome_mod, ...) {
  coefs <- summary(outcome_mod)$coefficients["z", ]
  estimate <- coefs[["Estimate"]]
  std_err <- coefs[["Std. Error"]]
  z <- estimate / std_err
  half_width <- qnorm(0.975) * std_err
  new_ipw(
    estimand = "ate",
    wt_mod = wt_mod,
    outcome_mod = outcome_mod,
    estimates = data.frame(
      effect = "rd",
      estimate = estimate,
      std.err = std_err,
      z = z,
      ci.lower = estimate - half_width,
      ci.upper = estimate + half_width,
      conf.level = 0.95,
      p.value = 2 * pnorm(-abs(z))
    ),
    se_method = "model-based",
    fit = NULL
  )
}

outcome_mod <- lm(y ~ z, data = dat, weights = wt_mod$wts)

# Dispatch on the weighting object finds the method.
ipw(wt_mod, outcome_mod)
#> Inverse Probability Weight Estimator
#> Estimand: ATE 
#> Effects: marginal (population-averaged) 
#> 
#> Weight Estimator:
#>   Call: <cg_toy_model> 
#> 
#> Outcome Model:
#>   Call: lm(formula = y ~ z, data = dat, weights = wt_mod$wts) 
#> 
#> Marginal estimates:
#>    estimate std.err      z ci.lower ci.upper conf.level p.value
#> rd  0.25000 0.22822 1.0954  -0.1973   0.6973       0.95  0.2733
```
