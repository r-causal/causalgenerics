# Model accessors for an inverse probability weighted result

The accessors a fitted model answers to, defined for the class
[`new_ipw()`](https://r-causal.github.io/causalgenerics/reference/new_ipw.md)
constructs.

- [`coef()`](https://rdrr.io/r/stats/coef.html) returns the effect
  estimates.

- [`vcov()`](https://rdrr.io/r/stats/vcov.html) returns the covariance
  of those estimates.

- [`confint()`](https://rdrr.io/r/stats/confint.html) returns their
  confidence limits.

- [`nobs()`](https://rdrr.io/r/stats/nobs.html) returns the number of
  observations the outcome model was fitted on.

- [`df.residual()`](https://rdrr.io/r/stats/df.residual.html) returns
  the residual degrees of freedom of the fitted variance object.

- [`weights()`](https://rdrr.io/r/stats/weights.html) returns the
  weights the outcome model was fitted with.

## Usage

``` r
# S3 method for class 'ipw'
coef(object, ..., effects = NULL)

# S3 method for class 'ipw'
vcov(object, ..., effects = NULL)

# S3 method for class 'ipw'
confint(object, parm, level = 0.95, ..., effects = NULL)

# S3 method for class 'ipw'
nobs(object, ...)

# S3 method for class 'ipw'
df.residual(object, ...)

# S3 method for class 'ipw'
weights(object, ...)
```

## Arguments

- object:

  An `ipw` object.

- ...:

  Further arguments. These methods ignore them.

- effects:

  The reading to report, either `"marginal"` or `"conditional"`. `NULL`,
  the default, reports the reading the result records; any other value
  overrides it for the one call and leaves the result as it is.

- parm:

  The rows to report an interval for, given either as the labels of the
  reported surface or as their positions in it. Missing means all of
  them.

- level:

  The confidence level. At the level the result stores, its own limits
  are returned; at any other level they are recomputed.

## Value

[`coef()`](https://rdrr.io/r/stats/coef.html) returns a named numeric
vector of effect estimates, one element per row of the `estimates`
frame, named by effect label. In the conditional reading it returns the
outcome model's coefficients as that model reports them, names included.

[`vcov()`](https://rdrr.io/r/stats/vcov.html) returns the square numeric
covariance matrix of those estimates, with the effect labels as dimnames
on both margins, exactly as the fitting package attached it. A result
that records no such matrix raises an error of class
`causalgenerics_no_vcov` rather than returning one built from the
standard errors, which would report the effects as uncorrelated. In the
conditional reading it returns the corrected covariance of the outcome
model's coefficients, and an outcome model that carries none raises an
error of class `causalgenerics_no_conditional_vcov`.

[`confint()`](https://rdrr.io/r/stats/confint.html) returns a matrix
with one row per effect `parm` selects, in the order `parm` gives them,
and two columns holding the lower and upper limit. The rows are named by
effect label and the columns by the two tail probabilities as
percentages, the way the
[`confint()`](https://rdrr.io/r/stats/confint.html) methods in stats
name theirs. A character `parm` that names an effect the result does not
report raises an error of class `causalgenerics_invalid_argument`. In
the conditional reading the rows are the outcome model's coefficients,
which `parm` names and indexes in the same two ways, and the limits are
the normal ones built from the corrected covariance at every level,
since the limits the result stores belong to the effects the marginal
reading reports.

[`coef()`](https://rdrr.io/r/stats/coef.html),
[`vcov()`](https://rdrr.io/r/stats/vcov.html), and
[`confint()`](https://rdrr.io/r/stats/confint.html) raise an error of
class `causalgenerics_invalid_argument_effects` when `effects` names
neither reading, and when the result's own field does.

[`nobs()`](https://rdrr.io/r/stats/nobs.html) returns a single integer,
the number of observations the outcome model was fitted on, which is the
number the estimates were computed from and not necessarily the number
the weighting model saw.

[`df.residual()`](https://rdrr.io/r/stats/df.residual.html) returns a
single integer, or `NA_integer_` when the result records no fitted
variance object or that object reports no residual degrees of freedom.

[`weights()`](https://rdrr.io/r/stats/weights.html) returns the outcome
model's weights as its model frame stores them, so a concrete weight
class such as `psw` comes back as itself, or `NULL` when the outcome
model was fitted unweighted.

## Details

These methods live here for the reason
[`print()`](https://rdrr.io/r/base/print.html) does. Two packages each
registering `coef.ipw()` would collide in the shared S3 method table,
and a caller writing against a result would then get whichever package
was installed last rather than the contract.

Everything they read is part of the
[`new_ipw()`](https://r-causal.github.io/causalgenerics/reference/new_ipw.md)
contract: the `estimates` frame, the `ipw_vcov` attribute attached to
it, `outcome_mod`, `fit`, and the `effects` field the section below
describes. They never branch on `se_method`, since the fields they read
already hold the result of whatever computation that names, and they
reach into `fit` only through ordinary S3 dispatch. A package whose
variance object is a bare list rather than a fitted model therefore uses
them unchanged.

Rows are named by the effect labels the
[`new_ipw()`](https://r-causal.github.io/causalgenerics/reference/new_ipw.md)
contract defines, so the name
[`coef()`](https://rdrr.io/r/stats/coef.html) gives an estimate is the
one [`vcov()`](https://rdrr.io/r/stats/vcov.html) and
[`confint()`](https://rdrr.io/r/stats/confint.html) use for it and the
one [`print()`](https://rdrr.io/r/base/print.html) labels its row with.

## The reading these methods report

[`coef()`](https://rdrr.io/r/stats/coef.html),
[`vcov()`](https://rdrr.io/r/stats/vcov.html), and
[`confint()`](https://rdrr.io/r/stats/confint.html) report the surface
the result's `effects` field names, and take an `effects` argument that
names one for a single call. The marginal reading is the causal contrast
estimates, which is what these methods reported before the field existed
and what a result that records no mode still reports. The conditional
reading is the outcome model's coefficient surface:
[`coef()`](https://rdrr.io/r/stats/coef.html) gives the model's
coefficients, [`vcov()`](https://rdrr.io/r/stats/vcov.html) gives their
covariance, and [`confint()`](https://rdrr.io/r/stats/confint.html)
bounds them.

The covariance the conditional reading reports is the block of the joint
estimation that the fitting package attached with
[`new_ipw_model()`](https://r-causal.github.io/causalgenerics/reference/new_ipw_model.md),
never the one the outcome model computed for itself. A model fitted as
though its weights were fixed understates its uncertainty, because the
weights were estimated from the same data, so there is nothing to fall
back on when no corrected block is there.
[`vcov()`](https://rdrr.io/r/stats/vcov.html) and
[`confint()`](https://rdrr.io/r/stats/confint.html) raise an error of
class `causalgenerics_no_conditional_vcov` instead, which the package
that produced the result answers by wrapping the outcome model with
[`new_ipw_model()`](https://r-causal.github.io/causalgenerics/reference/new_ipw_model.md)
before it builds the result.
[`coef()`](https://rdrr.io/r/stats/coef.html) needs no such block and
reports the coefficients either way.

[`nobs()`](https://rdrr.io/r/stats/nobs.html),
[`df.residual()`](https://rdrr.io/r/stats/df.residual.html), and
[`weights()`](https://rdrr.io/r/stats/weights.html) describe the fit
rather than a surface of it, so they answer the same way in either
reading and take no `effects` argument.

## See also

[`new_ipw()`](https://r-causal.github.io/causalgenerics/reference/new_ipw.md)
for the result class and the fields these methods read.

## Examples

``` r
dat <- data.frame(
  x = rep(c(-1.5, -0.5, 0.5, 1.5), each = 5),
  z = rep(c(0, 1), 10),
  y = rep(c(0, 1, 1, 0, 1), 4)
)

ps <- fitted(glm(z ~ x, family = binomial(), data = dat))
wts <- ifelse(dat$z == 1, 1 / ps, 1 / (1 - ps))

# Written out literally, in the shape the `ipw()` return contract documents.
estimates <- data.frame(
  effect = c("rd", "log(rr)"),
  estimate = c(0.199882, 0.560414),
  std.err = c(0.092425, 0.273519),
  z = c(2.1626, 2.0489),
  ci.lower = c(0.018732, 0.024326),
  ci.upper = c(0.381032, 1.096502),
  conf.level = 0.95,
  p.value = c(0.030570, 0.040470)
)

# The covariance of the two effects, which a method attaches to the estimates
# it returns. Both are computed from the same weighted means, so the
# off-diagonal entry is far from zero.
attr(estimates, "ipw_vcov") <- matrix(
  c(0.008542, 0.022753, 0.022753, 0.074813),
  nrow = 2,
  dimnames = list(estimates$effect, estimates$effect)
)

res <- new_ipw(
  estimand = "ate",
  wt_mod = glm(z ~ x, family = binomial(), data = dat),
  outcome_mod = glm(y ~ z, family = quasibinomial(), data = dat, weights = wts),
  estimates = estimates,
  se_method = "linearization",
  fit = NULL
)

coef(res)
#>       rd  log(rr) 
#> 0.199882 0.560414 
vcov(res)
#>               rd  log(rr)
#> rd      0.008542 0.022753
#> log(rr) 0.022753 0.074813

# At the stored level the stored limits come back.
confint(res)
#>            2.5 %   97.5 %
#> rd      0.018732 0.381032
#> log(rr) 0.024326 1.096502

# At any other level they are recomputed.
confint(res, parm = "rd", level = 0.9)
#>          5 %      95 %
#> rd 0.0478564 0.3519076

nobs(res)
#> [1] 20

# The linearization path records no fitted variance object.
df.residual(res)
#> [1] NA

head(weights(res))
#> [1] 1.785796 2.272594 1.785796 2.272594 1.785796 2.083669

# The conditional reading reports the outcome model's coefficient surface
# rather than the effects.
coef(res, effects = "conditional")
#>  (Intercept)            z 
#> 3.919028e-01 6.859992e-16 

# Its covariance is the corrected block a fitting package attaches to the
# outcome model. This one carries none, and the covariance the model computed
# for itself is not a substitute: it treats the estimated weights as fixed.
try(vcov(res, effects = "conditional"))
#> Error in vcov.ipw(res, effects = "conditional") : 
#>   The conditional reading reports the covariance the joint estimation of the weights and the outcome implies, and this result's outcome model records none; the package that produced the result attaches one by wrapping the model with `new_ipw_model()`.

outcome_mod <- glm(y ~ z, family = quasibinomial(), data = dat, weights = wts)

# Standing in for the outcome block of a stacked sandwich, which a method
# computes from the estimating equations of both steps.
corrected <- vcov(outcome_mod) * 1.4

conditional <- as_conditional(new_ipw(
  estimand = "ate",
  wt_mod = glm(z ~ x, family = binomial(), data = dat),
  outcome_mod = new_ipw_model(outcome_mod, corrected),
  estimates = estimates,
  se_method = "linearization",
  fit = NULL
))

# A result that records the conditional reading answers in it with nothing
# named at the call site.
coef(conditional)
#>  (Intercept)            z 
#> 3.919028e-01 6.859992e-16 
vcov(conditional)
#>             (Intercept)         z
#> (Intercept)    0.646421 -0.646421
#> z             -0.646421  1.292842
confint(conditional)
#>                 2.5 %   97.5 %
#> (Intercept) -1.183914 1.967720
#> z           -2.228542 2.228542
```
