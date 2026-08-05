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
coef(object, ...)

# S3 method for class 'ipw'
vcov(object, ...)

# S3 method for class 'ipw'
confint(object, parm, level = 0.95, ...)

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

- parm:

  The effects to report an interval for, given either as effect labels
  or as their positions in the result. Missing means all of them.

- level:

  The confidence level. At the level the result stores, its own limits
  are returned; at any other level they are recomputed.

## Value

[`coef()`](https://rdrr.io/r/stats/coef.html) returns a named numeric
vector of effect estimates, one element per row of the `estimates`
frame, named by effect label.

[`vcov()`](https://rdrr.io/r/stats/vcov.html) returns the square numeric
covariance matrix of those estimates, with the effect labels as dimnames
on both margins, exactly as the fitting package attached it. A result
that records no such matrix raises an error of class
`causalgenerics_no_vcov` rather than returning one built from the
standard errors, which would report the effects as uncorrelated.

[`confint()`](https://rdrr.io/r/stats/confint.html) returns a matrix
with one row per effect `parm` selects, in the order `parm` gives them,
and two columns holding the lower and upper limit. The rows are named by
effect label and the columns by the two tail probabilities as
percentages, the way the
[`confint()`](https://rdrr.io/r/stats/confint.html) methods in stats
name theirs. A character `parm` that names an effect the result does not
report raises an error of class `causalgenerics_invalid_argument`.

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
it, `outcome_mod`, and `fit`. They never branch on `se_method`, since
the fields they read already hold the result of whatever computation
that names, and they reach into `fit` only through ordinary S3 dispatch.
A package whose variance object is a bare list rather than a fitted
model therefore uses them unchanged.

Rows are named by the effect labels the
[`new_ipw()`](https://r-causal.github.io/causalgenerics/reference/new_ipw.md)
contract defines, so the name
[`coef()`](https://rdrr.io/r/stats/coef.html) gives an estimate is the
one [`vcov()`](https://rdrr.io/r/stats/vcov.html) and
[`confint()`](https://rdrr.io/r/stats/confint.html) use for it and the
one [`print()`](https://rdrr.io/r/base/print.html) labels its row with.

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
```
