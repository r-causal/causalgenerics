# Construct an inverse probability weighted result

`new_ipw()` is the low-level constructor for the object that every
[`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.md)
method returns. It is intended for package developers writing an
[`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.md)
method, not for end users, and it assumes its arguments are already
validated.

## Usage

``` r
new_ipw(estimand, wt_mod, outcome_mod, estimates, se_method, fit)

# S3 method for class 'ipw'
print(x, ...)

# S3 method for class 'ipw'
as.data.frame(x, row.names = NULL, optional = NULL, exponentiate = FALSE, ...)
```

## Arguments

- estimand:

  The causal estimand the method targeted, such as `"ate"` or `"att"`.

- wt_mod:

  The weighting object: the fitted model that produced the weights.

- outcome_mod:

  The fitted weighted outcome model.

- estimates:

  A data frame of effect estimates, in the shape the return value
  describes.

- se_method:

  The standard error method that ran, such as `"mestimation"` or
  `"linearization"`.

- fit:

  The fitted variance object, or `NULL` when the method has none.

- x:

  An `ipw` object.

- ...:

  Further arguments. [`print()`](https://rdrr.io/r/base/print.html)
  ignores them;
  [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) passes
  them to
  [`base::as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html).

- row.names, optional:

  Passed to
  [`base::as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html).

- exponentiate:

  If `TRUE`, exponentiate the log risk ratio and log odds ratio to
  produce risk ratios and odds ratios on their natural scale. The
  confidence interval bounds are also exponentiated. Standard errors, z
  statistics, and p-values remain on the log scale. Default is `FALSE`.

## Value

`new_ipw()` returns an S3 object of class `ipw`: a list of the following
six components, in this order.

- `estimand`:

  The causal estimand, such as `"ate"` or `"att"`.

- `wt_mod`:

  The weighting object: the fitted model that produced the weights.

- `outcome_mod`:

  The fitted outcome model.

- `estimates`:

  A data frame with one row per effect measure and the following
  columns: `effect` (the measure name), `estimate` (point estimate),
  `std.err` (standard error), `z` (z-statistic), `ci.lower` and
  `ci.upper` (confidence interval bounds), `conf.level`, and `p.value`.
  For a categorical exposure the data frame also has a `comparison`
  column, placed after `effect`, naming the non-reference level and
  reference level of each contrast.

- `se_method`:

  The standard error method used, such as `"mestimation"` or
  `"linearization"`.

- `fit`:

  The fitted object the variance estimator produced, or `NULL`. A method
  that stacks estimating equations records the M-estimator here; the
  linearization path has no such object and records `NULL`.

[`print()`](https://rdrr.io/r/base/print.html) returns its input
invisibly.
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) returns
the `estimates` component as a data frame.

## Details

The result layer is shared so that an IPW estimate reads the same way
whichever package produced it. A package supplying an
[`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.md)
method builds its return here and inherits the
[`print()`](https://rdrr.io/r/base/print.html) and
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) methods
registered against the class rather than writing its own. Two packages
each defining `print.ipw()` would collide in the shared S3 method table,
which is the situation this package exists to prevent.

The field names and their order are part of the contract, since callers
read fields by name and print the object positionally. `fit` is present
on every path, including the ones that have no fitted variance object to
report.

[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) returns
the `estimates` component. With `exponentiate = TRUE` it moves the
`log(rr)` and `log(or)` rows to their natural scale, exponentiating the
point estimate and the confidence limits and relabelling the two effects
`"rr"` and `"or"`. Standard errors, z statistics, and p-values stay on
the log scale, where the inference is done.

## The effect labels

Every row of `estimates` has a label, and it is the label rather than
the position that [`print()`](https://rdrr.io/r/base/print.html) writes
down the side of its table and that
[`coef()`](https://r-causal.github.io/causalgenerics/reference/ipw-accessors.md),
[`vcov()`](https://r-causal.github.io/causalgenerics/reference/ipw-accessors.md),
and
[`confint()`](https://r-causal.github.io/causalgenerics/reference/ipw-accessors.md)
name their results with. The label is the `effect` column on its own
when there is no `comparison` column, and `effect` and `comparison`
pasted together, such as `"rd b vs a"`, when there is. A categorical
exposure repeats each effect measure across its contrasts, so `effect`
alone would name several rows the same thing.

## The covariance of the effects

A method that can compute the covariance of the effects it reports
attaches it to the `estimates` data frame as an attribute named
`ipw_vcov`. The value is a square numeric matrix whose row order is the
row order of `estimates` and whose dimnames on both margins are the
effect labels above.
[`vcov()`](https://r-causal.github.io/causalgenerics/reference/ipw-accessors.md)
reads that attribute and raises an error when it is absent, so a method
that has no covariance to report attaches none rather than a substitute:
the standard errors in `estimates` give the diagonal of the matrix and
say nothing about the off-diagonal entries, which are not zero for
effects estimated from the same weighted means.

## See also

[`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.md),
the generic these results come from.

## Examples

``` r
dat <- data.frame(
  x = rep(c(-1.5, -0.5, 0.5, 1.5), each = 5),
  z = rep(c(0, 1), 10),
  y = rep(c(0, 1, 1, 0, 1), 4)
)

# Written out literally, in the shape the return contract documents. These
# stand in for what an `ipw()` method would compute from the models below.
estimates <- data.frame(
  effect = c("rd", "log(rr)", "log(or)"),
  estimate = c(0.199882, 0.560414, 0.878313),
  std.err = c(0.092425, 0.273519, 0.418661),
  z = c(2.1626, 2.0489, 2.0979),
  ci.lower = c(0.018732, 0.024326, 0.057753),
  ci.upper = c(0.381032, 1.096502, 1.698873),
  conf.level = 0.95,
  p.value = c(0.030570, 0.040470, 0.035910)
)

res <- new_ipw(
  estimand = "ate",
  wt_mod = glm(z ~ x, family = binomial(), data = dat),
  outcome_mod = glm(y ~ z, family = quasibinomial(), data = dat),
  estimates = estimates,
  se_method = "linearization",
  fit = NULL
)

res
#> Inverse Probability Weight Estimator
#> Estimand: ATE 
#> 
#> Weight Estimator:
#>   Call: glm(formula = z ~ x, family = binomial(), data = dat) 
#> 
#> Outcome Model:
#>   Call: glm(formula = y ~ z, family = quasibinomial(), data = dat) 
#> 
#> Estimates:
#>         estimate  std.err      z ci.lower ci.upper conf.level p.value  
#> rd      0.199882 0.092425 2.1626 0.018732  0.38103       0.95 0.03057 *
#> log(rr) 0.560414 0.273519 2.0489 0.024326  1.09650       0.95 0.04047 *
#> log(or) 0.878313 0.418661 2.0979 0.057753  1.69887       0.95 0.03591 *
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# The ratios on their natural scale.
as.data.frame(res, exponentiate = TRUE)
#>   effect estimate  std.err      z ci.lower ci.upper conf.level p.value
#> 1     rd 0.199882 0.092425 2.1626 0.018732 0.381032       0.95 0.03057
#> 2     rr 1.751397 0.273519 2.0489 1.024624 2.993676       0.95 0.04047
#> 3     or 2.406836 0.418661 2.0979 1.059453 5.467782       0.95 0.03591
```
