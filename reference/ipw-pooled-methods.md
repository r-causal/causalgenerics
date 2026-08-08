# Methods for a pooled inverse probability weighted result

The methods the class
[`pool_ipw()`](https://r-causal.github.io/causalgenerics/reference/pool_ipw.md)
returns carries.

- [`print()`](https://rdrr.io/r/base/print.html) summarizes the pooling
  and tabulates the pooled effects.

- [`coef()`](https://rdrr.io/r/stats/coef.html) returns the pooled
  estimates.

- [`vcov()`](https://rdrr.io/r/stats/vcov.html) returns their pooled
  covariance.

- [`confint()`](https://rdrr.io/r/stats/confint.html) returns their
  confidence limits.

- [`nobs()`](https://rdrr.io/r/stats/nobs.html) returns the smallest
  number of observations any pooled analysis was estimated from.

- [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) reports
  the pooled effects as a tidier-shaped table.

## Usage

``` r
# S3 method for class 'ipw_pooled'
print(x, ...)

# S3 method for class 'ipw_pooled'
coef(object, ...)

# S3 method for class 'ipw_pooled'
vcov(object, ...)

# S3 method for class 'ipw_pooled'
confint(object, parm, level = 0.95, ...)

# S3 method for class 'ipw_pooled'
nobs(object, ...)

# S3 method for class 'ipw_pooled'
as.data.frame(
  x,
  row.names = NULL,
  optional = FALSE,
  ...,
  conf.int = FALSE,
  conf.level = NULL,
  exponentiate = FALSE
)
```

## Arguments

- x:

  An `ipw_pooled` object.

- ...:

  Further arguments. These methods ignore them.

- object:

  An `ipw_pooled` object.

- parm:

  The rows to report an interval for, given either as the effect labels
  or as their positions. Missing means all of them.

- level:

  The confidence level. At the level a row stores, that row's own limits
  are returned; at any other level they are rebuilt from t on the row's
  degrees of freedom.

- row.names:

  A character vector of row names for the returned table, or `NULL` for
  the automatic ones.

- optional:

  Accepted for the
  [`base::as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html)
  generic. Every column of the table is named, so there is nothing for
  it to make optional.

- conf.int:

  If `TRUE`, append `conf.low` and `conf.high` columns after the rest of
  the table. Default is `FALSE`.

- conf.level:

  The confidence level the bounds report. `NULL`, the default, uses the
  level the result records.

- exponentiate:

  If `TRUE`, move the estimates that are on a log scale to their natural
  scale, as the section above describes. Default is `FALSE`.

## Value

[`print()`](https://rdrr.io/r/base/print.html) returns its input
invisibly.

[`coef()`](https://rdrr.io/r/stats/coef.html) returns a named numeric
vector of pooled estimates, named by effect label in the marginal
reading and by coefficient name in the conditional one.

[`vcov()`](https://rdrr.io/r/stats/vcov.html) returns the pooled
covariance of those estimates, with the same labels as dimnames on both
margins. A result whose analyses did not all carry a covariance records
none, and raises an error of class `causalgenerics_no_vcov_ipw_pooled`,
and of the general class `causalgenerics_no_vcov`, rather than returning
one built from the standard errors, which would report the effects as
uncorrelated.

[`confint()`](https://rdrr.io/r/stats/confint.html) returns a matrix
with one row per effect `parm` selects, in the order `parm` gives them,
and two columns holding the lower and upper limit, named by effect label
and by the two tail probabilities as percentages. A character `parm`
naming an effect the result does not report raises an error of class
`causalgenerics_invalid_argument`.

[`nobs()`](https://rdrr.io/r/stats/nobs.html) returns a single integer.

[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) returns a
plain data frame with the columns `term`, then `contrast` when the
result names contrasts, then `estimate`, `std.error`, `statistic`, `df`,
and `p.value`, with `conf.low` and `conf.high` appended when they are
asked for. The pooled covariance travels on it under the `ipw_vcov`
attribute unless the table was exponentiated.

## Details

These live here for the reason the `ipw` methods do. Two packages each
registering `print.ipw_pooled()` would collide in the shared S3 method
table, and a caller writing against a pooled result would then get
whichever package was installed last rather than the contract.

Everything they read is a field
[`pool_ipw()`](https://r-causal.github.io/causalgenerics/reference/pool_ipw.md)
wrote. They recompute no pooling, so a result read back out of a saved
file answers exactly as one just built does.

`ipw_pooled` deliberately does not inherit from `ipw`, and these methods
are deliberately separate from the ones registered against that class.
The two share column names and differ in what the numbers under them
mean: a pooled result reports several analyses rather than one fit and
refers its inference to t rather than to z, so a pooled result reaching
[`confint.ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw-accessors.md)
would come back with normal limits and nothing would say so.

## The degrees of freedom

Every pooled effect carries its own degrees of freedom, in the `df`
column of the `estimates` frame. They are the Barnard-Rubin adjusted
count, which depends on how much of that effect's variance came from
between the imputations, so two effects pooled from the same analyses
ordinarily differ in them.
[`confint()`](https://rdrr.io/r/stats/confint.html) and
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) refer
each row to t on its own count rather than to the normal or to a single
count for the table; with few imputations those counts can be in single
figures, where the difference is large.

There is deliberately no
[`df.residual()`](https://rdrr.io/r/stats/df.residual.html) method.
Residual degrees of freedom are a property of one fit, and a pooled
result has as many as it had imputations. The per-effect count above is
not a residual count and differs from row to row, so it is reported in
the column beside the interval and the p-value it was used for, and
[`df.residual()`](https://rdrr.io/r/stats/df.residual.html) on a pooled
result finds no method and no field and gives `NULL`.

## The confidence limits

[`confint()`](https://rdrr.io/r/stats/confint.html) returns the limits
the result stores for any row reported at the level asked for, and
rebuilds the rest. That is the rule
[`confint()`](https://rdrr.io/r/stats/confint.html) on an unpooled
result keeps, row by row: a stored pair need not be the one recomputing
gives, since it may have been rounded on its way into the frame.

[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) applies
the same rule to the frame as a whole, which is what
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) on an
unpooled result does. The stored pair comes back only when every row
records the level asked for, since a table mixing stored limits with
rebuilt ones would report two intervals under one pair of column
headings. `conf.level = NULL`, its default, names the level the frame
records, which is the level
[`pool_ipw()`](https://r-causal.github.io/causalgenerics/reference/pool_ipw.md)
built its limits at; a frame whose rows disagree records none, and the
limits are rebuilt at `0.95`.

## Exponentiating

On a marginal table `exponentiate = TRUE` means what it means for an
unpooled result. The rows labelled `log(rr)` and `log(or)` are matched
exactly, their estimate and their bounds move to the natural scale, and
the two terms are relabelled `rr` and `or`. The interval is settled
before the scale is, so a rebuilt bound is a t half width on the log
scale added to an estimate on the log scale and exponentiated
afterwards. The standard error, the statistic, and the p-value describe
the log scale and stay there, and the `ipw_vcov` attribute is dropped
rather than carried, since it would describe neither the table it sits
on nor anything else.

A conditional table has no rows labelled as ratios to pick out, so the
link the outcome models were fitted with settles the question for the
whole table: a `logit` link puts every coefficient on the log odds scale
and a `log` link puts every coefficient on the log risk scale, and both
are scales an exponential undoes. Every estimate moves and no term is
relabelled, since the terms are coefficient names and a coefficient does
not change its name with the scale its estimate is reported on. Every
other link raises an error of class `causalgenerics_exponentiate_link`,
and of the classes `causalgenerics_invalid_argument_exponentiate` and
`causalgenerics_invalid_argument`, rather than exponentiating
coefficients that describe nothing once exponentiated.

## See also

[`pool_ipw()`](https://r-causal.github.io/causalgenerics/reference/pool_ipw.md),
which produces these results, and
[`new_ipw()`](https://r-causal.github.io/causalgenerics/reference/new_ipw.md)
for the unpooled result they are pooled from.

## Examples

``` r
dat <- data.frame(
  x = rep(c(-1.5, -0.5, 0.5, 1.5), each = 5),
  z = rep(c(0, 1), 10),
  y = rep(c(0, 1, 1, 0, 1), 4)
)

imputed <- lapply(1:3, function(i) {
  completed <- dat
  completed$y[i] <- 1 - completed$y[i]
  completed
})

estimate <- list(c(0.20, 0.55), c(0.30, 0.60), c(0.40, 0.65))
std_err <- list(c(0.10, 0.25), c(0.20, 0.30), c(0.30, 0.35))

fits <- Map(
  function(completed, estimate, std.err) {
    new_ipw(
      estimand = "ate",
      wt_mod = glm(z ~ x, family = binomial(), data = completed),
      outcome_mod = glm(y ~ z, family = quasibinomial(), data = completed),
      estimates = data.frame(
        effect = c("rd", "log(rr)"),
        estimate = estimate,
        std.err = std.err,
        z = estimate / std.err,
        ci.lower = estimate - 1.96 * std.err,
        ci.upper = estimate + 1.96 * std.err,
        conf.level = 0.95,
        p.value = 2 * pnorm(-abs(estimate / std.err))
      ),
      se_method = "linearization",
      fit = NULL
    )
  },
  imputed,
  estimate,
  std_err
)

pooled <- pool_ipw(fits)

pooled
#> Pooled Inverse Probability Weight Estimator
#> Estimand: ATE 
#> Effects: marginal (population-averaged) 
#> Imputations: 3 
#> Complete-data df: 18 
#> 
#> Pooled marginal estimates:
#>         estimate std.err      t      df ci.lower ci.upper conf.level p.value  
#> rd       0.30000 0.24495 1.2247  9.6489 -0.24848  0.84848       0.95 0.24973  
#> log(rr)  0.60000 0.30822 1.9467 15.5637 -0.05489  1.25489       0.95 0.06986 .
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

coef(pooled)
#>      rd log(rr) 
#>     0.3     0.6 

# Wider than the normal limits would be, since the degrees of freedom are
# what few imputations leave.
confint(pooled)
#>              2.5 %    97.5 %
#> rd      -0.2484843 0.8484843
#> log(rr) -0.0548904 1.2548904

nobs(pooled)
#> [1] 20

# The tidier-shaped table, with the degrees of freedom the statistic beside
# it is referred to.
as.data.frame(pooled)
#>      term estimate std.error statistic        df    p.value
#> 1      rd      0.3 0.2449490  1.224745  9.648903 0.24972994
#> 2 log(rr)      0.6 0.3082207  1.946657 15.563733 0.06985665

# With an interval, and the ratio on its natural scale.
as.data.frame(pooled, conf.int = TRUE, exponentiate = TRUE)
#>   term estimate std.error statistic        df    p.value   conf.low conf.high
#> 1   rd 0.300000 0.2449490  1.224745  9.648903 0.24972994 -0.2484843 0.8484843
#> 2   rr 1.822119 0.3082207  1.946657 15.563733 0.06985665  0.9465889 3.5074540

# Residual degrees of freedom belong to one fit, so a pooled result has no
# method for them; the per-effect count is in the table above.
df.residual(pooled)
#> NULL
```
