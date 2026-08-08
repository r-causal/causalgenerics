# Pool inverse probability weighted results over multiply imputed data

`pool_ipw()` combines the
[`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.md)
results fitted to each of a set of multiply imputed datasets into a
single result, by Rubin's rules. It takes a plain list of results or the
`mira` object `mice::with()` returns, and gives back an object of class
`ipw_pooled`.

## Usage

``` r
pool_ipw(fits, ..., effects = NULL, dfcom = NULL, conf_level = NULL)
```

## Arguments

- fits:

  The results to pool: a list of `ipw` objects, one per imputed dataset,
  or a `mira`. At least two are needed, since the between-imputation
  variance is estimated from the spread across them.

- ...:

  These dots exist so that every argument after `fits` is matched by
  name. Passing anything through them is an error.

- effects:

  The reading to pool, either `"marginal"` or `"conditional"`. `NULL`,
  the default, pools the reading the results record. The marginal
  reading pools the causal contrast estimates; the conditional reading
  pools the outcome models' coefficients, with the standard errors
  implied by the corrected covariance each one carries.

- dfcom:

  The complete-data degrees of freedom. `NULL`, the default, reads them
  off the results as the section above describes.

- conf_level:

  The level the pooled bounds report. `NULL`, the default, uses the
  level the results stored, or `0.95` when they store none.

## Value

An S3 object of class `ipw_pooled`: a list of the following nine
components, in this order.

- `estimand`:

  The causal estimand every pooled result targeted.

- `estimates`:

  A data frame with one row per effect and the following columns:
  `effect` (the measure name), `contrast` after it when the results name
  contrasts, `estimate` (the pooled point estimate), `std.err` (the
  pooled standard error), `t` (the test statistic), `df` (the pooled
  degrees of freedom), `ci.lower` and `ci.upper`, `conf.level`, and
  `p.value`. The statistic and the p-value are referred to t on `df`
  rather than to the normal. When every pooled result carried the
  covariance of its effects, the pooled covariance is attached to this
  frame as the `ipw_vcov` attribute, with the effect labels as dimnames
  on both margins; when any of them carried none, no attribute is
  attached, since a matrix built from a subset of the imputations would
  sit beside estimates built from all of them.

- `pooling`:

  A data frame keyed by the same `effect` and `contrast` columns,
  holding `ubar` (the within-imputation variance), `b` (the
  between-imputation variance), `riv` (the relative increase in
  variance), `lambda` (the proportion of the total variance due to
  missingness), and `fmi` (the fraction of missing information).

- `se_method`:

  The standard error method every pooled result used.

- `effects`:

  The reading that was pooled, either `"marginal"` or `"conditional"`.

- `m`:

  The number of results pooled.

- `dfcom`:

  The complete-data degrees of freedom the adjustment used.

- `nobs`:

  The smallest number of observations any pooled result was estimated
  from.

- `outcome_link`:

  The link every pooled result's outcome model was fitted with, which is
  the scale the effects are reported on.

## Details

Each result is estimated on one completed dataset and carries a standard
error that accounts for the weights having been estimated, and pooling
those by Rubin's rules adds the uncertainty the imputation itself
contributed: the order matters, since a standard error that treated the
weights as fixed would be too small before the pooling ever saw it, and
pooling the imputed datasets rather than the estimates would understate
the uncertainty whatever the standard errors were.

For each effect the pooled estimate is the mean of the per-imputation
estimates, and the pooled variance is the mean of their squared standard
errors plus the between-imputation variance inflated by `1 + 1/m`
(Rubin, 1987). The degrees of freedom carry the Barnard-Rubin
small-sample adjustment (Barnard and Rubin, 1999), which is why the
complete-data count matters: with few imputations the pooled degrees of
freedom can be far below what a single analysis reports, and the
interval and the p-value are referred to t on them rather than to the
normal.

`mice` is not a dependency, and nothing here is imported from it. A
`mira` is recognized by its class and read through its `analyses`
element, which is all the object is, so a package that produces one by
another route is answered the same way.

## What the results have to agree on

The pooled estimate of an effect is an average of the per-imputation
ones, so the results have to be estimating the same things the same way.
A set that disagrees about its `estimand`, its `se_method`, the
presentation mode it records, the effects it reports, the confidence
level it stores, or the link its outcome model was fitted with is
refused with an error of class `causalgenerics_pool_mismatch`, and of a
second class naming which of those it was. The differing values travel
on the condition under `values`.

The effects have to agree as an ordered vector rather than as a set. The
labels are what say which row is which, so two results reporting the
same contrasts in different orders would otherwise have the `b vs a`
rows of one averaged with the `c vs a` rows of the other.

Two of those agreements are only required when the argument that would
settle the question is left at `NULL`. Naming `effects` says which
surface to pool and the stored mode is not read at all; naming
`conf_level` says what the bounds report and the stored level is not
read at all.

## The complete-data degrees of freedom

`dfcom` is the number of observations the analysis had minus the number
of parameters it fitted, before any data went missing. It is looked for
in three places, in order:

1.  The `dfcom` argument, when it is not `NULL`.

2.  [`df.residual()`](https://rdrr.io/r/stats/df.residual.html) on each
    result, which reports what the fitted variance object records, when
    at least one of them reports a number.

3.  [`df.residual()`](https://rdrr.io/r/stats/df.residual.html) on each
    result's outcome model, when at least one of those reports a number.

The smallest count is taken rather than the first, since the pooled
inference is no stronger than the weakest analysis supports, and results
that report nothing are passed over rather than making the minimum
missing. Wherever the count comes from it is at least 1, since the
adjustment has nothing to work with at zero.

When nothing reports a count, `dfcom` is `Inf` and a warning of class
`causalgenerics_pool_large_sample` says so. That is the honest answer
for a large sample and the widest degrees of freedom the adjustment can
give, which makes it the narrowest intervals, so it is said out loud
rather than assumed quietly.

## The scale the effects are pooled on

The estimates are pooled as they are stored. A result reports its ratio
effects on the log scale, under the labels `log(rr)` and `log(or)`,
which is the scale the inference is done on and the only one on which
averaging the per-imputation estimates means anything. Nothing here
exponentiates, and the pooled frame carries the same labels the results
did. Moving the ratios to their natural scale is presentation rather
than pooling, and belongs where the same choice is made for a single
result.

## References

Barnard, J. and Rubin, D. B. (1999). Small sample degrees of freedom
with multiple imputation. *Biometrika*, 86(4), 948-955.

Rubin, D. B. (1987). *Multiple Imputation for Nonresponse in Surveys*.
New York: John Wiley and Sons.

## See also

[`new_ipw()`](https://r-causal.github.io/causalgenerics/reference/new_ipw.md)
for the results this pools and the fields it reads.

## Examples

``` r
dat <- data.frame(
  x = rep(c(-1.5, -0.5, 0.5, 1.5), each = 5),
  z = rep(c(0, 1), 10),
  y = rep(c(0, 1, 1, 0, 1), 4)
)

# Three imputations of that dataset, differing in one filled-in cell each.
imputed <- lapply(1:3, function(i) {
  completed <- dat
  completed$y[i] <- 1 - completed$y[i]
  completed
})

# The estimates each analysis produced, written out literally in the shape
# the `ipw()` return contract documents. What is being shown is the pooling,
# so these stand in for what a method would compute from the models beside
# them.
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

# The pooled effects. The standard errors are wider than the average of the
# per-imputation ones, which is the uncertainty the imputation contributed.
pooled$estimates
#>    effect estimate   std.err        t        df   ci.lower  ci.upper conf.level
#> 1      rd      0.3 0.2449490 1.224745  9.648903 -0.2484843 0.8484843       0.95
#> 2 log(rr)      0.6 0.3082207 1.946657 15.563733 -0.0548904 1.2548904       0.95
#>      p.value
#> 1 0.24972994
#> 2 0.06985665

# How much of that uncertainty came from the imputation rather than the data.
pooled$pooling
#>    effect       ubar      b        riv     lambda       fmi
#> 1      rd 0.04666667 0.0100 0.28571429 0.22222222 0.3452017
#> 2 log(rr) 0.09166667 0.0025 0.03636364 0.03508772 0.1390444

# The complete-data count the adjustment used, read off the outcome models.
pooled$dfcom
#> [1] 18

# Naming one settles it directly, and the pooled degrees of freedom move.
pool_ipw(fits, dfcom = 500)$estimates$df
#> [1]  36.66623 370.84070

# The ratio effect stays on the log scale it was estimated on.
pooled$estimates$effect
#> [1] "rd"      "log(rr)"
```
