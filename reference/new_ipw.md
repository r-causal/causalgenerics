# Construct an inverse probability weighted result

`new_ipw()` is the low-level constructor for the object that every
[`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.md)
method returns. It is intended for package developers writing an
[`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.md)
method, not for end users, and it assumes its arguments are already
validated.

## Usage

``` r
new_ipw(
  estimand,
  wt_mod,
  outcome_mod,
  estimates,
  se_method,
  fit,
  effects = "marginal"
)

# S3 method for class 'ipw'
print(x, ...)

# S3 method for class 'ipw'
as.data.frame(
  x,
  row.names = NULL,
  optional = FALSE,
  ...,
  conf.int = FALSE,
  conf.level = 0.95,
  exponentiate = FALSE
)
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

- effects:

  The presentation mode the result reports its effects in, either
  `"marginal"` or `"conditional"`. A method that names no mode reports
  marginal effects.

- x:

  An `ipw` object.

- ...:

  Further arguments. Neither method reads them.

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

  The confidence level the bounds report. At the level every row of
  `estimates` records, the stored bounds are returned; at any other
  level they are the normal approximation built from the estimate and
  its standard error. Default is `0.95`.

- exponentiate:

  If `TRUE`, exponentiate the log risk ratio and log odds ratio to
  produce risk ratios and odds ratios on their natural scale,
  relabelling the two terms `"rr"` and `"or"`. The confidence bounds
  move with them. Standard errors, statistics, and p-values remain on
  the log scale, and the `ipw_vcov` attribute is dropped rather than
  carried, since it describes the estimates on the scale they were
  estimated on. Default is `FALSE`.

## Value

`new_ipw()` returns an S3 object of class `ipw`: a list of the following
seven components, in this order.

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

- `effects`:

  The presentation mode, either `"marginal"` or `"conditional"`. The
  marginal reading shows the causal contrast estimates and the
  conditional reading presents the outcome model's coefficient surface;
  both surfaces exist on every result. See
  [`as_marginal()`](https://r-causal.github.io/causalgenerics/reference/ipw-modes.md)
  and
  [`as_conditional()`](https://r-causal.github.io/causalgenerics/reference/ipw-modes.md).

[`print()`](https://rdrr.io/r/base/print.html) returns its input
invisibly.
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) returns a
plain data frame of the effect estimates under the tidier column names
described above, with the confidence bounds appended when they are asked
for.

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
report, and `effects` is present whether or not the method that built
the result named a mode.

[`print()`](https://rdrr.io/r/base/print.html) writes the estimand and
the call of each component model, then the table of the surface the
result's presentation mode names. The section below describes the two
modes and what each one tabulates.

[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) reports
the effect estimates as a tidier-shaped table rather than as a copy of
the `estimates` component. Its columns are `term`, then `comparison`
when the result names contrasts, then `estimate`, `std.error`,
`statistic`, and `p.value`. Those are the names the tidier convention
uses, so a fitting package's `tidy()` method is this table read as a
tibble and nothing more. The `estimates` component itself is unchanged
by any of what follows.

`conf.int = TRUE` appends `conf.low` and `conf.high` after the other
columns, and `conf.level` names the level they report. The level is an
argument rather than a column, since a column would repeat one number
down every row and be read as part of the table rather than as the level
the two bounds beside it were built at. The bounds `estimates` stores
are returned when every row of the frame records the level asked for.
They need not be the normal pair: a bootstrap or profile interval is
asymmetric about the estimate, and even a normal one rounded on its way
into the frame is not the number recomputing gives. At any other level,
and for a frame whose rows disagree about the level or record none, the
bounds are the normal approximation built from the estimate and its
standard error.

With `exponentiate = TRUE` the `log(rr)` and `log(or)` rows move to
their natural scale, exponentiating the point estimate and the
confidence bounds and relabelling the two terms `"rr"` and `"or"`.
Standard errors, statistics, and p-values stay on the log scale, where
the inference is done, and the interval is settled before the scale is:
bounds recomputed at another level are built on the log scale and
exponentiated afterwards. The covariance described below travels on the
returned table under the same `ipw_vcov` attribute while the rows are on
the scale they were estimated on, and is dropped when
`exponentiate = TRUE`, since a matrix left attached there would describe
neither the table it sits on nor anything else.

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

## The presentation mode

A result reports its effects in one of two readings, recorded in the
`effects` field. The `"marginal"` reading shows the causal contrast
estimates the method targeted; the `"conditional"` reading presents the
outcome model's coefficient surface. Both surfaces always exist on the
object, so the field says which one the result presents rather than
which one it holds.
[`as_marginal()`](https://r-causal.github.io/causalgenerics/reference/ipw-modes.md)
and
[`as_conditional()`](https://r-causal.github.io/causalgenerics/reference/ipw-modes.md)
are how a caller moves a result between them.

A printed result names its mode twice, since the two readings are
different tables of different numbers: once on an `Effects:` line beside
the estimand, and once in the heading of the table itself. The marginal
reading tabulates the effect estimates the result stores, under
`Marginal estimates:`. The conditional reading tabulates the outcome
model's coefficients, under `Conditional estimates (outcome model):`,
with the standard errors implied by the corrected covariance a fitting
package attaches through
[`new_ipw_model()`](https://r-causal.github.io/causalgenerics/reference/new_ipw_model.md).
Which coefficient each entry of that block belongs to is what its labels
say rather than its row order, so a block attached in another order
still prints each standard error beside the coefficient it belongs to.

An outcome model a fitting package never wrapped is still printed: the
coefficients are written on their own, followed by a note saying that no
covariance from the joint estimation is recorded, rather than beside the
standard errors the model computed for itself. An outcome model that
carries the
[`new_ipw_model()`](https://r-causal.github.io/causalgenerics/reference/new_ipw_model.md)
class with no covariance behind it is refused instead, with an error of
class `causalgenerics_no_vcov_ipw_model`, and one whose block cannot be
paired with its coefficients with an error of class
`causalgenerics_conditional_vcov_mismatch`. The note answers a package
that has not adopted the contract by telling the reader to wrap the
model, and that advice has already been taken in both of those cases:
the object is what is wrong, not the package that produced it.

An outcome model that reports no coefficients has no rows to tabulate
under that heading, and the printed form says so in place of the table.

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
the generic these results come from, and
[`as_marginal()`](https://r-causal.github.io/causalgenerics/reference/ipw-modes.md)
and
[`as_conditional()`](https://r-causal.github.io/causalgenerics/reference/ipw-modes.md)
for the presentation mode.

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
#> Effects: marginal (population-averaged) 
#> 
#> Weight Estimator:
#>   Call: glm(formula = z ~ x, family = binomial(), data = dat) 
#> 
#> Outcome Model:
#>   Call: glm(formula = y ~ z, family = quasibinomial(), data = dat) 
#> 
#> Marginal estimates:
#>         estimate  std.err      z ci.lower ci.upper conf.level p.value  
#> rd      0.199882 0.092425 2.1626 0.018732  0.38103       0.95 0.03057 *
#> log(rr) 0.560414 0.273519 2.0489 0.024326  1.09650       0.95 0.04047 *
#> log(or) 0.878313 0.418661 2.0979 0.057753  1.69887       0.95 0.03591 *
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# The tidier-shaped table the result reports.
as.data.frame(res)
#>      term estimate std.error statistic p.value
#> 1      rd 0.199882  0.092425    2.1626 0.03057
#> 2 log(rr) 0.560414  0.273519    2.0489 0.04047
#> 3 log(or) 0.878313  0.418661    2.0979 0.03591

# With an interval, and the ratios on their natural scale.
as.data.frame(res, conf.int = TRUE, exponentiate = TRUE)
#>   term estimate std.error statistic p.value conf.low conf.high
#> 1   rd 0.199882  0.092425    2.1626 0.03057 0.018732  0.381032
#> 2   rr 1.751397  0.273519    2.0489 0.04047 1.024624  2.993676
#> 3   or 2.406836  0.418661    2.0979 0.03591 1.059453  5.467782
```
