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
  effects = "marginal",
  readings = c("marginal", "conditional")
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
  exponentiate = FALSE,
  effects = NULL,
  call = sys.call()
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
  describes. It has to be a data frame, and it has to name each row's
  effect measure in an `effect` column of strings or of factor levels,
  with a measure named in every row, since every label the result
  reports is read off that column. Anything else is refused with an
  error of class `causalgenerics_invalid_argument_estimates`, and of the
  general class `causalgenerics_invalid_argument`.

- se_method:

  The standard error method that ran, such as `"mestimation"` or
  `"linearization"`.

- fit:

  The fitted variance object, or `NULL` when the method has none.

- effects:

  For `new_ipw()`, the presentation mode the result reports its effects
  in, either `"marginal"` or `"conditional"`; a method that names no
  mode reports marginal effects. For
  [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html), the
  reading to report: `NULL`, the default, reports the reading the result
  records, and any other value overrides it for the one call and leaves
  the result as it is. A reading the result does not support is refused
  with an error of class `causalgenerics_unsupported_reading`.

- readings:

  The readings the result supports, one or both of `"marginal"` and
  `"conditional"`. The default is both, which is what every result
  supported before the field existed. A method whose analysis has no
  meaning under one of them names the other, and the mode the result
  records has to be one of the readings named.

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

  If `TRUE`, move the estimates that are on a log scale to their natural
  scale. In the marginal reading those are the log risk ratio and the
  log odds ratio, and the two terms are relabeled `"rr"` and `"or"`; in
  the conditional reading the outcome model's link settles it for the
  whole table, as the section above describes. The confidence bounds
  move with the estimates. Standard errors, statistics, and p-values
  remain on the log scale, and the `ipw_vcov` attribute is dropped
  rather than carried, since it describes the estimates on the scale
  they were estimated on. Default is `FALSE`.

- call:

  The call to report a refusal against. A tidier that builds its table
  by calling this method passes the call a user wrote, so the refusal
  names the function they typed rather than the delegation behind it.
  The default is the method's own call, which is what a direct call
  reports.

## Value

`new_ipw()` returns an S3 object of class `ipw`: a list of the following
eight components, in this order.

- `estimand`:

  The causal estimand, such as `"ate"` or `"att"`.

- `wt_mod`:

  The weighting object: the fitted model that produced the weights.

- `outcome_mod`:

  The fitted outcome model.

- `estimates`:

  A data frame with one row per effect measure and the following
  columns: `effect` (the measure name, as a character or factor column
  naming a measure in every row, which the constructor requires),
  `estimate` (point estimate), `std.err` (standard error), `z`
  (z-statistic), `ci.lower` and `ci.upper` (confidence interval bounds),
  `conf.level`, and `p.value`. For a categorical exposure the data frame
  also has a `contrast` column, placed after `effect`, naming the
  non-reference level and reference level of each contrast. A frame
  stored against an earlier version of this contract names that column
  `comparison`. The older name is read as an alias for the canonical one
  wherever the column is read, so a result holding such a frame labels
  its rows and reports its table exactly as one holding a `contrast`
  column does. A method written now writes `contrast`. A result reported
  once per level of a grouping variable also has a `group` column,
  placed after the contrast column when one is present and after
  `effect` when the result names no contrasts, naming the subgroup each
  row was estimated in as a `"var = value"` string such as `"sex = 0"`.
  That column has one spelling and no alias. Both optional columns are
  absent rather than constant when the result reports one contrast or
  one group, since a column repeating a single value down the table
  would read as a contrast or a subgroup that was named. `group` must be
  character and must name a subgroup in every row; the constructor
  refuses anything else with an error of class
  `causalgenerics_invalid_argument_estimates`, and of the general class
  `causalgenerics_invalid_argument`, since a label pasted from such a
  column would relabel every row it keys without failing anywhere.

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
  conditional reading presents the outcome model's coefficient surface.
  See
  [`as_marginal()`](https://r-causal.github.io/causalgenerics/reference/ipw-modes.md)
  and
  [`as_conditional()`](https://r-causal.github.io/causalgenerics/reference/ipw-modes.md).

- `readings`:

  The readings the result supports, one or both of `"marginal"` and
  `"conditional"`, and always including the mode above. A result stored
  before the field existed carries fewer components and is read as
  supporting both readings.

[`print()`](https://rdrr.io/r/base/print.html) returns its input
invisibly.
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) returns a
plain data frame under the tidier column names described above, with the
confidence bounds appended when they are asked for: the effect estimates
in the marginal reading, and the outcome model's coefficients, one row
each, in the conditional one.

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
report, and `effects` and `readings` are present whether or not the
method that built the result named a mode or a set of readings.

[`print()`](https://rdrr.io/r/base/print.html) writes the estimand and
the call of each component model, then the table of the surface the
result's presentation mode names. The section below describes the two
modes and what each one tabulates.

[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) reports
the result as a tidier-shaped table rather than as a copy of the
`estimates` component. In the marginal reading its columns are `term`,
then `contrast` when the result names contrasts, then `group` when it
names subgroups, then `estimate`, `std.error`, `statistic`, and
`p.value`. Those are the names the tidier convention uses, so a fitting
package's `tidy()` method is this table read as a tibble and nothing
more. The `estimates` component itself is unchanged by any of what
follows.

The conditional reading reports the outcome model's coefficients under
those same headings, one row per coefficient, with `term` naming the
coefficient and `std.error` the standard error the corrected covariance
implies. It names no contrasts and no subgroups: those key the rows of
the effects table, and a coefficient is named by itself. It is the table
[`as.data.frame()`](https://r-causal.github.io/causalgenerics/reference/ipw-pooled-methods.md)
reports for a pooled result in the same reading, without the `df` column
a pooled result has and an unpooled one does not.

A result whose outcome model carries no corrected block has no such
standard error to report, and the reading is refused with an error of
class `causalgenerics_no_conditional_vcov` rather than reported from the
standard errors the model computed for itself. That is where the table
parts company with [`print()`](https://rdrr.io/r/base/print.html), which
writes those coefficients under a note saying that no covariance from
the joint estimation is recorded: a printed table is read by someone who
reads the note with it, and a data frame is read by code that would take
the column for the corrected one.

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
standard error. The conditional reading has no stored pair to prefer,
since the bounds the frame holds belong to the effects the other reading
reports, so its bounds are the normal approximation at every level.

With `exponentiate = TRUE` the `log(rr)` and `log(or)` rows of the
marginal table move to their natural scale, exponentiating the point
estimate and the confidence bounds and relabeling the two terms `"rr"`
and `"or"`. Standard errors, statistics, and p-values stay on the log
scale, where the inference is done, and the interval is settled before
the scale is: bounds recomputed at another level are built on the log
scale and exponentiated afterwards. The covariance described below
travels on the returned table under the same `ipw_vcov` attribute while
the rows are on the scale they were estimated on, and is dropped when
`exponentiate = TRUE`, since a matrix left attached there would describe
neither the table it sits on nor anything else.

A conditional table has no rows labeled as ratios to pick out, so the
link the outcome model was fitted with settles the question for the
whole table: a `logit` link puts every coefficient on the log odds scale
and a `log` link puts every coefficient on the log risk scale, and both
are scales an exponential undoes. Every estimate moves and no term is
relabeled, since a coefficient name names the term rather than the scale
its estimate is reported on. Every other link raises an error of class
`causalgenerics_exponentiate_link`, and of the classes
`causalgenerics_invalid_argument_exponentiate` and
`causalgenerics_invalid_argument`, rather than exponentiating
coefficients that describe nothing once exponentiated. That is the rule
[`as.data.frame()`](https://r-causal.github.io/causalgenerics/reference/ipw-pooled-methods.md)
keeps for a pooled result, in the same words.

## The effect labels

Every row of `estimates` has a label, and it is the label rather than
the position that [`print()`](https://rdrr.io/r/base/print.html) writes
down the side of its table and that
[`coef()`](https://r-causal.github.io/causalgenerics/reference/ipw-accessors.md),
[`vcov()`](https://r-causal.github.io/causalgenerics/reference/ipw-accessors.md),
and
[`confint()`](https://r-causal.github.io/causalgenerics/reference/ipw-accessors.md)
name their results with. The label is the columns that name a row,
pasted together in the order they appear in the frame's contract: the
`effect` column, then `contrast` when the result names contrasts, then
`group` when it names subgroups. A binary exposure estimated for the
whole sample is labeled `"rd"`, a categorical one `"rd b vs a"`, and one
reported by subgroup `"rd b vs a sex = 0"`. A categorical exposure
repeats each effect measure across its contrasts and a subgroup analysis
repeats each one across its groups, so `effect` alone would name several
rows the same thing.

## The presentation mode

A result reports its effects in one of two readings, recorded in the
`effects` field. The `"marginal"` reading shows the causal contrast
estimates the method targeted; the `"conditional"` reading presents the
outcome model's coefficient surface. The field says which one the result
presents, and
[`as_marginal()`](https://r-causal.github.io/causalgenerics/reference/ipw-modes.md)
and
[`as_conditional()`](https://r-causal.github.io/causalgenerics/reference/ipw-modes.md)
are how a caller moves a result between them.

Which readings a result can present at all is the separate fact the
`readings` field records. Both surfaces exist on most results and not on
all of them: an exposure entering the outcome model through several
columns has no single coefficient to read as the conditional effect, and
a package that builds such a result records the reading it can answer
for. A result records both unless the method that built it said
otherwise, and the mode it records has to be one of them, so a result
never presents a reading it does not support. Asking one for the reading
it does not support raises an error of class
`causalgenerics_unsupported_reading`, whether it is asked through the
mode generics, through
[`coef()`](https://r-causal.github.io/causalgenerics/reference/ipw-accessors.md),
[`vcov()`](https://r-causal.github.io/causalgenerics/reference/ipw-accessors.md),
or
[`confint()`](https://r-causal.github.io/causalgenerics/reference/ipw-accessors.md),
or through the `effects` argument of
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html). A result
stored before the field existed records no set and is read as supporting
both.

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

# A result whose analysis has no conditional reading records the one it
# supports, and the other is refused rather than reported.
marginal_only <- new_ipw(
  estimand = "ate",
  wt_mod = glm(z ~ x, family = binomial(), data = dat),
  outcome_mod = glm(y ~ z, family = quasibinomial(), data = dat),
  estimates = estimates,
  se_method = "linearization",
  fit = NULL,
  readings = "marginal"
)

try(as.data.frame(marginal_only, effects = "conditional"))
#> Error in as.data.frame.ipw(marginal_only, effects = "conditional") : 
#>   This result supports the marginal reading only, so there is no conditional reading of it to report; the package that produced it records the readings it supports when it builds the result.
```
