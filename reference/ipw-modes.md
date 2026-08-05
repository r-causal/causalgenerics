# The presentation mode of an inverse probability weighted result

The two readings of a result, defined for the class
[`new_ipw()`](https://r-causal.github.io/causalgenerics/reference/new_ipw.md)
constructs.

- `as_marginal()` returns the result reporting the causal contrast
  estimates.

- `as_conditional()` returns the result presenting the outcome model's
  coefficient surface.

## Usage

``` r
as_marginal(x, ...)

as_conditional(x, ...)
```

## Arguments

- x:

  An `ipw` object. These generics dispatch on this argument.

- ...:

  Arguments passed to methods.

## Value

`x` with its presentation mode set to the one asked for. The methods on
`ipw` change the `effects` field and nothing else, so every other field
comes back as it went in, the covariance attached to `estimates`
included.

## Details

Both surfaces exist on every result, so these generics record which one
the result presents rather than computing anything. They set the
`effects` field of the
[`new_ipw()`](https://r-causal.github.io/causalgenerics/reference/new_ipw.md)
contract and read nothing else, which makes them the supported way to
move a result between the two readings: a caller writes
`as_conditional(res)` rather than assigning to the field.

The methods on `ipw` are total. Every result has one of the two modes,
so asking for either is always answerable: they never error, asking
twice says what asking once said, and a result that goes out to the
other reading and back is the result that went in. A result built before
the field existed carries six fields rather than seven and reads as
marginal, which is the mode every method produced then.

The generics live here for the reason
[`print()`](https://rdrr.io/r/base/print.html) does. Two packages each
registering `as_conditional.ipw()` would collide in the shared S3 method
table, and a caller writing against a result would then get whichever
package was installed last rather than the contract. There is no
marginal or conditional reading of an object that is not an IPW result,
so the default method signals an error rather than inventing one.

## See also

[`new_ipw()`](https://r-causal.github.io/causalgenerics/reference/new_ipw.md)
for the result class and the field these generics set.

## Examples

``` r
dat <- data.frame(
  x = rep(c(-1.5, -0.5, 0.5, 1.5), each = 5),
  z = rep(c(0, 1), 10),
  y = rep(c(0, 1, 1, 0, 1), 4)
)

# Written out literally, in the shape the `ipw()` return contract documents.
estimates <- data.frame(
  effect = "rd",
  estimate = 0.199882,
  std.err = 0.092425,
  z = 2.1626,
  ci.lower = 0.018732,
  ci.upper = 0.381032,
  conf.level = 0.95,
  p.value = 0.030570
)

res <- new_ipw(
  estimand = "ate",
  wt_mod = glm(z ~ x, family = binomial(), data = dat),
  outcome_mod = glm(y ~ z, family = quasibinomial(), data = dat),
  estimates = estimates,
  se_method = "linearization",
  fit = NULL
)

# A method that names no mode reports marginal effects.
res$effects
#> [1] "marginal"

as_conditional(res)$effects
#> [1] "conditional"

# Asking for the reading a result already reports gives the result back, and
# the two readings round-trip.
identical(as_marginal(res), res)
#> [1] TRUE
identical(as_marginal(as_conditional(res)), res)
#> [1] TRUE

# Neither reading exists for an object that is not an IPW result.
try(as_marginal(1:3))
#> Error in as_marginal.default(1:3) : 
#>   No `as_marginal()` method for an object of class <integer>.
```
