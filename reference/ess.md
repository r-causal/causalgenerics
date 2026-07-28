# Effective sample size

`ess()` is the generic for the effective sample size, a summary of how
much information a set of weights retains relative to an unweighted
sample. When weights vary substantially, the effective sample size can
be much smaller than the number of observations.

## Usage

``` r
ess(x, ...)

# Default S3 method
ess(x, na.rm = FALSE, ...)
```

## Arguments

- x:

  An object to compute the effective sample size for, such as a weight
  vector or a fitted model. `ess()` dispatches on this argument.

- ...:

  Arguments passed to methods.

- na.rm:

  Should missing weights be dropped before computing? Missing weights
  otherwise make the result `NA`.

## Value

The default method returns a single number. Other methods return
whatever suits the object they are written for; a fitted-model method
may return a tibble.

## Details

The default method computes the Kish effective sample size, \\(\sum w)^2
/ \sum w^2\\, for any numeric vector. Every weight class in this
ecosystem is a double underneath, and
[`is.numeric()`](https://rdrr.io/r/base/numeric.html) is `TRUE` for one,
so the default is the whole implementation for weights and a concrete
weight class needs no method of its own. A method would be worth
registering only for an object that needs something other than the Kish
formula, such as a fitted model, where it might summarize by group and
return a tibble with one row per group.

The default errors when `x` is not numeric, which includes `NULL`,
rather than returning a value computed from nothing. Numeric input whose
quotient is `0 / 0` still returns `NaN`, which covers a zero-length
vector and weights that are all zero: both make the sum and the sum of
squares zero, and neither has an effective sample size to report.

The generic stays minimal, taking the object and `...`, so that a method
declares whatever further arguments it needs. A fitted-model method, for
example, takes the variable to group by that way.

## See also

`halfmoon::check_ess()`, which reports the effective sample size of one
or more weighting columns in a data frame, optionally by exposure group.
It is an ordinary function rather than an `ess()` method, so it does not
go through this generic.

## Examples

``` r
# Equal weights carry as much information as an unweighted sample.
ess(rep(0.5, 10))
#> [1] 10

# Weights that vary carry much less.
ess(c(rep(0.5, 9), 10))
#> [1] 2.056235
```
