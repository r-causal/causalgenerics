# causalgenerics

causalgenerics is the shared generics package for the r-causal
ecosystem. It provides a near-zero-dependency home for the S3 generics
that packages such as propensity, halfmoon, positively, and balancing
register methods on. Owning the generic definitions in one place means
that attaching several r-causal packages at once produces no masking:
each package contributes methods to a common generic instead of
redefining the function.

The package is modeled on the [generics](https://generics.r-lib.org)
package, which provides commonly used S3 generics for the same reason:
so that packages can share a definition instead of each defining their
own.

## Installation

You can install the released version of causalgenerics from
[CRAN](https://cran.r-project.org/) with:

``` r

install.packages("causalgenerics")
```

You can install the development version from
[GitHub](https://github.com/r-causal/causalgenerics) with:

``` r

# install.packages("pak")
pak::pak("r-causal/causalgenerics")
```

Most users will get causalgenerics as a dependency of another r-causal
package rather than installing it directly.

## Generics

causalgenerics owns the following generics:

- [`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.md):
  bring-your-own-model inverse probability weighted estimation of causal
  effects from a weighting model and a weighted outcome model.
- [`ess()`](https://r-causal.github.io/causalgenerics/reference/ess.md):
  the effective sample size of a set of weights or a fitted model.
- [`is_causal_wt()`](https://r-causal.github.io/causalgenerics/reference/causal-weights.md),
  [`estimand()`](https://r-causal.github.io/causalgenerics/reference/causal-weights.md),
  and `estimand<-()`: accessors for the metadata carried by causal
  weight vectors.

The generics are intentionally minimal. Method-specific arguments are
passed through `...`, and the concrete weight classes live in the
packages that own them.

## How ecosystem packages depend on it

Packages in the r-causal ecosystem import causalgenerics and register
their methods against these generics. Because the generic is defined
once, a user can attach any combination of those packages without one
masking another.
