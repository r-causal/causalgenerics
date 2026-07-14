# causalgenerics

<!-- badges: start -->
<!-- badges: end -->

causalgenerics is the shared generics package for the r-causal ecosystem. It
provides a near-zero-dependency home for the S3 generics that packages such as
propensity, halfmoon, positively, and balancing register methods on. Owning the
generic definitions in one place means that attaching several r-causal packages
at once produces no masking: each package contributes methods to a common
generic instead of redefining the function.

The package is modeled on the tidymodels generics package, which plays the same
role for the modeling ecosystem.

## Generics

causalgenerics owns the following generics:

- `ipw()`: bring-your-own-model inverse probability weighted estimation of
  causal effects from a weighting model and a weighted outcome model.
- `ess()`: the effective sample size of a set of weights or a fitted model.
- `is_causal_wt()`, `estimand()`, and `estimand<-()`: accessors for the
  metadata carried by causal weight vectors.

The generics are intentionally minimal. Method-specific arguments are passed
through `...`, and the classes themselves live in the packages that own them.

## How ecosystem packages depend on it

Packages in the r-causal ecosystem import causalgenerics and register their
methods against these generics. Because the generic is defined once, a user can
attach any combination of those packages without one masking another.

## Status

causalgenerics is under active development.
