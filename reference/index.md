# Package index

## Generics

The S3 generics the r-causal ecosystem shares. Packages register methods
on these rather than defining their own, so attaching several at once
produces no masking.

- [`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.md)
  : Inverse probability weighted estimation
- [`as_marginal()`](https://r-causal.github.io/causalgenerics/reference/ipw-modes.md)
  [`as_conditional()`](https://r-causal.github.io/causalgenerics/reference/ipw-modes.md)
  : The presentation mode of an inverse probability weighted result
- [`ess()`](https://r-causal.github.io/causalgenerics/reference/ess.md)
  : Effective sample size
- [`is_causal_wt()`](https://r-causal.github.io/causalgenerics/reference/causal-weights.md)
  [`estimand()`](https://r-causal.github.io/causalgenerics/reference/causal-weights.md)
  [`` `estimand<-`() ``](https://r-causal.github.io/causalgenerics/reference/causal-weights.md)
  : Causal weight accessors

## Shared classes

Constructors and helpers for the classes those generics return. Method
authors build on these instead of writing a result or weight class of
their own.

- [`new_ipw()`](https://r-causal.github.io/causalgenerics/reference/new_ipw.md)
  [`print(`*`<ipw>`*`)`](https://r-causal.github.io/causalgenerics/reference/new_ipw.md)
  [`as.data.frame(`*`<ipw>`*`)`](https://r-causal.github.io/causalgenerics/reference/new_ipw.md)
  : Construct an inverse probability weighted result
- [`new_ipw_model()`](https://r-causal.github.io/causalgenerics/reference/new_ipw_model.md)
  [`vcov(`*`<ipw_model>`*`)`](https://r-causal.github.io/causalgenerics/reference/new_ipw_model.md)
  : Attach a corrected covariance to a component model of an IPW fit
- [`coef(`*`<ipw>`*`)`](https://r-causal.github.io/causalgenerics/reference/ipw-accessors.md)
  [`vcov(`*`<ipw>`*`)`](https://r-causal.github.io/causalgenerics/reference/ipw-accessors.md)
  [`confint(`*`<ipw>`*`)`](https://r-causal.github.io/causalgenerics/reference/ipw-accessors.md)
  [`nobs(`*`<ipw>`*`)`](https://r-causal.github.io/causalgenerics/reference/ipw-accessors.md)
  [`df.residual(`*`<ipw>`*`)`](https://r-causal.github.io/causalgenerics/reference/ipw-accessors.md)
  [`weights(`*`<ipw>`*`)`](https://r-causal.github.io/causalgenerics/reference/ipw-accessors.md)
  : Model accessors for an inverse probability weighted result
- [`new_causal_wts()`](https://r-causal.github.io/causalgenerics/reference/new_causal_wts.md)
  : Construct a causal weight vector
- [`causal_wts_ptype2()`](https://r-causal.github.io/causalgenerics/reference/causal_wts_ptype2.md)
  : Common type of two causal weight vectors of the same class
