# Package index

## Generics

The S3 generics the r-causal ecosystem shares. Packages register methods
on these rather than defining their own, so attaching several at once
produces no masking.

- [`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.md)
  : Inverse probability weighted estimation
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
- [`new_causal_wts()`](https://r-causal.github.io/causalgenerics/reference/new_causal_wts.md)
  : Construct a causal weight vector
- [`causal_wts_ptype2()`](https://r-causal.github.io/causalgenerics/reference/causal_wts_ptype2.md)
  : Common type of two causal weight vectors of the same class
