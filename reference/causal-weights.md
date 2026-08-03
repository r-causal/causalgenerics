# Causal weight accessors

Generics for inspecting and setting the metadata carried by causal
weight vectors.

- `is_causal_wt()` tests whether an object is a causal weight vector.

- `estimand()` returns the causal estimand the weights target, such as
  `"ate"` or `"att"`.

- `estimand<-()` sets the causal estimand.

## Usage

``` r
is_causal_wt(x, ...)

estimand(x, ...)

estimand(x, ...) <- value
```

## Arguments

- x:

  An object to inspect or modify. These generics dispatch on this
  argument.

- ...:

  Arguments passed to methods.

- value:

  The causal estimand to assign.

## Value

`is_causal_wt()` returns a single logical value. `estimand()` returns
the estimand, typically a character string, or `NULL` when none is
recorded. `estimand<-()` returns `x` with the estimand updated.

## Details

This package defines the generics and the abstract weight class they are
written against. It supplies a method for each generic on `causal_wts`,
so a concrete class built with
[`new_causal_wts()`](https://r-causal.github.io/causalgenerics/reference/new_causal_wts.md)
inherits all three. The concrete classes themselves live in the packages
that own them, such as `propensity`, as do methods for any object that
is not a `causal_wts` vector. `is_causal_wt()` defaults to `FALSE` so
that any object that is not a recognized causal weight vector reports as
much; `estimand()` and `estimand<-()` have no meaningful default and
signal an error when no method is registered for the object.

## See also

[`new_causal_wts()`](https://r-causal.github.io/causalgenerics/reference/new_causal_wts.md)
for the class these generics have methods for, and the `propensity`
package for concrete weight classes.

## Examples

``` r
wts <- new_causal_wts(
  c(1.2, 0.8, 1.5),
  subclass = "cg_toy_wts",
  estimand = "ate"
)

is_causal_wt(wts)
#> [1] TRUE
estimand(wts)
#> [1] "ate"

# The estimand is metadata on the vector, so it can be replaced in place.
estimand(wts) <- "att"
estimand(wts)
#> [1] "att"

# An object that is not a causal weight vector reports `FALSE` rather than
# erroring.
is_causal_wt(1:3)
#> [1] FALSE

# `estimand()` has no meaningful default, so it errors instead.
try(estimand(1:3))
#> Error in estimand.default(1:3) : 
#>   No `estimand()` method for an object of class <integer>.
```
