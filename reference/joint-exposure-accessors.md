# Read the declaration a joint exposure carries

These three functions read back what
[`joint_exposure()`](https://r-causal.github.io/causalgenerics/reference/joint_exposure.md)
recorded, so that a package holding a column of exposures can act on the
crossing without rebuilding it.

`is_joint_exposure()` reports whether an object declares a crossing.
`joint_components()` reports which two treatments were crossed and the
levels each of them takes. `joint_reference()` reports the cell the
effects are reported against.

## Usage

``` r
is_joint_exposure(x)

joint_components(x)

joint_reference(x)
```

## Arguments

- x:

  A joint exposure, as built by
  [`joint_exposure()`](https://r-causal.github.io/causalgenerics/reference/joint_exposure.md).
  `is_joint_exposure()` accepts any object.

## Value

`is_joint_exposure()` returns a single `TRUE` or `FALSE`.
`joint_components()` returns a named list of two character vectors, one
per component, named for the treatments. `joint_reference()` returns the
reference cell's label as a single string.

## Details

They are plain functions rather than generics. The declaration lives in
the attributes of one class, and a package that meets a joint exposure
needs the same answer this package would give rather than a place to
register a different one.

`joint_components()` gives the component names in the order the crossing
was declared in, and each component's levels in the order the crossing
varies them, which is the order that settles the cell labels.
`joint_reference()` gives the first level of the vector, which is the
cell crossing the two components' reference levels.

## See also

[`joint_exposure()`](https://r-causal.github.io/causalgenerics/reference/joint_exposure.md),
which records what these read back.

## Examples

``` r
x <- joint_exposure(
  qsmk = c(0, 1, 0, 1),
  exercise = c("no", "no", "yes", "yes")
)

is_joint_exposure(x)
#> [1] TRUE
joint_components(x)
#> $qsmk
#> [1] "0" "1"
#> 
#> $exercise
#> [1] "no"  "yes"
#> 
joint_reference(x)
#> [1] "qsmk = 0, exercise = no"
```
