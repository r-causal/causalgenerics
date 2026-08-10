# Declare an exposure that crosses two treatments

`joint_exposure()` crosses two discrete treatments into one categorical
exposure and records the declaration on the result: which two treatments
were crossed, in what order, which levels each of them takes, and which
cell the effects are reported against. The package that builds the
weights, the package that checks balance over the cells, and the package
that reports the effects all read that one declaration, so they agree on
which cells exist without any of them owning it.

## Usage

``` r
joint_exposure(...)
```

## Arguments

- ...:

  Exactly two named vectors, one per treatment, recycled to a common
  length. Each may be a factor, a character vector, or a numeric or
  logical vector. The names are the treatment names the cell labels are
  written from, and the two must be distinct.

## Value

A vector of class
`c("joint_exposure", "factor", "vctrs_vctr", "integer")` with one level
per cell of the crossing, carrying the component levels and the cell
labels the accessors read back.

## Details

Conditioning on a variable and jointly intervening on two treatments are
different causal questions, and this class serves the second. Effect
modification asks how the effect of one treatment differs across the
levels of something else, which need not be a treatment and need not be
anything an analysis could set; what comes back is one effect of that
one treatment per subgroup, and a result reporting those carries the
`group` column
[`new_ipw()`](https://r-causal.github.io/causalgenerics/reference/new_ipw.md)
documents. Interaction asks what would happen under an intervention that
sets both treatments at once; what comes back is one exposure with a
cell per combination, reported against the reference cell. The packages
that implement
[`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.md)
methods develop the distinction, in the weights they fit and in the
effects they report. What is settled here is the declaration all of them
read.

The result is a factor. Its class vector is
`c("joint_exposure", "factor", "vctrs_vctr", "integer")`, with
`"factor"` ahead of `"vctrs_vctr"` so that the formula machinery treats
a joint exposure natively:
[`model.matrix()`](https://rdrr.io/r/stats/model.matrix.html) gives
exactly what a plain factor over the same cells gives, `assign` and
`contrasts` attributes included.

The cells are the crossing of the two components' levels, labeled
`"name1 = level1, name2 = level2"` and ordered with the first component
varying fastest, which is the order
[`interaction()`](https://rdrr.io/r/base/interaction.html) produces. The
reference cell crosses the two components' own reference levels, and
that ordering is what puts it first, where the contrasts a model builds
are read against it. A factor component's reference level is the first
of its [`levels()`](https://rdrr.io/r/base/levels.html); any other
component's is the smallest of its values, which makes `0` the reference
of a 0/1 treatment and `FALSE` the reference of a logical one.

Every cell of the crossing must be populated. A cell nobody falls in is
a positivity violation rather than a small sample: the joint effect of
the two treatments is not identified in those data, and no weight,
contrast, or balance check over the crossing can be computed.
Construction refuses it, as it refuses a component that never varies and
a component with a missing value.

The two components must also be named distinctly. The names are what
tell the treatments apart in the cells, so one name used twice would
label every cell with the same variable on either side of it and leave
nothing downstream able to tell which component a cell varies.
Construction refuses that as well.

The declaration survives every operation that keeps the full set of
cells, including `[`,
[`vctrs::vec_slice()`](https://vctrs.r-lib.org/reference/vec_slice.html),
[`sort()`](https://rdrr.io/r/base/sort.html), and a round trip through a
data frame, and it survives combining two joint exposures that declare
the same crossing. Operations that narrow or rewrite the level set give
it up and say so:
[`droplevels()`](https://rdrr.io/r/base/droplevels.html),
`x[i, drop = TRUE]`, and `levels<-` each warn and return a plain factor.
So does combining a joint exposure with a factor, with a character
vector, or with a joint exposure that declares a different crossing.
Casting the other way is refused outright, because a crossing cannot be
recovered from labels that merely look like one.

Two ways of combining escape all of that and give a bare factor without
a word, because neither reaches a method this package can register.
[`unlist()`](https://rdrr.io/r/base/unlist.html) combines the underlying
codes in base C code without dispatching at all, and
[`c()`](https://rdrr.io/r/base/c.html) dispatches on its first argument,
so a combine that begins with a plain factor never reaches this class.
Reach for
[`vctrs::vec_c()`](https://vctrs.r-lib.org/reference/vec_c.html) where
either could apply.

## See also

[`is_joint_exposure()`](https://r-causal.github.io/causalgenerics/reference/joint-exposure-accessors.md),
[`joint_components()`](https://r-causal.github.io/causalgenerics/reference/joint-exposure-accessors.md),
and
[`joint_reference()`](https://r-causal.github.io/causalgenerics/reference/joint-exposure-accessors.md),
which read the declaration back off the vector.

## Examples

``` r
x <- joint_exposure(
  qsmk = c(0, 1, 0, 1),
  exercise = c("no", "no", "yes", "yes")
)
x
#> <joint_exposure[4]: qsmk, exercise>
#> [1] qsmk = 0, exercise = no  qsmk = 1, exercise = no  qsmk = 0, exercise = yes
#> [4] qsmk = 1, exercise = yes
#> Reference: qsmk = 0, exercise = no

levels(x)
#> [1] "qsmk = 0, exercise = no"  "qsmk = 1, exercise = no" 
#> [3] "qsmk = 0, exercise = yes" "qsmk = 1, exercise = yes"
joint_components(x)
#> $qsmk
#> [1] "0" "1"
#> 
#> $exercise
#> [1] "no"  "yes"
#> 
joint_reference(x)
#> [1] "qsmk = 0, exercise = no"

# The formula machinery sees the factor it is, and reports every cell
# against the reference one.
colnames(model.matrix(~x, data.frame(x = x)))
#> [1] "(Intercept)"               "xqsmk = 1, exercise = no" 
#> [3] "xqsmk = 0, exercise = yes" "xqsmk = 1, exercise = yes"
```
