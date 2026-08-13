# Read the values an exposure takes

The three readings
[`detect_exposure_type()`](https://r-causal.github.io/causalgenerics/reference/detect_exposure_type.md)
classifies from, exported so that a package taking its own decision from
the same facts reads them the same way.

`observed_values()` returns the distinct values a vector takes.
`has_two_levels()` reports whether there are exactly two of them.
`is_categorical()` applies the unique-value heuristic that separates a
categorical numeric exposure from a continuous one.

## Usage

``` r
observed_values(.x)

has_two_levels(.x)

is_categorical(.exposure)
```

## Arguments

- .x:

  A vector.

- .exposure:

  An exposure vector.

## Value

`observed_values()` returns the distinct observed values, as a vector of
the same type as its input. `has_two_levels()` and `is_categorical()`
each return a single logical value.

## Details

A missing value is not a level. `observed_values()` drops the missing
values before taking the distinct ones, and that is what the other two
count, so an exposure recorded with missing values is read as the
exposure it is. Counting the missing values as a value of their own
would give a binary exposure three levels, which takes it off the binary
path entirely, and would make a categorical one look as though it had
one category more than it has.

`is_categorical()` treats a vector as categorical when the number of
values it takes is less than 20 percent of the number of observations
that are not missing. The threshold is strict, so a ratio of exactly 0.2
is not categorical. An exposure with nothing observed has no ratio to
take and is not categorical, which leaves it on the continuous branch of
detection rather than returning a missing value into the classification.

The distinct values are returned in the order they first appear rather
than sorted, since they are only ever counted and sorting them would
impose an ordering on a categorical exposure that nothing else here
honors. The values keep the type they came in as, so the observed values
of a factor are a factor over the same declared levels.

## See also

[`detect_exposure_type()`](https://r-causal.github.io/causalgenerics/reference/detect_exposure_type.md),
which combines these into a type.

## Examples

``` r
observed_values(c(0, 1, NA, 1))
#> [1] 0 1
has_two_levels(c(0, 1, NA, 1))
#> [1] TRUE

# 15 values over 300 observations is categorical; 300 over 300 is not.
is_categorical(rep_len(1:15, 300))
#> [1] TRUE
is_categorical(seq_len(300) / 300)
#> [1] FALSE
```
