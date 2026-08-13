# Refuse an exposure type the data cannot carry

`check_forced_type()` is the structural check on a declared exposure
type. A declaration wins over the detection heuristics wherever the data
can carry it, so this refuses only the declarations that no data of this
shape could support.

## Usage

``` r
check_forced_type(
  forced,
  .exposure,
  arg = ".exposure",
  call = rlang::caller_env()
)
```

## Arguments

- forced:

  The type the caller declared: `"binary"`, `"categorical"`, or
  `"continuous"`. `"auto"` is a request for detection rather than a
  declaration, and is refused here along with anything else that is not
  one of the three.

- .exposure:

  The exposure vector to classify.

- arg:

  The name the refusal gives the exposure argument.

- call:

  The environment the refusal reports as its calling context.

## Value

`NULL`, invisibly. The function is called for the refusal.

## Details

Two declarations can be contradicted. A binary exposure takes exactly
two observed values, and a continuous one is numeric. Every vector can
be read as the set of levels it takes, so a categorical declaration is
always possible and always passes.

Everything else stands. A numeric dose taking fifteen values declared
continuous is fitted as a dose even though detection would read it as
categorical, and a numeric exposure taking two values declared
continuous is the caller's decision. A check that refused those would
put the unique-value heuristic in charge of a reading the caller had
already made.

The refusal quotes what detection would have read, since dropping the
declaration is the route out of it. That reading is computed and not
announced: the caller declared a type, so an alert naming a different
one would read as the answer rather than as the alternative.

The refusal names `exposure_type` as the argument the declaration came
from, which is what the ecosystem calls it.

## See also

[`match_exposure_type()`](https://r-causal.github.io/causalgenerics/reference/match_exposure_type.md),
which resolves the declaration this checks, and
[`detect_exposure_type()`](https://r-causal.github.io/causalgenerics/reference/detect_exposure_type.md),
which supplies the reading it quotes.

## Examples

``` r
# A declaration the data can carry passes silently.
check_forced_type("binary", c(0, 1, NA, 1))

# So does one detection would have read the other way.
check_forced_type("continuous", rep_len(1:15, 300))

# A binary exposure takes two values, and this one takes three.
try(check_forced_type("binary", factor(c("low", "medium", "high"))))
#> Error in eval(expr, envir) : 
#>   `exposure_type` was set to "binary", but `.exposure` cannot be treated
#> that way.
#> ✖ A "binary" exposure takes exactly two observed values, and `.exposure` takes
#>   3.
#> ℹ Drop `exposure_type` to detect the type from the data, which reads
#>   `.exposure` as "categorical".

# A continuous exposure is numeric.
try(check_forced_type("continuous", c("a", "b", "c")))
#> Error in eval(expr, envir) : 
#>   `exposure_type` was set to "continuous", but `.exposure` cannot be
#> treated that way.
#> ✖ A "continuous" exposure is numeric, and `.exposure` is a character vector.
#> ℹ Drop `exposure_type` to detect the type from the data, which reads
#>   `.exposure` as "categorical".
```
