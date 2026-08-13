# Resolve the exposure type a function will work on

`match_exposure_type()` turns the `exposure_type` argument a function
exposes into the type it will actually work on. A caller either declares
one of the three types or asks for `"auto"`, and a function declares
through `valid_types` which of them it can answer.

## Usage

``` r
match_exposure_type(
  exposure_type = c("auto", "binary", "categorical", "continuous"),
  .exposure,
  valid_types = c("auto", "binary", "categorical", "continuous"),
  arg = ".exposure",
  announce = TRUE,
  call = rlang::caller_env()
)
```

## Arguments

- exposure_type:

  The type the caller declared, or `"auto"` to detect it from
  `.exposure`. A calling function passes its own argument through.

- .exposure:

  The exposure vector to classify.

- valid_types:

  The types the calling function can answer, including `"auto"` where it
  offers detection. The unsupported refusal names this set without
  `"auto"`, which is a request for detection rather than a type an
  exposure can have.

- arg:

  The name the announcement gives the exposure argument. A function
  whose exposure argument is named something else passes its own name,
  so the announcement never points at an argument the caller does not
  have.

- announce:

  Should the detected type be announced? The announcement reports which
  branch the analysis took, so a caller passes `FALSE` where that is not
  worth reporting, such as one detecting several exposures at once. It
  is suppressed for every call, whatever this argument says, when
  `options(causalgenerics.quiet = TRUE)`.

- call:

  The environment both refusals report as their calling context.

## Value

A single string: `"binary"`, `"categorical"`, or `"continuous"`.

## Details

Resolution has two tiers, and which one a request falls into is what the
refusal says. A type this package knows but the calling function does
not support is refused as unsupported, with an error of class
`causalgenerics_unsupported_exposure_type`, and the message names the
types the function does support. A string that is no exposure type at
all is refused as an unrecognized value instead, by
[`rlang::arg_match()`](https://rlang.r-lib.org/reference/arg_match.html),
which names the choices and suggests the nearest one. Both refusals
report the environment given in `call`, so a consumer package names its
own function rather than this one.

A detected type outside `valid_types` reaches the same unsupported
refusal as a declared one. Naming a type and letting it be detected
therefore give the same answer, and a caller cannot get further into a
function by leaving the type off.

A declared type that the function supports is returned as it stands. The
data is not consulted on that path and nothing is announced: the
declaration is how the caller models the exposure, and a heuristic over
unique values does not overrule it. A function that also wants a
declaration the data cannot carry refused calls
[`check_forced_type()`](https://r-causal.github.io/causalgenerics/reference/check_forced_type.md),
which is a separate decision.

A function that narrows `valid_types` narrows the default of its own
`exposure_type` argument to match. The two are matched against each
other, so a function offering all four choices while supporting three of
them reports the mismatch on the first call that leaves the argument at
its default.

## See also

[`detect_exposure_type()`](https://r-causal.github.io/causalgenerics/reference/detect_exposure_type.md),
which supplies the reading under `"auto"`, and
[`check_forced_type()`](https://r-causal.github.io/causalgenerics/reference/check_forced_type.md),
which refuses a declaration the data cannot carry.

## Examples

``` r
# A declared type is returned as it stands, without consulting the data.
match_exposure_type("continuous", c(0, 1, 2, 3))
#> [1] "continuous"

# `"auto"` reads the exposure and announces what it read.
match_exposure_type("auto", c(0, 1, 1, 0))
#> ℹ Treating `.exposure` as binary
#> [1] "binary"

# A function that works on two of the types says so.
try(match_exposure_type(
  "continuous",
  c(0, 1),
  valid_types = c("auto", "binary", "categorical")
))
#> Error in eval(expr, envir) : 
#>   Exposure type "continuous" is not supported.
#> ℹ Supported exposure types: "binary" and "categorical".
```
