# Attach a corrected covariance to a component model of an IPW fit

`new_ipw_model()` wraps one of the models an
[`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.md)
method was given so that the model reports the covariance the two-step
estimation implies rather than the one it computed for itself. It is
intended for package developers writing an
[`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.md)
method, not for end users, and beyond the shape of `vcov` it assumes its
arguments are already validated.

## Usage

``` r
new_ipw_model(model, vcov)

# S3 method for class 'ipw_model'
vcov(object, ...)
```

## Arguments

- model:

  The fitted component model to wrap, such as the outcome model or the
  weighting model of an
  [`ipw()`](https://r-causal.github.io/causalgenerics/reference/ipw.md)
  result.

- vcov:

  The corrected covariance of that model's parameters: a square numeric
  matrix with dimnames on both margins.

- object:

  An `ipw_model` object.

- ...:

  Further arguments. [`vcov()`](https://rdrr.io/r/stats/vcov.html)
  ignores them.

## Value

`new_ipw_model()` returns `model` with `"ipw_model"` prepended to its
class vector and `vcov` attached as the `ipw_vcov` attribute. Everything
else about the model is untouched.
[`vcov()`](https://rdrr.io/r/stats/vcov.html) returns that matrix.

## Details

A weighted outcome model fitted as though its weights were fixed
understates its own uncertainty, because the weights were estimated from
the same data. A method that solves the estimating equations of both
steps as one stacked system has the corrected block on hand, and
wrapping the model is how that block reaches a caller who writes
`vcov(result$outcome_mod)`.

The wrapper is two changes and no more: `"ipw_model"` is prepended to
the model's class vector, and the matrix is attached as the `ipw_vcov`
attribute. Nothing is registered for the class beyond
[`vcov()`](https://rdrr.io/r/stats/vcov.html), so
[`predict()`](https://rdrr.io/r/stats/predict.html),
[`coef()`](https://rdrr.io/r/stats/coef.html),
[`model.frame()`](https://rdrr.io/r/stats/model.frame.html), and every
other generic walk past it to the model's own methods. That is what lets
the corrected covariance travel with the model rather than beside it.

The dimnames of `vcov` are the fitting package's to choose. They are
ordinarily `names(coef(model))`, but a method whose stacked system is
parameterized some other way labels the block with its own parameter
names. The constructor checks that both margins are labelled rather than
what the labels say, since only the fitting package knows which block of
the sandwich belongs to which model.

## See also

[`new_ipw()`](https://r-causal.github.io/causalgenerics/reference/new_ipw.md)
for the result these models are components of, and
[ipw-accessors](https://r-causal.github.io/causalgenerics/reference/ipw-accessors.md)
for the accessors on the result itself.

## Examples

``` r
dat <- data.frame(
  x = rep(c(-1.5, -0.5, 0.5, 1.5), each = 5),
  z = rep(c(0, 1), 10),
  y = rep(c(0, 1, 1, 0, 1), 4)
)

ps <- fitted(glm(z ~ x, family = binomial(), data = dat))
wts <- ifelse(dat$z == 1, 1 / ps, 1 / (1 - ps))
outcome_mod <- glm(y ~ z, family = quasibinomial(), data = dat, weights = wts)

# Standing in for the outcome block of a stacked sandwich. A real method
# computes this from the estimating equations of both steps; treating the
# weights as estimated rather than fixed ordinarily widens it.
corrected <- vcov(outcome_mod) * 1.4

wrapped <- new_ipw_model(outcome_mod, corrected)

class(wrapped)
#> [1] "ipw_model" "glm"       "lm"       

# The corrected covariance rather than the model's own.
vcov(wrapped)
#>             (Intercept)         z
#> (Intercept)    0.646421 -0.646421
#> z             -0.646421  1.292842

# Everything else about the model reaches its own methods.
coef(wrapped)
#>  (Intercept)            z 
#> 3.919028e-01 1.012353e-15 
```
