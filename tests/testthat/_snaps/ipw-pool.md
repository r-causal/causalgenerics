# pool_ipw() refuses a `fits` argument that is not a list

    Code
      pool_ipw(1:3)
    Condition
      Error in `pool_ipw()`:
      ! `fits` must be a list of `ipw` results.

# pool_ipw() refuses elements that are not results

    Code
      pool_ipw(fits[[1]])
    Condition
      Error in `pool_ipw()`:
      ! `fits` must be a list of `ipw` results rather than a single result.

---

    Code
      pool_ipw(one_bad)
    Condition
      Error in `pool_ipw()`:
      ! `fits` must hold `ipw` results throughout, but the element at position 2 is not one.

---

    Code
      pool_ipw(two_bad)
    Condition
      Error in `pool_ipw()`:
      ! `fits` must hold `ipw` results throughout, but the elements at positions 2 and 3 are not.

# pool_ipw() refuses fewer than two results

    Code
      pool_ipw(list())
    Condition
      Error in `pool_ipw()`:
      ! `fits` must hold at least two results, since the between-imputation variance is estimated from the spread across them, but it holds 0.

---

    Code
      pool_ipw(fits[1])
    Condition
      Error in `pool_ipw()`:
      ! `fits` must hold at least two results, since the between-imputation variance is estimated from the spread across them, but it holds 1.

# pool_ipw() refuses results targeting different estimands

    Code
      pool_ipw(fits)
    Condition
      Error in `pool_ipw()`:
      ! The estimand must be the same in every result, but they report "ate" and "att".

# pool_ipw() refuses results whose standard errors differ in kind

    Code
      pool_ipw(fits)
    Condition
      Error in `pool_ipw()`:
      ! The standard error method must be the same in every result, but they report "linearization" and "mestimation".

# pool_ipw() refuses results recording different presentation modes

    Code
      pool_ipw(fits)
    Condition
      Error in `pool_ipw()`:
      ! The presentation mode must be the same in every result, but they report "marginal" and "conditional".

# pool_ipw() refuses results reporting different effects

    Code
      pool_ipw(fits)
    Condition
      Error in `pool_ipw()`:
      ! The effects reported must be the same in every result, but they report ("rd b vs a", "log(rr) b vs a", "log(or) b vs a", "rd c vs a", "log(rr) c vs a", "log(or) c vs a") and ("rd b vs a", "log(rr) b vs a", "log(or) b vs a", "rd d vs a", "log(rr) d vs a", "log(or) d vs a").

# pool_ipw() refuses stored levels that disagree with no level named

    Code
      pool_ipw(fits)
    Condition
      Error in `pool_ipw()`:
      ! The confidence level must be the same in every result, but they report 0.95 and 0.9.

# pool_ipw() refuses outcome models fitted on different links

    Code
      pool_ipw(fits)
    Condition
      Error in `pool_ipw()`:
      ! The outcome model's link must be the same in every result, but they report "logit" and "identity".

# pool_ipw() refuses an `effects` value that names no reading

    Code
      pool_ipw(fits, effects = "fixed")
    Condition
      Error in `pool_ipw()`:
      ! `effects` must be a single string, either "marginal" or "conditional".

# pool_ipw() refuses a `conf_level` that is not a probability

    Code
      pool_ipw(fits, conf_level = 2)
    Condition
      Error in `pool_ipw()`:
      ! `conf_level` must be a single number greater than 0 and less than 1.

# pool_ipw() refuses a dfcom that is not a count

    Code
      pool_ipw(fits, dfcom = "18")
    Condition
      Error in `pool_ipw()`:
      ! `dfcom` must be a single number, or `Inf` to assume a large sample.

# pool_ipw() refuses arguments it has no name for

    Code
      pool_ipw(fits, "conditional")
    Condition
      Error in `pool_ipw()`:
      ! `...` must be empty, since every argument after `fits` is matched by name, but an unnamed argument was passed.

---

    Code
      pool_ipw(fits, effect = "conditional")
    Condition
      Error in `pool_ipw()`:
      ! `...` must be empty, since every argument after `fits` is matched by name, but `effect` was passed.

# pool_ipw() assumes a large sample when nothing reports a df

    Code
      pooled <- pool_ipw(fits)
    Condition
      Warning in `pool_ipw()`:
      No result reports the residual degrees of freedom of its complete-data analysis, so a large sample is assumed; pass `dfcom` to name the count the analyses were fitted with.

# pool_ipw() refuses the conditional mode without corrected blocks

    Code
      pool_ipw(fits, effects = "conditional")
    Condition
      Error in `pool_ipw()`:
      ! The conditional reading reports the covariance the joint estimation of the weights and the outcome implies, and this result's outcome model records none; the package that produced the result attaches one by wrapping the model with `new_ipw_model()`.

