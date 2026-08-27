# as_conditional() refuses a result that supports marginal only

    Code
      as_conditional(res)
    Condition
      Error in `as_conditional.ipw()`:
      ! This result supports the marginal reading only, so there is no conditional reading of it to report; the package that produced it records the readings it supports when it builds the result.

# as_marginal() refuses a result that supports conditional only

    Code
      as_marginal(res)
    Condition
      Error in `as_marginal.ipw()`:
      ! This result supports the conditional reading only, so there is no marginal reading of it to report; the package that produced it records the readings it supports when it builds the result.

# a pooled result refuses the reading it could not pool

    Code
      as_conditional(res)
    Condition
      Error in `as_conditional.ipw_pooled()`:
      ! This pooled result carries no conditional reading, since pooling that reading over the same results was refused. The conditional reading reports the covariance the joint estimation of the weights and the outcome implies, and this result's outcome model records none; the package that produced the result attaches one by wrapping the model with `new_ipw_model()`.

# a pooled result with no other reading stored refuses its own way

    Code
      as_conditional(legacy)
    Condition
      Error in `as_conditional.ipw_pooled()`:
      ! This pooled result carries no conditional reading: it was pooled before both readings were kept, so there is nothing recorded to present; pooling the results again gives a result that carries both.

# the no-method error names the generic and the class

    Code
      as_marginal(x)
    Condition
      Error in `as_marginal.default()`:
      ! No `as_marginal()` method for an object of class <cg_unregistered>.

---

    Code
      as_conditional(x)
    Condition
      Error in `as_conditional.default()`:
      ! No `as_conditional()` method for an object of class <cg_unregistered>.

