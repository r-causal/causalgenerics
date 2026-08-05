# vcov() errors when no covariance was attached

    Code
      vcov(res)
    Condition
      Error in `vcov.ipw()`:
      ! This `ipw` object records no covariance to report; the package that produced it attaches one when it supports the `ipw_vcov` contract.

# confint() errors on labels the result does not have

    Code
      confint(res, parm = "rr")
    Condition
      Error in `confint.ipw()`:
      ! `parm` must name effects the result reports, which are "rd", "log(rr)", "log(or)".

---

    Code
      confint(res, parm = c("rd", "rr"))
    Condition
      Error in `confint.ipw()`:
      ! `parm` must name effects the result reports, which are "rd", "log(rr)", "log(or)".

# estimand<-() has no method for a result

    Code
      estimand(res) <- "att"
    Condition
      Error in `estimand<-.default`:
      ! No `estimand<-()` method for an object of class <ipw>.

# vcov() refuses the conditional mode without a corrected block

    Code
      vcov(res)
    Condition
      Error in `vcov.ipw()`:
      ! The conditional reading reports the covariance the joint estimation of the weights and the outcome implies, and this result's outcome model records none; the package that produced the result attaches one by wrapping the model with `new_ipw_model()`.

# confint() refuses a conditional parm the outcome model lacks

    Code
      confint(res, parm = "rd")
    Condition
      Error in `confint.ipw()`:
      ! `parm` must name coefficients the outcome model reports, which are "(Intercept)", "z".

# confint() refuses the conditional mode without a corrected block

    Code
      confint(res)
    Condition
      Error in `confint.ipw()`:
      ! The conditional reading reports the covariance the joint estimation of the weights and the outcome implies, and this result's outcome model records none; the package that produced the result attaches one by wrapping the model with `new_ipw_model()`.

# the accessors refuse an effects argument that is not a reading

    Code
      coef(res, effects = "banana")
    Condition
      Error in `coef.ipw()`:
      ! `effects` must be a single string, either "marginal" or "conditional".

---

    Code
      vcov(res, effects = 1)
    Condition
      Error in `vcov.ipw()`:
      ! `effects` must be a single string, either "marginal" or "conditional".

---

    Code
      confint(res, effects = c("marginal", "conditional"))
    Condition
      Error in `confint.ipw()`:
      ! `effects` must be a single string, either "marginal" or "conditional".

# the accessors refuse a stored mode that is not a reading

    Code
      coef(res)
    Condition
      Error in `coef.ipw()`:
      ! `effects` must be a single string, either "marginal" or "conditional".

