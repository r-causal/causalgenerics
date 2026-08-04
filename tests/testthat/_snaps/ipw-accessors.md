# vcov() errors when no covariance was attached

    Code
      vcov(res)
    Condition
      Error in `vcov.ipw()`:
      ! This `ipw` result records no covariance of the effects it reports; refit it with a current version of the package that produced it.

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

