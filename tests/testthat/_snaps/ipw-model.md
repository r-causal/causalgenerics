# vcov() errors when a wrapped model carries no covariance

    Code
      vcov(stripped)
    Condition
      Error in `vcov.ipw_model()`:
      ! This `ipw_model` object records no covariance to report; the package that produced it attaches one when it supports the `ipw_vcov` contract.

# the covariance error states the contract

    Code
      new_ipw_model(mod, c(theta1 = 1, theta2 = 1))
    Condition
      Error in `new_ipw_model()`:
      ! `vcov` must be a square numeric matrix.

---

    Code
      new_ipw_model(mod, matrix(c(1, 0, 0, 1), nrow = 2))
    Condition
      Error in `new_ipw_model()`:
      ! `vcov` must have dimnames on both margins.

