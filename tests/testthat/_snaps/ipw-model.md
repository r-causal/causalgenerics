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

