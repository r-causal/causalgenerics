# the effects error states the contract

    Code
      ipw_with_effects("everything")
    Condition
      Error in `new_ipw()`:
      ! `effects` must be a single string, either "marginal" or "conditional".

---

    Code
      ipw_with_effects(c("marginal", "conditional"))
    Condition
      Error in `new_ipw()`:
      ! `effects` must be a single string, either "marginal" or "conditional".

# print() summarizes a binary-exposure result

    Code
      print(res)
    Output
      Inverse Probability Weight Estimator
      Estimand: ATE 
      Effects: marginal (population-averaged) 
      
      Weight Estimator:
        Call: glm(formula = z ~ x, family = binomial(), data = dat) 
      
      Outcome Model:
        Call: glm(formula = y ~ z, family = quasibinomial(), data = dat) 
      
      Marginal estimates:
              estimate  std.err      z ci.lower ci.upper conf.level p.value  
      rd      0.199882 0.092425 2.1626 0.018732  0.38103       0.95 0.03057 *
      log(rr) 0.560414 0.273519 2.0489 0.024326  1.09650       0.95 0.04047 *
      log(or) 0.878313 0.418661 2.0979 0.057753  1.69887       0.95 0.03591 *
      ---
      Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# print() keys rows by effect and comparison for a categorical result

    Code
      print(res)
    Output
      Inverse Probability Weight Estimator
      Estimand: ATE 
      Effects: marginal (population-averaged) 
      
      Weight Estimator:
        Call: glm(formula = z ~ x, family = binomial(), data = dat) 
      
      Outcome Model:
        Call: glm(formula = y ~ z, family = quasibinomial(), data = dat) 
      
      Marginal estimates:
                     estimate  std.err      z  ci.lower ci.upper conf.level   p.value
      rd b vs a      0.081945 0.050387 1.6263 -0.016811  0.18070       0.95 0.1038822
      log(rr) b vs a 0.168870 0.104633 1.6139 -0.036207  0.37395       0.95 0.1065433
      log(or) b vs a 0.328762 0.203058 1.6191 -0.069225  0.72675       0.95 0.1054363
      rd c vs a      0.166939 0.045182 3.6948  0.078384  0.25549       0.95 0.0002200
      log(rr) c vs a 0.318293 0.091898 3.4635  0.138176  0.49841       0.95 0.0005331
      log(or) c vs a 0.676435 0.185786 3.6409  0.312300  1.04057       0.95 0.0002717
                        
      rd b vs a         
      log(rr) b vs a    
      log(or) b vs a    
      rd c vs a      ***
      log(rr) c vs a ***
      log(or) c vs a ***
      ---
      Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# print() summarizes a continuous-outcome result

    Code
      print(res)
    Output
      Inverse Probability Weight Estimator
      Estimand: ATE 
      Effects: marginal (population-averaged) 
      
      Weight Estimator:
        Call: glm(formula = z ~ x, family = binomial(), data = dat) 
      
      Outcome Model:
        Call: glm(formula = y ~ z, family = quasibinomial(), data = dat) 
      
      Marginal estimates:
           estimate std.err      z ci.lower ci.upper conf.level   p.value    
      diff  2.25255 0.17524 12.854   1.9091    2.596       0.95 < 2.2e-16 ***
      ---
      Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# print() names the conditional reading in its metadata and its table

    Code
      print(res)
    Output
      Inverse Probability Weight Estimator
      Estimand: ATE 
      Effects: conditional (outcome model) 
      
      Weight Estimator:
        Call: glm(formula = z ~ x, family = binomial(), data = dat) 
      
      Outcome Model:
        Call: glm(formula = y ~ z, family = quasibinomial(), data = dat) 
      
      Conditional estimates (outcome model):
                  Estimate Std. Error z value Pr(>|z|)
      (Intercept) -0.84730    0.86066 -0.9845   0.3249
      z            1.69460    1.21716  1.3923   0.1638

# print() pairs each coefficient with its own standard error

    Code
      print(res)
    Output
      Inverse Probability Weight Estimator
      Estimand: ATE 
      Effects: conditional (outcome model) 
      
      Weight Estimator:
        Call: glm(formula = z ~ x, family = binomial(), data = dat) 
      
      Outcome Model:
        Call: glm(formula = y ~ z, family = quasibinomial(), data = dat) 
      
      Conditional estimates (outcome model):
                  Estimate Std. Error z value  Pr(>|z|)    
      (Intercept)  -0.8473     0.5000 -1.6946   0.09015 .  
      z             1.6946     0.2500  6.7784 1.215e-11 ***
      ---
      Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# print() reports a conditional result with no corrected covariance

    Code
      print(res)
    Output
      Inverse Probability Weight Estimator
      Estimand: ATE 
      Effects: conditional (outcome model) 
      
      Weight Estimator:
        Call: glm(formula = z ~ x, family = binomial(), data = dat) 
      
      Outcome Model:
        Call: glm(formula = y ~ z, family = quasibinomial(), data = dat) 
      
      Conditional estimates (outcome model):
                  Estimate
      (Intercept)  -0.8473
      z             1.6946
      
      Standard errors are not reported: this result's outcome model records
      no covariance from the joint estimation of the weights and the outcome.
      The package that produced it attaches one by wrapping the model with
      `new_ipw_model()`.

# print() says a conditional result has no coefficients to tabulate

    Code
      print(res)
    Output
      Inverse Probability Weight Estimator
      Estimand: ATE 
      Effects: conditional (outcome model) 
      
      Weight Estimator:
        Call: glm(formula = z ~ x, family = binomial(), data = dat) 
      
      Outcome Model:
        Call: <cg_no_coef> 
      
      Conditional estimates (outcome model):
      The outcome model reports no coefficients, so there is no
      conditional table to print.

# print() refuses a stored mode that is not a reading

    Code
      print(res)
    Condition
      Error in `print.ipw()`:
      ! `effects` must be a single string, either "marginal" or "conditional".

# the as.data.frame() argument errors state the contract

    Code
      as.data.frame(res, conf.level = 95)
    Condition
      Error in `as.data.frame.ipw()`:
      ! `conf.level` must be a single number greater than 0 and less than 1.

---

    Code
      as.data.frame(res, conf.int = NA)
    Condition
      Error in `as.data.frame.ipw()`:
      ! `conf.int` must be a single logical value, either `TRUE` or `FALSE`.

---

    Code
      as.data.frame(res, exponentiate = 1)
    Condition
      Error in `as.data.frame.ipw()`:
      ! `exponentiate` must be a single logical value, either `TRUE` or `FALSE`.

