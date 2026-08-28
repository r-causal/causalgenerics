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

# the readings error states the contract

    Code
      ipw_with_readings("everything")
    Condition
      Error in `new_ipw()`:
      ! `readings` must be a character vector of one or both of "marginal" and "conditional", each named at most once.

---

    Code
      ipw_with_readings(c("marginal", "marginal"))
    Condition
      Error in `new_ipw()`:
      ! `readings` must be a character vector of one or both of "marginal" and "conditional", each named at most once.

---

    Code
      ipw_with_readings("marginal", effects = "conditional")
    Condition
      Error in `new_ipw()`:
      ! `readings` must include the reading the result records, since a result cannot record a reading it does not support, and this one records the "conditional" reading.

# new_ipw() refuses estimates that are not a data frame

    Code
      ipw_result(as.list(binary_estimates()))
    Condition
      Error in `new_ipw()`:
      ! `estimates` must be a data frame of effect estimates, since every surface of a result reports from its rows and its columns, but it is <list>.

# new_ipw() refuses an estimates frame that names no effects

    Code
      ipw_result(effectless_estimates())
    Condition
      Error in `new_ipw()`:
      ! `estimates` must name each row's effect measure in an `effect` column, since every label a result reports begins with it, but this frame carries none.

# new_ipw() refuses an effect column that names no measures

    Code
      ipw_result(numeric_effect_estimates())
    Condition
      Error in `new_ipw()`:
      ! `estimates` must name its effects in a character or factor `effect` column, since the labels a row is keyed by are read off it as strings, but the column this frame carries is <integer>.

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

# print() keys rows by effect and contrast for a categorical result

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

# as.data.frame() refuses to exponentiate an unexponentiable link

    Code
      as.data.frame(res, exponentiate = TRUE)
    Condition
      Error in `as.data.frame.ipw()`:
      ! `exponentiate` needs coefficients on a scale an exponential undoes, and the outcome models were fitted with the "identity" link, whose coefficients are not on one; only the "logit" and "log" links exponentiate in the conditional reading.

# as.data.frame() refuses a conditional reading with no block

    Code
      as.data.frame(res)
    Condition
      Error in `as.data.frame.ipw()`:
      ! The conditional reading reports the covariance the joint estimation of the weights and the outcome implies, and this result's outcome model records none; the package that produced the result attaches one by wrapping the model with `new_ipw_model()`.

# as.data.frame() refuses a reading the result does not support

    Code
      as.data.frame(res, effects = "conditional")
    Condition
      Error in `as.data.frame.ipw()`:
      ! This result supports the marginal reading only, so there is no conditional reading of it to report; the package that produced it records the readings it supports when it builds the result.

# print() keys rows by the effect and the group together

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
                      estimate  std.err      z  ci.lower ci.upper conf.level p.value
      rd sex = 0      0.151234 0.081422 1.8574 -0.008350  0.31082       0.95 0.06325
      log(rr) sex = 0 0.421887 0.240118 1.7570 -0.048736  0.89251       0.95 0.07892
      log(or) sex = 0 0.664215 0.371244 1.7892 -0.063410  1.39184       0.95 0.07359
      rd sex = 1      0.248531 0.104663 2.3746  0.043395  0.45367       0.95 0.01757
      log(rr) sex = 1 0.698742 0.308951 2.2617  0.093209  1.30428       0.95 0.02372
      log(or) sex = 1 1.092408 0.472183 2.3135  0.166946  2.01787       0.95 0.02069
                       
      rd sex = 0      .
      log(rr) sex = 0 .
      log(or) sex = 0 .
      rd sex = 1      *
      log(rr) sex = 1 *
      log(or) sex = 1 *
      ---
      Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# print() keys rows by effect, contrast, and group together

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
                             estimate  std.err      z  ci.lower ci.upper conf.level
      rd b vs a sex = 0      0.081945 0.050387 1.6263 -0.016812  0.18070       0.95
      log(rr) b vs a sex = 0 0.168870 0.104633 1.6139 -0.036207  0.37395       0.95
      rd c vs a sex = 0      0.166939 0.045182 3.6948  0.078384  0.25549       0.95
      log(rr) c vs a sex = 0 0.318293 0.091898 3.4635  0.138176  0.49841       0.95
      rd b vs a sex = 1      0.062318 0.041205 1.5124 -0.018442  0.14308       0.95
      log(rr) b vs a sex = 1 0.129441 0.085734 1.5098 -0.038595  0.29748       0.95
      rd c vs a sex = 1      0.128507 0.037164 3.4578  0.055667  0.20135       0.95
      log(rr) c vs a sex = 1 0.245106 0.075611 3.2417  0.096911  0.39330       0.95
                               p.value    
      rd b vs a sex = 0      0.1038832    
      log(rr) b vs a sex = 0 0.1065433    
      rd c vs a sex = 0      0.0002200 ***
      log(rr) c vs a sex = 0 0.0005331 ***
      rd b vs a sex = 1      0.1304349    
      log(rr) b vs a sex = 1 0.1310950    
      rd c vs a sex = 1      0.0005445 ***
      log(rr) c vs a sex = 1 0.0011883 ** 
      ---
      Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# new_ipw() refuses a group column that is not character

    Code
      ipw_result(factor_group_estimates())
    Condition
      Error in `new_ipw()`:
      ! `estimates` must name its subgroups in a character `group` column, since a group is written as a "var = value" string, but the column this frame carries is <factor>.

# new_ipw() refuses a group column that leaves a row unnamed

    Code
      ipw_result(unnamed_group_estimates(1L))
    Condition
      Error in `new_ipw()`:
      ! `estimates` must name a subgroup in every row of its `group` column, but 1 row records none.

---

    Code
      ipw_result(unnamed_group_estimates(2L))
    Condition
      Error in `new_ipw()`:
      ! `estimates` must name a subgroup in every row of its `group` column, but 2 rows record none.

