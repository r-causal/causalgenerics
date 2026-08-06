# print() summarizes a pooled binary result

    Code
      print(res)
    Output
      Pooled Inverse Probability Weight Estimator
      Estimand: ATE 
      Effects: marginal (population-averaged) 
      Imputations: 3 
      Complete-data df: 17 
      
      Pooled marginal estimates:
              estimate std.err      t      df  ci.lower ci.upper conf.level p.value  
      rd       0.30000 0.24495 1.2247  9.1975 -0.252304   0.8523       0.95 0.25111  
      log(rr)  0.60000 0.30822 1.9467 14.6302 -0.058406   1.2584       0.95 0.07104 .
      log(or)  0.95000 0.46637 2.0370 13.9847 -0.050365   1.9504       0.95 0.06104 .
      ---
      Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# print() keys pooled rows by effect and contrast

    Code
      print(res)
    Output
      Pooled Inverse Probability Weight Estimator
      Estimand: ATE 
      Effects: marginal (population-averaged) 
      Imputations: 3 
      Complete-data df: 17 
      
      Pooled marginal estimates:
                     estimate  std.err      t      df  ci.lower ci.upper conf.level
      rd b vs a      0.100000 0.064807 1.5430 12.0585 -0.041127  0.24113       0.95
      log(rr) b vs a 0.200000 0.115614 1.7299 13.1864 -0.049411  0.44941       0.95
      log(or) b vs a 0.380000 0.217945 1.7436 13.7448 -0.088261  0.84826       0.95
      rd c vs a      0.200000 0.068557 2.9173  4.4904  0.017581  0.38242       0.95
      log(rr) c vs a 0.350000 0.110454 3.1688 10.5821  0.105717  0.59428       0.95
      log(or) c vs a 0.750000 0.211424 3.5474 10.1008  0.279555  1.22045       0.95
                      p.value   
      rd b vs a      0.148648   
      log(rr) b vs a 0.106971   
      log(or) b vs a 0.103548   
      rd c vs a      0.037681 * 
      log(rr) c vs a 0.009355 **
      log(or) c vs a 0.005210 **
      ---
      Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# print() names the conditional reading of a pooled result

    Code
      print(res)
    Output
      Pooled Inverse Probability Weight Estimator
      Estimand: ATE 
      Effects: conditional (outcome model) 
      Imputations: 3 
      Complete-data df: 18 
      
      Pooled conditional estimates (outcome model):
                  estimate  std.err       t     df ci.lower ci.upper conf.level
      (Intercept)  0.41759  0.90068  0.4636 7.6522  -1.6759   2.5111       0.95
      z           -0.14728  1.19411 -0.1233 9.8813  -2.8123   2.5177       0.95
                  p.value
      (Intercept)  0.6558
      z            0.9043

# as.data.frame() refuses to exponentiate an identity-link table

    Code
      as.data.frame(res, exponentiate = TRUE)
    Condition
      Error in `as.data.frame.ipw_pooled()`:
      ! `exponentiate` needs coefficients on a scale an exponential undoes, and the outcome models were fitted with the "identity" link, whose coefficients are not on one; only the "logit" and "log" links exponentiate in the conditional reading.

