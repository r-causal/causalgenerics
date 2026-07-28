# print() summarizes a binary-exposure result

    Code
      print(res)
    Output
      Inverse Probability Weight Estimator
      Estimand: ATE 
      
      Weight Estimator:
        Call: glm(formula = z ~ x, family = binomial(), data = dat) 
      
      Outcome Model:
        Call: glm(formula = y ~ z, family = quasibinomial(), data = dat) 
      
      Estimates:
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
      
      Weight Estimator:
        Call: glm(formula = z ~ x, family = binomial(), data = dat) 
      
      Outcome Model:
        Call: glm(formula = y ~ z, family = quasibinomial(), data = dat) 
      
      Estimates:
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
      
      Weight Estimator:
        Call: glm(formula = z ~ x, family = binomial(), data = dat) 
      
      Outcome Model:
        Call: glm(formula = y ~ z, family = quasibinomial(), data = dat) 
      
      Estimates:
           estimate std.err      z ci.lower ci.upper conf.level   p.value    
      diff  2.25255 0.17524 12.854   1.9091    2.596       0.95 < 2.2e-16 ***
      ---
      Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

