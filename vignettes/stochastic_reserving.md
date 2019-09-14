---
title: "Example Stochastic Reserving"
author: "Roger Hayne"
date: "9/6/2019"
output: 
   - rmarkdown::pdf_document
vignette: >
  %\VignetteEngine{knitr::knitr}
  %\VignetteIndexEntry{Stochastic Reserving}
  %\usepackage[UTF-8]{inputenc}
---




```r
library(mvtnorm)
library(MASS)
library(abind)
library(stochasticreserver)
```

## Initialize Triangle

Input (B0) is a development array of cumulative averages with a the
 exposures (claims) used in the denominator appended as the last column.
 Assumption is for the same development increments as exposure increments
 and that all development lags with no development have # been removed.
Data elements that are not available are indicated as such.  This should
work (but not tested for) just about any subset of an upper triangular
data matrix.  
Another requirement of this code is that the matrix contain
no columns that are all zero.


```r
B0 = matrix(c(670.25868,1480.24821,1938.53579,2466.25469,2837.84888,3003.52391,
            3055.38674,3132.93838,3141.18638,3159.72524,
            767.98833,1592.50266,2463.79447,3019.71976,3374.72689,3553.61387,3602.27898,
            3627.28386,3645.5656,NA,
            740.57952,1615.79681,2345.85028,2910.52511,3201.5226,3417.71335,3506.58672,
            3529.00243,NA,NA,
            862.11956,1754.90405,2534.77727,3270.85361,3739.88962,4003.00219,4125.30694,
            NA,NA,NA,
            840.94172,1859.02531,2804.54535,3445.34665,3950.47098,4185.95298,NA,NA,NA,NA,
            848.00496,2052.922,3076.13789,3861.03111,4351.57694,NA,NA,NA,NA,NA,
            901.77403,1927.88718,3003.58919,3881.41744,NA,NA,NA,NA,NA,NA,
            935.19866,2103.97736,3181.75054,NA,NA,NA,NA,NA,NA,NA,
            759.32467,1584.91057,NA,NA,NA,NA,NA,NA,NA,NA,
            723.30282,NA,NA,NA,NA,NA,NA,NA,NA,NA),10,10,byrow = TRUE)
dnom = c(39161.,38672.4628,41801.048,42263.2794,41480.8768,40214.3872,43598.5056,
       42118.324,43479.4248,49492.4106)

size <- nrow(B0)
# Identify model to be used
#   Berquist for the Berquist-Sherman Incremental Severity
#   CapeCod for the Cape Cod
#   Hoerl for the Generalized Hoerl Curve Model with trend
#   Wright for the Generalized Hoerl Curve with individual accident year levels
#   Chain for the Chain Ladder model
#model = "Berquist"
#model = "CapeCod"
model = "Hoerl"
#model = "Wright"
#model = "Chain"
# Toggle graphs off if desired
graphs = TRUE

# Toggle simulations off if desired
simulation = TRUE
```

```r
# Set tau to have columns with entries 1 through size
tau = t(array((1:size), c(size, size)))

# Calculate incremental average matrix
A0 = cbind(B0[, 1], (B0[, (2:size)] + 0 * B0[, (1:(size - 1))]) -
             (B0[, (1:(size - 1))] + 0 * B0[, (2:size)]))

# Generate a matrix to reflect exposure count in the variance structure
logd = log(matrix(dnom, size, size))

# Set up matrix of rows and columns, makes later calculations simpler
rowNum = row(A0)
colNum = col(A0)

#. msk is a mask matrix of allowable data, upper triangular assuming same
#' development increments as exposure increments
#' msn is a mask matrix that picks off the first forecast diagonal
#' msd is a mask matrix that picks off the to date diagonal
msk = (size - rowNum) >= colNum - 1
msn = (size - rowNum) == colNum - 2
msd = (size - rowNum) == colNum - 1

# Amount paid to date
paid_to_date = rowSums(B0 * msd, na.rm = TRUE)
```

## START OF MODEL SPECIFIC CODE 


```r
  if (model == "Berquist") {
    model_lst <- berquist(tau, B0, paid_to_date, msk)
  } else if (model == "CapeCod") {
    model_lst <- capecod(tau, B0, paid_to_date, msk)
  } else if (model == "Hoerl") {
    model_lst <- hoerl(tau, B0, paid_to_date, msk)
  } else if (model == "Wright") {
    model_lst <- wright(tau, B0, paid_to_date, msk)
  } else if (model == "Chain") {
    model_lst <- chain(tau, B0, paid_to_date, msk)
  }
g_obj <- model_lst$g_obj
g_grad <- model_lst$g_grad
g_hess <- model_lst$g_hess
a0 <- model_lst$a0
```

## Negative Loglikelihood Function to be Minimized

Note that the general
form of the model has parameters in addition to those in the loss model,
namely the power for the variance and the constant of proprtionality that
varies by column.  So if the original model has k parameters with size
columns of data, the total objective function has k + size + 1 parameters



```r
l.obj = function(a, A, g_obj) {
  npar = length(a) - 2
  e = g_obj(a[1:npar])
  v = exp(-outer(logd[, 1], rep(a[npar + 1], size), "-")) * (e^2)^a[npar + 2]
  t1 = log(2 * pi * v) / 2
  t2 = (A - e) ^ 2 / (2 * v)
  sum(t1 + t2, na.rm = TRUE)
}
# Gradient of the objective function
l.grad = function(a, A, g_obj, g_grad) {
  npar = length(a) - 2
  p = a[npar + 2]
  Av = aperm(array(A, c(size, size, npar)), c(3, 1, 2))
  e = g_obj(a[1:npar])
  ev = aperm(array(e, c(size, size, npar)), c(3, 1, 2))
  v = exp(-outer(logd[, 1], rep(a[npar + 1], size), "-")) * (e^2)^p
  vv = aperm(array(v, c(size, size, npar)), c(3, 1, 2))
  dt = rowSums(g_grad(a[1:npar]) * ((p / ev) + (ev - Av) / vv - p * 
                                      (Av - ev)^2 / (vv * ev)),
               na.rm = TRUE,
               dims = 1)
  yy = 1 - (A - e) ^ 2 / v
  dk = sum(yy / 2, na.rm = TRUE)
  dp = sum(yy * log(e ^ 2) / 2, na.rm = TRUE)
  c(dt, dk, dp)
}
```

## Hessian of the objective function

-   e is the expectated value matrix
-   v is the matrix of variances
-   A, e, v all have shape c(size, size)
-   The variables _v are copies of the originals to shape c(npar,size,size), paralleling
    the gradient of g.
-   The variables _m are copies of the originals to shape c(npar,npar,size,size),
    paralleling the hessian of g


```r
l.hess = function(a, A, g_obj, g_grad, g_hess) {
  npar = length(a) - 2
  p = a[npar + 2]
  Av = aperm(array(A, c(size, size, npar)), c(3, 1, 2))
  Am = aperm(array(A, c(size, size, npar, npar)), c(3, 4, 1, 2))
  e = g_obj(a[1:npar])
  ev = aperm(array(e, c(size, size, npar)), c(3, 1, 2))
  em = aperm(array(e, c(size, size, npar, npar)), c(3, 4, 1, 2))
  v = exp(-outer(logd[, 1], rep(a[npar + 1], size), "-")) * (e ^ 2) ^ p
  vv = aperm(array(v, c(size, size, npar)), c(3, 1, 2))
  vm = aperm(array(v, c(size, size, npar, npar)), c(3, 4, 1, 2))
  g1 = g_grad(a[1:npar])
  gg = aperm(array(g1, c(npar, size, size, npar)), c(4, 1, 2, 3))
  gg = gg * aperm(gg, c(2, 1, 3, 4))
  gh = g_hess(a[1:npar])
  dtt = rowSums(
    gh * (p / em + (em - Am) / vm - p * (Am - em) ^ 2 / (vm * em)) +
      gg * (
        1 / vm + 4 * p * (Am - em) / (vm * em) + p * (2 * p + 1) * (Am - em) ^ 2 /
          (vm * em ^ 2) - p / em ^ 2
      ),
    dims = 2,
    na.rm = TRUE
  )
  dkt = rowSums((g1 * (Av - ev) + p * g1 * (Av - ev) ^ 2 / ev) / vv, na.rm = TRUE)
  dtp = rowSums(g1 * (1 / ev + (
    log(ev ^ 2) * (Av - ev) + (p * log(ev ^ 2) - 1) * (Av - ev) ^ 2 / ev
  ) / vv),
  na.rm = TRUE)
  dkk = sum((A - e) ^ 2 / (2 * v), na.rm = TRUE)
  dpk = sum(log(e ^ 2) * (A - e) ^ 2 / (2 * v), na.rm = TRUE)
  dpp = sum(log(e ^ 2) ^ 2 * (A - e) ^ 2 / (2 * v), na.rm = TRUE)
  m1 = rbind(array(dkt), c(dtp))
  rbind(cbind(dtt, t(m1)), cbind(m1, rbind(cbind(dkk, c(
    dpk
  )), c(dpk, dpp))))
}
```


End of funciton specificaitons now on to the minimization

## Minimization

### Get starting values for kappa and p parameters, default 10 and 1

```r
ttt = c(10, 1)
```

For starting values use fitted objective function and assume variance for a
cell is estimated by the square of the difference between actual and expected
averages.  Note since log(0) is -Inf we need to go through some machinations
to prep the y values for the fit


```r
E = g_obj(a0)
yyy = (A0 - E)^2
yyy = logd + log(((yyy != 0) * yyy) - (yyy == 0))
sss = na.omit(data.frame(x = c(log(E^2)), y = c(yyy)))
ttt = array(coef(lm(sss$y ~ sss$x)))[1:2]
a0 = c(a0, ttt)

set.seed(1) # to check reproducibility with original code
max = list(iter.max = 10000, eval.max = 10000)
```

### Actual minimization































