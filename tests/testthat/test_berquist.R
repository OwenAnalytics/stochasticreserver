context("berquist")
berquist_yyy <- structure(c(18.5903477820371, 17.0857629741298, 16.4587546640197,
                            19.5555511727509, 18.3251047554907, 17.7498111021868, 19.2022607846289,
                            19.5371976228424, 20.0849340151416, 21.0357586325281, 18.2296970960545,
                            18.4280430363886, 16.6735693713869, 16.956730057397, 19.4191436997024,
                            21.6144754167453, 18.3151895803851, 20.8674197297801, 21.282675587742,
                            NA, 22.0279418095115, 19.5319868588964, 19.1534144589673, 17.98379475847,
                            19.9928755730403, 20.8288904829514, 21.2701214231012, 21.0456346056213,
                            NA, NA, 19.2817602076264, 18.8619631743866, 19.1048634416956,
                            19.6065402998731, 16.8270074708046, 19.9430372932894, 21.114109011366,
                            NA, NA, NA, 15.9342800754823, 17.9417192825488, 20.0970835134835,
                            18.6908252201051, 19.4533495293854, 18.7231044654091, NA, NA,
                            NA, NA, 18.005245792937, 17.5305627296255, 6.58288886618252,
                            18.1188262737508, 15.0543948069065, NA, NA, NA, NA, NA, 17.1805821794974,
                            17.5064005751957, 14.2769112609376, 17.9106964171237, NA, NA,
                            NA, NA, NA, NA, 17.568301918838, 16.6106667961334, 17.0138589184648,
                            NA, NA, NA, NA, NA, NA, NA, 14.2075512576392, 13.1001374511748,
                            NA, NA, NA, NA, NA, NA, NA, NA, 12.8989360791041, NA, NA, NA,
                            NA, NA, NA, NA, NA, NA), .Dim = c(10L, 10L))

berquist_E <- structure(c(725.265411834231, 741.901255755787, 758.918685919382,
                          776.326454995502, 794.13351641984, 812.349028998367, 830.982361618031,
                          850.043098065509, 869.541041956491, 889.486221778029, 855.920068918206,
                          875.552816383496, 895.635891849066, 916.179624755979, 937.194581477828,
                          958.691570755399, 980.68164925599, 1003.17612726025, 1026.18657447947,
                          1049.72482600626, 765.104920510847, 782.654586927508, 800.606800476046,
                          818.970794619338, 837.756014613689, 856.972122366866, 876.629001407553,
                          896.736761968812, 917.30574618813, 938.346533426747, 605.442721043376,
                          619.330120671669, 633.536063841295, 648.067857174836, 662.932974891221,
                          678.139062649979, 693.69394148368, 709.605611820565, 725.882257599456,
                          742.532250479044, 386.170851487497, 395.028681886652, 404.089689607652,
                          413.358535049518, 422.839985509613, 432.538917635635, 442.460319933846,
                          452.609295334845, 462.991063818184, 473.610965097186, 206.729699468256,
                          211.471581485774, 216.322230870176, 221.28414248322, 226.359868412788,
                          231.552019285523, 236.863265609559, 242.296339148066, 247.854034324288,
                          253.539209658816, 79.0453134069634, 80.8584227530207, 82.7131204660421,
                          84.6103604831249, 86.5511186223889, 88.5363930848739, 90.5672049679506,
                          92.6445987905072, 94.7696430301826, 96.9434306729222, 44.5541122125271,
                          45.5760763717672, 46.6214819304856, 47.6908665779979, 48.7847803369194,
                          49.903785846061, 51.0484586498147, 52.219387494176, 53.4171746295578,
                          54.6424361205483, 14.3955726480275, 14.7257724560238, 15.0635462533199,
                          15.4090677690111, 15.7625147171203, 16.1240688880025, 16.4939162418462,
                          16.8722470043192, 17.2592557644087, 17.6551415745052, 21.7343796624638,
                          22.2329140498713, 22.7428836169014, 23.2645506590675, 23.7981834883107,
                          24.3440565710031, 24.9024506691147, 25.4736529846201, 26.057957307216,
                          26.6556641654282), .Dim = c(10L, 10L))
berquist_a0 <- c(709.002597747623, 836.727551631865, 747.948774807513, 591.866722172126,
                 377.511642512502, 202.094145898798, 77.2728598811514, 43.5550638201962,
                 14.0727769958082, 21.247024040686, 0.0226784826650085, 10.1469982236692,
                 0.674215930439166)
B0 <- stochasticreserver::B0
A0 <- stochasticreserver::A0
dnom <- stochasticreserver::dnom
size <- nrow(B0)
# Set up matrix of rows and columns, makes later calculations simpler
rowNum = row(A0)
colNum = col(A0)
# Generate a matrix to reflect exposure count in the variance structure
logd = log(matrix(dnom, size, size))

#. upper_triangle_mask is a mask matrix of allowable data, upper triangular assuming same
#' development increments as exposure increments
#' msn is a mask matrix that picks off the first forecast diagonal
#' msd is a mask matrix that picks off the to date diagonal
upper_triangle_mask = (size - rowNum) >= colNum - 1
msd = (size - rowNum) == colNum - 1

# Amount paid to date
paid_to_date = rowSums(B0 * msd, na.rm = TRUE)
model_lst <- berquist(B0, paid_to_date, upper_triangle_mask)
g_obj <- model_lst$g_obj
g_grad <- model_lst$g_grad
g_hess <- model_lst$g_hess
a0 <- model_lst$a0

## Negative Loglikelihood Function to be Minimized

# Note that the general
# form of the model has parameters in addition to those in the loss model,
# namely the power for the variance and the constant of proprtionality that
# varies by column.  So if the original model has k parameters with size
# columns of data, the total objective function has k + size + 1 parameters


l.obj = function(a, A) {
  npar = length(a) - 2
  e = g_obj(a[1:npar])
  v = exp(-outer(logd[, 1], rep(a[npar + 1], size), "-")) * (e^2)^a[npar + 2]
  t1 = log(2 * pi * v) / 2
  t2 = (A - e) ^ 2 / (2 * v)
  sum(t1 + t2, na.rm = TRUE)
}
# Gradient of the objective function
l.grad = function(a, A) {
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

## Hessian of the objective function

# -   e is the expectated value matrix
# -   v is the matrix of variances
# -   A, e, v all have shape c(size, size)
# -   The variables _v are copies of the originals to shape c(npar,size,size), paralleling
# the gradient of g.
# -   The variables _m are copies of the originals to shape c(npar,npar,size,size),
# paralleling the hessian of g


l.hess = function(a, A) {
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


# End of funciton specificaitons now on to the minimization

## Minimization

### Get starting values for kappa and p parameters, default 10 and 1

ttt = c(10, 1)

# For starting values use fitted objective function and assume variance for a
# cell is estimated by the square of the difference between actual and expected
# averages.  Note since log(0) is -Inf we need to go through some machinations
# to prep the y values for the fit

# ```{r minimization-setup}
E = g_obj(a0)
yyy = (A0 - E)^2
yyy = logd + log(((yyy != 0) * yyy) - (yyy == 0))
sss = na.omit(data.frame(x = c(log(E^2)), y = c(yyy)))
ttt = array(coef(lm(sss$y ~ sss$x)))[1:2]
a0 = c(a0, ttt)

test_that("berquist returns the correct vector", {
  expect_equal(berquist_a0, a0)
})
test_that("berquist functions are correctly formed", {
  expect_equal(yyy, berquist_yyy)
  expect_equal(E, berquist_E)
})

