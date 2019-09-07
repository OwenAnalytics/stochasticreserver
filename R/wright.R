#if (model=="Wright") {
#
# g itself
# Basic design is for g to be a function of a single parameter vector, however
# in the simulations it is necessary to work on a matrix of parameters, one
# row for each simulated parameter, so g.obj must be flexible enough to handle
# both.
# Here g.obj is Wright's operational time model with separate level by
# accident year
#' @export
wright <- function(theta, tau, B0, ptd, msk) {
  g.obj = function(theta) {
    if (is.vector(theta))
    {
      exp(array(theta[1:10], c(10, 10)) + theta[11] * tau + theta[12] * tau ^
            2 +
            theta[13] * log(tau))
    }
    else
    {
      exp(
        array(theta[, 1:10], c(nrow(theta), 10, 10)) +
          array(theta[, 11], c(nrow(theta), 10, 10)) *
          aperm(array(tau, c(
            10, 10, nrow(theta)
          )), c(3, 1, 2)) +
          array(theta[, 12], c(nrow(theta), 10, 10)) *
          aperm(array(tau ^ 2, c(
            10, 10, nrow(theta)
          )), c(3, 1, 2)) +
          array(theta[, 13], c(nrow(theta), 10, 10)) *
          aperm(array(log(tau), c(
            10, 10, nrow(theta)
          )), c(3, 1, 2))
      )
    }
  }

  # Gradient of g
  # Note the gradient is a 3-dimensional function of the parameters theta
  # with dimensions 13, 10, 10.  The first dimension
  # represents the parameters involved in the derivatives
  g.grad = function(theta) {
    abind(array(outer(1:10, 1:10, "=="), c(10, 10, 10)),
          abind(tau, abind(tau ^ 2, log(tau),
                           along = 0.5),
                along = 1),
          along = 1) * aperm(array(g.obj(theta), c(10, 10, 13)), c(3, 1, 2))
  }

  # Hessian of g
  # Note the Hessian is a 4-dimensional function of the parameters theta
  # with dimensions 13, 13, 10, 10.  First two dimensions
  # represent the parameters involved in the partial derivatives
  g.hess = function(theta)  {
    aa = array(abind(array(outer(1:10, 1:10, "=="), c(10, 10, 10)),
                     abind(
                       tau, abind(tau ^ 2, log(tau),
                                  along = 0.5),
                       along = 1
                     ),
                     along = 1),
               c(13, 10, 10, 13))
    aperm(aa, c(4, 1, 2, 3)) * aperm(aa, c(1, 4, 2, 3)) *
      aperm(array(g.obj(theta), c(10, 10, 13, 13)), c(3, 4, 1, 2))
  }

  # Base starting values on classic chain ladder forecasts
  ptd = ((!ptd == 0) * ptd) + (ptd == 0) * mean(ptd)
  tmp = c((
    colSums(B0[, 2:10] + 0 * B0[, 1:9], na.rm = TRUE) /
      colSums(B0[, 1:9] + 0 * B0[, 2:10], na.rm = TRUE)
  ),
  1)
  yy = 1 / (cumprod(tmp[11 - (1:10)]))[11 - (1:10)]
  xx = yy - c(0, yy[1:9])
  ww = t(array(xx, c(10, 10)))
  uv = ptd / ((10 == rowSums(msk)) + (10 > rowSums(msk)) * rowSums(msk *
                                                                     ww))
  tmp = na.omit(data.frame(
    x1 = c(tau),
    x2 = c(tau ^ 2),
    x3 = log(c(tau)),
    y = c(log(outer(uv, xx)))
  ))
  ccs = array(coef(lm(tmp$y ~ tmp$x1 + tmp$x2 + tmp$x3)))[1:4]
  a0 = c(log(uv / rowSums(exp(
    ccs[2] * tau + ccs[3] * tau ^ 2 + ccs[4] * log(tau)
  ))), ccs[2:4])
  return(list(
    g.obj = g.obj,
    g.grad = g.grad,
    g.hess = g.hess,
    a0 = a0
  ))
}
