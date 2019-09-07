#' Create list for Berquist-Sherman incremental severity model
#'
#' g - Assumed loss emergence model, a function of the parameters a.
#' Note g must be matrix-valued with 10 rows and 10 columns
#'
#' g itself
#' Basic design is for g to be a function of a single parameter vector, however
#' in the simulations it is necessary to work on a matrix of parameters, one
#' row for each simulated parameter, so g.obj must be flexible enough to handle
#' both.
#' Here g.obj is the Berquist-Sherman incremental severity model
#' @param theta do not know
#' @param tau do not know
#' @param B0 development triangle
#' @param ptd do not know
#' @param msk mask for triangle
#'
#' @importFrom stats coef lm na.omit
#' @import abind
#' @export
berquist <- function(tau, B0, ptd, msk) {
  g.obj = function(theta) {
    if (is.vector(theta))
    {
      outer(exp(theta[11] * (1:10)), theta[1:10])
    }
    else
    {
      array(exp(outer(theta[, 11], (1:10))), c(nrow(theta), 10, 10)) *
        aperm(array(theta[, (1:10)], c(nrow(theta), 10, 10)), c(1, 3, 2))
    }
  }

  # Gradient of g
  # Note the gradient is a 3-dimensional function of the parameters theta
  # with dimensions 11 (=length(theta)), 10, 10.  The first dimension
  # represents the parameters involved in the derivatives
  g.grad = function(theta) {
    abind(aperm(array(rep(
      exp(theta[11] * (1:10)), 10 * 10 * 10
    ), c(10, 10, 10)), c(2, 1, 3)) *
      outer((1:10), outer(rep(1, 10), (1:10)), "=="),
    outer((1:10) * exp(theta[11] * (1:10)), theta[1:10]),
    along = 1)
  }

  # Hessian of g
  # Note the Hessian is a 4-dimensional function of the parameters theta
  # with dimensions 11 (=length(theta)), 11, 10, 10.  First two dimensions
  # represent the parameters involved in the partial derivatives
  g.hess = function(theta)  {
    aa = aperm(outer(diag(rep(1, 10)),
                     array((1:10) * exp((
                       1:10
                     ) * theta[11]), c(10, 1))), c(1, 4, 3, 2))
    abind(abind(array(0, c(10, 10, 10, 10)), aa, along = 2),
          abind(aperm(aa, c(2, 1, 3, 4)),
                array(outer((1:10) ^ 2 * exp((1:10) * theta[11]), theta[(1:10)]
                ), c(1, 1, 10, 10)),
                along = 2),
          along = 1)
  }

  # Set up starting values.  Essentially start with classical chain ladder
  # ultimate estimates and estimate trend from that and incremental average
  # start values based on incrementals from classic chain ladder
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
  tmp = na.omit(data.frame(x = 1:10, y = log(uv)))
  trd = 0.01
  trd = array(coef(lm(tmp$y ~ tmp$x))[2])[1]
  a0 = c((xx * mean(uv / (exp(trd)^(1:10)))), trd)
  return(list(g.obj = g.obj, g.grad = g.grad, g.hess = g.hess, a0 = a0))
}
