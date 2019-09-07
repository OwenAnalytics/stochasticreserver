#' Create list for Generalized Hoerl Curve Model with trend
#
#' g itself
#' Basic design is for g to be a function of a single parameter vector, however
#' in the simulations it is necessary to work on a matrix of parameters, one
#' row for each simulated parameter, so g.obj must be flexible enough to handle
#' both.
#' Here g.obj is Wright's operational time model with trend added
#' @param tau do not know
#' @param B0 development triangle
#' @param ptd do not know
#' @param msk mask for triangle
#'
#' @importFrom stats coef lm na.omit
#' @import abind
#' @export
hoerl <- function(tau, B0, ptd, msk) {
  size <- nrow(B0)
  g.obj = function(theta) {
    if (is.vector(theta))
    {
      exp(theta[1] +
            colSums(abind(
              tau, abind(tau ^ 2, log(tau), along = 0.5), along = 1
            ) *
              array(theta[c(2, 3, 4)], c(3, size, size))) +
            theta[5] * array((1:size), c(size, size)))
    }
    else
    {
      exp(
        array(theta[, 1], c(nrow(theta), size, size)) +
          colSums(abind(
            aperm(array(tau, c(
              size, size, nrow(theta)
            )), c(3, 1, 2)),
            abind(aperm(array(
              tau ^ 2, c(size, size, nrow(theta))
            ), c(3, 1, 2)),
            aperm(array(
              log(tau), c(size, size, nrow(theta))
            ), c(3, 1, 2)), along = 0.5),
            along = 1
          ) *
            aperm(array(
              theta[, c(2, 3, 4)], c(nrow(theta), 3, size, size)
            ), c(2, 1, 3, 4))) +
          array(theta[, 5], c(nrow(theta), size, size)) *
          aperm(array((1:size), c(
            size, nrow(theta), size
          )), c(2, 1, 3))
      )
    }
  }

  # Gradient of g
  # Note the gradient is a 3-dimensional function of the parameters theta
  # with dimensions 5 (=length(theta)), size, size.  The first dimension
  # represents the parameters involved in the derivatives
  g.grad = function(theta) {
    abind(array(1, c(size, size)),
          abind(tau, abind(
            tau ^ 2,
            abind(log(tau),
                  array((1:size), c(size, size)), along = 0.5),
            along = 1
          ),
          along = 1),
          along = 1) * aperm(array(g.obj(theta), c(size, size, 5)), c(3, 1, 2))
  }

  # Hessian of g
  # Note the Hessian is a 4-dimensional function of the parameters theta
  # with dimensions 5 (=length(theta)), 5, size, size.  First two dimensions
  # represent the parameters involved in the partial derivatives
  g.hess = function(theta)  {
    if (length(theta) != 5)
      stop("theta does not equal 5 in hoerl()")
    aa = aperm(array(abind(
      array(1, c(size, size)),
      abind(tau, abind(
        tau ^ 2,
        abind(log(tau),
              array((1:size), c(size, size)), along =
                0.5),
        along = 1
      ),
      along = 1),
      along = 1
    ), c(5, size, size, 5)),
    c(4, 1, 2, 3))
    aa * aperm(aa, c(2, 1, 3, 4)) * aperm(array(g.obj(theta), c(size, size, 5, 5)), c(3, 4, 1, 2))
  }

  # Base starting values on classic chain ladder forecasts and inherent trend
  ptd = ((!ptd == 0) * ptd) + (ptd == 0) * mean(ptd)
  tmp = c((
    colSums(B0[, 2:size] + 0 * B0[, 1:(size - 1)], na.rm = TRUE) /
      colSums(B0[, 1:(size - 1)] + 0 * B0[, 2:size], na.rm = TRUE)
  ),
  1)
  yy = 1 / (cumprod(tmp[(size + 1) - (1:size)]))[(size + 1) - (1:size)]
  xx = yy - c(0, yy[1:(size - 1)])
  ww = t(array(xx, c(size, size)))
  uv = ptd / ((size == rowSums(msk)) + (size > rowSums(msk)) * rowSums(msk *
                                                                     ww))
  tmp = na.omit(data.frame(x = 1:size, y = log(uv)))
  trd = 0.01
  trd = array(coef(lm(tmp$y ~ tmp$x))[2])[1]
  tmp = na.omit(data.frame(
    x1 = c(tau),
    x2 = c(tau ^ 2),
    x3 = log(c(tau)),
    y = c(log(outer(uv, xx)))
  ))
  ccs = array(coef(lm(tmp$y ~ tmp$x1 + tmp$x2 + tmp$x3)))[1:4]
  a0 = c(c(ccs), trd)
  return(list(
    g.obj = g.obj,
    g.grad = g.grad,
    g.hess = g.hess,
    a0 = a0
  ))
}
