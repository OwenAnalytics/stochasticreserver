#' Create list for Kramer Chain Ladder parmaterization model
#' g - Assumed loss emergence model, a function of the parameters a.
#' Note g must be matrix-valued with 10 rows and 10 columns

#' g itself
#' Basic design is for g to be a function of a single parameter vector, however
#' in the simulations it is necessary to work on a matrix of parameters, one
#' row for each simulated parameter, so g.obj must be flexible enough to handle
#' both.
#' Here g.obj is nonlinear and based on the Kramer Chain Ladder parmaterization
#' @param tau do not know
#' @param B0 development triangle
#' @param ptd do not know
#' @param msk mask for triangle
#'
#' @importFrom stats coef lm na.omit
#' @import abind
#' @export
capecod <- function(tau, B0, ptd, msk) {
  g.obj = function(theta) {
    if (is.vector(theta))
    {
      theta[1] * outer((c(1, theta[2:10])), c(1, theta[11:19]))
    }
    else
    {
      array(theta[, 1], c(nrow(theta), 10, 10)) *
        array(cbind(1, theta[, 2:10]), c(nrow(theta), 10, 10)) *
        aperm(array(cbind(1, theta[, 11:19]), c(nrow(theta), 10, 10)), c(1, 3, 2))
    }
  }

  # Gradient of g
  # Note the gradient is a 3-dimensional function of the parameters theta
  # with dimensions 19 (=length(theta)), 10, 10.  The first dimension
  # represents the parameters involved in the derivatives
  g.grad = function(theta) {
    abind(outer(c(1, theta[2:10]), c(1, theta[11:19])),
          theta[1] *
            abind(
              aperm(array(t(
                outer((2:10), (1:10), "==")
              ), c(10, 9, 10)), c(2, 1, 3)) *
                aperm(array(c(1, theta[11:19]), c(10, 10, 9)), c(3, 2, 1)),
              aperm(array(t(
                outer((2:10), (1:10), "==")
              ), c(10, 9, 10)), c(2, 3, 1)) *
                aperm(array(c(1, theta[2:10]), c(10, 10, 9)), c(3, 1, 2)),
              along = 1
            ),
          along = 1)
  }

  # Hessian of g
  # Note the Hessian is a 4-dimensional function of the parameters theta
  # with dimensions 19 (=length(theta)), 19, 10, 10.  First two dimensions
  # represent the parameters involved in the partial derivatives
  g.hess = function(theta)  {
    a1 = abind(aperm(array(t(
      outer((2:10), (1:10), "==")
    ), c(10, 9, 10)), c(2, 1, 3)) *
      aperm(array(c(1, theta[11:19]), c(10, 10, 9)), c(3, 2, 1)),
    aperm(array(t(
      outer((2:10), (1:10), "==")
    ), c(10, 9, 10)), c(2, 3, 1)) *
      aperm(array(c(1, theta[2:10]), c(10, 10, 9)), c(3, 1, 2)),
    along = 1)
    a2 = theta[1] * (aperm(array(1:10, c(10, 10, 9, 9)), c(3, 4, 1, 2)) ==
                       array(2:10, c(9, 9, 10, 10))) *
      (aperm(array(1:10, c(10, 10, 9, 9)), c(4, 3, 2, 1)) ==
         aperm(array(2:10, c(9, 9, 10, 10)), c(2, 1, 4, 3)))
    abind(abind(array(0, c(10, 10)), a1, along = 1),
          abind(a1,
                abind(
                  abind(array(0, c(9, 9, 10, 10)), a2, along = 2),
                  abind(aperm(a2, c(2, 1, 3, 4)), array(0, c(9, 9, 10, 10)), along =
                          2),
                  along = 1
                ), along = 1),
          along = 2)
  }

  # Set up starting values based on development factors for columns and
  # relative sizes for rows
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
  a0 = c((uv[1] * xx[1]), (uv[2:10] / uv[1]), (xx[2:10] / xx[1]))
  return(list(
    g.obj = g.obj,
    g.grad = g.grad,
    g.hess = g.hess,
    a0 = a0
  ))
}
