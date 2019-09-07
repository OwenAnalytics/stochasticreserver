#' Create list for Kramer Chain Ladder parmaterization model
#' g - Assumed loss emergence model, a function of the parameters a.
#' Note g must be matrix-valued with size rows and size columns

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
  size <- nrow(B0)
  g.obj = function(theta) {
    if (is.vector(theta))
    {
      theta[1] * outer((c(1, theta[2:size])),
                       c(1, theta[(size + 1):((2 * size) - 1)]))
    }
    else
    {
      array(theta[, 1], c(nrow(theta), size, size)) *
        array(cbind(1, theta[, 2:size]), c(nrow(theta), size, size)) *
        aperm(array(cbind(1, theta[, (size + 1):((2 * size) - 1)]),
                    c(nrow(theta), size, size)), c(1, 3, 2))
    }
  }

  # Gradient of g
  # Note the gradient is a 3-dimensional function of the parameters theta
  # with dimensions ((size * 2) - 1) (=length(theta)), size, size.  The first dimension
  # represents the parameters involved in the derivatives
  g.grad = function(theta) {
    abind(outer(c(1, theta[2:size]), c(1, theta[(size + 1):((2 * size) - 1)])),
          theta[1] *
            abind(
              aperm(array(t(
                outer((2:size), (1:size), "==")
              ), c(size, size - 1, size)), c(2, 1, 3)) *
                aperm(array(c(1, theta[(size + 1):((2 * size) - 1)]), c(size, size, size - 1)), c(3, 2, 1)),
              aperm(array(t(
                outer((2:size), (1:size), "==")
              ), c(size, size - 1, size)), c(2, 3, 1)) *
                aperm(array(c(1, theta[2:size]), c(size, size, size - 1)), c(3, 1, 2)),
              along = 1
            ),
          along = 1)
  }

  # Hessian of g
  # Note the Hessian is a 4-dimensional function of the parameters theta
  # with dimensions ((size * 2) - 1) (=length(theta)), ((size * 2) - 1),
  # size, size.  First two dimensions
  # represent the parameters involved in the partial derivatives
  g.hess = function(theta)  {
    a1 = abind(aperm(array(t(
      outer((2:size), (1:size), "==")
    ), c(size, size - 1, size)), c(2, 1, 3)) *
      aperm(array(c(1, theta[(size + 1):((2 * size) - 1)]), c(size, size, size - 1)), c(3, 2, 1)),
    aperm(array(t(
      outer((2:size), (1:size), "==")
    ), c(size, size - 1, size)), c(2, 3, 1)) *
      aperm(array(c(1, theta[2:size]), c(size, size, size - 1)), c(3, 1, 2)),
    along = 1)
    a2 = theta[1] * (aperm(array(1:size, c(size, size, size - 1, size - 1)), c(3, 4, 1, 2)) ==
                       array(2:size, c(size - 1, size - 1, size, size))) *
      (aperm(array(1:size, c(size, size, size - 1, size - 1)), c(4, 3, 2, 1)) ==
         aperm(array(2:size, c(size - 1, size - 1, size, size)), c(2, 1, 4, 3)))
    abind(abind(array(0, c(size, size)), a1, along = 1),
          abind(a1,
                abind(
                  abind(array(0, c(size - 1, size - 1, size, size)), a2, along = 2),
                  abind(aperm(a2, c(2, 1, 3, 4)), array(0, c(size - 1, size - 1, size, size)), along =
                          2),
                  along = 1
                ), along = 1),
          along = 2)
  }

  # Set up starting values based on development factors for columns and
  # relative sizes for rows
  ptd = ((!ptd == 0) * ptd) + (ptd == 0) * mean(ptd)
  tmp = c((
    colSums(B0[, 2:size] + 0 * B0[, 1:(size - 1)], na.rm = TRUE) /
      colSums(B0[, 1:(size - 1)] + 0 * B0[, 2:size], na.rm = TRUE)
  ),
  1)
  yy = 1 / (cumprod(tmp[11 - (1:size)]))[11 - (1:size)]
  xx = yy - c(0, yy[1:(size - 1)])
  ww = t(array(xx, c(size, size)))
  uv = ptd / ((size == rowSums(msk)) + (size > rowSums(msk)) * rowSums(msk *
                                                                     ww))
  a0 = c((uv[1] * xx[1]), (uv[2:size] / uv[1]), (xx[2:size] / xx[1]))
  return(list(
    g.obj = g.obj,
    g.grad = g.grad,
    g.hess = g.hess,
    a0 = a0
  ))
}
