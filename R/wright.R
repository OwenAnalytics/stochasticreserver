#' Create list for Generalized Hoerl Curve with individual accident year levels
#' (Wright's)
#'
#' g itself
#' Basic design is for g to be a function of a single parameter vector, however
#' in the simulations it is necessary to work on a matrix of parameters, one
#' row for each simulated parameter, so g_obj must be flexible enough to handle
#' both.
#' Here g_obj is Wright's operational time model with separate level by
#' accident year
#' @param tau do not know
#' @param B0 development triangle
#' @param paid_to_date numeric vector of length \code{size}. It is the lower diagnal of
#' the development triangle in row order. It represents the amount paid to date.
#' @param msk is a mask matrix of allowable data, upper triangular assuming same
#' development increments as exposure increments
#'
#' @importFrom stats coef lm na.omit
#' @import abind
#' @export
wright <- function(tau, B0, paid_to_date, msk) {
  size <- nrow(B0)
  g_obj = function(theta) {
    if (is.vector(theta))
    {
      exp(array(theta[1:size], c(size, size)) + theta[(size + 1)] * tau +
            theta[(size + 2)] * tau^2 +
            theta[(size + 3)] * log(tau))
    }
    else
    {
      exp(
        array(theta[, 1:size], c(nrow(theta), size, size)) +
          array(theta[, (size + 1)], c(nrow(theta), size, size)) *
          aperm(array(tau, c(
            size, size, nrow(theta)
          )), c(3, 1, 2)) +
          array(theta[, (size + 2)], c(nrow(theta), size, size)) *
          aperm(array(tau ^ 2, c(
            size, size, nrow(theta)
          )), c(3, 1, 2)) +
          array(theta[, (size + 3)], c(nrow(theta), size, size)) *
          aperm(array(log(tau), c(
            size, size, nrow(theta)
          )), c(3, 1, 2))
      )
    }
  }

  # Gradient of g
  # Note the gradient is a 3-dimensional function of the parameters theta
  # with dimensions (size + 3), size, size.  The first dimension
  # represents the parameters involved in the derivatives
  g_grad = function(theta) {
    abind(array(outer(1:size, 1:size, "=="), c(size, size, size)),
          abind(tau, abind(tau ^ 2, log(tau),
                           along = 0.5),
                along = 1),
          along = 1) * aperm(array(g_obj(theta), c(size, size, (size + 3))), c(3, 1, 2))
  }

  # Hessian of g
  # Note the Hessian is a 4-dimensional function of the parameters theta
  # with dimensions (size + 3), (size + 3), size, size.  First two dimensions
  # represent the parameters involved in the partial derivatives
  g_hess = function(theta)  {
    aa = array(abind(array(outer(1:size, 1:size, "=="), c(size, size, size)),
                     abind(
                       tau, abind(tau ^ 2, log(tau),
                                  along = 0.5),
                       along = 1
                     ),
                     along = 1),
               c((size + 3), size, size, (size + 3)))
    aperm(aa, c(4, 1, 2, 3)) * aperm(aa, c(1, 4, 2, 3)) *
      aperm(array(g_obj(theta), c(size, size, (size + 3), (size + 3))), c(3, 4, 1, 2))
  }

  # Base starting values on classic chain ladder forecasts
  paid_to_date = ((!paid_to_date == 0) * paid_to_date) + (paid_to_date == 0) * mean(paid_to_date)
  tmp = c((
    colSums(B0[, 2:size] + 0 * B0[, 1:(size - 1)], na.rm = TRUE) /
      colSums(B0[, 1:(size - 1)] + 0 * B0[, 2:size], na.rm = TRUE)
  ),
  1)
  yy = 1 / (cumprod(tmp[(size + 1) - (1:size)]))[(size + 1) - (1:size)]
  xx = yy - c(0, yy[1:(size - 1)])
  ww = t(array(xx, c(size, size)))
  uv = paid_to_date / ((size == rowSums(msk)) + (size > rowSums(msk)) * rowSums(msk *
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
    g_obj = g_obj,
    g_grad = g_grad,
    g_hess = g_hess,
    a0 = a0
  ))
}
