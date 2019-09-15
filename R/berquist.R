#' Create list for Berquist-Sherman incremental severity model
#'
#' g - Assumed loss emergence model, a function of the parameters a.
#' Note g must be matrix-valued with size rows and size columns
#'
#' g itself
#' Basic design is for g to be a function of a single parameter vector, however
#' in the simulations it is necessary to work on a matrix of parameters, one
#' row for each simulated parameter, so g_obj must be flexible enough to handle
#' both.
#' Here g_obj is the Berquist-Sherman incremental severity model
#' @param B0 development triangle
#' @param paid_to_date numeric vector of length \code{size}. It is the lower diagnal of
#' the development triangle in row order. It represents the amount paid to date.
#' @param msk is a mask matrix of allowable data, upper triangular assuming same
#' development increments as exposure increments
#'
#' @importFrom stats coef lm na.omit
#' @import abind
#' @export
berquist <- function(B0, paid_to_date, msk) {
  size <- nrow(B0)
  g_obj = function(theta) {
    if (is.vector(theta))
    {
      outer(exp(theta[(size + 1)] * (1:size)), theta[1:size])
    }
    else
    {
      array(exp(outer(theta[, (size + 1)], (1:size))), c(nrow(theta), size, size)) *
        aperm(array(theta[, (1:size)], c(nrow(theta), size, size)), c(1, 3, 2))
    }
  }

  # Gradient of g
  # Note the gradient is a 3-dimensional function of the parameters theta
  # with dimensions (size + 1) (=length(theta)), size, size.  The first dimension
  # represents the parameters involved in the derivatives
  g_grad = function(theta) {
    if (length(theta) != (size + 1))
      stop("theta is not equal to (size - 1) in berquist()")
    abind(aperm(array(rep(
      exp(theta[(size + 1)] * (1:size)), size * size * size
    ), c(size, size, size)), c(2, 1, 3)) *
      outer((1:size), outer(rep(1, size), (1:size)), "=="),
    outer((1:size) * exp(theta[(size + 1)] * (1:size)), theta[1:size]),
    along = 1)
  }

  # Hessian of g
  # Note the Hessian is a 4-dimensional function of the parameters theta
  # with dimensions (size + 1) (=length(theta)), (size + 1), size, size.  First two dimensions
  # represent the parameters involved in the partial derivatives
  g_hess = function(theta)  {
    if (length(theta) != (size + 1))
      stop("theta is not equal to (size - 1) in berquist()")
    aa = aperm(outer(diag(rep(1, size)),
                     array((1:size) * exp((
                       1:size
                     ) * theta[(size + 1)]), c(size, 1))), c(1, 4, 3, 2))
    abind(abind(array(0, c(size, size, size, size)), aa, along = 2),
          abind(aperm(aa, c(2, 1, 3, 4)),
                array(outer((1:size) ^ 2 * exp((1:size) * theta[(size + 1)]), theta[(1:size)]
                ), c(1, 1, size, size)),
                along = 2),
          along = 1)
  }

  # Set up starting values.  Essentially start with classical chain ladder
  # ultimate estimates and estimate trend from that and incremental average
  # start values based on incrementals from classic chain ladder
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
  tmp = na.omit(data.frame(x = 1:size, y = log(uv)))
  trd = 0.01
  trd = array(coef(lm(tmp$y ~ tmp$x))[2])[1]
  a0 = c((xx * mean(uv / (exp(trd)^(1:size)))), trd)
  return(list(g_obj = g_obj, g_grad = g_grad, g_hess = g_hess, a0 = a0))
}
