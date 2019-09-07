#' Create list for Cape Cod model
#'
#' g itself
#' Basic design is for g to be a function of a single parameter vector, however
#' in the simulations it is necessary to work on a matrix of parameters, one
#' row for each simulated parameter, so g.obj must be flexible enough to handle
#' both.
#' Here g.obj is a version of the Cape Cod model but with the restriction
#' that the expected cumulative averge payments to date equal the actual
#' average paid to date.  Because of this restriction the incremental averages
#' are expressed as a percentage times the expected ultimate by row.
#' Formulae all assume a full, square development triangle.
#' @param tau do not know
#' @param B0 development triangle
#' @param ptd do not know
#' @param msk mask for triangle
#'
#' @importFrom stats coef lm na.omit
#' @import abind
#' @export
chain <- function(tau, B0, ptd, msk) {
  size <- nrow(B0)
  g.obj = function(theta) {
    if (is.vector(theta))
    {
      th = t(array(c(theta[1:(size - 1)], (
        1 - sum(theta[1:(size - 1)])
      )), c(size, size)))
      uv = ptd / ((size == rowSums(msk)) + (size > rowSums(msk)) * rowSums(msk *
                                                                         th))
      th * array(uv, c(size, size))
    }
    else
    {
      th = aperm(array(cbind(theta[, 1:(size - 1)], (
        1 - rowSums(theta[, 1:(size - 1)])
      )),
      c(nrow(theta), size, size)), c(1, 3, 2))
      mska = aperm(array(msk, c(size, size, nrow(theta))), c(3, 1, 2))
      ptda = t(array(ptd, c(size, nrow(theta))))
      uva = ptda / ((size == rowSums(mska, dims = 2))
                    + (size > rowSums(mska, dims = 2)) * rowSums(mska * th, dims =
                                                                 2))
      th * array(uva, c(nrow(theta), size, size))
    }
  }

  v1 = aperm(array((1:size), c(size, size, (size - 1))),
             c(3, 2, 1)) ==
    array((1:(size - 1)), c((size - 1), size, size))
  v2 = aperm(array(msk[, 1:(size - 1)], c(size, (size - 1), size)),
             c(2, 1, 3)) &
    aperm(array(msk, c(size, size, (size - 1))),
          c(3, 1, 2))
  v2[, 1, ] = FALSE
  rsm = rowSums(msk)

  # Gradient of g
  # Note the gradient is a 3-dimensional function of the parameters theta
  # with dimensions (size - 1) (=length(theta)), size, size.  The first dimension
  # represents the parameters involved in the derivatives
  g.grad = function(theta) {
    if (length(theta) != (size - 1))
      stop("theta is not equal to (size - 1) in chain()")
    th = t(array(c(theta, (1 - sum(
      theta
    ))), c(size, size)))
    psm = rowSums(msk * th)
    psc = rowSums(th[, 1:(size - 1)] * (1 - msk[, 1:(size - 1)]))
    uv = ptd / ((size == rsm) + (size > rsm) * psm)
    uva = aperm(array(uv, c(size, size, (size - 1))), c(3, 1, 2))
    thj = aperm(array(outer((1 / psm), c(
      theta, 1 - sum(theta)
    )),
    c(size, size, (size - 1))),
    c(3, 1, 2))
    xx = uva * (v1 - v2 * thj)
    xx[, , size] = -uva[, , 1] * ((1 - v2[, , 1]) + v2[, , 1] * t(array((1 -
                                                                         psc) / psm, c(size, (size - 1)))))
    xx[, 1, size] = -uv[1]
    xx
  }

  d1 = array((1:(size - 1)), c((size - 1), (size - 1), size, size))
  d2 = aperm(array((1:(size - 1)), c((size - 1), (size - 1), size, size)),
             c(2, 1, 3, 4))
  d3 = aperm(array((1:size), c(size, (size - 1), (size - 1), size)),
             c(2, 3, 1, 4))
  d4 = aperm(array((1:size), c(size, (size - 1), (size - 1), size)),
             c(2, 3, 4, 1))
  rsma = aperm(array(rsm, c(size, (size - 1), (size - 1), size)),
               c(2, 3, 1, 4))
  mm1 = (d1 <= rsma) & (d2 <= rsma)
  mm2 = ((d4 == d1) & (d2 <= rsma)) | ((d4 == d2) & (d1 <= rsma))
  mm3 = (d1 == d4) & (d2 == d4) & (d1 <= rsma) & (d2 <= rsma)
  mm5 = (((d1 > rsma) &
            (d2 <= rsma)) | ((d2 > rsma) & (d1 <= rsma))) & !((d1 > rsma) &
                                                                (d2 > rsma))

  # Hessian of g
  # Note the Hessian is a 4-dimensional function of the parameters theta
  # with dimensions (size - 1) (=length(theta)), (size - 1), size, size.  First two dimensions
  # represent the parameters involved in the partial derivatives
  g.hess = function(theta)  {
    if (length(theta) != (size - 1))
      stop("theta is not equal to (size - 1) in chain()")
    th = t(array(c(theta, (1 - sum(theta))), c(size, size)))
    psm = rowSums(msk * th)
    psc = rowSums(th[, 1:(size - 1)] * (1 - msk[, 1:(size - 1)]))
    uv = ptd / (((size == rowSums(msk)) + (size > rowSums(msk)) * psm) ^ 2)
    psc = rowSums(th[, 1:(size - 1)] * (1 - msk[, 1:(size - 1)]))
    uva = aperm(array(uv, c(size, size, (size - 1), (size - 1))), c(3, 4, 1, 2))
    thj = aperm(array(outer((1 / psm), 2 * c(
      theta, 1 - sum(theta)
    )),
    c(size, size, (size - 1), (size - 1))),
    c(3, 4, 1, 2))
    xx = uva * (thj * mm1 - mm3 - mm2)
    xx[, , , size] = uva[, , , 1] *
      (mm5[, , , 1] + mm1[, , , 1] * aperm(array(2 * (1 - psc) / psm,
                                                 c(size, (size - 1), (size - 1))),
                                           c(2, 3, 1)))
    xx[, , 1, ] = 0
    xx
  }

  # Starting values essentially classical chain ladder

  tmp = c((
    colSums(B0[, 2:size] + 0 * B0[, 1:(size - 1)], na.rm = TRUE) /
      colSums(B0[, 1:(size - 1)] + 0 * B0[, 2:size], na.rm = TRUE)
  ),
  1)
  yy = 1 / (cumprod(tmp[(size + 1) - (1:size)]))[(size + 1) - (1:size)]
  a0 = (yy - c(0, yy[1:(size - 1)]))[1:(size - 1)]
  return(list(
    g.obj = g.obj,
    g.grad = g.grad,
    g.hess = g.hess,
    a0 = a0
  ))
}
