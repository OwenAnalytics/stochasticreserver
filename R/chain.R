#if(model=="Chain") {
#
# g itself
# Basic design is for g to be a function of a single parameter vector, however
# in the simulations it is necessary to work on a matrix of parameters, one
# row for each simulated parameter, so g.obj must be flexible enough to handle
# both.
# Here g.obj is a version of the Cape Cod model but with the restriction
# that the expected cumulative averge payments to date equal the actual
# average paid to date.  Because of this restriction the incremental averages
# are expressed as a percentage times the expected ultimate by row.
# Formulae all assume a full, square development triangle.
#' @export
chain <- function(tau, B0, ptd, msk) {
  g.obj = function(theta) {
    if (is.vector(theta))
    {
      th = t(array(c(theta[1:9], (
        1 - sum(theta[1:9])
      )), c(10, 10)))
      uv = ptd / ((10 == rowSums(msk)) + (10 > rowSums(msk)) * rowSums(msk *
                                                                         th))
      th * array(uv, c(10, 10))
    }
    else
    {
      th = aperm(array(cbind(theta[, 1:9], (
        1 - rowSums(theta[, 1:9])
      )),
      c(nrow(theta), 10, 10)), c(1, 3, 2))
      mska = aperm(array(msk, c(10, 10, nrow(theta))), c(3, 1, 2))
      ptda = t(array(ptd, c(10, nrow(theta))))
      uva = ptda / ((10 == rowSums(mska, dims = 2))
                    + (10 > rowSums(mska, dims = 2)) * rowSums(mska * th, dims =
                                                                 2))
      th * array(uva, c(nrow(theta), 10, 10))
    }
  }

  v1 = aperm(array((1:10), c(10, 10, 9)),
             c(3, 2, 1)) ==
    array((1:9), c(9, 10, 10))
  v2 = aperm(array(msk[, 1:9], c(10, 9, 10)),
             c(2, 1, 3)) &
    aperm(array(msk, c(10, 10, 9)),
          c(3, 1, 2))
  v2[, 1, ] = FALSE
  rsm = rowSums(msk)

  # Gradient of g
  # Note the gradient is a 3-dimensional function of the parameters theta
  # with dimensions 9 (=length(theta)), 10, 10.  The first dimension
  # represents the parameters involved in the derivatives
  g.grad = function(theta) {
    th = t(array(c(theta, (1 - sum(
      theta
    ))), c(10, 10)))
    psm = rowSums(msk * th)
    psc = rowSums(th[, 1:9] * (1 - msk[, 1:9]))
    uv = ptd / ((10 == rsm) + (10 > rsm) * psm)
    uva = aperm(array(uv, c(10, 10, 9)), c(3, 1, 2))
    thj = aperm(array(outer((1 / psm), c(
      theta, 1 - sum(theta)
    )),
    c(10, 10, 9)),
    c(3, 1, 2))
    xx = uva * (v1 - v2 * thj)
    xx[, , 10] = -uva[, , 1] * ((1 - v2[, , 1]) + v2[, , 1] * t(array((1 -
                                                                         psc) / psm, c(10, 9))))
    xx[, 1, 10] = -uv[1]
    xx
  }

  d1 = array((1:9), c(9, 9, 10, 10))
  d2 = aperm(array((1:9), c(9, 9, 10, 10)),
             c(2, 1, 3, 4))
  d3 = aperm(array((1:10), c(10, 9, 9, 10)),
             c(2, 3, 1, 4))
  d4 = aperm(array((1:10), c(10, 9, 9, 10)),
             c(2, 3, 4, 1))
  rsma = aperm(array(rsm, c(10, 9, 9, 10)),
               c(2, 3, 1, 4))
  mm1 = (d1 <= rsma) & (d2 <= rsma)
  mm2 = ((d4 == d1) & (d2 <= rsma)) | ((d4 == d2) & (d1 <= rsma))
  mm3 = (d1 == d4) & (d2 == d4) & (d1 <= rsma) & (d2 <= rsma)
  mm5 = (((d1 > rsma) &
            (d2 <= rsma)) | ((d2 > rsma) & (d1 <= rsma))) & !((d1 > rsma) &
                                                                (d2 > rsma))

  # Hessian of g
  # Note the Hessian is a 4-dimensional function of the parameters theta
  # with dimensions 9 (=length(theta)), 9, 10, 10.  First two dimensions
  # represent the parameters involved in the partial derivatives
  g.hess = function(theta)  {
    th = t(array(c(theta, (1 - sum(
      theta
    ))), c(10, 10)))
    psm = rowSums(msk * th)
    psc = rowSums(th[, 1:9] * (1 - msk[, 1:9]))
    uv = ptd / (((10 == rowSums(msk)) + (10 > rowSums(msk)) * psm) ^ 2)
    psc = rowSums(th[, 1:9] * (1 - msk[, 1:9]))
    uva = aperm(array(uv, c(10, 10, 9, 9)), c(3, 4, 1, 2))
    thj = aperm(array(outer((1 / psm), 2 * c(
      theta, 1 - sum(theta)
    )),
    c(10, 10, 9, 9)),
    c(3, 4, 1, 2))
    xx = uva * (thj * mm1 - mm3 - mm2)
    xx[, , , 10] = uva[, , , 1] * (mm5[, , , 1] + mm1[, , , 1] *
                                     aperm(array(2 * (1 - psc) / psm, c(10, 9, 9)),
                                           c(2, 3, 1)))
    xx[, , 1, ] = 0
    xx
  }

  # Starting values essentially classical chain ladder

  tmp = c((
    colSums(B0[, 2:10] + 0 * B0[, 1:9], na.rm = TRUE) /
      colSums(B0[, 1:9] + 0 * B0[, 2:10], na.rm = TRUE)
  ),
  1)
  yy = 1 / (cumprod(tmp[11 - (1:10)]))[11 - (1:10)]
  a0 = (yy - c(0, yy[1:9]))[1:9]
  return(list(
    g.obj = g.obj,
    g.grad = g.grad,
    g.hess = g.hess,
    a0 = a0
  ))
}
