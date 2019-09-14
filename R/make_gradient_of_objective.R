# Make Gradient of the objective function
#' @param a do not know
#' @param A do not know
#' @param g_obj objective function
#' @param g_grad gradient function
#' @export
make_gradient_of_objective <- function(a, A, g_obj, g_grad) {
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
