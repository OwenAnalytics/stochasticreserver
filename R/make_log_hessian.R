#' Make Hessian of the objective function
#'
#' -   m_expected_value is the expectated value matrix
#' -   m_variance is the matrix of variances
#' -   A, m_expected_value, m_variance all have shape c(size, size)
#' -   The variables _v are copies of the originals to shape c(npar,size,size),
#' paralleling the gradient of g.
#' -   The variables _m are copies of the originals to shape
#' c(npar,npar,size,size), paralleling the hessian of g
#' hessian-of-objective-function
#' @param a do not know
#' @param A do not know
#' @param g_obj objective function
#' @param g_grad gradient function
#' @param g_hess hessian function
#' @export
make_log_hessian <- function(a, A, g_obj, g_grad, g_hess) {
  npar = length(a) - 2
  p = a[npar + 2]
  Av = aperm(array(A, c(size, size, npar)), c(3, 1, 2))
  Am = aperm(array(A, c(size, size, npar, npar)), c(3, 4, 1, 2))
  m_expected_value = g_obj(a[1:npar])
  ev = aperm(array(m_expected_value, c(size, size, npar)), c(3, 1, 2))
  em = aperm(array(m_expected_value, c(size, size, npar, npar)), c(3, 4, 1, 2))
  m_variance = exp(-outer(logd[, 1], rep(a[npar + 1], size), "-")) * (m_expected_value ^ 2) ^ p
  vv = aperm(array(m_variance, c(size, size, npar)), c(3, 1, 2))
  vm = aperm(array(m_variance, c(size, size, npar, npar)), c(3, 4, 1, 2))
  g1 = g_grad(a[1:npar])
  gg = aperm(array(g1, c(npar, size, size, npar)), c(4, 1, 2, 3))
  gg = gg * aperm(gg, c(2, 1, 3, 4))
  gh = g_hess(a[1:npar])
  dtt = rowSums(
    gh * (p / em + (em - Am) / vm - p * (Am - em) ^ 2 / (vm * em)) +
      gg * (
        1 / vm + 4 * p * (Am - em) / (vm * em) + p * (2 * p + 1) * (Am - em) ^ 2 /
          (vm * em ^ 2) - p / em ^ 2
      ),
    dims = 2,
    na.rm = TRUE
  )
  dkt = rowSums((g1 * (Av - ev) + p * g1 * (Av - ev) ^ 2 / ev) / vv, na.rm = TRUE)
  dtp = rowSums(g1 * (1 / ev + (
    log(ev ^ 2) * (Av - ev) + (p * log(ev ^ 2) - 1) * (Av - ev) ^ 2 / ev
  ) / vv),
  na.rm = TRUE)
  dkk = sum((A - m_expected_value) ^ 2 / (2 * m_variance), na.rm = TRUE)
  dpk = sum(log(m_expected_value ^ 2) * (A - m_expected_value) ^ 2 / (2 * m_variance), na.rm = TRUE)
  dpp = sum(log(m_expected_value ^ 2) ^ 2 * (A - m_expected_value) ^ 2 / (2 * m_variance), na.rm = TRUE)
  m1 = rbind(array(dkt), c(dtp))
  rbind(cbind(dtt, t(m1)), cbind(m1, rbind(cbind(dkk, c(
    dpk
  )), c(dpk, dpp))))
}
