#' Make Negative Loglikelihood Function to be Minimized
#'
#' Note that the general
#' form of the model has parameters in addition to those in the loss model,
#' namely the power for the variance and the constant of proprtionality that
#' varies by column.  So if the original model has k parameters with size
#' columns of data, the total objective function has k + size + 1 parameters
#' @param a do not know
#' @param A do not know
#' @param g_obj objective function
#' @export
make_negative_log_likelihood <- function(a, A, g_obj) {
  npar = length(a) - 2
  e = g_obj(a[1:npar])
  v = exp(-outer(logd[, 1], rep(a[npar + 1], size), "-")) * (e^2)^a[npar + 2]
  t1 = log(2 * pi * v) / 2
  t2 = (A - e) ^ 2 / (2 * v)
  sum(t1 + t2, na.rm = TRUE)
}
