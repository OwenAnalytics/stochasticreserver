# Hessian of g
# Note the Hessian is a 4-dimensional function of the parameters theta
# with dimensions 11 (=length(theta)), 11, 10, 10.  First two dimensions
# represent the parameters involved in the partial derivatives
#' @export
g_hess <- function(theta)  {
  aa <- aperm(outer(diag(rep(1,10)),
                 array((1:10)*exp((1:10)*theta[11]),c(10,1))),c(1,4,3,2))
  abind(abind(array(0,c(10,10,10,10)),aa,along=2),
        abind(aperm(aa,c(2,1,3,4)),
              array(outer((1:10)^2*exp((1:10)*theta[11]),theta[(1:10)]),c(1,1,10,10)),
              along=2),
        along=1)
}
