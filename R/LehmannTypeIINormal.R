##' The Lehmann Type II Normal Distribution
##' @description Density, distribution function, quantile function and random
##' generation for the Lehmann Type II Normal distribution with location
##' parameter equal to xi, scale parameter equal to eta and power parameter
##' equal to alpha.

##' @usage dl2norm(x, xi = 0, eta = 1, alpha, log = FALSE)

##' @param x vector of quantiles.
##' @param xi vector of location parameters.
##' @param eta vector of scale parameters.
##' @param alpha vector of power parameters.

##' @source Silva, K.T. A. (2021) Distribuições Logística Exponencial-Inversa de Lehmann Tipo I e II.
##' @author Silva, K.T.A.; MONTALVO, G.S.A.

##' @export
dL2norm <- function(x, xi=0, eta = 1, alpha, log = FALSE){
  d <- alpha/eta * (1 - pnorm((x - xi)/eta))^(alpha-1) * dnorm((x-xi)/eta)
  if(log == TRUE) return(log(d))
  d
}
