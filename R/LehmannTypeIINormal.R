##' @title The Lehmann Type II Normal Distribution
##' @name LehmannTypeIINormal
##' @aliases dL2norm rL2norm
##' @description Density, distribution function, quantile function and random
##' generation for the Lehmann Type II Normal distribution with location
##' parameter equal to xi, scale parameter equal to eta and power parameter
##' equal to alpha.

##' @usage dL2norm(x, xi = 0, eta = 1, alpha, log = FALSE)
##' rL2norm(n, xi = 0, eta = 1, alpha)

##' @param x vector of quantiles.
##' @param n number of observations.
##' @param xi vector of location parameters.
##' @param eta vector of scale parameters.
##' @param alpha vector of power parameters.
##' @param log logical; if TRUE, probabilities p are givem log(p)
##' @return dL2norm gives the density, pL2norm gives the distribution function,
##' qL2norm gives the quantile function, and rL2norm generates random sample.

##' @encoding UTF-8
##' @references  Silva, Keyliane Travassos Almeida da. Distribuições Logística
##' Exponencial-Inversa de Lehmann Tipo I e II.–2021. 46 f. Trabalho de
##' Conclusão de Curso (graduação) – Universidade Federal do Ceará, Centro de
##' Ciências, Curso de Estatística, Fortaleza, 2021. Orientação:
##' Prof. Dr. Gualberto Segundo Agamez Montalvo.

##' @author Silva, K.T.A.; MONTALVO, G.S.A.

##' @import stats

##' @export
# -------------- DENSIDADE --------------------
dL2norm <- function(x, xi=0, eta = 1, alpha, log = FALSE){
  d <- alpha/eta * (1 - pnorm((x - xi)/eta))^(alpha-1) * dnorm((x-xi)/eta)
  if(log == TRUE) return(log(d))
  d
}
# --------------------------------------------

##' @export
# -- FUNCAO GERADORA DE NUMEROS ALEATORIOS ---
rL2norm <- function(n, xi=0, eta = 1, alpha){
  q <- runif(n)
  x <- eta*(qnorm( 1 - (1-q)^(1/alpha)))+xi
  return(x)
}
# --------------------------------------------
