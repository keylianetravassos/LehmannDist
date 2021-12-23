##' @title The Lehmann Type I Logistic Inverse Exponential
##' @name LehmannType1LIE
##' @aliases dllie pllie
##' @description Density, distribution function, quantile function and random
##' generation for the Lehmann Type I Logistic Inverse Exponential distribution with location
##' parameter equal to xi, scale parameter equal to eta and power parameter
##' equal to alpha.

##' @usage dllie(x, beta, lambda, alpha)
##' pllie(x, beta, lambda, alpha)

##' @param x vector of quantiles.
##' @param beta vector of positive shape parameter.
##' @param lambda vector of positive scale parameter.
##' @param alpha vector of positive power parameters.
##' @return dllie gives the density, pllie gives the distribution function,
##' qllie gives the quantile function, and rllie generates random sample.

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
dllie <- function(x, beta, lambda, alpha){
  if(alpha < 0 || beta < 0 || lambda < 0) return(NaN)
  (alpha*beta*lambda)/(x^2)*(exp(lambda/x)*(exp(lambda/x)-1)^(beta-1))/((1 + (exp(lambda/x)-1)^beta)^(1+alpha))
}
# --------------------------------------------

##' @export
# ------------- DISTRIBUICAO ------------------
pllie <- function(x, beta, lambda, alpha){
  if(alpha < 0 || beta < 0 || lambda < 0) return(NaN)
  (1/((1 + (exp(lambda/x) -1)^beta)))^alpha
}

# --------------------------------------------


##' @export
# -- FUNCAO GERADORA DE NUMEROS ALEATORIOS ---

# --------------------------------------------

##' @export
# -------- FUNCAO DE LOG-VEROSSIMILHANCA ----------
