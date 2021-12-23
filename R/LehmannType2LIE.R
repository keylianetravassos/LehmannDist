##' @title The Lehmann Type II Logistic Inverse Exponential
##' @name LehmannType2LIE
##' @aliases dl2lie pl2lie
##' @description Density, distribution function, quantile function and random
##' generation for the Lehmann Type II Logistic Exponential Inverse distribution with location
##' parameter equal to xi, scale parameter equal to eta and power parameter
##' equal to alpha.

##' @usage dl2lie(x, beta, lambda, alpha)
##' pl2lie(x, beta, lambda, alpha)

##' @param x vector of quantiles.
##' @param beta vector of positive shape parameter.
##' @param lambda vector of positive scale parameter.
##' @param alpha vector of positive power parameters.
##' @return dl2lie gives the density, pl2lie gives the distribution function,
##' ql2lie gives the quantile function, and rl2lie generates random sample.

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
dl2lie <- function(x, beta, lambda, alpha){
  if(alpha < 0 || beta < 0 || lambda < 0) return(NaN)
  (alpha*beta*lambda)/(x^2)*(exp(lambda/x)*(exp(lambda/x)-1)^(beta-1))/((1 + (exp(lambda/x)-1)^beta)^2)*(1 - 1/(1 + (exp(lambda/x)-1)^beta))^(alpha-1)
}
# --------------------------------------------

##' @export
# ------------- DISTRIBUICAO ------------------
pl2lie <- function(x, beta, lambda, alpha){
  if(alpha < 0 || beta < 0 || lambda < 0) return(NaN)
  1 - (1 - 1/(1 + (exp(lambda/x)-1)^beta))^alpha
}

# --------------------------------------------


##' @export
# -- FUNCAO GERADORA DE NUMEROS ALEATORIOS ---

# --------------------------------------------

##' @export
# -------- FUNCAO DE LOG-VEROSSIMILHANCA ----------

