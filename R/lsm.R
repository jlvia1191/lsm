# lsm.R

#' @title Estimation of the log likelihood of the saturated model
#' @description This package calculates the estimation of the log likelihood of the saturated model, when the values of the outcome variable are either 0 or 1.
#' @param formula An expression of the form y ~ model, where y is the outcome variable (binary or dichotomous: its values are 0 or 1).
#' @param data an optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model. If not found in data, the variables are taken from environment(formula), typically the environment from which dsm is called.
#' @return Value of the estimation.
#' @details The saturated model is characterized by the assumptions 1 and 2 presented in section 5 by Llinas [1].
#' @references [1] Humberto Jesus Llinas. (2006). Accuracies in the theory of the logistic models.Revista Colombiana De Estadistica,29(2), 242-244.
#' @references [2] Hosmer, D. (2013). Wiley Series in Probability and Statistics Ser. : Applied Logistic Regression (3). New York: John Wiley &amp; Sons, Incorporated.
#' @author Humberto Llinas [aut], Universidad del Norte, Barranquilla-Colombia \ Jorge Villalba [cre, aut], Unicolombo, Cartagena-Colombia \ Omar Fabregas [aut], Universidad del Norte, Barranquilla-Colombia.
#' @examples  x1 <- c(68, 72, 68, 76, 69, 71, 68, 61, 69, 68)
#'  x2 <- c(0.00, 55.90, 0.00, 20.00, 55.90, 0.00, 27.20, 24.00, 0.00, 27.20)
#'  y <- c (0, 1, 0, 0, 1, 0, 0, 1, 0, 1)
#'  data <- data.frame (y, x1, x2)
#'  lsm(y~x1+x2, data)
#' @export
#' @import stats
#' @import reshape2

lsm <- function(formula,data){

  mf <- model.frame(formula = formula, data=data)
  z <- dcast(data = mf, formula = formula, sum, value.var = "y")
  n <- dcast(data=mf, formula = formula, fun.aggregate = length, value.var = "y")
  zj <- as.matrix(z[2, ])
  nj <- as.matrix(n[1, ] + n[2, ])
  pj <- zj / nj

  sat <- sum(ifelse ( (pj) == 0 | (pj) == 1, 0, zj * log (pj) + (nj-zj) *log (1-pj) ))

  return(sat)

}





