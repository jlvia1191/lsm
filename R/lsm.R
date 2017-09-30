# lsm.R

#' @title Estimation of the log Likelihood of the Saturated Model
#' @description When the values of the outcome variable Y are either 0 or 1, the function lsm() calculates the estimation of the log likelihood in the saturated model. This model is characterized by Llinas (2006, ISSN:2389-8976) in section 2.3 through the assumptions 1 and 2. The function LogLik() works (almost perfectly) when the number of independent variables K is high, but for small K it calculates wrong values in some cases. For this reason, when Y is dichotomous and the data are grouped or ungrouped, it is recommended to use function lsm() because it works very well for all K.
#' @param formula An expression of the form y ~ model, where y is the outcome variable (binary or dichotomous: its values are 0 or 1).
#' @param data an optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model. If not found in data, the variables are taken from environment(formula), typically the environment from which lsm() is called.
#' @return  Value of the estimation and  the total of the population.
#' @details The saturated model is characterized by the assumptions 1 and 2 presented in section 2.3 by Llinas (2006, ISSN:2389-8976).
#' @references [1] Humberto Jesus Llinas. (2006). Accuracies in the theory of the logistic models.Revista Colombiana De Estadistica,29(2), 242-244.
#' @references [2] Hosmer, D. (2013). Wiley Series in Probability and Statistics Ser. : Applied Logistic Regression (3). New York: John Wiley &amp; Sons, Incorporated.
#' @author Humberto Llinas Solano [aut], Universidad del Norte, Barranquilla-Colombia \\ Omar Fabregas Cera [aut], Universidad del Norte, Barranquilla-Colombia \\ Jorge Villalba Acevedo [cre, aut], Unicolombo, Cartagena-Colombia.
#' @examples  x1 <- c(68, 72, 68, 76, 69, 71, 68, 61, 69, 68)
#'  x2 <- c(0.00, 55.90, 0.00, 20.00, 55.90, 0.00, 27.20, 24.00, 0.00, 27.20)
#'  y <- c (0, 1, 0, 0, 1, 0, 0, 1, 0, 1)
#'  data <- data.frame (y, x1, x2)
#'  lsm(y ~ x1 + x2, data)
#'
#'
#' # Other example.
#'   y <- c(1, 0, 1, 0, 1, 1, 1, 1, 0, 0, 1, 1)
#'  x1 <-	c(2, 2, 2, 5, 5, 5, 5, 8, 8, 11, 11, 11)
#'  x2 <-	c(3, 3, 3, 6, 6, 6, 6, 9, 9, 12, 12, 12)
#'  x3 <-	c(4, 4, 4, 7, 7, 7, 7, 10, 10, 13, 13, 13)
#'  x5 <-	c(1, 1, 1, 9, 9, 9, 9, 4, 4, 2, 2, 2)
#'  x4 <-	c(5, 5, 5, 6, 6, 6, 6, 7, 7, 8, 8, 8)
#'  data <- data.frame (y, x1, x2, x3, x4, x5) ;data
#'  lsm(y ~ x1 + x2 + x3 + x4 + x5, data)
#'
#' ## For more ease, use the following notation.
#'  lsm(y~., data)
#'
#' @export
#' @import stats


lsm <- function(formula,data){
  L <- as.formula(formula)
  mf <- model.frame(formula = L, data = data)
  res <- do.call(rbind, tapply(mf[, 1], t(apply(t(t(mf[, -1])), 1, paste0,      collapse = "")), function(x) c(p = mean(x), z = sum(x), n = length(x))))
  sj <- (res[, 2]*log(res[, 1]) + (res[, 3] - res[, 2])*log(1 - res[, 1]))
  sat <- round(sum(ifelse ( (res[, 1]) == 0 | (res[, 1]) == 1, 0, sj)), 4)
  na <- list(log_Likelihood = sat, populations = length(res) / 3)
  x <- na
  class(x) <- "lsm"
  x

  print.lsm <- function(x, ...){
    cat("\nlog_Likelihood:\n")
    print(x$"log_Likelihood")
    cat("\npopulations:\n")
    print(x$populations)
  }

}







