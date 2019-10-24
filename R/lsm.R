# lsm.R

#' @title Estimation of the log Likelihood of the Saturated Model
#' @description When the values of the outcome variable \code{Y} are either 0 or 1, the function \code{lsm()} calculates the estimation of the log likelihood in the saturated model. This model is characterized by Llinas (2006, ISSN:2389-8976) in section 2.3 through the assumptions 1 and 2. If \code{Y} is dichotomous and the data are grouped in \code{J} populations, it is recommended to use the function \code{lsm()} because it works very well for all \code{K}.

#' @param formula An expression of the form y ~ model, where y is the outcome variable (binary or dichotomous: its values are 0 or 1).
#' @param family an optional funtion for example binomial.
#' @param data an optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model. If not found in data, the variables are taken from environment(formula), typically the environment from which \code{lsm()} is called.
#' @return  \code{lsm} returns an object of class "\code{lsm}".
#'
#' An object of class "\code{lsm}" is a list containing at least the
#'  following components:
#'
#' \item{log_Likelihood}{Estimation of the log likelihood.}
#'
#' \item{populations}{Total number \code{J} of populations in the model.}
#'
#' \item{z_j}{Value of \code{Zj} (the sum of the observations in the \code{jth} population).}
#'
#' \item{n_j}{Number of the observations in the \code{jth} population.}
#'
#' \item{p_j}{Estimation of \code{pj} in the \code{jth} population.}
#'
#' \item{fitted.values}{Value of the log_Likelihood in the \code{jth} population.}
#'
#' \item{v_j}{Variance of the Bernoulli variables in the \code{jth} population.}
#'
#' \item{m_j}{Expected value of \code{Zj}.}
#'
#' \item{V_j}{Variance of \code{Zj}.}
#'
#' \item{V}{Variance and covariance matrix of \code{Z}, the vector that contains all the \code{Zj}.}
#'
#' \item{S_p}{Score vector of the model.}
#'
#' \item{I_p}{Information matrix of the model.}
#'
#' \item{Zast_j}{Standardized variable of \code{Zj}.}
#'
#' @details The saturated model is characterized by the assumptions 1 and 2 presented in section 2.3 by Llinas (2006, ISSN:2389-8976).
#' @references [1] Humberto Jesus Llinas. (2006). Accuracies in the theory of the logistic models. Revista Colombiana De Estadistica,29(2), 242-244.
#' @references [2] Hosmer, D. (2013). Wiley Series in Probability and Statistics Ser. : Applied Logistic Regression (3). New York: John Wiley & Sons, Incorporated.
#' @author Humberto Llinas Solano [aut], Universidad del Norte, Barranquilla-Colombia \\ Omar Fabregas Cera [aut], Universidad del Norte, Barranquilla-Colombia \\ Jorge Villalba Acevedo [cre, aut], Cartagena-Colombia.
#'@examples
#' # Hosmer, D. (2013) page 3: Age and coranary Heart Disease (CHD) Status of 20 subjects:
#'
#' AGE <- c(20, 23, 24, 25, 25, 26, 26, 28, 28, 29, 30, 30, 30, 30, 30, 30, 30, 32, 33, 33)
#' CHD <- c(0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0)
#'
#'  data <- data.frame (CHD, AGE)
#'  lsm(CHD ~ AGE , family = binomial,  data)
#'
#' # Other case.
#'
#'  y <- c(0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1)
#'  x1 <- c(2, 2, 2, 2, 2, 5, 5, 5, 5, 6, 6, 6, 8, 8, 11, 11, 11, 1)
#'  x2 <- c(3, 3, 3, 3, 3, 6, 6, 6, 6, 8, 8, 8, 9, 9, 12, 12, 12, 12)
#'  x3 <- c(4, 4, 4, 4, 4, 7, 7, 7, 7, 9, 9, 9, 10, 10, 13, 13, 13, 13)
#'  x4 <- c(1, 1, 1, 1, 1, 9, 9, 9, 9, 10, 10, 10, 4, 4, 2, 2, 2, 2)
#'  x5 <- c(32, 32, 32, 32, 32, 20, 20, 20, 20, 21, 21, 21, 19, 19, 16, 16, 16, 16)
#'  x6 <- c(15, 15, 15, 15, 15, 18, 18, 18, 18, 16, 16, 16, 25, 25, 20, 20, 20, 20)
#'  x7 <- c(28, 28, 28, 28, 28, 23, 23, 23, 23, 32, 32, 32, 24, 24, 32, 32, 32, 32)
#'  x8 <- c(0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0)
#'  x9 <- c(6, 6, 6, 6, 6, 10, 10, 10, 10, 11, 11, 11, 7, 7, 21, 21, 21, 21)
#'  x10 <- c(5, 5, 5, 5, 5, 6, 6, 6, 6, 7, 7, 7, 7, 7, 8, 8, 8, 8)
#'
#'  data <- data.frame (y, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10)
#'  lsm(y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10, family = binomial, data)
#'
#' ## For more ease, use the following notation.
#'  lsm(y~., family = binomial, data)
#'
#' ## Other case.
#'
#'   y <- as.factor(c(1, 0, 1, 0, 1, 1, 1, 1, 0, 0, 1, 1))
#'  x1 <- as.factor(c(2, 2, 2, 5, 5, 5, 5, 8, 8, 11, 11, 11))
#'  x2 <- as.factor(c(3, 3, 3, 6, 6, 6, 6, 9, 9, 12, 12, 12))
#'  x3 <- as.factor(c(4, 4, 4, 7, 7, 7, 7, 10, 10, 13, 13, 13))
#'  x4 <- as.factor(c(1, 1, 1, 9, 9, 9, 9, 4, 4, 2, 2, 2))
#'  x5 <- as.factor(c(5, 5, 5, 6, 6, 6, 6, 7, 7, 8, 8, 8))
#'
#'  data <- data.frame (y, x1, x2, x3, x4, x5)
#'  lsm(y ~ x1 + x2 + x3 + x4 + x5, family = binomial, data)
#'
#' ## For more ease, use the following notation.
#'  lsm(y~. , family = binomial, data)
#'
#' @export

lsm <- function(formula, family, data)
{
  mf <- model.frame(formula = formula, data = data)
  res <-do.call(rbind, (tapply(as.vector(mf[, 1]), t(apply((mf[, -1,drop =FALSE]), 1, paste0,collapse = "-")),function(x) c(z = sum(as.numeric(x)), n = length(as.numeric(x)),p = mean(as.numeric(x))))))
  zj<- res[, 1]
  nj <- res[, 2]
  pj <- res[, 3]
  vj <- pj * (1 - pj)
  mj <- nj * pj
  Vj <- nj * vj
  V <- diag(vj)
  sp <- as.matrix((zj - nj * pj)/ vj)
  ip <- diag(nj / vj)
  Zj <- (zj - nj * pj) / sqrt(nj * vj)
  sj <- (res[, 1] * log(res[, 3]) + (res[, 2] - res[, 1]) * log(1 - res[, 3]))
  Lj <-ifelse((res[, 3]) == 0 | (res[, 3]) == 1, 0, sj)
  Sat <- sum (Lj)
  #Lc <- ifelse(Math.factor(mf[, 1]) == 0,log(1-Math.factor(mf[, 1])),log(Math.factor(mf[, 1])))
  Com  <- 0 
  Media <- mean(mf[, 1])
  n <- length(mf[, 1])
  Nul <- n*(Media*log(Media)+(1-Media)*log(1-Media))
  
  coef <- coefficients(glm(formula, family , data))
  B <-as.matrix(coef)
  X <- as.matrix(cbind(1,mf[,-1]))
  Bt <- as.matrix(coef)
  g_t <- X%*%Bt
  p_ <- function(g_){exp(g_t)/(1+exp(g_t))}
  p_t <- p_(g_t)
  q_i <- 1 - p_t
  q_i
  
  Logi <- sum((mf[, 1]*log(p_t)+(1-mf[, 1])*log(q_i)))
  
  I <- p_t * q_i 
  h <-as.vector(I)
  V <- diag(h,length(h))
  
  o <-(t(X)%*% V %*%X)
  varB <-solve(o)
  SEBj <-(diag(varB))^(1/2) 
  W <- t(B)%*%o%*%B
  z <- coef/SEBj
  p_vz =2*pnorm(abs(z), lower.tail=FALSE)
  exb = exp(coef)
  OR = exp(coef[2])
  
  
  #Comparativo de Modelos#
  
  #Logístico vs Completo#
  k <-1
  Dvu <- 2*(Com-Logi)
  gu <- (n-(k+1))
  p_vu <- pchisq(c(Dvu), df=gu, lower.tail=FALSE)
  #Decisi1<-ifelse(p_val1<0.05,"Se rechaza H_0","No se rechaza H_0")
  
  #Nulo vs Logístico#
  Dvd <- 2*(Logi-Nul)
  p_vd <- pchisq(c(Dvd), df=k, lower.tail=FALSE)
  #Decisi<-ifelse(p_val<0.05,"Se rechaza H_0","No se rechaza H_0")
  
  #Logítico vs Saturado#
  J <- length(res) / 3
  Dvt <- 2*(Sat - Logi)
  gt <- (J-(k+1))
  p_vt <- pchisq(c(Dvt), df=gt, lower.tail=FALSE)
  #Decisi2<-ifelse(p_val2<0.05,"Se rechaza H_0","No se rechaza H_0")
  
 
  Ela <-list(coefficients = coef,
             z=z,
             mcov = varB,
             mcor= cor(varB),
             Std_Error = SEBj,
             Wald = W,
             p.valor = p_vz,
             ExpB = exb, 
             obs = n,
             Df = length(coef)-1,
             B = B,
             X = X,
             Bt = Bt,
             populations = J,
             log_lik_LOGT = Logi,
             log_Lik_SAT = Sat,
             log_Lik_COM = Com,
             log_Lik_NUL = Nul,
             Deviuno = Dvu,
             Devidos = Dvd, 
             Devitres = Dvt,
             glu = gu,
             k = k,
             glt = gt,
             p_vuno = p_vu,
             p_vdos = p_vd,
             p_vtres = p_vt,
             logit = g_t,
             z_j = as.matrix(zj), 
             n_j = nj,
             p_j = pj, 
             fitted.values = Lj, 
             v_j = vj, 
             m_j = as.matrix(mj),
             V_j = Vj, 
             V = V,
             S_p = sp, 
             I_p = ip, 
             Zast_j = as.matrix(Zj))
  
  Ela$call <- match.call()
  class(Ela) <- "lsm"
  return(Ela)
}

#' @export
print.lsm <- function(x, ...)
{
 
  TB <- cbind(x$coefficients, x$Std_Error, x$ExpB, x$z, x$p.valor)
  colnames(TB) <- c("coef", "Std. Error", "Exp(B)", "z by Wald", "P.value(>|z|)")
  
  cat("\nCall:\n")
  print(x$call)
 
  cat("\nEstimated Regression Model (Maximum Likelihood) \n")
  printCoefmat(TB, P.values=TRUE,has.Pvalue=TRUE)
  
  cat("\nLikelihood Ratio Test=", x$Devidos, " on ", x$Df, " Df, p.value=", 1 -
        pchisq(x$Devidos, x$Df), ", n=", x$obs, "\n", sep = "")
  
  
  cat("\nLog_Likelihood: \n")
  LL <- cbind(0, x$log_Lik_NUL, x$log_lik_LOGT, x$log_Lik_SAT)
  dimnames(LL) <- list("Estimation",c("Full","Null","Logit","Saturate"))
  print(t(LL))
  
  cat("\nPopulations by saturate model:" ,x$populations,"\n\n", sep = "")
 
}

#' @export
summary.lsm<- function(object, ...)
{
  TAB <- cbind(Deviance = c(object$Devidos,object$Deviuno,object$Devitres),
               Df = c(object$k,object$glu,object$glt),
               P.value = c(object$p_vdos,object$p_vuno,object$p_vtres))
  row.names(TAB) <-c("Null_vs_Logit","Logit_vs_Full", "Logit_vs_Saturate")
  
  res <- list(Call=object$call,
              anova=TAB)
  class(res) <- "summary.lsm"
  return(res)
}

#' @export
print.summary.lsm <- function(x, ...)
{
  cat("\nCall:\n")
  print(x$Call)
  cat("\nAnalysis of Deviance: \n")
  printCoefmat(x$anova, P.values=TRUE,has.Pvalue=TRUE)
}

#' @export
confint.lsm <-  function(object, parm, level =0.95, ...)
{ 
  
  if ((length(level) != 1L) || is.na(level) || (level <= 
  0) || (level >= 1)){
    stop("'conf.level' must be a single number between 0 and 1")
    }
  
  alpha <- 1 - level
  z <- qnorm(1-alpha/2)
  li = object$coefficients - z*object$Std_Error
  ls = object$coefficients + z*object$Std_Error
  ret<-cbind(li, ls)
  colnames(ret)<-c("lower","upper")
  odds<-cbind(exp(li[-1]), exp(ls[-1]))
  colnames(odds)<-c("lower","upper")
  sal <- list(confint=ret,ratios=odds,level*100)
  class(sal) <- "confint.lsm"
  print(sal)
}

#' @export
print.confint.lsm<- function(x, ...)
{
  cat("\n95,0% confidence intervals for odds ratios \n")
  print(x$confint)
  
  cat("\n95,0% confidence intervals for odds ratios \n")
  print(x$ratios)
}


