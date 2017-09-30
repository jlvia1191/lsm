lsm
---

When the values of the outcome variable Y are either 0 or 1, the function lsm() calculates the estimation of the log likelihood in the saturated model. This model is characterized by Llinas (2006, ISSN:2389-8976) in section 2.3 through the assumptions 1 and 2. The function LogLik() works (almost perfectly) when the number of independent variables K is high, but for small K it calculates wrong values in some cases. For this reason, when Y is dichotomous and the data are grouped or ungrouped, it is recommended the function lsm() because it works very well for all K.

Details
-------

The saturated model is characterized by the assumptions 1 and 2 presented in section 2.3 by Llinas (2006, ISSN:2389-8976).

References
----------

\[1\] Humberto Jesus Llinas. (2006). Accuracies in the theory of the logistic models. Revista Colombiana De Estadistica,29(2), 242-244.

\[2\] Hosmer, D. (2013). Wiley Series in Probability and Statistics Ser. : Applied Logistic Regression (3). New York: John Wiley & Sons, Incorporated.

Author(s)
---------

Humberto Llinas Solano \[aut\], Universidad del Norte, Barranquilla-Colombia \\ Omar Fabregas Cera \[aut\], Universidad del Norte, Barranquilla-Colombia \\ Jorge Villalba Acevedo \[cre, aut\], Unicolombo, Cartagena-Colombia.

Installation
------------

``` r
library(devtools)
install_github("jlvia1191/lsm")

```

Example Usage
-------------

``` r
 library(lsm)

 x1 <- c(68, 72, 68, 76, 69, 71, 68, 61, 69, 68)
 x2 <- c(0.00, 55.90, 0.00, 20.00, 55.90, 0.00, 27.20, 24.00, 0.00, 27.20)
 y <- c (0, 1, 0, 0, 1, 0, 0, 1, 0, 1)
 data <- data.frame (y, x1, x2)

 lsm(y ~ x1 + x2, data)
 $log_Likelihood
 [1] -1.3863
 
 $populations
 [1] 8
 
 attr(,"class")
 [1] "lsm"
```

\# Other example.

``` r
   y <- c(1,    0, 1,   0,  1,  1,  1,  1,  0,  0,  1,  1)
  x1 <- c(2, 2, 2,  5,  5,  5,  5,  8,  8,  11, 11, 11)
  x2 <- c(3, 3, 3,  6,  6,  6,  6,  9,  9,  12, 12, 12)
  x3 <- c(4, 4, 4,  7,  7,  7,  7,  10, 10, 13, 13, 13)
  x5 <- c(1, 1, 1,  9,  9,  9,  9,  4,  4,  2,  2,  2)
  x4 <- c(5, 5, 5,  6,  6,  6,  6,  7,  7,  8,  8,  8)
  data <- data.frame (y, x1, x2, x3, x4, x5) 
  lsm(y ~ x1 + x2 + x3 + x4 + x5, data)
 $log_Likelihood
 [1] -7.4547
 
 $populations
 [1] 4
 
 attr(,"class")
 [1] "lsm"
  
## For more ease, use the following notation.
  
  lsm(y~., data)
 $log_Likelihood
 [1] -7.4547
 
 $populations
 [1] 4
 
 attr(,"class")
 [1] "lsm"
```
