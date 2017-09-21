
lsm
---

An implementation of the Estimation of the log likelihood of the saturated model.

This package calculates the estimation of the log likelihood of the saturated model, when the values of the outcome variable are either 0 or 1.

Details
-------

The saturated model is characterized by the assumptions 1 and 2 presented in section 5 by Llinas \[1\].

References
----------

\[1\] Humberto Jesús Llinás. (2006). Accuracies in the theory of the logistic models.Revista Colombiana De Estadistica,29(2), 242-244.

\[2\] Hosmer, D. (2013). Wiley Series in Probability and Statistics Ser. : Applied Logistic Regression (3). New York: John Wiley & Sons, Incorporated.

Author(s)
------

Humberto Llinas [aut], Universidad del Norte, Barranquilla-Colombia \\ Jorge Villalba [aut,cre], Unicolombo, Cartagena-Colombia \\ Omar Fabregas [aut], Universidad del Norte, Barranquilla-Colombia.

Installation
------------

``` r
library(devtools)
install_github("jlvia1191/lsm")

```

Example Usage
-------------

``` r
 x1 <- c(68, 72, 68, 76, 69, 71, 68, 61, 69, 68)
 x2 <- c(0.00, 55.90, 0.00, 20.00, 55.90, 0.00, 27.20, 24.00, 0.00, 27.20)
 y <- c (0, 1, 0, 0, 1, 0, 0, 1, 0, 1)
 data <- data.frame (y, x1, x2)
 library(lsm)
 lsm(y~x1+x2, data)
#> [1] -1.386294
```
