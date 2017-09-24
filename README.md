lsm
---

When the values of the outcome variable Y are either 0 or 1, the function lsm calculates the estimation of the log likelihood in the saturated model. This model is characterized by Llinas \[1\] in section 2.3 through the assumptions 1 and 2. The function LogLik works (almost perfectly) when the mumber of independent variables K is high, but for small K it calculates wrong values. For this reason, when Y is dichotomous and the data are grouped or ungrouped, it is recommended the function lsm because it works very well for all K.

Details
-------

The saturated model is characterized by the assumptions 1 and 2 presented in section 5 by Llinas \[1\].

References
----------

\[1\] Humberto Jesus Llinas. (2006). Accuracies in the theory of the logistic models.Revista Colombiana De Estadistica,29(2), 242-244.

\[2\] Hosmer, D. (2013). Wiley Series in Probability and Statistics Ser. : Applied Logistic Regression (3). New York: John Wiley & Sons, Incorporated.

Author(s)
---------

Humberto Llinas \[aut\], Universidad del Norte, Barranquilla-Colombia \\ Jorge Villalba \[aut, cre\], Unicolombo, Cartagena-Colombia \\ Omar Fabregas \[aut\], Universidad del Norte, Barranquilla-Colombia.

Installation
------------

``` r
library(devtools)
install_github("jlvia1191/lsm")
#> Downloading GitHub repo jlvia1191/lsm@master
#> from URL https://api.github.com/repos/jlvia1191/lsm/zipball/master
#> Installing lsm
#> "C:/PROGRA~1/R/R-33~1.3/bin/x64/R" --no-site-file --no-environ --no-save  \
#>   --no-restore --quiet CMD INSTALL  \
#>   "C:/Users/jorgeR/AppData/Local/Temp/Rtmp6jrM7A/devtoolsa6c62ea2ee2/jlvia1191-lsm-1f29ada"  \
#>   --library="C:/Users/jorgeR/Documents/R/win-library/3.3"  \
#>   --install-tests
#> 
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

``` r
y <- c(1,   0,  1,  0,  1,  1,  1,  1,  0,  0,  1,  1)
  x1 <- c(2,    2,  2,  5,  5,  5,  5,  8,  8,  11, 11, 11)
  x2 <- c(3,    3,  3,  6,  6,  6,  6,  9,  9,  12, 12, 12)
  x3 <- c(4,    4,  4,  7,  7,  7,  7,  10, 10, 13, 13, 13)
  x5 <- c(1,    1,  1,  9,  9,  9,  9,  4,  4,  2,  2,  2)
  x4 <- c(5,    5,  5,  6,  6,  6,  6,  7,  7,  8,  8,  8)
  data <- data.frame (y, x1, x2, x3, x4, x5) ;data
#>    y x1 x2 x3 x4 x5
#> 1  1  2  3  4  5  1
#> 2  0  2  3  4  5  1
#> 3  1  2  3  4  5  1
#> 4  0  5  6  7  6  9
#> 5  1  5  6  7  6  9
#> 6  1  5  6  7  6  9
#> 7  1  5  6  7  6  9
#> 8  1  8  9 10  7  4
#> 9  0  8  9 10  7  4
#> 10 0 11 12 13  8  2
#> 11 1 11 12 13  8  2
#> 12 1 11 12 13  8  2
   library(lsm)
  lsm(y ~ x1+x2+x3+x4+x5, data)
#> [1] -7.45472
```
